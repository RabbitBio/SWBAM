#include "swbam_mpi.h"
#include "markdup-pending-index.h"
#include "markdup-boundary-progress.h"
#include "markdup-halo-plan.h"
#include "bounded-cache.h"
#include "ordered-body.h"
#include "swbam/mpi_runtime.h"
#include "swbam/operators/bam_read_batch.h"
#include "swbam/cpe_write_pipeline.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <new>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#endif

#include <mpi.h>

extern "C" {
    void slave_mpi_decompress_bam2bam_passthrough();
    void slave_mpi_markdup_extract();
    void slave_mpi_compressfunc();
    void slave_mpi_sort_compress_payload();
}

int MpiCommonInitMemWriter(MemWriter &w, size_t cap);
int MpiCommonBuildBamHeaderMemory(
    sam_hdr_t *hdr, int compress_level,
    char **data, size_t *size);
int MpiMarkdupFindDuplicatesStreamingHash(
    std::vector<MpiMarkdupCandidateShared> *owner_candidates,
    const std::vector<unsigned char> &owner_qnames,
    int comm_size,
    std::vector<std::vector<uint64_t> > *duplicates_by_source,
    MpiMarkdupStats *stats);

namespace {

const int kMarkdupNB = 64;
const int kMarkdupTagCount = 8;
const int kMarkdupExchangeChunk = INT_MAX / 2;

enum MdOnePassRecordState {
    kMdOnePassKeep = 0,
    kMdOnePassUnknown = 1,
    kMdOnePassDuplicate = 2
};

const unsigned char kMarkdupBgzfEofBlock[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

static int MdEnvFlagEnabled(const char *name) {
    const char *value = getenv(name);
    return value && value[0] && value[0] != '0';
}

struct MdBlockSet {
    bam_block *blocks;
    unsigned char *data;
    int n;
};

struct MdRecordSet {
    bam1_t *records;
    unsigned char *data;
    bam1_t **ptrs;
    uint32_t *bam_lens;
    int records_per_block;
    int n_blocks;
    size_t arena_stride;
};





struct MdCandidateWorkspace {
    MpiMarkdupCandidateShared *candidates;
    unsigned char *qnames;
    int candidates_per_block;
    size_t qname_stride;
};

struct MdExtractBatchView {
    MdBlockSet *uncompressed;
    Bam2BamPara *decomp;
    int n_blocks;
    uint64_t ordinal_base;
    uint64_t ordinal_end;
    size_t fixed_workspace_bytes;
};

struct MdMarkdupCandidateBufferDeleter {
    void operator()(MpiMarkdupCandidateShared *p) const {
        if (p) aligned_free_custom((unsigned char *)p);
    }
};



struct MdByteBufferDeleter {
    void operator()(unsigned char *p) const {
        if (p) aligned_free_custom(p);
    }
};



struct MdCoordBoundary {
    int has_records;
    int32_t first_ref;
    int64_t first_coord;
    int32_t last_ref;
    int64_t last_coord;
};

struct MdRemoteCandidateRef {
    size_t index;
    int owner;
};

struct MdBoundaryCandidateId {
    uint64_t global_order;
    uint8_t kind;

    bool operator==(const MdBoundaryCandidateId &other) const {
        return global_order == other.global_order && kind == other.kind;
    }
};

struct MdBoundaryCandidateIdHash {
    size_t operator()(const MdBoundaryCandidateId &id) const {
        uint64_t value = id.global_order ^
            ((uint64_t)id.kind * 0x9e3779b97f4a7c15ull);
        value ^= value >> 30;
        value *= 0xbf58476d1ce4e5b9ull;
        value ^= value >> 27;
        value *= 0x94d049bb133111ebull;
        value ^= value >> 31;
        return (size_t)value;
    }
};

struct MdOnePassBoundaryPlan {
    int enabled;
    bool partial;
    bool conservative;
    uint64_t first_probe_end_block;
    uint64_t last_probe_begin_block;
    std::vector<MdCoordBoundary> boundaries;
    std::vector<swbam_markdup::HaloInterval> halos;
    std::vector<MpiMarkdupKeyShared> closed_owner_keys;
    std::vector<swbam_markdup::HaloInterval> closed_owner_ranges;
    std::unordered_set<MdBoundaryCandidateId,
                       MdBoundaryCandidateIdHash> candidates;
    std::unordered_map<uint64_t, uint8_t> decisions;

    MdOnePassBoundaryPlan()
        : enabled(0), partial(false), conservative(false),
          first_probe_end_block(0), last_probe_begin_block(0) {}

    bool MayContainBatch(
            const std::vector<MpiMarkdupCandidateShared> &batch) const {
        if (!enabled || candidates.empty() || batch.empty()) return false;
        // Extraction preserves file order. Only the two probed block ranges
        // can contain precomputed IDs; interior batches need no hash lookups.
        return (batch.front().global_order >> 32) < first_probe_end_block ||
               (batch.back().global_order >> 32) >= last_probe_begin_block;
    }
};

typedef int (*MdRemoteCandidateObserver)(
    void *opaque, const MpiMarkdupCandidateShared &candidate,
    int owner);



static int MdAllRanksOk(int local_ok) {
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    return global_ok;
}

static int MdAllRanksOkTimed(int local_ok, MpiMarkdupStats *stats) {
    double t0 = GetTime();
    const int global_ok = MdAllRanksOk(local_ok);
    stats->t_rank_sync += GetTime() - t0;
    return global_ok;
}

static double MdReduceMax(double local) {
    double result = 0.0;
    MPI_Reduce(&local, &result, 1, MPI_DOUBLE, MPI_MAX, 0,
               MPI_COMM_WORLD);
    return result;
}

static int MdParseMemory(const std::string &text, size_t *value) {
    *value = 0;
    if (text.empty()) return 0;
    errno = 0;
    char *end = nullptr;
    unsigned long long base = strtoull(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str()) return -1;
    while (*end == ' ' || *end == '\t') ++end;
    unsigned long long multiplier = 1;
    if (*end) {
        char suffix = *end++;
        if (suffix == 'K' || suffix == 'k') multiplier = 1024ull;
        else if (suffix == 'M' || suffix == 'm') multiplier = 1024ull * 1024ull;
        else if (suffix == 'G' || suffix == 'g') multiplier = 1024ull * 1024ull * 1024ull;
        else if (suffix == 'T' || suffix == 't') multiplier = 1024ull * 1024ull * 1024ull * 1024ull;
        else return -1;
        if (*end == 'B' || *end == 'b') ++end;
        while (*end == ' ' || *end == '\t') ++end;
        if (*end) return -1;
    }
    if (base > (unsigned long long)SIZE_MAX / multiplier) return -1;
    *value = (size_t)(base * multiplier);
    return 0;
}

static size_t MdVectorBytes(
        const std::vector<MpiMarkdupCandidateShared> &candidates,
        const std::vector<unsigned char> &qnames,
        const std::vector<uint8_t> &bitmap) {
    size_t total = candidates.capacity() * sizeof(MpiMarkdupCandidateShared);
    if (qnames.capacity() > SIZE_MAX - total) return SIZE_MAX;
    total += qnames.capacity();
    if (bitmap.capacity() > SIZE_MAX - total) return SIZE_MAX;
    return total + bitmap.capacity();
}


static int MdCheckMemory(size_t bytes, size_t limit,
                         MpiMarkdupStats *stats, int rank,
                         const char *stage) {
    if ((long long)bytes > stats->tracked_peak_bytes) {
        stats->tracked_peak_bytes =
            bytes > (size_t)LLONG_MAX ? LLONG_MAX : (long long)bytes;
    }
    if (limit != 0 && bytes > limit) {
        fprintf(stderr,
                "[rank %d] ERROR: markdup memory limit exceeded at %s: "
                "required=%zu limit=%zu. Increase -m.\n",
                rank, stage, bytes, limit);
        return -1;
    }
    return 0;
}

static int MdAllocateBlockSet(MdBlockSet *set, int n) {
    memset(set, 0, sizeof(*set));
    set->n = n;
    set->blocks = (bam_block *)aligned_alloc_custom(
        64, (size_t)n * sizeof(bam_block));
    set->data = aligned_alloc_custom(
        64, (size_t)n * BGZF_MAX_BLOCK_SIZE);
    if (!set->blocks || !set->data) return -1;
    memset(set->blocks, 0, (size_t)n * sizeof(bam_block));
    for (int i = 0; i < n; ++i) {
        set->blocks[i].data =
            set->data + (size_t)i * BGZF_MAX_BLOCK_SIZE;
        set->blocks[i].block_id = i;
    }
    return 0;
}

static void MdFreeBlockSet(MdBlockSet *set) {
    if (set->blocks) aligned_free_custom((unsigned char *)set->blocks);
    if (set->data) aligned_free_custom(set->data);
    memset(set, 0, sizeof(*set));
}

static int MdAllocateRecordSet(MdRecordSet *set, int n_blocks,
                               int records_per_block,
                               size_t arena_stride) {
    memset(set, 0, sizeof(*set));
    const size_t total = (size_t)n_blocks * records_per_block;
    set->records_per_block = records_per_block;
    set->n_blocks = n_blocks;
    set->arena_stride = arena_stride;
    set->records = (bam1_t *)aligned_alloc_custom(
        64, total * sizeof(bam1_t));
    set->data = aligned_alloc_custom(
        64, (size_t)n_blocks * arena_stride);
    set->ptrs = (bam1_t **)aligned_alloc_custom(
        64, total * sizeof(bam1_t *));
    set->bam_lens = (uint32_t *)aligned_alloc_custom(
        64, total * sizeof(uint32_t));
    if (!set->records || !set->data || !set->ptrs || !set->bam_lens) {
        return -1;
    }
    memset(set->records, 0, total * sizeof(bam1_t));
    memset(set->bam_lens, 0, total * sizeof(uint32_t));
    for (size_t i = 0; i < total; ++i) {
        set->records[i].mempolicy = BAM_USER_OWNS_DATA;
        set->ptrs[i] = set->records + i;
    }
    return 0;
}

static void MdFreeRecordSet(MdRecordSet *set) {
    if (set->records) aligned_free_custom((unsigned char *)set->records);
    if (set->data) aligned_free_custom(set->data);
    if (set->ptrs) aligned_free_custom((unsigned char *)set->ptrs);
    if (set->bam_lens) aligned_free_custom(
        (unsigned char *)set->bam_lens);
    memset(set, 0, sizeof(*set));
}

static int MdAllocateCandidateWorkspace(MdCandidateWorkspace *workspace,
                                        int records_per_block) {
    memset(workspace, 0, sizeof(*workspace));
    workspace->candidates_per_block = records_per_block * 2;
    workspace->qname_stride = BGZF_MAX_BLOCK_SIZE;
    workspace->candidates =
        (MpiMarkdupCandidateShared *)aligned_alloc_custom(
            64, (size_t)kMarkdupNB *
                    workspace->candidates_per_block *
                    sizeof(MpiMarkdupCandidateShared));
    workspace->qnames = aligned_alloc_custom(
        64, (size_t)kMarkdupNB * workspace->qname_stride);
    return workspace->candidates && workspace->qnames ? 0 : -1;
}

static void MdFreeCandidateWorkspace(MdCandidateWorkspace *workspace) {
    if (workspace->candidates) {
        aligned_free_custom((unsigned char *)workspace->candidates);
    }
    if (workspace->qnames) aligned_free_custom(workspace->qnames);
    memset(workspace, 0, sizeof(*workspace));
}

static int MdMemReadBlock(MemReader &reader, bam_block *block) {
    if (reader.pos >= reader.size ||
        reader.size - reader.pos < BLOCK_HEADER_LENGTH) {
        return -1;
    }
    const unsigned char *base =
        (const unsigned char *)reader.base + reader.pos;
    int length = (int)base[16] | ((int)base[17] << 8);
    length++;
    if (length <= 0 || (size_t)length > reader.size - reader.pos ||
        length > BGZF_MAX_BLOCK_SIZE) {
        return -1;
    }
    memcpy(block->data, base, (size_t)length);
    block->length = (unsigned int)length;
    block->pos = 0;
    block->errcode = 0;
    block->block_address = (int64_t)reader.pos;
    reader.pos += (size_t)length;
    return 0;
}

static void MdInitDecompPara(Bam2BamPara *para, int block_id,
                             MdBlockSet *input, MdBlockSet *uncompressed,
                             MdRecordSet *records, int active) {
    const int records_per_block = records->records_per_block;
    memset(para, 0, sizeof(*para));
    para->block_id = block_id;
    para->output_records =
        records->ptrs + (size_t)block_id * records_per_block;
    para->record_base =
        records->records + (size_t)block_id * records_per_block;
    para->bam_lens =
        records->bam_lens + (size_t)block_id * records_per_block;
    para->record_capacity = records_per_block;
    para->data_arena =
        records->data + (size_t)block_id * records->arena_stride;
    para->data_arena_capacity = records->arena_stride;
    para->filter.min_mapq = -1;
    para->filter.max_mapq = -1;
    para->filter.require_flag = 0;
    para->filter.exclude_flag = 0;
    para->filter.ref_tid = -1;
    para->filter.min_read_len = -1;
    para->filter.max_read_len = -1;
    if (active) {
        para->input_block = input->blocks + block_id;
        para->un_comp_block = uncompressed->blocks + block_id;
        para->status = 0;
    } else {
        para->input_block = nullptr;
        para->un_comp_block = nullptr;
        para->status = -1;
    }
}

static int MdReadGroup(MemReader &reader, MdBlockSet *input,
                       int *n_blocks, MpiMarkdupStats *stats) {
    double t0 = GetTime();
    int count = 0;
    for (int i = 0; i < kMarkdupNB; ++i) {
        if (MdMemReadBlock(reader, input->blocks + i) != 0) break;
        input->blocks[i].block_id = i;
        count++;
    }
    *n_blocks = count;
    stats->t_read += GetTime() - t0;
    return 0;
}

class MdInputPass {
public:
    explicit MdInputPass(MemReader *reader)
        : memory_reader_(reader), backend_(nullptr), spans_(nullptr),
          span_count_(0), position_(0) {}

    MdInputPass(const swbam::BamInputBackend *backend,
                const swbam::BgzfBlockSpan *spans,
                size_t span_count)
        : memory_reader_(nullptr), backend_(backend), spans_(spans),
          span_count_(span_count), position_(0) {}

    bool uses_backend() const { return backend_ != nullptr; }

    void Reset() {
        if (memory_reader_) memory_reader_->pos = 0;
        position_ = 0;
    }

    int Read(MdBlockSet *input, swbam::BgzfBlockBatch *view,
             int *n_blocks, MpiMarkdupStats *stats) {
        if (!input || !n_blocks || !stats) return -1;
        if (!backend_) {
            if (!memory_reader_) return -1;
            return MdReadGroup(*memory_reader_, input, n_blocks, stats);
        }
        if (!view || (span_count_ > 0 && !spans_) ||
            position_ > span_count_) {
            return -1;
        }
        double t0 = GetTime();
        size_t count = span_count_ - position_;
        if (count > view->capacity()) count = view->capacity();
        const int ret = count > 0
            ? backend_->ReadBatch(spans_ + position_, count, view) : 0;
        if (ret == 0) position_ += count;
        *n_blocks = ret == 0 ? (int)count : 0;
        stats->t_read += GetTime() - t0;
        return ret;
    }

private:
    MemReader *memory_reader_;
    const swbam::BamInputBackend *backend_;
    const swbam::BgzfBlockSpan *spans_;
    size_t span_count_;
    size_t position_;
};

static const char *MdExtractErrorText(int status) {
    switch (status) {
    case -10:
        return "missing MC:Z tag; run fixmate -m first";
    case -11:
        return "MC tag has wrong type; run fixmate -m first";
    case -12:
        return "missing ms integer tag; run fixmate -m first";
    case -13:
        return "ms tag has wrong type; run fixmate -m first";
    case -20:
        return "input is not coordinate sorted";
    default:
        return "candidate extraction failed";
    }
}

typedef int (*MdCandidateBatchPostProcessor)(
    void *context,
    std::vector<MpiMarkdupCandidateShared> *candidates,
    std::vector<unsigned char> *qnames,
    int progress_tid, int progress_pos,
    uint64_t records_so_far,
    MdExtractBatchView *batch);

static int MdExtractCandidatesPass(
        MdInputPass &source, long long global_block_begin,
        int rank, int include_fails, size_t memory_limit,
        std::vector<MpiMarkdupCandidateShared> *local_candidates,
        std::vector<unsigned char> *local_qnames,
        uint64_t *local_records,
        int *range_has_records,
        int *range_first_tid, int *range_first_pos,
        int *range_last_tid, int *range_last_pos,
        MpiMarkdupStats *stats,
        MdCandidateBatchPostProcessor batch_post_processor = nullptr,
        void *batch_context = nullptr,
        long long target_batches = 0) {

    // 工作区分配
    const int records_per_block = (int)MPI_RECORDS_PER_BLOCK;
    struct MdExtractSlot {
        MdBlockSet input; // Borrowed runtime descriptors; never freed by the slot.
        MdBlockSet uncompressed;
        MdRecordSet records;
        MdCandidateWorkspace candidate_workspace;
        Bam2BamPara decomp[kMarkdupNB];
        MpiMarkdupExtractPara extract[kMarkdupNB];
        int n_blocks;
        long long block_group_base;
        uint64_t ordinal_base;
        uint64_t ordinal_end;
        int progress_tid;
        int progress_pos;
    };

    MdExtractSlot slots[2];
    memset(slots, 0, sizeof(slots));

    auto free_slot = [](MdExtractSlot *slot) {
        MdFreeBlockSet(&slot->uncompressed);
        MdFreeRecordSet(&slot->records);
        MdFreeCandidateWorkspace(&slot->candidate_workspace);
        slot->n_blocks = 0;
        slot->progress_tid = -1;
        slot->progress_pos = -1;
    };
    auto free_slots = [&]() {
        free_slot(slots + 0);
        free_slot(slots + 1);
    };

    for (int s = 0; s < 2; ++s) {
        if (MdAllocateBlockSet(&slots[s].uncompressed, kMarkdupNB) != 0 ||
            MdAllocateRecordSet(&slots[s].records, kMarkdupNB,
                                records_per_block,
                                MPI_BAM_BLOCK_ARENA_SIZE) != 0 ||
            MdAllocateCandidateWorkspace(
                &slots[s].candidate_workspace,
                records_per_block) != 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: failed to allocate markdup "
                    "candidate pipeline workspace.\n",
                    rank);
            free_slots();
            return -1;
        }
    }

    // 固定内存估算：candidate extract 使用双槽流水线。
    const MdCandidateWorkspace &cw0 = slots[0].candidate_workspace;
    const size_t fixed_workspace_per_slot =
        (size_t)kMarkdupNB * BGZF_MAX_BLOCK_SIZE * 2 +
        (size_t)kMarkdupNB * MPI_BAM_BLOCK_ARENA_SIZE +
        (size_t)kMarkdupNB * records_per_block *
            (sizeof(bam1_t) + sizeof(bam1_t *) + sizeof(uint32_t)) +
        (size_t)kMarkdupNB * cw0.candidates_per_block *
            sizeof(MpiMarkdupCandidateShared) +
        (size_t)kMarkdupNB * cw0.qname_stride;
    size_t fixed_workspace = fixed_workspace_per_slot;
    if (fixed_workspace > SIZE_MAX - fixed_workspace_per_slot) {
        free_slots();
        return -1;
    }
    fixed_workspace += fixed_workspace_per_slot;
    // The runtime also owns decoded double buffers; this operator keeps its
    // own transferable decoded arenas for the one-pass pending window.
    fixed_workspace += (size_t)2 * kMarkdupNB * BGZF_MAX_BLOCK_SIZE;
    if (MdCheckMemory(fixed_workspace, memory_limit, stats, rank,
                      "candidate pipeline workspace") != 0) {
        free_slots();
        return -1;
    }

    // ordinal 是本 rank 内 record 编号，作为索引
    uint64_t ordinal = 0;
    int previous_tid = -1;
    int previous_pos = -1;
    int have_previous = 0;
    long long processed_batches = 0;
    *range_has_records = 0;
    *range_first_tid = -1;
    *range_first_pos = -1;
    *range_last_tid = -1;
    *range_last_pos = -1;

    auto prepare_extract = [&](MdExtractSlot *slot) {
        uint64_t block_ordinal = slot->ordinal_base;
        for (int b = 0; b < kMarkdupNB; ++b) {
            memset(slot->extract + b, 0, sizeof(slot->extract[b]));
            slot->extract[b].block_id = b;
            slot->extract[b].status = b < slot->n_blocks ? 0 : -1;
            if (b >= slot->n_blocks) continue;
            slot->extract[b].records = slot->decomp[b].output_records;
            slot->extract[b].n_records = slot->decomp[b].n_total_records;
            slot->extract[b].ordinal_base = block_ordinal;
            slot->extract[b].global_block_index =
                (uint64_t)(slot->block_group_base + b);
            slot->extract[b].source_rank = rank;
            slot->extract[b].include_fails = include_fails;
            slot->extract[b].candidates =
                slot->candidate_workspace.candidates +
                (size_t)b *
                    slot->candidate_workspace.candidates_per_block;
            slot->extract[b].candidate_capacity =
                slot->candidate_workspace.candidates_per_block;
            slot->extract[b].qname_arena =
                slot->candidate_workspace.qnames +
                (size_t)b * slot->candidate_workspace.qname_stride;
            slot->extract[b].qname_capacity =
                slot->candidate_workspace.qname_stride;
            block_ordinal +=
                (uint64_t)slot->decomp[b].n_total_records;
        }
        slot->ordinal_end = block_ordinal;
    };

    auto launch_extract = [&](MdExtractSlot *slot) {
        __real_athread_spawn((void *)slave_mpi_markdup_extract,
                             slot->extract, 1);
    };

    auto join_extract = [&](MdExtractSlot *slot) {
        (void)slot;
        double cpe_sync_t0 = GetTime();
        athread_join();
        const double sync_dt = GetTime() - cpe_sync_t0;
        stats->t_cpe_sync += sync_dt;
        stats->t_candidate_extract += sync_dt;
    };

    auto process_slot = [&](MdExtractSlot *slot) -> int {
        double status_total = 0.0;
        double merge_total = 0.0;
        int batch_last_tid = previous_tid;
        int batch_last_pos = previous_pos;
        for (int b = 0; b < slot->n_blocks; ++b) {
            double status_t0 = GetTime();
            // 检查从核提取是否成功
            if (slot->extract[b].status != 0) {
                stats->t_extract_status +=
                    status_total + (GetTime() - status_t0);
                fprintf(stderr,
                        "[rank %d] ERROR: %s at global block %lld "
                        "record=%d status=%d.\n",
                        rank,
                        MdExtractErrorText(slot->extract[b].status),
                        slot->block_group_base + b,
                        slot->extract[b].record_index,
                        slot->extract[b].status);
                return -1;
            }

            // 检查 coordinate sorted
            // MdExtractCandidates 内部检查 rank 内 block 顺序
            // MdValidateRankBoundaries 检查 rank 与 rank 之间的坐标顺序
            if (slot->extract[b].has_records) {
                if (!*range_has_records) {
                    *range_first_tid = slot->extract[b].first_tid;
                    *range_first_pos = slot->extract[b].first_pos;
                    *range_has_records = 1;
                }
                if (slot->extract[b].first_tid >= 0 && have_previous &&
                    (slot->extract[b].first_tid < previous_tid ||
                     (slot->extract[b].first_tid == previous_tid &&
                      slot->extract[b].first_pos < previous_pos))) {
                    fprintf(stderr,
                            "[rank %d] ERROR: input is not coordinate sorted "
                            "near global block %lld.\n",
                            rank, slot->block_group_base + b);
                    stats->t_extract_status +=
                        status_total + (GetTime() - status_t0);
                    return -1;
                }
                previous_tid = slot->extract[b].last_tid;
                previous_pos = slot->extract[b].last_pos;
                *range_last_tid = slot->extract[b].last_tid;
                *range_last_pos = slot->extract[b].last_pos;
                if (batch_post_processor) {
                    batch_last_tid = slot->extract[b].last_tid;
                    batch_last_pos = slot->extract[b].last_pos;
                }
                have_previous = 1;
            }
            status_total += GetTime() - status_t0;

            // Append contiguous CPE output, then rebase only pair QNAME offsets.
            double merge_t0 = GetTime();
            const uint64_t qname_base =
                (uint64_t)local_qnames->size();
            if (slot->extract[b].qname_used > 0) {
                local_qnames->insert(
                    local_qnames->end(), slot->extract[b].qname_arena,
                    slot->extract[b].qname_arena +
                        slot->extract[b].qname_used);
            }
            const size_t candidate_base = local_candidates->size();
            const int candidate_count = slot->extract[b].n_candidates;
            if (candidate_count > 0) {
                local_candidates->insert(local_candidates->end(),
                    slot->extract[b].candidates,
                    slot->extract[b].candidates + candidate_count);
            }
            for (size_t i = candidate_base; i < local_candidates->size(); ++i) {
                MpiMarkdupCandidateShared &candidate = (*local_candidates)[i];
                if (!candidate.key.single) {
                    candidate.qname_offset += qname_base;
                }
            }
            stats->total_records += slot->decomp[b].n_total_records;
            stats->examined_records += slot->extract[b].examined;
            stats->excluded_records += slot->extract[b].excluded;
            stats->pair_candidates +=
                slot->extract[b].pair_candidates;
            stats->single_candidates +=
                slot->extract[b].single_candidates;
            merge_total += GetTime() - merge_t0;
        }
        stats->t_extract_status += status_total;
        stats->t_extract_merge += merge_total;
        ordinal = slot->ordinal_end;
        stats->group_count++;

        // 动态内存检查
        double memcheck_t0 = GetTime();
        size_t dynamic = MdVectorBytes(
            *local_candidates, *local_qnames,
            std::vector<uint8_t>());
        if (dynamic == SIZE_MAX ||
            fixed_workspace > SIZE_MAX - dynamic ||
            MdCheckMemory(fixed_workspace + dynamic, memory_limit,
                          stats, rank,
                          "local candidate accumulation") != 0) {
            stats->t_extract_memcheck += GetTime() - memcheck_t0;
            return -1;
        }
        stats->t_extract_memcheck += GetTime() - memcheck_t0;
        slot->progress_tid = batch_last_tid;
        slot->progress_pos = batch_last_pos;
        return 0;
    };

    auto post_process_slot = [&](MdExtractSlot *slot) -> int {
        if (!batch_post_processor) return 0;
        MdExtractBatchView batch = {};
        batch.uncompressed = &slot->uncompressed;
        batch.decomp = slot->decomp;
        batch.n_blocks = slot->n_blocks;
        batch.ordinal_base = slot->ordinal_base;
        batch.ordinal_end = slot->ordinal_end;
        batch.fixed_workspace_bytes = fixed_workspace;
        if (batch_post_processor(batch_context, local_candidates,
                           local_qnames, slot->progress_tid,
                           slot->progress_pos, slot->ordinal_end,
                           &batch) != 0) {
            return -1;
        }
        local_candidates->clear();
        local_qnames->clear();
        processed_batches++;
        return 0;
    };


    class MdCandidateReadSource : public swbam::cpe::CpeReadBatchSource {
    public:
        MdCandidateReadSource(MdInputPass &pass, MpiMarkdupStats *stats)
            : pass_(pass), stats_(stats) {}
        int ReadNext(swbam::BgzfBlockBatch *batch, size_t *count) {
            MdBlockSet view = {batch->blocks(), nullptr, (int)batch->capacity()};
            int n = 0;
            const int ret = pass_.Read(&view, batch, &n, stats_);
            *count = (size_t)n;
            return ret;
        }
    private:
        MdInputPass &pass_;
        MpiMarkdupStats *stats_;
    };

    // Delay the previous batch's MPE merge until the next decode is running.
    // Its record/candidate pools are distinct; the MPI post-processing runs
    // only after that decode joins, since it may itself launch compression.
    class MdCandidateReadOperator : public swbam::operators::BamReadBatchOperator {
    public:
        MdCandidateReadOperator(MdExtractSlot *slots, long long block_begin, int rank,
                     MpiMarkdupStats *stats,
                     std::function<int(MdExtractSlot *)> merge,
                     std::function<int(MdExtractSlot *)> post,
                     std::function<void(MdExtractSlot *)> extract)
            : BamReadBatchOperator({"markdup-candidates",
                  (void *)slave_mpi_decompress_bam2bam_passthrough, kMarkdupNB}),
              slots_(slots), next_block_(block_begin), rank_(rank), stats_(stats),
              merge_(merge), post_(post), extract_(extract), current_(1),
              have_previous_(false), merge_wall_(0.0) {}
        int Initialize() { return 0; }
        void Shutdown() {}
        int Prepare(const swbam::BgzfBlockBatch &input,
                    swbam::BgzfBlockBatch *, size_t count) {
            const uint64_t ordinal_base = have_previous_
                ? slots_[current_].ordinal_end : 0;
            current_ ^= 1;
            MdExtractSlot &slot = slots_[current_];
            slot.input = {const_cast<bam_block *>(input.blocks()), nullptr,
                          (int)input.capacity()};
            slot.n_blocks = (int)count;
            slot.block_group_base = next_block_;
            slot.ordinal_base = slot.ordinal_end = ordinal_base;
            next_block_ += (long long)count;
            stats_->input_blocks += (long long)count;
            for (int b = 0; b < kMarkdupNB; ++b)
                MdInitDecompPara(slot.decomp + b, b, &slot.input,
                    &slot.uncompressed, &slot.records, b < slot.n_blocks);
            return 0;
        }
        void *kernel_arguments() { return slots_[current_].decomp; }
        int DuringKernel() {
            const double t0 = GetTime();
            const int ret = have_previous_ ? merge_(&slots_[current_ ^ 1]) : 0;
            merge_wall_ = GetTime() - t0;
            return ret;
        }
        void ObserveKernel(double wall, size_t) {
            // Like the old loop, exclude overlapped MPE candidate merging.
            const double exposed = std::max(0.0, wall - merge_wall_);
            stats_->t_candidate_decomp += exposed;
            stats_->t_cpe_sync += exposed;
        }
        int Validate(size_t count) const {
            const MdExtractSlot &slot = slots_[current_];
            for (size_t b = 0; b < count; ++b) {
                if (slot.decomp[b].status != 0) {
                    fprintf(stderr,
                        "[rank %d] ERROR: markdup candidate BAM decode failed "
                        "at global block %lld status=%d record=%d.\n",
                        rank_, slot.block_group_base + (long long)b,
                        slot.decomp[b].status, slot.decomp[b].record_index);
                    return -1;
                }
            }
            return 0;
        }
        int PostProcessBatch(size_t, long long *records) {
            if (have_previous_ && post_(&slots_[current_ ^ 1]) != 0) return -1;
            MdExtractSlot &slot = slots_[current_];
            extract_(&slot);
            *records = (long long)(slot.ordinal_end - slot.ordinal_base);
            have_previous_ = true;
            return 0;
        }
        int Finish() {
            if (!have_previous_) return 0;
            MdExtractSlot *last = &slots_[current_];
            return merge_(last) == 0 ? post_(last) : -1;
        }
    private:
        MdExtractSlot *slots_;
        long long next_block_;
        int rank_;
        MpiMarkdupStats *stats_;
        std::function<int(MdExtractSlot *)> merge_, post_;
        std::function<void(MdExtractSlot *)> extract_;
        int current_;
        bool have_previous_;
        double merge_wall_;
    };

    MdCandidateReadSource input(source, stats);
    MdCandidateReadOperator op(slots, global_block_begin, rank, stats,
        process_slot, post_process_slot, [&](MdExtractSlot *slot) {
            prepare_extract(slot);
            launch_extract(slot);
            join_extract(slot);
        });
    swbam::cpe::CpeReadPipelineOptions options;
    options.overlap_input = false;
    int ret = swbam::cpe::RunCpeReadPipeline(&input, &op, nullptr, options);
    while (ret == 0 && batch_post_processor &&
           processed_batches < target_batches) {
        ret = batch_post_processor(batch_context, local_candidates,
            local_qnames, previous_tid, previous_pos, ordinal, nullptr);
        ++processed_batches;
    }
    *local_records = ordinal;
    free_slots();
    return ret;
}


static bool MdKeyLess(const MpiMarkdupKeyShared &a,
                      const MpiMarkdupKeyShared &b) {
    if (a.single != b.single) return a.single < b.single;
    if (a.this_ref != b.this_ref) return a.this_ref < b.this_ref;
    if (a.this_coord != b.this_coord) {
        return a.this_coord < b.this_coord;
    }
    if (a.orientation != b.orientation) {
        return a.orientation < b.orientation;
    }
    if (!a.single) {
        if (a.other_ref != b.other_ref) {
            return a.other_ref < b.other_ref;
        }
        if (a.other_coord != b.other_coord) {
            return a.other_coord < b.other_coord;
        }
        if (a.leftmost != b.leftmost) {
            return a.leftmost < b.leftmost;
        }
    }
    return false;
}

static bool MdKeyEqual(const MpiMarkdupKeyShared &a,
                       const MpiMarkdupKeyShared &b) {
    return !MdKeyLess(a, b) && !MdKeyLess(b, a);
}


static int MdSendrecvBytes(const void *send_data,
                           unsigned long long send_bytes, int dest,
                           void *recv_data,
                           unsigned long long recv_bytes, int src,
                           int tag) {
    unsigned long long sent = 0;
    unsigned long long received = 0;
    while (sent < send_bytes || received < recv_bytes) {
        int send_count = (int)std::min(
            (unsigned long long)kMarkdupExchangeChunk,
            send_bytes - sent);
        int recv_count = (int)std::min(
            (unsigned long long)kMarkdupExchangeChunk,
            recv_bytes - received);
        const unsigned char *send_ptr =
            send_count ? (const unsigned char *)send_data + sent
                       : nullptr;
        unsigned char *recv_ptr =
            recv_count ? (unsigned char *)recv_data + received
                       : nullptr;
        const int send_dest = send_count ? dest : MPI_PROC_NULL;
        const int recv_src = recv_count ? src : MPI_PROC_NULL;
        if (MPI_Sendrecv((void *)send_ptr, send_count, MPI_BYTE,
                         send_dest, tag,
                         recv_ptr, recv_count, MPI_BYTE,
                         recv_src, tag, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE) != MPI_SUCCESS) {
            return -1;
        }
        sent += (unsigned long long)send_count;
        received += (unsigned long long)recv_count;
    }
    return 0;
}

/*
把本 rank 提取出来的 local_candidates
按照 duplicate key 的 hash 分配给 owner rank
通过 MPI 交换
让每个 owner rank 拿到自己负责判断 duplicate 的候选集合
*/


static int MdCompareCoord(int32_t a_ref, int64_t a_coord,
                          int32_t b_ref, int64_t b_coord) {
    if (a_ref != b_ref) return a_ref < b_ref ? -1 : 1;
    if (a_coord != b_coord) return a_coord < b_coord ? -1 : 1;
    return 0;
}

static int MdBuildCoordBoundaries(
        int range_has_records, int range_first_tid, int range_first_pos,
        int range_last_tid, int range_last_pos, int comm_size,
        std::vector<MdCoordBoundary> *boundaries) {
    long long local[5] = {
        (long long)range_has_records,
        range_has_records ? (long long)range_first_tid + 1 : 0,
        range_has_records ? (long long)range_first_pos + 1 : 0,
        range_has_records ? (long long)range_last_tid + 1 : 0,
        range_has_records ? (long long)range_last_pos + 1 : 0
    };
    std::vector<long long> gathered((size_t)comm_size * 5, 0);
    MPI_Allgather(local, 5, MPI_LONG_LONG, gathered.data(), 5,
                  MPI_LONG_LONG, MPI_COMM_WORLD);
    boundaries->assign((size_t)comm_size, MdCoordBoundary());
    for (int r = 0; r < comm_size; ++r) {
        const long long *src = gathered.data() + (size_t)r * 5;
        MdCoordBoundary &b = (*boundaries)[(size_t)r];
        b.has_records = src[0] != 0;
        b.first_ref = (int32_t)src[1];
        b.first_coord = (int64_t)src[2];
        b.last_ref = (int32_t)src[3];
        b.last_coord = (int64_t)src[4];
    }
    return 0;
}

static int MdBuildFirstCoordBoundaries(
        int range_has_records, int range_first_tid, int range_first_pos,
        int comm_size, std::vector<MdCoordBoundary> *boundaries) {
    // Unplaced tail records are preserved, but define no coordinate owner range.
    const int has_coordinate = range_has_records && range_first_tid >= 0 &&
        range_first_pos >= 0;
    long long local[3] = {
        (long long)has_coordinate,
        has_coordinate ? (long long)range_first_tid + 1 : 0,
        has_coordinate ? (long long)range_first_pos + 1 : 0
    };
    std::vector<long long> gathered((size_t)comm_size * 3, 0);
    MPI_Allgather(local, 3, MPI_LONG_LONG, gathered.data(), 3,
                  MPI_LONG_LONG, MPI_COMM_WORLD);
    boundaries->assign((size_t)comm_size, MdCoordBoundary());
    int have_previous = 0;
    int32_t previous_ref = 0;
    int64_t previous_coord = 0;
    int local_ok = 1;
    for (int r = 0; r < comm_size; ++r) {
        const long long *src = gathered.data() + (size_t)r * 3;
        MdCoordBoundary &b = (*boundaries)[(size_t)r];
        b.has_records = src[0] != 0;
        b.first_ref = (int32_t)src[1];
        b.first_coord = (int64_t)src[2];
        b.last_ref = 0;
        b.last_coord = 0;
        if (!b.has_records) continue;
        if (have_previous &&
            MdCompareCoord(b.first_ref, b.first_coord,
                           previous_ref, previous_coord) < 0) {
            local_ok = 0;
        }
        previous_ref = b.first_ref;
        previous_coord = b.first_coord;
        have_previous = 1;
    }
    return MdAllRanksOk(local_ok) ? 0 : -1;
}

static int MdCoordinateOwner(const MpiMarkdupKeyShared &key,
                             const std::vector<MdCoordBoundary> &boundaries,
                             int fallback_rank) {
    int last_nonempty = fallback_rank;
    const int32_t ref = key.this_ref;
    const int64_t coord = key.this_coord;
    for (size_t r = 0; r < boundaries.size(); ++r) {
        const MdCoordBoundary &b = boundaries[r];
        if (!b.has_records) continue;
        last_nonempty = (int)r;
        if (MdCompareCoord(ref, coord, b.last_ref, b.last_coord) <= 0) {
            return (int)r;
        }
    }
    return last_nonempty;
}

static int MdCoordinateOwnerFast(const MpiMarkdupKeyShared &key,
                                 const std::vector<MdCoordBoundary> &boundaries,
                                 int rank) {
    if (rank < 0 || rank >= (int)boundaries.size()) {
        return MdCoordinateOwner(key, boundaries, rank);
    }

    const MdCoordBoundary &local = boundaries[(size_t)rank];
    const int32_t ref = key.this_ref;
    const int64_t coord = key.this_coord;
    if (local.has_records) {
        const int cmp_first =
            MdCompareCoord(ref, coord, local.first_ref,
                           local.first_coord);
        const int cmp_last =
            MdCompareCoord(ref, coord, local.last_ref,
                           local.last_coord);
        if (cmp_first >= 0 && cmp_last <= 0) {
            return rank;
        }

        if (cmp_last > 0) {
            int last_nonempty = rank;
            for (size_t r = (size_t)rank + 1;
                 r < boundaries.size(); ++r) {
                const MdCoordBoundary &b = boundaries[r];
                if (!b.has_records) continue;
                last_nonempty = (int)r;
                if (MdCompareCoord(ref, coord, b.last_ref,
                                   b.last_coord) <= 0) {
                    return (int)r;
                }
            }
            return last_nonempty;
        }
    }

    return MdCoordinateOwner(key, boundaries, rank);
}

static int MdCoordinateOwnerByFirst(
        const MpiMarkdupKeyShared &key,
        const std::vector<MdCoordBoundary> &boundaries,
        int fallback_rank, int next_rank = -2) {
    // First-coordinate cuts assign every key to one stable rank without
    // needing a full-file pass to discover each rank's final coordinate.
    // Almost every candidate belongs to its source rank. Check that interval
    // directly before falling back to the full boundary walk.
    if (fallback_rank >= 0 &&
        fallback_rank < (int)boundaries.size()) {
        const MdCoordBoundary &local =
            boundaries[(size_t)fallback_rank];
        if (local.has_records &&
            MdCompareCoord(key.this_ref, key.this_coord,
                           local.first_ref, local.first_coord) >= 0) {
            if (next_rank == -2) {
                next_rank = -1;
                for (size_t r = (size_t)fallback_rank + 1;
                     r < boundaries.size(); ++r) {
                    if (boundaries[r].has_records) {
                        next_rank = (int)r;
                        break;
                    }
                }
            }
            if (next_rank < 0 ||
                MdCompareCoord(
                    key.this_ref, key.this_coord,
                    boundaries[(size_t)next_rank].first_ref,
                    boundaries[(size_t)next_rank].first_coord) < 0) {
                return fallback_rank;
            }
        }
    }
    int owner = -1;
    for (size_t r = 0; r < boundaries.size(); ++r) {
        const MdCoordBoundary &b = boundaries[r];
        if (!b.has_records) continue;
        if (owner < 0) owner = (int)r;
        if (MdCompareCoord(key.this_ref, key.this_coord,
                           b.first_ref, b.first_coord) < 0) {
            break;
        }
        owner = (int)r;
    }
    return owner >= 0 ? owner : fallback_rank;
}


static int MdExchangeCandidatesByCoordinate(
        std::vector<MpiMarkdupCandidateShared> *local_candidates_in,
        std::vector<unsigned char> *local_qnames_in,
        int range_has_records, int range_first_tid, int range_first_pos,
        int range_last_tid, int range_last_pos,
        int rank, int comm_size, size_t memory_limit,
        std::vector<MpiMarkdupCandidateShared> *owner_candidates,
        std::vector<unsigned char> *owner_qnames,
        MpiMarkdupStats *stats,
        const std::vector<MdCoordBoundary> *fixed_boundaries = nullptr,
        int owner_by_first = 0,
        MdRemoteCandidateObserver remote_observer = nullptr,
        void *remote_observer_opaque = nullptr) {
    double exchange_total_t0 = GetTime();
    std::vector<MpiMarkdupCandidateShared> &local_candidates =
        *local_candidates_in;
    std::vector<unsigned char> &local_qnames = *local_qnames_in;
    std::vector<MdCoordBoundary> dynamic_boundaries;
    if (!fixed_boundaries) {
        MdBuildCoordBoundaries(range_has_records, range_first_tid,
                               range_first_pos, range_last_tid,
                               range_last_pos, comm_size,
                               &dynamic_boundaries);
    }
    const std::vector<MdCoordBoundary> &boundaries = fixed_boundaries
        ? *fixed_boundaries : dynamic_boundaries;

    std::vector<unsigned long long> send_candidate_counts(
        (size_t)comm_size, 0);
    std::vector<unsigned long long> recv_candidate_counts(
        (size_t)comm_size, 0);
    std::vector<unsigned long long> send_qname_counts(
        (size_t)comm_size, 0);
    std::vector<unsigned long long> recv_qname_counts(
        (size_t)comm_size, 0);
    std::vector<size_t> send_candidate_offsets((size_t)comm_size + 1, 0);
    std::vector<size_t> send_qname_offsets((size_t)comm_size + 1, 0);
    std::unique_ptr<MpiMarkdupCandidateShared,
                    MdMarkdupCandidateBufferDeleter> send_candidates_flat;
    std::unique_ptr<unsigned char, MdByteBufferDeleter> send_qnames_flat;
    std::vector<MdRemoteCandidateRef> remote_candidates;
    size_t self_candidate_count = 0;
    size_t self_qname_count = 0;
    int local_ok = 1;

    double pack_t0 = GetTime();
    int next_owner_rank = -1;
    if (owner_by_first && rank >= 0 && rank < (int)boundaries.size()) {
        for (size_t r = (size_t)rank + 1; r < boundaries.size(); ++r) {
            if (boundaries[r].has_records) {
                next_owner_rank = (int)r;
                break;
            }
        }
    }
    for (size_t i = 0; i < local_candidates.size(); ++i) {
        const MpiMarkdupCandidateShared &candidate = local_candidates[i];
        const int owner = owner_by_first
            ? MdCoordinateOwnerByFirst(candidate.key, boundaries, rank, next_owner_rank)
            : MdCoordinateOwnerFast(candidate.key, boundaries, rank);
        if (owner < 0 || owner >= comm_size) {
            local_ok = 0;
            break;
        }
        if (owner == rank) {
            self_candidate_count++;
            if (!candidate.key.single) {
                if (candidate.qname_offset > local_qnames.size() ||
                    candidate.qname_len >
                        local_qnames.size() -
                            (size_t)candidate.qname_offset) {
                    fprintf(stderr,
                            "[rank %d] ERROR: invalid local markdup QNAME range.\n",
                            rank);
                    local_ok = 0;
                    break;
                }
                self_qname_count += (size_t)candidate.qname_len;
            }
            continue;
        }
        if (remote_observer &&
            remote_observer(remote_observer_opaque,
                            candidate, owner) != 0) {
            local_ok = 0;
            break;
        }
        MdRemoteCandidateRef ref;
        ref.index = i;
        ref.owner = owner;
        remote_candidates.push_back(ref);
        if (send_candidate_counts[(size_t)owner] == ULLONG_MAX) {
            local_ok = 0;
            break;
        }
        send_candidate_counts[(size_t)owner]++;
        if (!candidate.key.single) {
            if (candidate.qname_offset > local_qnames.size() ||
                candidate.qname_len >
                    local_qnames.size() -
                        (size_t)candidate.qname_offset) {
                fprintf(stderr,
                        "[rank %d] ERROR: invalid local markdup QNAME range.\n",
                        rank);
                local_ok = 0;
                break;
            }
            const unsigned long long qname_len =
                (unsigned long long)candidate.qname_len;
            if (send_qname_counts[(size_t)owner] >
                ULLONG_MAX - qname_len) {
                local_ok = 0;
                break;
            }
            send_qname_counts[(size_t)owner] += qname_len;
        }
    }
    if (!MdAllRanksOk(local_ok)) return -1;

    size_t total_send_candidates = 0;
    size_t total_send_qnames = 0;
    for (int i = 0; i < comm_size; ++i) {
        send_candidate_offsets[(size_t)i] = total_send_candidates;
        send_qname_offsets[(size_t)i] = total_send_qnames;
        if (send_candidate_counts[(size_t)i] >
                (unsigned long long)(SIZE_MAX - total_send_candidates) ||
            send_qname_counts[(size_t)i] >
                (unsigned long long)(SIZE_MAX - total_send_qnames)) {
            local_ok = 0;
            break;
        }
        total_send_candidates +=
            (size_t)send_candidate_counts[(size_t)i];
        total_send_qnames += (size_t)send_qname_counts[(size_t)i];
    }
    send_candidate_offsets[(size_t)comm_size] = total_send_candidates;
    send_qname_offsets[(size_t)comm_size] = total_send_qnames;
    if (!MdAllRanksOk(local_ok)) return -1;

    if (total_send_candidates > 0) {
        if (total_send_candidates >
            SIZE_MAX / sizeof(MpiMarkdupCandidateShared)) {
            local_ok = 0;
        } else {
            send_candidates_flat.reset(
                (MpiMarkdupCandidateShared *)aligned_alloc_custom(
                    64, total_send_candidates *
                            sizeof(MpiMarkdupCandidateShared)));
            if (!send_candidates_flat) local_ok = 0;
        }
    }
    if (total_send_qnames > 0) {
        send_qnames_flat.reset(
            aligned_alloc_custom(64, total_send_qnames));
        if (!send_qnames_flat) local_ok = 0;
    }
    if (!MdAllRanksOk(local_ok)) return -1;

    double setup_t0 = GetTime();
    MPI_Alltoall(send_candidate_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG,
                 recv_candidate_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
    MPI_Alltoall(send_qname_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG,
                 recv_qname_counts.data(), 1,
                 MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

    unsigned long long recv_candidate_total = 0;
    unsigned long long recv_qname_total = 0;
    for (int i = 0; i < comm_size; ++i) {
        if (recv_candidate_total >
                ULLONG_MAX - recv_candidate_counts[(size_t)i] ||
            recv_qname_total >
                ULLONG_MAX - recv_qname_counts[(size_t)i]) {
            local_ok = 0;
            break;
        }
        recv_candidate_total += recv_candidate_counts[(size_t)i];
        recv_qname_total += recv_qname_counts[(size_t)i];
    }
    if (!MdAllRanksOk(local_ok)) return -1;

    if (recv_candidate_total >
            (unsigned long long)(SIZE_MAX /
                sizeof(MpiMarkdupCandidateShared)) ||
        recv_qname_total > (unsigned long long)SIZE_MAX ||
        self_candidate_count >
            SIZE_MAX - (size_t)recv_candidate_total ||
        self_qname_count > SIZE_MAX - (size_t)recv_qname_total) {
        local_ok = 0;
    }

    size_t send_candidate_bytes_total = 0;
    size_t exchange_memory = 0;
    if (local_ok) {
        send_candidate_bytes_total =
            total_send_candidates * sizeof(MpiMarkdupCandidateShared);
        exchange_memory = send_candidate_bytes_total;
        if (total_send_qnames > SIZE_MAX - exchange_memory) {
            local_ok = 0;
        } else {
            exchange_memory += total_send_qnames;
        }
    }
    size_t local_memory =
        local_candidates.capacity() *
            sizeof(MpiMarkdupCandidateShared) +
        local_qnames.capacity();
    if (remote_candidates.capacity() >
        SIZE_MAX / sizeof(remote_candidates[0])) {
        local_ok = 0;
    } else if (local_ok) {
        const size_t remote_ref_memory =
            remote_candidates.capacity() *
            sizeof(remote_candidates[0]);
        if (local_memory > SIZE_MAX - remote_ref_memory) {
            local_ok = 0;
        } else {
            local_memory += remote_ref_memory;
        }
    }
    if (local_ok) {
        size_t recv_owner_bytes = 0;
        if ((size_t)recv_candidate_total >
            SIZE_MAX / sizeof(MpiMarkdupCandidateShared)) {
            local_ok = 0;
        } else {
            recv_owner_bytes =
                (size_t)recv_candidate_total *
                sizeof(MpiMarkdupCandidateShared);
        }
        if (local_ok &&
            (size_t)recv_qname_total > SIZE_MAX - recv_owner_bytes) {
            local_ok = 0;
        } else if (local_ok) {
            recv_owner_bytes += (size_t)recv_qname_total;
        }
        if (local_ok &&
            (local_memory > SIZE_MAX - exchange_memory ||
             local_memory + exchange_memory >
                 SIZE_MAX - recv_owner_bytes ||
             MdCheckMemory(local_memory + exchange_memory +
                               recv_owner_bytes,
                           memory_limit, stats, rank,
                           "coordinate candidate owner exchange") != 0)) {
            local_ok = 0;
        }
    }
    if (!MdAllRanksOk(local_ok)) return -1;
    stats->t_candidate_exchange_setup += GetTime() - setup_t0;

    std::vector<size_t> candidate_write = send_candidate_offsets;
    std::vector<size_t> qname_write = send_qname_offsets;
    for (size_t i = 0; i < remote_candidates.size(); ++i) {
        const MdRemoteCandidateRef &ref = remote_candidates[i];
        MpiMarkdupCandidateShared candidate =
            local_candidates[ref.index];
        const int owner = ref.owner;
        if (!candidate.key.single) {
            const size_t qname_src = (size_t)candidate.qname_offset;
            const size_t qname_len = (size_t)candidate.qname_len;
            const size_t qname_dst = qname_write[(size_t)owner];
            candidate.qname_offset =
                (uint64_t)(qname_dst -
                           send_qname_offsets[(size_t)owner]);
            if (qname_len > 0) {
                memcpy(send_qnames_flat.get() + qname_dst,
                       local_qnames.data() + qname_src, qname_len);
            }
            qname_write[(size_t)owner] += qname_len;
        }
        send_candidates_flat.get()[candidate_write[(size_t)owner]++] =
            candidate;
    }
    for (size_t i = remote_candidates.size(); i > 0; --i) {
        const size_t remove_index = remote_candidates[i - 1].index;
        const size_t last_index = local_candidates.size() - 1;
        if (remove_index != last_index) {
            local_candidates[remove_index] =
                local_candidates[last_index];
        }
        local_candidates.pop_back();
    }
    remote_candidates.clear();
    owner_candidates->swap(local_candidates);
    owner_qnames->swap(local_qnames);
    const size_t needed_candidates =
        owner_candidates->size() + (size_t)recv_candidate_total;
    const size_t needed_qnames =
        owner_qnames->size() + (size_t)recv_qname_total;
    if (owner_candidates->capacity() < needed_candidates) {
        owner_candidates->reserve(needed_candidates);
    }
    if (owner_qnames->capacity() < needed_qnames) {
        owner_qnames->reserve(needed_qnames);
    }
    stats->t_candidate_pack += GetTime() - pack_t0;
    if (!MdAllRanksOk(local_ok)) return -1;

    double mpi_t0 = GetTime();
    for (int step = 1; step < comm_size; ++step) {
        int dest = (rank + step) % comm_size;
        int src = (rank - step + comm_size) % comm_size;
        unsigned long long send_candidate_count =
            send_candidate_counts[(size_t)dest];
        unsigned long long send_qname_count =
            send_qname_counts[(size_t)dest];
        unsigned long long recv_candidate_count =
            recv_candidate_counts[(size_t)src];
        unsigned long long recv_qname_count =
            recv_qname_counts[(size_t)src];
        if (recv_candidate_count >
            (unsigned long long)(SIZE_MAX /
                sizeof(MpiMarkdupCandidateShared)) ||
            recv_qname_count > (unsigned long long)SIZE_MAX) {
            return -1;
        }
        size_t candidate_base = owner_candidates->size();
        size_t qname_base = owner_qnames->size();
        owner_candidates->resize(
            candidate_base + (size_t)recv_candidate_count);
        owner_qnames->resize(qname_base + (size_t)recv_qname_count);

        unsigned long long send_candidate_bytes =
            send_candidate_count * sizeof(MpiMarkdupCandidateShared);
        unsigned long long recv_candidate_bytes =
            recv_candidate_count * sizeof(MpiMarkdupCandidateShared);
        if (MdSendrecvBytes(
                send_candidate_count
                    ? send_candidates_flat.get() +
                          send_candidate_offsets[(size_t)dest]
                    : nullptr,
                send_candidate_bytes, dest,
                recv_candidate_count
                    ? owner_candidates->data() + candidate_base
                    : nullptr,
                recv_candidate_bytes, src, 4111) != 0 ||
            MdSendrecvBytes(
                send_qname_count
                    ? send_qnames_flat.get() +
                          send_qname_offsets[(size_t)dest]
                    : nullptr,
                send_qname_count, dest,
                recv_qname_count
                    ? owner_qnames->data() + qname_base
                    : nullptr,
                recv_qname_count, src, 4112) != 0) {
            return -1;
        }
        for (size_t i = candidate_base;
             i < owner_candidates->size(); ++i) {
            if (!(*owner_candidates)[i].key.single) {
                (*owner_candidates)[i].qname_offset +=
                    (uint64_t)qname_base;
            }
        }
        stats->mpi_candidate_bytes +=
            (long long)send_candidate_bytes +
            (long long)send_qname_count;
    }

    stats->t_candidate_exchange_mpi += GetTime() - mpi_t0;
    stats->t_candidate_exchange += GetTime() - exchange_total_t0;
    stats->owner_candidates = (long long)owner_candidates->size();
    stats->qname_bytes = (long long)owner_qnames->size();
    return 0;
}

static uint8_t MdBoundaryCandidateKind(
        const MpiMarkdupCandidateShared &candidate) {
    if (!candidate.key.single) return 0;
    return candidate.paired_marker ? 2 : 1;
}

static MdBoundaryCandidateId MdMakeBoundaryCandidateId(
        const MpiMarkdupCandidateShared &candidate) {
    MdBoundaryCandidateId id;
    id.global_order = candidate.global_order;
    id.kind = MdBoundaryCandidateKind(candidate);
    return id;
}

static bool MdBoundaryDecisionCandidate(
        const MpiMarkdupCandidateShared &candidate) {
    return !candidate.key.single || !candidate.paired_marker;
}

static int MdRouteBoundaryDuplicateOrders(
        std::vector<std::vector<uint64_t> > *duplicates_by_source,
        int rank, int comm_size, MdOnePassBoundaryPlan *plan,
        MpiMarkdupStats *stats) {
    if (!duplicates_by_source || !plan || !stats ||
        duplicates_by_source->size() != (size_t)comm_size) {
        return -1;
    }
    std::vector<unsigned long long> send_counts((size_t)comm_size, 0);
    std::vector<unsigned long long> recv_counts((size_t)comm_size, 0);
    for (int r = 0; r < comm_size; ++r) {
        std::vector<uint64_t> &values =
            (*duplicates_by_source)[(size_t)r];
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()),
                     values.end());
        send_counts[(size_t)r] =
            (unsigned long long)values.size();
    }
    MPI_Alltoall(send_counts.data(), 1, MPI_UNSIGNED_LONG_LONG,
                 recv_counts.data(), 1, MPI_UNSIGNED_LONG_LONG,
                 MPI_COMM_WORLD);

    int local_ok = 1;
    for (int step = 0; step < comm_size; ++step) {
        const int dest = (rank + step) % comm_size;
        const int src = (rank - step + comm_size) % comm_size;
        const unsigned long long send_count =
            send_counts[(size_t)dest];
        const unsigned long long recv_count =
            recv_counts[(size_t)src];
        if (recv_count >
            (unsigned long long)(SIZE_MAX / sizeof(uint64_t))) {
            local_ok = 0;
        }
        if (!MdAllRanksOk(local_ok)) return -1;

        std::vector<uint64_t> received;
        if (step == 0) {
            received = (*duplicates_by_source)[(size_t)rank];
        } else {
            try {
                received.resize((size_t)recv_count);
            } catch (...) {
                local_ok = 0;
            }
            if (!MdAllRanksOk(local_ok)) return -1;
            if (MdSendrecvBytes(
                    send_count
                        ? (*duplicates_by_source)[(size_t)dest].data()
                        : nullptr,
                    send_count * sizeof(uint64_t), dest,
                    recv_count ? received.data() : nullptr,
                    recv_count * sizeof(uint64_t), src, 4141) != 0) {
                return -1;
            }
            stats->mpi_result_bytes +=
                (long long)(send_count * sizeof(uint64_t));
        }
        for (size_t i = 0; i < received.size(); ++i) {
            std::unordered_map<uint64_t, uint8_t>::iterator found =
                plan->decisions.find(received[i]);
            if (found == plan->decisions.end()) {
                fprintf(stderr,
                        "[rank %d] ERROR: boundary duplicate order=%llu "
                        "was not selected by the source rank.\n",
                        rank, (unsigned long long)received[i]);
                local_ok = 0;
                continue;
            }
            found->second = (uint8_t)kMdOnePassDuplicate;
        }
    }
    return MdAllRanksOk(local_ok) ? 0 : -1;
}

static int MdAppendBoundaryCandidates(
        const std::vector<MpiMarkdupCandidateShared> &input_candidates,
        const std::vector<unsigned char> &input_qnames,
        const std::vector<int64_t> &positions,
        int max_read_length, MdOnePassBoundaryPlan *plan,
        std::vector<MpiMarkdupCandidateShared> *selected_candidates,
        std::vector<unsigned char> *selected_qnames) {
    if (positions.size() != input_candidates.size()) return -1;
    std::unordered_set<MdBoundaryCandidateId, MdBoundaryCandidateIdHash> seen;
    for (size_t i = 0; i < input_candidates.size(); ++i) {
        const MpiMarkdupCandidateShared &input = input_candidates[i];
        if (!swbam_markdup::HaloContainsPosition(
                plan->halos, input.key.this_ref, positions[i])) {
            continue;
        }
        const MdBoundaryCandidateId id =
            MdMakeBoundaryCandidateId(input);
        if (!seen.insert(id).second) continue;

        MpiMarkdupCandidateShared candidate = input;
        if (!candidate.key.single) {
            const size_t qname_offset =
                (size_t)candidate.qname_offset;
            const size_t qname_len = (size_t)candidate.qname_len;
            if (qname_offset > input_qnames.size() ||
                qname_len > input_qnames.size() - qname_offset) {
                return -1;
            }
            candidate.qname_offset =
                (uint64_t)selected_qnames->size();
            selected_qnames->insert(
                selected_qnames->end(),
                input_qnames.data() + qname_offset,
                input_qnames.data() + qname_offset + qname_len);
        }
        // Boundary duplicate routing uses a globally unique record id.
        candidate.ordinal = candidate.global_order;
        selected_candidates->push_back(candidate);
        if (MdBoundaryDecisionCandidate(candidate)) {
            plan->decisions[candidate.global_order] =
                !plan->partial && swbam_markdup::HaloCoversKey(
                    plan->halos, candidate.key, max_read_length)
                    ? (uint8_t)kMdOnePassKeep : (uint8_t)kMdOnePassUnknown;
        }
    }
    return 0;
}

struct MdBoundaryProbeContext {
    size_t first_end;
    size_t last_begin;
    size_t processed_blocks;
    long long global_block_begin;
    size_t memory_limit;
    bool overflow;
    int first_last_tid;
    int first_last_pos;
    int last_first_tid;
    int last_first_pos;
    bool captured;
    std::vector<MpiMarkdupCandidateShared> candidates;
    std::vector<unsigned char> qnames;
    std::vector<int64_t> positions;

    MdBoundaryProbeContext()
        : first_end(0), last_begin(0), processed_blocks(0),
          global_block_begin(0), memory_limit(0), overflow(false), first_last_tid(-1),
          first_last_pos(-1), last_first_tid(-1),
          last_first_pos(-1), captured(false) {}
};

static int MdCaptureBoundaryProbe(
        void *opaque,
        std::vector<MpiMarkdupCandidateShared> *candidates,
        std::vector<unsigned char> *qnames, int, int, uint64_t,
        MdExtractBatchView *batch) {
    MdBoundaryProbeContext *probe =
        (MdBoundaryProbeContext *)opaque;
    if (!probe || !candidates || !qnames ||
        !batch || !batch->decomp || batch->n_blocks <= 0) {
        return -1;
    }
    for (size_t i = 0; i < (size_t)batch->n_blocks; ++i) {
        Bam2BamPara &para = batch->decomp[i];
        if (para.n_total_records <= 0 || !para.output_records) continue;
        const size_t index = probe->processed_blocks + i;
        if (index <= probe->first_end) {
            bam1_t *record = para.output_records[para.n_total_records - 1];
            if (!record) return -1;
            probe->first_last_tid = record->core.tid;
            probe->first_last_pos = record->core.pos;
        }
        if (index >= probe->last_begin && probe->last_first_tid == -1) {
            bam1_t *record = para.output_records[0];
            if (!record) return -1;
            probe->last_first_tid = record->core.tid;
            probe->last_first_pos = record->core.pos;
        }
    }
    // Retain metadata only. Decoded records remain owned by the extractor.
    const size_t retained = (probe->candidates.size() + candidates->size()) *
        (sizeof(MpiMarkdupCandidateShared) + sizeof(int64_t)) +
        probe->qnames.size() + qnames->size();
    if (probe->memory_limit &&
        (batch->fixed_workspace_bytes >= probe->memory_limit ||
         retained > (probe->memory_limit - batch->fixed_workspace_bytes) / 3)) {
        probe->overflow = true;
        probe->candidates.clear();
        probe->qnames.clear();
        probe->positions.clear();
    }
    if (!probe->overflow) {
        const size_t qname_base = probe->qnames.size();
        probe->qnames.insert(probe->qnames.end(), qnames->begin(), qnames->end());
        for (size_t i = 0; i < candidates->size(); ++i) {
            MpiMarkdupCandidateShared candidate = (*candidates)[i];
            const uint64_t index = (candidate.global_order >> 32) -
                (uint64_t)probe->global_block_begin - probe->processed_blocks;
            const uint64_t record_index = candidate.global_order & 0xffffffffull;
            if (index >= (uint64_t)batch->n_blocks ||
                record_index >= (uint64_t)batch->decomp[index].n_total_records)
                return -1;
            bam1_t *record = batch->decomp[index].output_records[record_index];
            if (!record) return -1;
            if (!candidate.key.single) candidate.qname_offset += qname_base;
            probe->candidates.push_back(candidate);
            probe->positions.push_back((int64_t)record->core.pos + 1);
        }
    }
    probe->processed_blocks += batch->n_blocks;
    probe->captured = true;
    return 0;
}

static int MdBuildOnePassBoundaryPlan(
        const swbam::BamInputBackend &backend,
        const swbam::BgzfBlockSpan *spans, size_t span_count,
        long long global_block_begin, int rank, int comm_size,
        int include_fails, int max_read_length, size_t memory_limit,
        MdOnePassBoundaryPlan *plan, MpiMarkdupStats *stats) {
    if (!plan || !stats || max_read_length <= 0 ||
        (span_count > 0 && !spans)) {
        return -1;
    }
    const double boundary_t0 = GetTime();
    plan->enabled = 0;
    plan->partial = false;
    plan->conservative = false;
    plan->boundaries.clear();
    plan->halos.clear();
    plan->closed_owner_keys.clear();
    plan->closed_owner_ranges.clear();
    plan->candidates.clear();
    plan->decisions.clear();

    // Keep the existing conservative window when rank partitions are empty.
    const int local_nonempty = span_count > 0 ? 1 : 0;
    int all_nonempty = 0;
    MPI_Allreduce(&local_nonempty, &all_nonempty, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (!all_nonempty) {
        stats->t_boundary_check += GetTime() - boundary_t0;
        return 0;
    }

    const size_t side_blocks = std::min(span_count, (size_t)kMarkdupNB / 2);
    int first_has = 0, first_tid = -1, first_pos = -1;
    int last_tid = -1, last_pos = -1;
    std::vector<swbam::BgzfBlockSpan> probe_spans;
    std::vector<size_t> probe_indices;
    MdBoundaryProbeContext probe;
    int local_ok = 1;
    {
        const size_t last_actual_begin = span_count - side_blocks;
        plan->first_probe_end_block = (uint64_t)global_block_begin + side_blocks;
        plan->last_probe_begin_block = (uint64_t)global_block_begin + last_actual_begin;
        if (last_actual_begin <= side_blocks) {
            probe_spans.assign(spans, spans + span_count);
            probe_indices.resize(span_count);
            for (size_t i = 0; i < span_count; ++i) probe_indices[i] = i;
        } else {
            probe_spans.reserve(side_blocks * 2);
            probe_indices.reserve(side_blocks * 2);
            for (size_t i = 0; i < side_blocks; ++i) {
                probe_spans.push_back(spans[i]);
                probe_indices.push_back(i);
            }
            for (size_t i = last_actual_begin; i < span_count; ++i) {
                probe_spans.push_back(spans[i]);
                probe_indices.push_back(i);
            }
        }
        size_t last_probe_begin = 0;
        while (last_probe_begin < probe_indices.size() &&
               probe_indices[last_probe_begin] < last_actual_begin) {
            ++last_probe_begin;
        }
        if (last_probe_begin >= probe_indices.size()) return -1;

        std::vector<MpiMarkdupCandidateShared> ignored_candidates;
        std::vector<unsigned char> ignored_qnames;
        uint64_t ignored_records = 0;
        MpiMarkdupStats probe_stats = {};
        const bool whole_rank = probe_indices.size() == span_count;
        probe.first_end = whole_rank ? span_count - 1 : side_blocks - 1;
        probe.last_begin = last_probe_begin;
        probe.global_block_begin = global_block_begin;
        probe.memory_limit = memory_limit;

        MdInputPass probe_source(
            &backend, probe_spans.data(), probe_spans.size());
        local_ok = MdExtractCandidatesPass(
            probe_source, global_block_begin, rank, include_fails,
            memory_limit, &ignored_candidates, &ignored_qnames,
            &ignored_records, &first_has, &first_tid, &first_pos,
            &last_tid, &last_pos, &probe_stats,
            MdCaptureBoundaryProbe, &probe) == 0;
        if (!MdAllRanksOk(local_ok)) return -1;

        local_ok = probe.captured && (probe.overflow ||
            probe.candidates.size() == (size_t)(
                probe_stats.pair_candidates + probe_stats.single_candidates));
        if (!local_ok) {
            fprintf(stderr,
                    "[rank %d] ERROR: incomplete markdup boundary probe "
                    "candidates: retained=%zu expected=%lld.\n",
                    rank, probe.candidates.size(),
                    probe_stats.pair_candidates + probe_stats.single_candidates);
        }
        if (!MdAllRanksOk(local_ok)) return -1;

        for (size_t i = 0; i < probe.candidates.size(); ++i) {
            MpiMarkdupCandidateShared &candidate = probe.candidates[i];
            const uint64_t fake_block = candidate.global_order >> 32;
            if (fake_block < (uint64_t)global_block_begin) {
                local_ok = 0;
                break;
            }
            const uint64_t probe_index64 =
                fake_block - (uint64_t)global_block_begin;
            if (probe_index64 >= probe_indices.size()) {
                local_ok = 0;
                break;
            }
            const uint64_t record_index = candidate.global_order &
                0xffffffffull;
            candidate.global_order =
                ((uint64_t)global_block_begin +
                 (uint64_t)probe_indices[(size_t)probe_index64]) << 32 |
                record_index;
        }
        if (!MdAllRanksOk(local_ok)) return -1;
        if (MdBuildFirstCoordBoundaries(
                first_has, first_tid, first_pos, comm_size,
                &plan->boundaries) != 0) {
            return -1;
        }

        int previous_rank = -1;
        int next_rank = -1;
        for (int r = rank - 1; r >= 0; --r) {
            if (plan->boundaries[(size_t)r].has_records) {
                previous_rank = r;
                break;
            }
        }
        for (int r = rank + 1; r < comm_size; ++r) {
            if (plan->boundaries[(size_t)r].has_records) {
                next_rank = r;
                break;
            }
        }

        // Invalid/unretained metadata cannot support even a partial prepass.
        // Incomplete coordinate coverage, in contrast, still allows DUP-only
        // decisions among the candidates that were actually observed.
        if (!MdAllRanksOk(first_has && first_tid >= 0 && !probe.overflow)) {
            plan->boundaries.clear();
            stats->t_boundary_check += GetTime() - boundary_t0;
            if (rank == 0)
                printf("MPI BAM markdup boundary halo=fallback pos_halo=1 blocks=%zu "
                       "reason=unavailable-probe\n", probe_spans.size());
            return 0;
        }
        if (local_ok && !whole_rank && previous_rank >= 0) {
            const MdCoordBoundary &cut =
                plan->boundaries[(size_t)rank];
            const int32_t halo_ref = probe.first_last_tid + 1;
            const int64_t halo_coord = (int64_t)probe.first_last_pos + 1;
            const int64_t needed = cut.first_coord >
                    INT64_MAX - max_read_length
                ? INT64_MAX : cut.first_coord + max_read_length;
            local_ok = probe.first_last_tid >= 0 &&
                MdCompareCoord(halo_ref, halo_coord,
                               cut.first_ref, needed) > 0;
        }
        if (local_ok && !whole_rank && next_rank >= 0) {
            const MdCoordBoundary &cut =
                plan->boundaries[(size_t)next_rank];
            const int32_t halo_ref = probe.last_first_tid + 1;
            const int64_t halo_coord = (int64_t)probe.last_first_pos + 1;
            const int64_t needed = cut.first_coord <
                    INT64_MIN + max_read_length
                ? INT64_MIN : cut.first_coord - max_read_length;
            // An unmapped file tail needs no right halo on the final rank.
            local_ok = probe.last_first_tid >= 0 &&
                MdCompareCoord(halo_ref, halo_coord,
                               cut.first_ref, needed) < 0;
        }
        int coverage_ok = 0;
        MPI_Allreduce(&local_ok, &coverage_ok, 1, MPI_INT, MPI_MIN,
                      MPI_COMM_WORLD);
        plan->partial = coverage_ok == 0;
    }

    bool seen_nonempty = false;
    for (size_t r = 0; r < plan->boundaries.size(); ++r) {
        const MdCoordBoundary &cut = plan->boundaries[r];
        if (!cut.has_records || cut.first_ref <= 0) continue;
        if (!seen_nonempty) { seen_nonempty = true; continue; }
        swbam_markdup::HaloInterval halo = {cut.first_ref,
            swbam_markdup::HaloOffset(cut.first_coord, -max_read_length),
            swbam_markdup::HaloOffset(cut.first_coord, max_read_length)};
        plan->halos.push_back(halo);
    }
    swbam_markdup::MergeHalos(&plan->halos);

    std::vector<MpiMarkdupCandidateShared> selected_candidates;
    std::vector<unsigned char> selected_qnames;
    local_ok = MdAppendBoundaryCandidates(
        probe.candidates, probe.qnames, probe.positions,
        max_read_length, plan, &selected_candidates,
        &selected_qnames) == 0;
    if (!MdAllRanksOk(local_ok)) return -1;

    std::vector<MpiMarkdupCandidateShared> owner_candidates;
    std::vector<unsigned char> owner_qnames;
    MpiMarkdupStats boundary_stats = {};
    if (MdExchangeCandidatesByCoordinate(
            &selected_candidates, &selected_qnames,
            first_has, first_tid, first_pos, last_tid, last_pos,
            rank, comm_size, memory_limit,
            &owner_candidates, &owner_qnames, &boundary_stats,
            &plan->boundaries, 1) != 0) {
        return -1;
    }
    int conservative = plan->partial ? 1 : 0;
    for (size_t i = 0; i < owner_candidates.size(); ++i) {
        const MpiMarkdupCandidateShared &candidate = owner_candidates[i];
        // Partial coverage never certifies a final KEEP or a closed marker.
        if (plan->partial) continue;
        if (swbam_markdup::HaloCoversKey(plan->halos, candidate.key, max_read_length)) {
            plan->closed_owner_keys.push_back(candidate.key);
            continue;
        }
        // A partial key can use the ordinary window only when its entire
        // position support belongs strictly to this rank's coordinate range.
        int64_t low = 0, high = 0;
        bool local_support = swbam_markdup::KeyPositionSupport(
            candidate.key, max_read_length, &low, &high);
        const MdCoordBoundary &first = plan->boundaries[(size_t)rank];
        if (rank > 0 && MdCompareCoord(candidate.key.this_ref, low,
                                      first.first_ref, first.first_coord) <= 0)
            local_support = false;
        for (int r = rank + 1; r < comm_size; ++r) {
            const MdCoordBoundary &next = plan->boundaries[(size_t)r];
            if (!next.has_records) continue;
            if (MdCompareCoord(candidate.key.this_ref, high,
                               next.first_ref, next.first_coord) >= 0)
                local_support = false;
            break;
        }
        if (!local_support || candidate.source_rank != rank) conservative = 1;
    }
    std::sort(plan->closed_owner_keys.begin(), plan->closed_owner_keys.end(), MdKeyLess);
    plan->closed_owner_keys.erase(std::unique(plan->closed_owner_keys.begin(),
        plan->closed_owner_keys.end(), MdKeyEqual), plan->closed_owner_keys.end());
    std::map<int32_t, swbam_markdup::HaloInterval> key_ranges;
    for (const auto &key : plan->closed_owner_keys) {
        auto found = key_ranges.find(key.this_ref);
        if (found == key_ranges.end()) {
            key_ranges.emplace(key.this_ref, swbam_markdup::HaloInterval{
                key.this_ref, key.this_coord, key.this_coord});
        } else {
            found->second.low = std::min(found->second.low, key.this_coord);
            found->second.high = std::max(found->second.high, key.this_coord);
        }
    }
    for (const auto &entry : key_ranges) plan->closed_owner_ranges.push_back(entry.second);
    int any_conservative = 0;
    MPI_Allreduce(&conservative, &any_conservative, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    plan->conservative = any_conservative != 0;
    std::vector<std::vector<uint64_t> > duplicates_by_source;
    if (MpiMarkdupFindDuplicatesStreamingHash(
            &owner_candidates, owner_qnames, comm_size,
            &duplicates_by_source, &boundary_stats) != 0 ||
        MdRouteBoundaryDuplicateOrders(
            &duplicates_by_source, rank, comm_size, plan,
            &boundary_stats) != 0) {
        return -1;
    }

    size_t final_keep = 0, final_dup = 0, pending = 0, marker_pending = 0;
    for (size_t i = 0; i < probe.candidates.size(); ++i) {
        const MpiMarkdupCandidateShared &candidate = probe.candidates[i];
        if (!swbam_markdup::HaloContainsPosition(
                plan->halos, candidate.key.this_ref, probe.positions[i])) continue;
        if (!MdBoundaryDecisionCandidate(candidate)) {
            if (!plan->partial && swbam_markdup::HaloCoversKey(
                    plan->halos, candidate.key, max_read_length))
                plan->candidates.insert(MdMakeBoundaryCandidateId(candidate));
            else ++marker_pending;
            continue;
        }
        const auto found = plan->decisions.find(candidate.global_order);
        if (found == plan->decisions.end()) return -1;
        if (found->second == kMdOnePassUnknown) {
            ++pending;
            continue;
        }
        plan->candidates.insert(MdMakeBoundaryCandidateId(candidate));
        if (found->second == kMdOnePassDuplicate) ++final_dup;
        else ++final_keep;
    }
    for (auto it = plan->decisions.begin(); it != plan->decisions.end();) {
        if (it->second == kMdOnePassUnknown) it = plan->decisions.erase(it);
        else ++it;
    }

    stats->pair_duplicates += boundary_stats.pair_duplicates;
    stats->single_duplicates += boundary_stats.single_duplicates;
    stats->owner_candidates += boundary_stats.owner_candidates;
    stats->qname_bytes += boundary_stats.qname_bytes;
    stats->mpi_candidate_bytes += boundary_stats.mpi_candidate_bytes;
    stats->mpi_result_bytes += boundary_stats.mpi_result_bytes;
    plan->enabled = 1;
    stats->t_boundary_check += GetTime() - boundary_t0;
    if (rank == 0) {
        printf("MPI BAM markdup boundary halo=%s pos_halo=1 blocks=%zu\n",
               plan->partial ? "partial" : "enabled", probe_spans.size());
    }
    fprintf(stderr, "[md-halo rank=%d] enabled=1 partial=%d conservative=%d probe_blocks=%zu "
            "final_keep=%zu final_dup=%zu pending=%zu marker_pending=%zu "
            "candidate_bytes=%lld result_bytes=%lld\n", rank, (int)plan->partial,
            (int)plan->conservative,
            probe_spans.size(), final_keep, final_dup, pending, marker_pending,
            boundary_stats.mpi_candidate_bytes, boundary_stats.mpi_result_bytes);
    return 0;
}





static void MdInitEmptyRawCompress(MpiSortRawCompressPara *para,
                                   int block_id) {
    memset(para, 0, sizeof(*para));
    para->block_id = block_id;
    para->status = -1;
    para->compress_level = 1;
}

static void MdSetupRawCompress(MpiSortRawCompressPara *para,
                               int block_id, bam_block *payload,
                               bam_block *output,
                               int compress_level) {
    memset(para, 0, sizeof(*para));
    para->block_id = block_id;
    para->un_comp_block = payload;
    para->un_comp_size = payload ? payload->pos : 0;
    para->output_block = output;
    para->output_size = 0;
    para->status = para->un_comp_size > 0 ? 0 : -1;
    para->compress_level = compress_level;
}

static void MdResetPayloadBlock(bam_block *block, int block_id) {
    block->pos = 0;
    block->length = 0;
    block->errcode = 0;
    block->block_id = block_id;
    block->block_address = 0;
}

static uint32_t MdLoadLe32(const unsigned char *p) {
    uint32_t v = 0;
    memcpy(&v, p, sizeof(v));
    return v;
}

static int MdAppendRawPayload(bam_block *payload,
                              const unsigned char *record,
                              uint32_t raw_len) {
    if (!payload || !record || raw_len == 0 ||
        raw_len > BGZF_BLOCK_SIZE) {
        return -1;
    }
    if (payload->pos > 0 &&
        payload->pos + (int)raw_len > BGZF_BLOCK_SIZE) {
        return 1;
    }
    if (payload->pos + (int)raw_len > BGZF_BLOCK_SIZE) {
        return -1;
    }
    memcpy(payload->data + payload->pos, record, raw_len);
    payload->pos += (int)raw_len;
    payload->length = payload->pos;
    return 0;
}


static int MdValidateRankBoundaries(
        int local_has_records, int first_tid, int first_pos,
        int last_tid, int last_pos, int rank, int comm_size) {
    long long local[5] = {
        local_has_records, first_tid, first_pos, last_tid, last_pos
    };
    std::vector<long long> all((size_t)comm_size * 5);
    MPI_Allgather(local, 5, MPI_LONG_LONG,
                  all.data(), 5, MPI_LONG_LONG,
                  MPI_COMM_WORLD);
    int previous_tid = -1;
    int previous_pos = -1;
    int have_previous = 0;
    for (int r = 0; r < comm_size; ++r) {
        const long long *entry = all.data() + (size_t)r * 5;
        if (!entry[0]) continue;
        int current_first_tid = (int)entry[1];
        int current_first_pos = (int)entry[2];
        if (current_first_tid >= 0 && have_previous &&
            (current_first_tid < previous_tid ||
             (current_first_tid == previous_tid &&
              current_first_pos < previous_pos))) {
            if (rank == 0) {
                fprintf(stderr,
                        "ERROR: markdup input is not coordinate sorted "
                        "across rank boundary before rank %d.\n", r);
            }
            return -1;
        }
        previous_tid = (int)entry[3];
        previous_pos = (int)entry[4];
        have_previous = 1;
    }
    return 0;
}


struct MdOnePassBlockInfo {
    uint64_t ordinal_base;
    uint64_t spool_offset;
    uint64_t block_order;
    int n_records;
    int unresolved;
    enum Phase { Waiting, Ready, Compressing, Stored } phase;
    bool dirty;
};

struct MdOnePassPendingBatch {
    MdBlockSet blocks;
    std::vector<MdOnePassBlockInfo> infos;
    std::vector<uint8_t> states;
    std::shared_ptr<swbam_usage::CachedTemporary> spool;
    uint64_t ordinal_base;
    uint64_t ordinal_end;
    int n_blocks;
    int remaining = 0;
    size_t diagnostic_spills = 0;

    MdOnePassPendingBatch()
        : ordinal_base(0), ordinal_end(0), n_blocks(0) {
        memset(&blocks, 0, sizeof(blocks));
    }

    ~MdOnePassPendingBatch() {
        FreeBlocks();
    }

    void FreeBlocks() {
        if (!blocks.data && blocks.blocks)
            for (int b = 0; b < blocks.n; ++b)
                if (blocks.blocks[b].data) aligned_free_custom(blocks.blocks[b].data);
        MdFreeBlockSet(&blocks);
    }

    MdOnePassPendingBatch(MdOnePassPendingBatch &&other)
        : blocks(other.blocks), infos(std::move(other.infos)),
          states(std::move(other.states)),
          spool(std::move(other.spool)),
          ordinal_base(other.ordinal_base), ordinal_end(other.ordinal_end),
          n_blocks(other.n_blocks), remaining(other.remaining),
          diagnostic_spills(other.diagnostic_spills) {
        memset(&other.blocks, 0, sizeof(other.blocks));
        other.n_blocks = 0;
    }

    MdOnePassPendingBatch &operator=(MdOnePassPendingBatch &&other) {
        if (this != &other) {
            FreeBlocks();
            blocks = other.blocks;
            infos = std::move(other.infos);
            states = std::move(other.states);
            spool = std::move(other.spool);
            ordinal_base = other.ordinal_base;
            ordinal_end = other.ordinal_end;
            n_blocks = other.n_blocks;
            remaining = other.remaining;
            diagnostic_spills = other.diagnostic_spills;
            memset(&other.blocks, 0, sizeof(other.blocks));
            other.n_blocks = 0;
        }
        return *this;
    }

private:
    MdOnePassPendingBatch(const MdOnePassPendingBatch &);
    MdOnePassPendingBatch &operator=(const MdOnePassPendingBatch &);
};

struct MdOnePassReadyBlock {
    MdOnePassPendingBatch *batch;
    int block_index;
};

class MdOnePassCompressOperator : public swbam::cpe::CpeWriteKernelOperator {
public:
    explicit MdOnePassCompressOperator(int level)
        : CpeWriteKernelOperator({"markdup-window-compress",
              (void *)slave_mpi_sort_compress_payload, kMarkdupNB, false}),
          level_(level) {}
    int Initialize() { return 0; }
    void Shutdown() {}
    int Prepare(const swbam::cpe::CpeWriteBatchInput &input,
                swbam::BgzfBlockBatch *, swbam::BgzfBlockBatch *output) {
        if (input.kind != swbam::cpe::CpeWriteBatchInput::kUncompressedBgzf)
            return -1;
        for (int i = 0; i < kMarkdupNB; ++i) {
            MdInitEmptyRawCompress(paras_ + i, i);
            if ((size_t)i < input.count) {
                bam_block *block = const_cast<bam_block *>(
                    input.uncompressed->blocks() + i);
                if (block->length == 0 || block->length > BGZF_BLOCK_SIZE)
                    return -1;
                MdSetupRawCompress(paras_ + i, i, block,
                                   output->blocks() + i, level_);
            }
        }
        return 0;
    }
    void *kernel_arguments() { return paras_; }
    void ObserveKernel(double, size_t) {}
    int Validate(size_t count) const {
        for (size_t i = 0; i < count; ++i) {
            if (paras_[i].status != 0 || !paras_[i].output_block ||
                paras_[i].output_block->length == 0 ||
                paras_[i].output_block->length > BGZF_MAX_BLOCK_SIZE)
                return -1;
        }
        return 0;
    }
    long long OutputBytes(size_t count) const {
        long long bytes = 0;
        for (size_t i = 0; i < count; ++i)
            bytes += paras_[i].output_block->length;
        return bytes;
    }
    int Finish() { return 0; }
private:
    int level_;
    MpiSortRawCompressPara paras_[kMarkdupNB];
};

class MdOnePassOutput : public swbam::cpe::CompressedBgzfBatchPostProcessor {
public:
    MdOnePassOutput(swbam_usage::OrderedBody *sink, MpiMarkdupStats *stats,
                    long long *written, const uint64_t *keys)
        : sink_(sink), stats_(stats), written_(written), keys_(keys) {}
    int PostProcessCompressedBatch(const bam_block *blocks, size_t count) {
        const double t0 = GetTime();
        int ret = 0;
        for (size_t i = 0; i < count; ++i) {
            if (sink_->AppendOrdered(keys_[i], blocks[i].data, blocks[i].length) != 0) {
                ret = -1;
                break;
            }
            ++stats_->bgzf_blocks;
            ++*written_;
        }
        stats_->t_write += GetTime() - t0;
        return ret;
    }
private:
    swbam_usage::OrderedBody *sink_;
    MpiMarkdupStats *stats_;
    long long *written_;
    const uint64_t *keys_;
};

struct MdOnePassContext {
    int rank;
    int comm_size;
    int max_read_length;
    size_t memory_limit;
    size_t extract_workspace_bytes;
    int *range_has_records;
    int *range_first_tid;
    int *range_first_pos;
    swbam::RankBodySink *sink;
    MpiMarkdupStats *stats;
    const MdOnePassBoundaryPlan *boundary_plan;
    MpiMarkdupStreamingWindow *window;
    std::vector<MdCoordBoundary> boundaries;
    std::vector<MpiMarkdupCandidateShared> owner_candidates;
    std::vector<unsigned char> owner_qnames;
    std::vector<std::vector<MpiMarkdupDecisionShared> > decisions_by_source;
    std::list<MdOnePassPendingBatch> pending;
    bool use_pending_index = false;
    size_t pending_index_peak = 0;
    swbam_markdup::PendingIndex<MdOnePassPendingBatch> pending_index;
    std::vector<long long> boundary_progress;
    bool release_first_boundary = false;
    bool release_last_boundary = false;
    std::map<uint64_t, MdOnePassReadyBlock> ready;
    std::vector<MdOnePassReadyBlock> dirty_blocks;
    swbam_usage::OrderedBody *ordered_body = nullptr;
    uint64_t ordered_next = 0;
    bool reordered = false;
    uint64_t prepared_keys[kMarkdupNB];
    size_t catchup_batches = 0;
    std::vector<MdBlockSet> free_block_sets;
    MdBlockSet fallback_payload;
    bam_block prepared_payload[kMarkdupNB];
    swbam::BgzfBlockBatch prepared_view;
    std::unique_ptr<MdOnePassCompressOperator> compress_op;
    std::unique_ptr<MdOnePassOutput> output_post;
    std::unique_ptr<swbam::cpe::CpeWritePipelineSession> write_session;
    std::vector<MdOnePassReadyBlock> prepared_inputs;
    int prepared_outputs;
    int compression_running;
    int compress_level;
    long long captured_blocks;
    long long retired_blocks;
    uint64_t captured_records;
    uint64_t retired_records;
    long long launched_output_blocks;
    long long written_output_blocks;
    long long processed_batches;
    long long initial_candidate_bytes = 0;
    long long initial_result_bytes = 0;
    unsigned long long remote_candidates = 0;
    unsigned long long remote_markers = 0;
    size_t resident_pending = 0;
    size_t spilled_batches = 0;
    std::weak_ptr<swbam_usage::CachedTemporary> spill_file;
    bool pending_diagnostics = MdEnvFlagEnabled("SWBAM_MARKDUP_PENDING_DIAG");
    size_t diagnostic_unique_spills = 0;
    size_t diagnostic_reloads = 0;
    size_t diagnostic_ready_spilled = 0;
    size_t diagnostic_waiting_spilled = 0;
    uint64_t diagnostic_spill_bytes = 0;
    uint64_t diagnostic_reload_bytes = 0;
    double diagnostic_spill_seconds = 0.0;
    double diagnostic_reload_seconds = 0.0;
    size_t diagnostic_last_reported_spills = 0;

    MdOnePassContext()
        : rank(0), comm_size(1), max_read_length(300), memory_limit(0),
          extract_workspace_bytes(0),
          range_has_records(nullptr), range_first_tid(nullptr),
          range_first_pos(nullptr), sink(nullptr), stats(nullptr),
          boundary_plan(nullptr), window(nullptr), prepared_outputs(0),
          compression_running(0), captured_blocks(0), retired_blocks(0),
          captured_records(0), retired_records(0),
          launched_output_blocks(0), written_output_blocks(0),
          processed_batches(0) {
        memset(&fallback_payload, 0, sizeof(fallback_payload));
        memset(prepared_payload, 0, sizeof(prepared_payload));
        compress_level = 1;
    }
};

static void MdOnePassCleanup(MdOnePassContext *context) {
    if (!context) return;
    if (context->stats) {
        fprintf(stderr, "[md-comm rank=%d] halo_enabled=%d conservative=%d "
                "main_candidates=%llu main_markers=%llu main_candidate_bytes=%lld "
                "main_result_bytes=%lld\n", context->rank,
                context->boundary_plan ? context->boundary_plan->enabled : 0,
                !context->boundary_plan || !context->boundary_plan->enabled ||
                    context->boundary_plan->conservative ? 1 : 0,
                context->remote_candidates, context->remote_markers,
                context->stats->mpi_candidate_bytes - context->initial_candidate_bytes,
                context->stats->mpi_result_bytes - context->initial_result_bytes);
    }
    if (context->spilled_batches)
        fprintf(stderr, "[rank %d] markdup deferred payload spill batches=%zu\n",
                context->rank, context->spilled_batches);
    if (context->pending_diagnostics) {
        fprintf(stderr, "[md-pending-summary rank=%d] batches=%lld spills=%zu "
                "unique_spills=%zu respills=%zu reload_blocks=%zu ready_spilled=%zu "
                "waiting_spilled=%zu spill_bytes=%llu reload_bytes=%llu "
                "spill_s=%.6f reload_s=%.6f\n", context->rank,
                context->processed_batches, context->spilled_batches,
                context->diagnostic_unique_spills,
                context->spilled_batches - context->diagnostic_unique_spills,
                context->diagnostic_reloads, context->diagnostic_ready_spilled,
                context->diagnostic_waiting_spilled,
                (unsigned long long)context->diagnostic_spill_bytes,
                (unsigned long long)context->diagnostic_reload_bytes,
                context->diagnostic_spill_seconds, context->diagnostic_reload_seconds);
        fflush(stderr);
    }
    fprintf(stderr, "[md-scheduler rank=%d] catchup_batches=%zu reordered=%d raw_spills=%zu "
            "raw_write_bytes=%llu raw_read_bytes=%llu lookup=%s index_peak=%zu\n",
            context->rank, context->catchup_batches, (int)context->reordered,
            context->spilled_batches,
            (unsigned long long)context->diagnostic_spill_bytes,
            (unsigned long long)context->diagnostic_reload_bytes,
            context->use_pending_index ? "indexed" : "linear", context->pending_index_peak);
    // Session destruction joins an in-flight CPE before any borrowed arena dies.
    context->write_session.reset();
    context->compression_running = 0;
    if (context->window) {
        MpiMarkdupStreamingWindowDestroy(context->window);
        context->window = nullptr;
    }
    context->pending_index.Clear();
    context->pending.clear();
    for (size_t i = 0; i < context->free_block_sets.size(); ++i) {
        MdFreeBlockSet(&context->free_block_sets[i]);
    }
    context->free_block_sets.clear();
    MdFreeBlockSet(&context->fallback_payload);
}

static void MdOnePassRecycleBlockSet(MdOnePassContext *context,
                                     MdBlockSet *blocks) {
    if (!context || !blocks || !blocks->blocks) return;
    if (context->free_block_sets.size() < 2) {
        try {
            context->free_block_sets.push_back(*blocks);
            memset(blocks, 0, sizeof(*blocks));
            return;
        } catch (...) {
        }
    }
    MdFreeBlockSet(blocks);
}

static int MdOnePassAcquireBlockSet(MdOnePassContext *context,
                                    MdBlockSet *blocks) {
    memset(blocks, 0, sizeof(*blocks));
    if (!context->free_block_sets.empty()) {
        *blocks = context->free_block_sets.back();
        context->free_block_sets.pop_back();
        return 0;
    }
    return MdAllocateBlockSet(blocks, kMarkdupNB);
}

static int MdOnePassSpillPending(MdOnePassContext *context) {
    const size_t limit = swbam_usage::CacheMemoryLimit();
    if (context->resident_pending <= limit) return 0;
    for (auto &batch : context->pending) {
        if (context->resident_pending <= limit / 2) break;
        if (!batch.blocks.data) continue;
        bool unresolved_only = true;
        for (const auto &info : batch.infos)
            if (info.phase != MdOnePassBlockInfo::Stored &&
                (info.phase != MdOnePassBlockInfo::Waiting || !info.unresolved))
                unresolved_only = false;
        if (!unresolved_only) continue;
        const double diagnostic_t0 = context->pending_diagnostics ? GetTime() : 0.0;
        std::shared_ptr<swbam_usage::CachedTemporary> spool = context->spill_file.lock();
        // Rotate files so retired batches do not pin every earlier spill byte.
        if (!spool || spool->size() >= (64ull << 20)) {
            spool = std::make_shared<swbam_usage::CachedTemporary>();
            context->spill_file = spool;
        }
        for (int b = 0; b < batch.n_blocks; ++b) {
            if (batch.infos[b].phase == MdOnePassBlockInfo::Stored) continue;
            batch.infos[b].spool_offset = spool->size();
            const bam_block &block = batch.blocks.blocks[b];
            if (spool->Append(block.data, block.length)) return -1;
            context->diagnostic_spill_bytes += block.length;
            if (context->pending_diagnostics) {
                if (batch.infos[b].unresolved)
                    ++context->diagnostic_waiting_spilled;
                else
                    ++context->diagnostic_ready_spilled;
            }
        }
        aligned_free_custom(batch.blocks.data);
        batch.blocks.data = nullptr;
        for (int b = 0; b < batch.blocks.n; ++b) batch.blocks.blocks[b].data = nullptr;
        context->resident_pending -= (size_t)batch.blocks.n * BGZF_MAX_BLOCK_SIZE;
        batch.spool = std::move(spool);
        ++context->spilled_batches;
        if (context->pending_diagnostics) {
            if (batch.diagnostic_spills++ == 0)
                ++context->diagnostic_unique_spills;
            context->diagnostic_spill_seconds += GetTime() - diagnostic_t0;
        }
    }
    return 0;
}

static int MdOnePassLoadBlock(MdOnePassContext *context,
                              MdOnePassPendingBatch *batch, int index) {
    bam_block &block = batch->blocks.blocks[index];
    if (block.data) return 0;
    if (!batch->spool) return -1;
    const double diagnostic_t0 = context->pending_diagnostics ? GetTime() : 0.0;
    block.data = aligned_alloc_custom(64, BGZF_MAX_BLOCK_SIZE);
    if (!block.data) return -1;
    context->resident_pending += BGZF_MAX_BLOCK_SIZE;
    if (batch->spool->ReadAt(batch->infos[index].spool_offset,
                            block.data, block.length)) return -1;
    context->diagnostic_reload_bytes += block.length;
    if (context->pending_diagnostics) {
        ++context->diagnostic_reloads;
        context->diagnostic_reload_seconds += GetTime() - diagnostic_t0;
    }
    return 0;
}

static void MdOnePassDiagnosePending(MdOnePassContext *context,
                                      int progress_tid, int progress_pos) {
    if (!context->pending_diagnostics) return;
    const long long batch_index = context->processed_batches + 1;
    if (batch_index != 1 && batch_index % 256 != 0 &&
        !(context->spilled_batches && !context->diagnostic_last_reported_spills)) return;
    size_t ready = 0, waiting = 0, ready_behind = 0, spooled = 0, unknown = 0;
    bool blocked = false;
    uint64_t first_unknown = UINT64_MAX;
    for (const auto &batch : context->pending) {
        if (batch.spool) ++spooled;
        for (int b = 0; b < batch.n_blocks; ++b) {
            const MdOnePassBlockInfo &info = batch.infos[b];
            if (info.phase == MdOnePassBlockInfo::Stored) continue;
            if (info.unresolved) {
                ++waiting;
                unknown += (size_t)info.unresolved;
                if (!blocked) {
                    const size_t begin = (size_t)(info.ordinal_base - batch.ordinal_base);
                    for (int r = 0; r < info.n_records; ++r)
                        if (batch.states[begin + (size_t)r] == kMdOnePassUnknown) {
                            first_unknown = info.ordinal_base + (uint64_t)r;
                            break;
                        }
                }
                blocked = true;
            } else {
                ++ready;
                if (blocked) ++ready_behind;
            }
        }
    }
    fprintf(stderr, "[md-pending rank=%d] batch=%lld progress=%d:%d "
            "pending_batches=%zu resident=%zu spooled_batches=%zu "
            "ready_blocks=%zu waiting_blocks=%zu ready_behind=%zu unknown=%zu "
            "first_unknown=%llu spills=%zu unique_spills=%zu reloads=%zu\n",
            context->rank, batch_index, progress_tid, progress_pos,
            context->pending.size(), context->resident_pending, spooled,
            ready, waiting, ready_behind, unknown,
            (unsigned long long)first_unknown, context->spilled_batches,
            context->diagnostic_unique_spills, context->diagnostic_reloads);
    fflush(stderr);
    context->diagnostic_last_reported_spills = context->spilled_batches;
}

static void MdOnePassTouchBlock(MdOnePassContext *context,
                                MdOnePassPendingBatch *batch,
                                MdOnePassBlockInfo *info) {
    if (info->dirty) return;
    info->dirty = true;
    context->dirty_blocks.push_back({batch, (int)(info - batch->infos.data())});
}

static int MdOnePassQueueReady(MdOnePassContext *context) {
    for (const auto &entry : context->dirty_blocks) {
        MdOnePassBlockInfo &info = entry.batch->infos[entry.block_index];
        info.dirty = false;
        if (info.phase != MdOnePassBlockInfo::Waiting) return -1;
        if (!info.unresolved) {
            info.phase = MdOnePassBlockInfo::Ready;
            if (!context->ready.emplace(info.block_order, entry).second) return -1;
        }
    }
    context->dirty_blocks.clear();
    return 0;
}

static int MdOnePassFindRecord(
        MdOnePassContext *context, uint64_t ordinal,
        MdOnePassPendingBatch **batch_out,
        MdOnePassBlockInfo **info_out, uint8_t **state_out) {
    if (context->use_pending_index) {
        MdOnePassPendingBatch *batch = context->pending_index.Find(ordinal);
        if (!batch) return -1;
        auto info = std::upper_bound(batch->infos.begin(), batch->infos.end(), ordinal,
            [](uint64_t value, const MdOnePassBlockInfo &block) {
                return value < block.ordinal_base;
            });
        if (info == batch->infos.begin()) return -1;
        --info;
        if (ordinal - info->ordinal_base >= (uint64_t)info->n_records) return -1;
        if (info->phase != MdOnePassBlockInfo::Waiting) {
            fprintf(stderr, "[rank %d] ERROR: late markdup decision for sealed block=%llu.\n",
                    context->rank, (unsigned long long)info->block_order);
            return -1;
        }
        *batch_out = batch;
        *info_out = &*info;
        *state_out = &batch->states[(size_t)(ordinal - batch->ordinal_base)];
        MdOnePassTouchBlock(context, batch, &*info);
        return 0;
    }
    for (std::list<MdOnePassPendingBatch>::iterator it =
             context->pending.begin();
         it != context->pending.end(); ++it) {
        if (ordinal < it->ordinal_base || ordinal >= it->ordinal_end) {
            continue;
        }
        const size_t state_index = (size_t)(ordinal - it->ordinal_base);
        for (int b = 0; b < it->n_blocks; ++b) {
            MdOnePassBlockInfo &info = it->infos[(size_t)b];
            if (ordinal >= info.ordinal_base &&
                ordinal < info.ordinal_base +
                    (uint64_t)info.n_records) {
                *batch_out = &*it;
                *info_out = &info;
                *state_out = &it->states[state_index];
                if (info.phase != MdOnePassBlockInfo::Waiting) {
                    fprintf(stderr, "[rank %d] ERROR: late markdup decision for sealed block=%llu.\n",
                            context->rank, (unsigned long long)info.block_order);
                    return -1;
                }
                MdOnePassTouchBlock(context, &*it, &info);
                return 0;
            }
        }
        return -1;
    }
    return -1;
}

static int MdOnePassMarkUnknown(MdOnePassContext *context,
                                uint64_t ordinal) {
    MdOnePassPendingBatch *batch = nullptr;
    MdOnePassBlockInfo *info = nullptr;
    uint8_t *state = nullptr;
    if (MdOnePassFindRecord(context, ordinal, &batch, &info,
                            &state) != 0) {
        fprintf(stderr,
                "[rank %d] ERROR: one-pass cannot find remote-owned "
                "ordinal=%llu pending_batches=%zu.\n",
                context->rank, (unsigned long long)ordinal,
                context->pending.size());
        return -1;
    }
    (void)batch;
    if (*state == kMdOnePassKeep) {
        *state = kMdOnePassUnknown;
        info->unresolved++;
        return 0;
    }
    return *state == kMdOnePassUnknown ? 0 : -1;
}

static int MdOnePassObserveRemoteCandidate(
        void *opaque,
        const MpiMarkdupCandidateShared &candidate,
        int owner) {
    MdOnePassContext *context = (MdOnePassContext *)opaque;
    if (!context || owner < 0 || owner >= context->comm_size) return -1;
    ++context->remote_candidates;
    if (candidate.key.single && candidate.paired_marker) {
        ++context->remote_markers;
        return 0;
    }
    if (MdOnePassMarkUnknown(context, candidate.ordinal) == 0) return 0;
    fprintf(stderr,
            "[rank %d] ERROR: one-pass cannot pend remote candidate "
            "ordinal=%llu owner=%d.\n",
            context->rank, (unsigned long long)candidate.ordinal,
            owner);
    return -1;
}

static int MdOnePassApplyDecision(
        MdOnePassContext *context,
        const MpiMarkdupDecisionShared &decision) {
    MdOnePassPendingBatch *batch = nullptr;
    MdOnePassBlockInfo *info = nullptr;
    uint8_t *state = nullptr;
    if (MdOnePassFindRecord(context, decision.ordinal, &batch, &info,
                            &state) != 0) {
        fprintf(stderr,
                "[rank %d] ERROR: one-pass remote decision cannot find "
                "ordinal=%llu duplicate=%u pending_batches=%zu.\n",
                context->rank,
                (unsigned long long)decision.ordinal,
                (unsigned)decision.duplicate,
                context->pending.size());
        return -1;
    }
    (void)batch;
    const uint8_t final_state = decision.duplicate
        ? kMdOnePassDuplicate : kMdOnePassKeep;
    if (*state == kMdOnePassUnknown) {
        *state = final_state;
        info->unresolved--;
        return info->unresolved >= 0 ? 0 : -1;
    }
    if (*state == final_state) return 0;
    fprintf(stderr,
            "[rank %d] ERROR: one-pass conflicting remote decision "
            "ordinal=%llu duplicate=%u state=%u.\n",
            context->rank, (unsigned long long)decision.ordinal,
            (unsigned)decision.duplicate, (unsigned)*state);
    return -1;
}

static int MdOnePassApplyLocalAction(
        void *opaque, uint64_t ordinal, int action) {
    MdOnePassContext *context = (MdOnePassContext *)opaque;
    if (!context) return -1;
    MdOnePassPendingBatch *batch = nullptr;
    MdOnePassBlockInfo *info = nullptr;
    uint8_t *state = nullptr;
    if (MdOnePassFindRecord(context, ordinal, &batch, &info,
                            &state) != 0) {
        fprintf(stderr,
                "[rank %d] ERROR: one-pass local action cannot find "
                "ordinal=%llu action=%d pending_batches=%zu.\n",
                context->rank, (unsigned long long)ordinal, action,
                context->pending.size());
        return -1;
    }
    (void)batch;
    if (action == RB_MARKDUP_DECISION_PENDING) {
        if (*state == kMdOnePassKeep) {
            *state = kMdOnePassUnknown;
            info->unresolved++;
            return 0;
        }
        if (*state == kMdOnePassUnknown) return 0;
        fprintf(stderr,
                "[rank %d] ERROR: one-pass cannot pend ordinal=%llu "
                "state=%u.\n", context->rank,
                (unsigned long long)ordinal, (unsigned)*state);
        return -1;
    }
    const uint8_t final_state =
        action == RB_MARKDUP_DECISION_DUPLICATE
            ? kMdOnePassDuplicate : kMdOnePassKeep;
    if (*state == kMdOnePassUnknown) {
        *state = final_state;
        info->unresolved--;
        return info->unresolved >= 0 ? 0 : -1;
    }
    if (action == RB_MARKDUP_DECISION_DUPLICATE &&
        *state == kMdOnePassKeep) {
        *state = kMdOnePassDuplicate;
        return 0;
    }
    if (*state == final_state) return 0;
    fprintf(stderr,
            "[rank %d] ERROR: one-pass conflicting local action "
            "ordinal=%llu action=%d state=%u.\n",
            context->rank, (unsigned long long)ordinal, action,
            (unsigned)*state);
    return -1;
}

static int MdOnePassCaptureBatch(
        MdOnePassContext *context, MdExtractBatchView *view,
        const std::vector<MpiMarkdupCandidateShared> &candidates) {
    if (!view || view->n_blocks <= 0) return 0;
    if (!view->uncompressed || !view->uncompressed->blocks ||
        !view->decomp || view->ordinal_end < view->ordinal_base) {
        return -1;
    }
    context->extract_workspace_bytes = view->fixed_workspace_bytes;
    // Metadata remains necessary for late decisions, even after payload spill.
    // Reject before allocating another batch when that irreducible state fills
    // its share of the budget; never silently drop candidates or records.
    size_t metadata = 0;
    for (const auto &entry : context->pending) {
        const size_t n = sizeof(entry) + entry.infos.capacity() * sizeof(MdOnePassBlockInfo) +
            entry.states.capacity() + (size_t)entry.blocks.n * sizeof(bam_block);
        if (n > SIZE_MAX - metadata) return -1;
        metadata += n;
    }
    if (context->use_pending_index) {
        const size_t index_bytes = context->pending_index.memory_bytes();
        if (index_bytes > SIZE_MAX - metadata) return -1;
        metadata += index_bytes;
    }
    const size_t next_metadata = sizeof(MdOnePassPendingBatch) +
        (size_t)view->n_blocks * (sizeof(MdOnePassBlockInfo) + sizeof(bam_block)) +
        (size_t)(view->ordinal_end - view->ordinal_base) +
        (context->use_pending_index
            ? swbam_markdup::PendingIndex<MdOnePassPendingBatch>::EntryBytes() : 0);
    if (context->memory_limit && (metadata > context->memory_limit / 4 ||
        next_metadata > context->memory_limit / 4 - metadata)) {
        fprintf(stderr, "[rank %d] ERROR: markdup deferred decision metadata exceeds work budget.\n", context->rank);
        return -1;
    }

    MdBlockSet replacement = {};
    if (MdOnePassAcquireBlockSet(context, &replacement) != 0) {
        return -1;
    }
    MdBlockSet retained = *view->uncompressed;
    *view->uncompressed = replacement;

    try {
        context->pending.emplace_back();
    } catch (...) {
        MdFreeBlockSet(view->uncompressed);
        *view->uncompressed = retained;
        return -1;
    }
    MdOnePassPendingBatch &batch = context->pending.back();
    batch.blocks = retained;
    context->resident_pending += (size_t)retained.n * BGZF_MAX_BLOCK_SIZE;
    batch.ordinal_base = view->ordinal_base;
    batch.ordinal_end = view->ordinal_end;
    batch.n_blocks = view->n_blocks;
    batch.remaining = view->n_blocks;
    try {
        batch.infos.resize((size_t)view->n_blocks);
        batch.states.assign(
            (size_t)(view->ordinal_end - view->ordinal_base),
            (uint8_t)kMdOnePassKeep);
    } catch (...) {
        context->pending.pop_back();
        return -1;
    }

    uint64_t ordinal = view->ordinal_base;
    for (int b = 0; b < view->n_blocks; ++b) {
        MdOnePassBlockInfo &info = batch.infos[(size_t)b];
        info.ordinal_base = ordinal;
        info.n_records = view->decomp[b].n_total_records;
        info.unresolved = 0;
        info.block_order = (uint64_t)context->captured_blocks + (uint64_t)b;
        info.phase = MdOnePassBlockInfo::Waiting;
        info.dirty = false;
        MdOnePassTouchBlock(context, &batch, &info);
        ordinal += (uint64_t)info.n_records;
    }
    if (ordinal != view->ordinal_end) return -1;
    if (context->use_pending_index) {
        try {
            if (!context->pending_index.Insert(&batch)) return -1;
            context->pending_index_peak = std::max(context->pending_index_peak,
                                                   context->pending_index.size());
        } catch (const std::bad_alloc &) {
            return -1;
        }
    }
    context->captured_blocks += view->n_blocks;
    context->captured_records += view->ordinal_end - view->ordinal_base;

    if (context->boundary_plan &&
        context->boundary_plan->MayContainBatch(candidates)) {
        for (size_t i = 0; i < candidates.size(); ++i) {
            const MpiMarkdupCandidateShared &candidate = candidates[i];
            if (!MdBoundaryDecisionCandidate(candidate)) {
                continue;
            }
            const MdBoundaryCandidateId id =
                MdMakeBoundaryCandidateId(candidate);
            if (context->boundary_plan->candidates.find(id) ==
                context->boundary_plan->candidates.end()) {
                continue;
            }
            std::unordered_map<uint64_t, uint8_t>::const_iterator decision =
                context->boundary_plan->decisions.find(
                    candidate.global_order);
            if (decision == context->boundary_plan->decisions.end() ||
                candidate.ordinal < batch.ordinal_base ||
                candidate.ordinal >= batch.ordinal_end) {
                return -1;
            }
            batch.states[(size_t)(candidate.ordinal -
                                  batch.ordinal_base)] =
                decision->second;
        }
    }

    return 0;
}

static int MdOnePassPatchBlock(
        MdOnePassPendingBatch *batch, int block_index,
        MpiMarkdupStats *stats) {
    if (!batch || block_index < 0 || block_index >= batch->n_blocks) {
        return -1;
    }
    MdOnePassBlockInfo &info = batch->infos[(size_t)block_index];
    if (info.unresolved != 0) return -1;
    bam_block *block = batch->blocks.blocks + block_index;
    size_t pos = 0;
    for (int r = 0; r < info.n_records; ++r) {
        if (pos + 36u > block->length) return -1;
        unsigned char *raw = block->data + pos;
        const uint32_t block_len = MdLoadLe32(raw);
        const uint64_t raw_len = (uint64_t)block_len + 4u;
        if (block_len < 32 || raw_len > block->length - pos) {
            return -1;
        }
        const size_t state_index = (size_t)(
            info.ordinal_base + (uint64_t)r - batch->ordinal_base);
        const uint8_t state = batch->states[state_index];
        if (state == kMdOnePassUnknown) return -1;
        if (state == kMdOnePassDuplicate) {
            uint32_t flag_nc = MdLoadLe32(raw + 16);
            flag_nc |= (uint32_t)BAM_FDUP << 16;
            memcpy(raw + 16, &flag_nc, sizeof(flag_nc));
            stats->marked_records++;
        }
        pos += (size_t)raw_len;
    }
    if (pos != block->length) return -1;
    block->pos = block->length;
    return 0;
}

static int MdOnePassPendingBytes(MdOnePassContext *context,
                                 size_t *bytes_out) {
    size_t total = 0;
    for (std::list<MdOnePassPendingBatch>::const_iterator it =
             context->pending.begin();
         it != context->pending.end(); ++it) {
        const size_t block_bytes =
            (size_t)it->blocks.n *
            (sizeof(bam_block) + (it->blocks.data ? BGZF_MAX_BLOCK_SIZE : 0));
        const size_t info_bytes =
            it->infos.capacity() * sizeof(MdOnePassBlockInfo);
        const size_t state_bytes = it->states.capacity();
        if (block_bytes > SIZE_MAX - total ||
            info_bytes > SIZE_MAX - total - block_bytes ||
            state_bytes > SIZE_MAX - total - block_bytes - info_bytes) {
            return -1;
        }
        total += block_bytes + info_bytes + state_bytes;
        if (!it->blocks.data) {
            for (int b = 0; b < it->blocks.n; ++b) {
                if (!it->blocks.blocks[b].data) continue;
                if (BGZF_MAX_BLOCK_SIZE > SIZE_MAX - total) return -1;
                total += BGZF_MAX_BLOCK_SIZE;
            }
        }
    }
    const size_t queue_bytes = context->ready.size() *
        (sizeof(std::pair<const uint64_t, MdOnePassReadyBlock>) + 4 * sizeof(void *)) +
        context->dirty_blocks.capacity() * sizeof(MdOnePassReadyBlock) +
        context->ordered_body->memory_bytes() + context->pending_index.memory_bytes();
    if (queue_bytes > SIZE_MAX - total) return -1;
    total += queue_bytes;
    const size_t pool_block_bytes =
        (sizeof(bam_block) + BGZF_MAX_BLOCK_SIZE) * kMarkdupNB;
    if (context->free_block_sets.size() >
            (SIZE_MAX - total) / pool_block_bytes) {
        return -1;
    }
    total += context->free_block_sets.size() * pool_block_bytes;
    *bytes_out = total;
    return 0;
}

// A legal BGZF block may inflate to 64 KiB, while the CPE compressor accepts
// at most BGZF_BLOCK_SIZE bytes. Repack only that uncommon case, preserving
// record boundaries and the original record order.
static int MdOnePassPrepareOversizedBlock(
        MdOnePassContext *context, bam_block *source, int n_records,
        int *output_count) {
    if (!context || !source || !output_count || n_records <= 0 ||
        *output_count < 0 || *output_count > kMarkdupNB - 2) {
        return -1;
    }
    int active = *output_count;
    MdResetPayloadBlock(context->fallback_payload.blocks + active, active);
    size_t pos = 0;
    for (int r = 0; r < n_records; ++r) {
        if (pos + 36u > source->length) return -1;
        const unsigned char *raw = source->data + pos;
        const uint32_t block_len = MdLoadLe32(raw);
        const uint64_t raw_len64 = (uint64_t)block_len + 4u;
        if (block_len < 32 || raw_len64 > BGZF_BLOCK_SIZE ||
            raw_len64 > (uint64_t)source->length - pos) {
            return -1;
        }
        int append_status = MdAppendRawPayload(
            context->fallback_payload.blocks + active, raw,
            (uint32_t)raw_len64);
        if (append_status == 1) {
            context->prepared_payload[active] =
                context->fallback_payload.blocks[active];
            ++active;
            if (active >= kMarkdupNB) return -1;
            MdResetPayloadBlock(
                context->fallback_payload.blocks + active, active);
            append_status = MdAppendRawPayload(
                context->fallback_payload.blocks + active, raw,
                (uint32_t)raw_len64);
        }
        if (append_status != 0) return -1;
        pos += (size_t)raw_len64;
    }
    if (pos != source->length ||
        context->fallback_payload.blocks[active].pos <= 0) {
        return -1;
    }
    context->prepared_payload[active] = context->fallback_payload.blocks[active];
    *output_count = active + 1;
    return 0;
}

static int MdOnePassPrepareOutput(MdOnePassContext *context) {
    if (context->prepared_outputs != 0 ||
        !context->prepared_inputs.empty()) {
        return -1;
    }
    int output_count = 0;
    double patch_t0 = GetTime();
    auto commit_prepared = [&]() -> int {
        context->prepared_outputs = output_count;
        context->stats->t_pack += GetTime() - patch_t0;
        return 0;
    };
    while (!context->ready.empty() && output_count < kMarkdupNB) {
        const MdOnePassReadyBlock ready = context->ready.begin()->second;
        MdOnePassPendingBatch *batch = ready.batch;
        const int b = ready.block_index;
        MdOnePassBlockInfo &info = batch->infos[(size_t)b];
        if (info.unresolved || info.phase != MdOnePassBlockInfo::Ready) return -1;
        bam_block *source = batch->blocks.blocks + b;
        if (source->length > BGZF_BLOCK_SIZE && output_count > kMarkdupNB - 2)
            return commit_prepared();
        if (!context->reordered && info.block_order != context->ordered_next) {
            if (context->ordered_body->Activate()) return -1;
            context->reordered = true;
        }
        context->ordered_next = info.block_order + 1;
        if (source->length == 0 && info.n_records == 0) {
            info.phase = MdOnePassBlockInfo::Compressing;
            context->prepared_inputs.push_back(ready);
            context->ready.erase(context->ready.begin());
            continue;
        }
        if (MdOnePassLoadBlock(context, batch, b)) return -1;
        if (MdOnePassPatchBlock(batch, b, context->stats) != 0) {
            context->stats->t_pack += GetTime() - patch_t0;
            return -1;
        }
        const int first_output = output_count;
        if (source->length > BGZF_BLOCK_SIZE) {
            if (MdOnePassPrepareOversizedBlock(
                    context, source, info.n_records,
                    &output_count) != 0) {
                fprintf(stderr,
                        "[rank %d] ERROR: one-pass markdup cannot "
                        "repack decoded BGZF payload=%u.\n",
                        context->rank, source->length);
                context->stats->t_pack += GetTime() - patch_t0;
                return -1;
            }
        } else {
            context->prepared_payload[output_count] = *source;
            output_count++;
        }
        for (int p = first_output; p < output_count; ++p)
            context->prepared_keys[p] = info.block_order * 2 + (uint64_t)(p - first_output);
        info.phase = MdOnePassBlockInfo::Compressing;
        context->prepared_inputs.push_back(ready);
        context->ready.erase(context->ready.begin());
    }
    return commit_prepared();
}

static int MdOnePassLaunchOutput(MdOnePassContext *context) {
    if (context->compression_running || context->prepared_outputs == 0) {
        return 0;
    }
    const double launch_t0 = GetTime();
    const int ret = context->write_session->Start(
        swbam::cpe::CpeWriteBatchInput(&context->prepared_view,
                                     (size_t)context->prepared_outputs));
    context->stats->t_cpe_launch += GetTime() - launch_t0;
    if (ret != 0) return -1;
    context->compression_running = 1;
    context->launched_output_blocks += context->prepared_outputs;
    return 0;
}

static int MdOnePassFinishOutput(MdOnePassContext *context) {
    if (!context->compression_running && context->prepared_outputs > 0 &&
        MdOnePassLaunchOutput(context) != 0) {
        return -1;
    }
    const bool had_compression = context->compression_running != 0;
    if (context->compression_running) {
        const double compress_t0 = GetTime();
        const double sync_t0 = GetTime();
        const int ret = context->write_session->Complete();
        context->stats->t_cpe_sync += GetTime() - sync_t0;
        context->stats->t_compress += GetTime() - compress_t0;
        context->compression_running = 0;
        if (ret != 0 || context->write_session->Flush() != 0) return -1;
    }

    for (size_t i = 0; i < context->prepared_inputs.size(); ++i) {
        if (context->pending.empty()) return -1;
        const MdOnePassReadyBlock &ready = context->prepared_inputs[i];
        MdOnePassPendingBatch &front = *ready.batch;
        MdOnePassBlockInfo &info = front.infos[(size_t)ready.block_index];
        if (info.phase != MdOnePassBlockInfo::Compressing) return -1;
        if (!had_compression &&
            front.infos[(size_t)ready.block_index].n_records != 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: retiring a non-empty one-pass "
                    "block without a compression launch: block=%d "
                    "records=%d raw=%u prepared_inputs=%zu "
                    "prepared_outputs=%d pending=%zu.\n",
                    context->rank, ready.block_index,
                    front.infos[(size_t)ready.block_index].n_records,
                    front.blocks.blocks[ready.block_index].length,
                    context->prepared_inputs.size(),
                    context->prepared_outputs,
                    context->pending.size());
            return -1;
        }
        context->retired_blocks++;
        context->retired_records +=
            (uint64_t)front.infos[(size_t)ready.block_index].n_records;
        info.phase = MdOnePassBlockInfo::Stored;
        --front.remaining;
        if (!front.blocks.data && front.blocks.blocks[ready.block_index].data) {
            aligned_free_custom(front.blocks.blocks[ready.block_index].data);
            front.blocks.blocks[ready.block_index].data = nullptr;
            context->resident_pending -= BGZF_MAX_BLOCK_SIZE;
        }
        if (!front.remaining) {
            if (front.blocks.data)
                context->resident_pending -= (size_t)front.blocks.n * BGZF_MAX_BLOCK_SIZE;
            if (front.blocks.data) {
                MdBlockSet recycled = front.blocks;
                memset(&front.blocks, 0, sizeof(front.blocks));
                MdOnePassRecycleBlockSet(context, &recycled);
            }
        }
    }
    context->prepared_inputs.clear();
    context->prepared_outputs = 0;
    context->pending.remove_if([context](const MdOnePassPendingBatch &batch) {
        if (batch.remaining != 0) return false;
        if (context->use_pending_index) context->pending_index.Erase(&batch);
        return true;
    });
    return 0;
}

static int MdOnePassRelievePressure(MdOnePassContext *context) {
    const size_t limit = swbam_usage::CacheMemoryLimit();
    const auto pressured = [&](size_t target) {
        return context->resident_pending > target ||
            context->ready.size() > target / BGZF_MAX_BLOCK_SIZE;
    };
    if (!pressured(limit)) return 0;
    // Drain available work before spilling. An unresolved window must still
    // allow input progress, so stop as soon as the ready queue is empty.
    while (pressured(limit / 2) && !context->ready.empty()) {
        if (MdOnePassPrepareOutput(context) || MdOnePassLaunchOutput(context) ||
            MdOnePassFinishOutput(context)) return -1;
        ++context->catchup_batches;
    }
    return MdOnePassSpillPending(context);
}

static int MdOnePassRouteDecisions(MdOnePassContext *context) {
    const double exchange_t0 = GetTime();
    std::vector<unsigned long long> send_counts(
        (size_t)context->comm_size, 0);
    std::vector<unsigned long long> recv_counts(
        (size_t)context->comm_size, 0);

    int local_ok = 1;
    for (int dest = 0; dest < context->comm_size; ++dest) {
        std::vector<MpiMarkdupDecisionShared> &values =
            context->decisions_by_source[(size_t)dest];
        std::sort(values.begin(), values.end(),
                  [](const MpiMarkdupDecisionShared &a,
                     const MpiMarkdupDecisionShared &b) {
                      if (a.ordinal != b.ordinal) {
                          return a.ordinal < b.ordinal;
                      }
                      return a.duplicate > b.duplicate;
                  });
        size_t write = 0;
        for (size_t i = 0; i < values.size(); ++i) {
            if (write > 0 &&
                values[write - 1].ordinal == values[i].ordinal) {
                values[write - 1].duplicate =
                    values[write - 1].duplicate || values[i].duplicate;
            } else {
                values[write++] = values[i];
            }
        }
        values.resize(write);
        if (dest == context->rank) {
            for (size_t i = 0; i < values.size(); ++i) {
                if (MdOnePassApplyDecision(context, values[i]) != 0) {
                    local_ok = 0;
                    break;
                }
            }
        } else {
            send_counts[(size_t)dest] =
                (unsigned long long)values.size();
        }
        if (!local_ok) break;
    }

    if (!MdAllRanksOk(local_ok)) return -1;

    MPI_Alltoall(send_counts.data(), 1, MPI_UNSIGNED_LONG_LONG,
                 recv_counts.data(), 1, MPI_UNSIGNED_LONG_LONG,
                 MPI_COMM_WORLD);
    for (int step = 1; step < context->comm_size; ++step) {
        const int dest = (context->rank + step) % context->comm_size;
        const int src = (context->rank - step + context->comm_size) %
                        context->comm_size;
        const unsigned long long send_count = send_counts[(size_t)dest];
        const unsigned long long recv_count = recv_counts[(size_t)src];
        if (send_count > ULLONG_MAX /
                sizeof(MpiMarkdupDecisionShared) ||
            recv_count > (unsigned long long)(SIZE_MAX /
                sizeof(MpiMarkdupDecisionShared))) {
            return -1;
        }
        std::vector<MpiMarkdupDecisionShared> received;
        try {
            received.resize((size_t)recv_count);
        } catch (...) {
            local_ok = 0;
        }
        if (!MdAllRanksOk(local_ok)) return -1;
        const unsigned long long send_bytes =
            send_count * sizeof(MpiMarkdupDecisionShared);
        const unsigned long long recv_bytes =
            recv_count * sizeof(MpiMarkdupDecisionShared);
        if (MdSendrecvBytes(
                send_count ? context->decisions_by_source[(size_t)dest].data()
                           : nullptr,
                send_bytes, dest,
                recv_count ? received.data() : nullptr,
                recv_bytes, src, 4131) != 0) {
            return -1;
        }
        context->stats->mpi_result_bytes += (long long)send_bytes;
        for (size_t i = 0; i < received.size(); ++i) {
            if (MdOnePassApplyDecision(context, received[i]) != 0) {
                local_ok = 0;
                break;
            }
        }
        if (!MdAllRanksOk(local_ok)) return -1;
    }
    for (size_t i = 0; i < context->decisions_by_source.size(); ++i) {
        context->decisions_by_source[i].clear();
    }
    context->stats->t_result_exchange += GetTime() - exchange_t0;
    return 0;
}

static int MdOnePassInitWindow(MdOnePassContext *context) {
    if (context->window) return 0;
    if (!context->range_has_records || !context->range_first_tid ||
        !context->range_first_pos) {
        return -1;
    }
    const bool have_boundary_plan =
        context->boundary_plan && context->boundary_plan->enabled;
    const bool conservative = !have_boundary_plan || context->boundary_plan->conservative;
    context->use_pending_index = conservative;
    if (have_boundary_plan) {
        context->boundaries = context->boundary_plan->boundaries;
    } else if (MdBuildFirstCoordBoundaries(
                   *context->range_has_records,
                   *context->range_first_tid,
                   *context->range_first_pos,
                   context->comm_size, &context->boundaries) != 0) {
        return -1;
    }
    int window_last_tid = -1;
    int window_last_pos = -1;
    if (conservative) {
        for (int r = context->rank + 1; r < context->comm_size; ++r) {
            const MdCoordBoundary &next = context->boundaries[(size_t)r];
            if (!next.has_records) continue;
            window_last_tid = next.first_ref - 1;
            window_last_pos = (int)(next.first_coord - 1);
            break;
        }
    }
    context->window = MpiMarkdupStreamingWindowCreate(
        context->comm_size, context->max_read_length,
        conservative ? *context->range_first_tid : -1,
        conservative ? *context->range_first_pos : -1,
        window_last_tid, window_last_pos);
    return context->window ? 0 : -1;
}

static int MdOnePassGatherBoundaryProgress(MdOnePassContext *context, int local_ok,
                                          bool active, int tid, int pos) {
    const long long local[4] = {local_ok, active ? 1ll : 0ll,
                               (long long)tid + 1, (long long)pos + 1};
    context->boundary_progress.resize((size_t)context->comm_size * 4);
    const double t0 = GetTime();
    const int ret = MPI_Allgather(local, 4, MPI_LONG_LONG,
        context->boundary_progress.data(), 4, MPI_LONG_LONG, MPI_COMM_WORLD);
    context->stats->t_rank_sync += GetTime() - t0;
    if (ret != MPI_SUCCESS) return 0;
    for (int r = 0; r < context->comm_size; ++r)
        if (!context->boundary_progress[(size_t)r * 4]) return 0;
    const MdCoordBoundary &first = context->boundaries[(size_t)context->rank];
    context->release_first_boundary = swbam_markdup::BoundaryIsComplete(
        first.first_ref, first.first_coord, context->max_read_length,
        context->boundaries, context->boundary_progress);
    for (int r = context->rank + 1; r < context->comm_size; ++r) {
        const MdCoordBoundary &last = context->boundaries[(size_t)r];
        if (!last.has_records) continue;
        context->release_last_boundary = swbam_markdup::BoundaryIsComplete(
            last.first_ref, last.first_coord, context->max_read_length,
            context->boundaries, context->boundary_progress);
        break;
    }
    return 1;
}

static int MdPostProcessOnePassBatch(
        void *opaque,
        std::vector<MpiMarkdupCandidateShared> *local_candidates,
        std::vector<unsigned char> *local_qnames,
        int progress_tid, int progress_pos,
        uint64_t records_so_far,
        MdExtractBatchView *batch) {
    MdOnePassContext *context = (MdOnePassContext *)opaque;
    if (!context || !context->stats || !context->sink) return -1;

    static const bool diagnostics = MdEnvFlagEnabled("SWBAM_MARKDUP_DIAG");
    const long long batch_index = context->processed_batches + 1;
    const bool trace = diagnostics &&
        (batch_index == 1 || batch_index % 512 == 0 || !batch);
    const auto mark = [&](const char *stage) {
        if (!trace) return;
        fprintf(stderr,
                "[rank %d] markdup_diag batch=%lld stage=%s progress=%d:%d "
                "records=%llu pending=%zu captured=%lld retired=%lld "
                "output_bytes=%llu\n",
                context->rank, batch_index, stage, progress_tid, progress_pos,
                (unsigned long long)records_so_far, context->pending.size(),
                context->captured_blocks, context->retired_blocks,
                (unsigned long long)context->sink->size());
        fflush(stderr);
    };
    mark("begin");

    if (MdOnePassLaunchOutput(context) != 0) return -1;
    int local_ok = MdOnePassInitWindow(context) == 0;
    if (!MdAllRanksOk(local_ok)) {
        MdOnePassFinishOutput(context);
        return -1;
    }
    mark("window-ready");
    double buffer_t0 = GetTime();
    local_ok = MdOnePassCaptureBatch(
        context, batch, *local_candidates) == 0;
    context->stats->t_stream_buffer += GetTime() - buffer_t0;
    // Replace this existing status collective only on fallback. Exchange and
    // process the current candidates before applying the certified watermark.
    const int capture_ok = context->use_pending_index
        ? MdOnePassGatherBoundaryProgress(context, local_ok, batch != nullptr, progress_tid, progress_pos)
        : MdAllRanksOk(local_ok);
    if (!capture_ok) {
        MdOnePassFinishOutput(context);
        return -1;
    }
    mark("capture-done");

    if (context->boundary_plan &&
        context->boundary_plan->MayContainBatch(*local_candidates)) {
        size_t write = 0;
        for (size_t i = 0; i < local_candidates->size(); ++i) {
            const MpiMarkdupCandidateShared &candidate =
                (*local_candidates)[i];
            const bool precomputed =
                context->boundary_plan->candidates.find(
                    MdMakeBoundaryCandidateId(candidate)) !=
                context->boundary_plan->candidates.end();
            if (!precomputed) {
                if (write != i) {
                    (*local_candidates)[write] = candidate;
                }
                ++write;
            }
        }
        local_candidates->resize(write);
    }

    const long long previous_owner_candidates =
        context->stats->owner_candidates;
    const long long previous_qname_bytes = context->stats->qname_bytes;
    mark("candidate-exchange.begin");
    if (MdExchangeCandidatesByCoordinate(
            local_candidates, local_qnames,
            *context->range_has_records,
            *context->range_first_tid,
            *context->range_first_pos,
            -1, -1, context->rank, context->comm_size,
            context->memory_limit,
            &context->owner_candidates, &context->owner_qnames,
            context->stats, &context->boundaries, 1,
            MdOnePassObserveRemoteCandidate, context) != 0) {
        MdOnePassFinishOutput(context);
        return -1;
    }
    context->stats->owner_candidates += previous_owner_candidates;
    context->stats->qname_bytes += previous_qname_bytes;
    mark("candidate-exchange.done");

    // A late contributor to a closed key invalidates the L-based coverage
    // certificate. Fail rather than silently keeping a precomputed winner.
    if (context->boundary_plan && context->boundary_plan->enabled &&
        !context->boundary_plan->closed_owner_keys.empty()) {
        const auto &closed = context->boundary_plan->closed_owner_keys;
        for (const auto &candidate : context->owner_candidates) {
            if (!swbam_markdup::HaloContainsPosition(
                    context->boundary_plan->closed_owner_ranges,
                    candidate.key.this_ref, candidate.key.this_coord)) continue;
            if (std::binary_search(closed.begin(), closed.end(), candidate.key, MdKeyLess)) {
                fprintf(stderr, "[rank %d] ERROR: late candidate for a finalized halo key "
                        "ref=%d coord=%lld source=%d; increase -l to cover "
                        "the alignment and clipping displacement.\n",
                        context->rank, candidate.key.this_ref,
                        (long long)candidate.key.this_coord, candidate.source_rank);
                local_ok = 0;
                break;
            }
        }
    }
    if (!MdAllRanksOk(local_ok)) {
        MdOnePassFinishOutput(context);
        return -1;
    }

    mark("window-process.begin");
    local_ok = MpiMarkdupStreamingWindowProcessOptimistic(
        context->window, context->owner_candidates,
        context->owner_qnames, progress_tid, progress_pos, context->rank,
        MdOnePassApplyLocalAction, context,
        &context->decisions_by_source, context->stats) == 0;
    if (local_ok && context->use_pending_index) {
        const double t0 = GetTime();
        local_ok = MpiMarkdupStreamingWindowReleaseBoundaryOptimistic(
            context->window, context->release_first_boundary, context->release_last_boundary,
            context->rank, MdOnePassApplyLocalAction, context, &context->decisions_by_source) == 0;
        context->stats->t_stream_finalize += GetTime() - t0;
    }
    if (!local_ok) {
        fprintf(stderr,
                "[rank %d] ERROR: one-pass optimistic window failed "
                "at progress=%d:%d.\n",
                context->rank, progress_tid, progress_pos);
    }
    if (!MdAllRanksOk(local_ok)) {
        MdOnePassFinishOutput(context);
        return -1;
    }
    mark("window-process.done");
    mark("decision-exchange.begin");
    if (MdOnePassRouteDecisions(context) != 0) {
        MdOnePassFinishOutput(context);
        return -1;
    }
    mark("decision-exchange.done");
    context->owner_candidates.clear();
    context->owner_qnames.clear();

    mark("compress.begin");
    local_ok = MdOnePassFinishOutput(context) == 0;
    mark("compress.done");
    try {
        if (local_ok) local_ok = MdOnePassQueueReady(context) == 0 &&
            MdOnePassRelievePressure(context) == 0 &&
            MdOnePassPrepareOutput(context) == 0;
    } catch (const std::bad_alloc &) {
        fprintf(stderr, "[rank %d] ERROR: markdup output scheduling allocation failed.\n",
                context->rank);
        local_ok = 0;
    }
    MdOnePassDiagnosePending(context, progress_tid, progress_pos);

    const double memcheck_t0 = GetTime();
    size_t pending_bytes = 0;
    const size_t window_bytes =
        MpiMarkdupStreamingWindowMemory(context->window);
    const size_t output_workspace_bytes =
        (size_t)3 * kMarkdupNB *
        (sizeof(bam_block) + BGZF_MAX_BLOCK_SIZE);
    size_t boundary_plan_bytes = context->boundary_progress.capacity() * sizeof(long long);
    if (context->boundary_plan && context->boundary_plan->enabled) {
        const MdOnePassBoundaryPlan &plan = *context->boundary_plan;
        boundary_plan_bytes +=
            plan.candidates.bucket_count() * sizeof(void *) +
            plan.decisions.bucket_count() * sizeof(void *) +
            plan.candidates.size() *
                (sizeof(MdBoundaryCandidateId) + 2 * sizeof(void *)) +
            plan.decisions.size() *
                (sizeof(std::pair<const uint64_t, uint8_t>) +
                 2 * sizeof(void *)) +
            plan.boundaries.capacity() * sizeof(MdCoordBoundary) +
            plan.halos.capacity() * sizeof(swbam_markdup::HaloInterval) +
            plan.closed_owner_ranges.capacity() * sizeof(swbam_markdup::HaloInterval) +
            plan.closed_owner_keys.capacity() * sizeof(MpiMarkdupKeyShared);
    }
    local_ok = local_ok && MdOnePassPendingBytes(context, &pending_bytes) == 0 &&
        window_bytes != SIZE_MAX &&
        context->extract_workspace_bytes <=
            SIZE_MAX - output_workspace_bytes &&
        pending_bytes <= SIZE_MAX - output_workspace_bytes -
            context->extract_workspace_bytes &&
        window_bytes <= SIZE_MAX - output_workspace_bytes -
            context->extract_workspace_bytes - pending_bytes &&
        boundary_plan_bytes <= SIZE_MAX - output_workspace_bytes -
            context->extract_workspace_bytes - pending_bytes -
            window_bytes;
    if (local_ok) {
        const size_t active_bytes = output_workspace_bytes +
            context->extract_workspace_bytes + pending_bytes +
            window_bytes + boundary_plan_bytes;
        if ((long long)pending_bytes >
            context->stats->stream_pending_peak_bytes) {
            context->stats->stream_pending_peak_bytes =
                pending_bytes > (size_t)LLONG_MAX
                    ? LLONG_MAX : (long long)pending_bytes;
        }
        local_ok = MdCheckMemory(
            active_bytes, context->memory_limit, context->stats,
            context->rank, "one-pass decoded block window") == 0;
    }
    context->stats->t_extract_memcheck += GetTime() - memcheck_t0;
    if (!MdAllRanksOk(local_ok)) return -1;
    ++context->processed_batches;
    mark("done");
    return 0;
}

static int MdDrainOnePassOutput(MdOnePassContext *context) {
    while (!context->pending.empty() ||
           context->prepared_outputs > 0 ||
           !context->prepared_inputs.empty()) {
        if (context->prepared_outputs == 0 &&
            context->prepared_inputs.empty() &&
            MdOnePassPrepareOutput(context) != 0) {
            return -1;
        }
        if (context->prepared_outputs == 0 &&
            !context->prepared_inputs.empty()) {
            if (MdOnePassFinishOutput(context) != 0) return -1;
            continue;
        }
        if (context->prepared_outputs == 0) {
            size_t unknown = 0;
            uint64_t first_unknown = UINT64_MAX;
            for (std::list<MdOnePassPendingBatch>::const_iterator it =
                     context->pending.begin();
                 it != context->pending.end(); ++it) {
                for (size_t i = 0; i < it->states.size(); ++i) {
                    if (it->states[i] == kMdOnePassUnknown) {
                        unknown++;
                        if (first_unknown == UINT64_MAX) {
                            first_unknown = it->ordinal_base + i;
                        }
                    }
                }
            }
            fprintf(stderr,
                    "[rank %d] ERROR: one-pass output blocked by %zu "
                    "unresolved records; first_ordinal=%llu "
                    "pending_batches=%zu.\n",
                    context->rank, unknown,
                    (unsigned long long)first_unknown,
                    context->pending.size());
            return -1;
        }
        if (MdOnePassLaunchOutput(context) != 0 ||
            MdOnePassFinishOutput(context) != 0) {
            return -1;
        }
    }
    return 0;
}

static int MdMarkdupOnePassToSink(
        MdInputPass &source, long long global_block_begin,
        long long local_block_count, int rank, int comm_size,
        int include_fails, int max_read_length, size_t memory_limit,
        int compress_level, swbam::RankBodySink *sink,
        const MdOnePassBoundaryPlan *boundary_plan,
        uint64_t *local_records,
        int *range_has_records,
        int *range_first_tid, int *range_first_pos,
        int *range_last_tid, int *range_last_pos,
        MpiMarkdupStats *stats) {
    if (!sink || (compress_level != 0 && compress_level != 1 &&
                  compress_level != 6)) {
        return -1;
    }
    MdOnePassContext context;
    context.rank = rank;
    context.comm_size = comm_size;
    context.max_read_length = max_read_length;
    context.memory_limit = memory_limit;
    context.range_has_records = range_has_records;
    context.range_first_tid = range_first_tid;
    context.range_first_pos = range_first_pos;
    context.sink = sink;
    context.ordered_body = dynamic_cast<swbam_usage::OrderedBody *>(sink);
    if (!context.ordered_body) return -1;
    context.stats = stats;
    context.initial_candidate_bytes = stats->mpi_candidate_bytes;
    context.initial_result_bytes = stats->mpi_result_bytes;
    context.boundary_plan = boundary_plan;
    context.compress_level = compress_level;
    context.decisions_by_source.resize((size_t)comm_size);
    context.compress_op.reset(new MdOnePassCompressOperator(compress_level));
    context.output_post.reset(new MdOnePassOutput(
        context.ordered_body, stats, &context.written_output_blocks, context.prepared_keys));
    context.write_session.reset(new swbam::cpe::CpeWritePipelineSession);
    if (MdAllocateBlockSet(&context.fallback_payload, kMarkdupNB) != 0 ||
        context.prepared_view.Attach(context.prepared_payload, kMarkdupNB) != 0 ||
        context.write_session->Initialize(context.output_post.get(),
                                          context.compress_op.get()) != 0) {
        MdOnePassCleanup(&context);
        return -1;
    }

    long long max_local_blocks = 0;
    MPI_Allreduce(&local_block_count, &max_local_blocks, 1,
                  MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    const long long target_batches = max_local_blocks > 0
        ? (max_local_blocks + kMarkdupNB - 1) / kMarkdupNB : 0;
    std::vector<MpiMarkdupCandidateShared> local_candidates;
    std::vector<unsigned char> local_qnames;
    source.Reset();
    const int extract_status = MdExtractCandidatesPass(
        source, global_block_begin, rank, include_fails,
        memory_limit, &local_candidates, &local_qnames,
        local_records, range_has_records,
        range_first_tid, range_first_pos,
        range_last_tid, range_last_pos, stats,
        MdPostProcessOnePassBatch, &context, target_batches);
    int local_ok = extract_status == 0;
    if (!MdAllRanksOkTimed(local_ok, stats)) {
        MdOnePassCleanup(&context);
        return -1;
    }

    double boundary_t0 = GetTime();
    if (MdValidateRankBoundaries(
            *range_has_records, *range_first_tid, *range_first_pos,
            *range_last_tid, *range_last_pos,
            rank, comm_size) != 0) {
        local_ok = 0;
    }
    stats->t_boundary_check += GetTime() - boundary_t0;
    if (!MdAllRanksOkTimed(local_ok, stats)) {
        MdOnePassCleanup(&context);
        return -1;
    }

    if (context.prepared_outputs > 0 &&
        MdOnePassLaunchOutput(&context) != 0) {
        local_ok = 0;
    }
    const double finalize_t0 = GetTime();
    if (local_ok && context.window &&
        MpiMarkdupStreamingWindowFinalizeOptimistic(
            context.window, context.rank,
            MdOnePassApplyLocalAction, &context,
            &context.decisions_by_source) != 0) {
        local_ok = 0;
    }
    stats->t_stream_finalize += GetTime() - finalize_t0;
    if (!MdAllRanksOk(local_ok)) {
        MdOnePassFinishOutput(&context);
        MdOnePassCleanup(&context);
        return -1;
    }
    if (context.window && MdOnePassRouteDecisions(&context) != 0) {
        local_ok = 0;
    }
    if (MdOnePassFinishOutput(&context) != 0) local_ok = 0;
    if (local_ok && MdOnePassQueueReady(&context) != 0) local_ok = 0;
    if (local_ok && MdDrainOnePassOutput(&context) != 0) {
        local_ok = 0;
    }
    if (local_ok && context.write_session->Finish() != 0) local_ok = 0;
    if (!MdAllRanksOkTimed(local_ok, stats)) {
        MdOnePassCleanup(&context);
        return -1;
    }
    if (!context.pending.empty()) local_ok = 0;
    if (context.captured_blocks != context.retired_blocks ||
        context.captured_records != context.retired_records ||
        context.captured_records != *local_records) {
        fprintf(stderr,
                "[rank %d] ERROR: one-pass output conservation failed "
                "blocks=%lld/%lld records=%llu/%llu input=%llu.\n",
                rank, context.retired_blocks, context.captured_blocks,
                (unsigned long long)context.retired_records,
                (unsigned long long)context.captured_records,
                (unsigned long long)*local_records);
        local_ok = 0;
    }
    if (context.launched_output_blocks !=
            context.written_output_blocks) {
        fprintf(stderr,
                "[rank %d] ERROR: one-pass compressed output "
                "conservation failed blocks=%lld/%lld.\n",
                rank, context.written_output_blocks,
                context.launched_output_blocks);
        local_ok = 0;
    }
    if (!MdAllRanksOkTimed(local_ok, stats)) {
        MdOnePassCleanup(&context);
        return -1;
    }
    MdOnePassCleanup(&context);
    return 0;
}

static int MdMarkdupOnePassToSink(
        MemReader &reader, long long global_block_begin,
        long long local_block_count, int rank, int comm_size,
        int include_fails, int max_read_length, size_t memory_limit,
        int compress_level, swbam::RankBodySink *sink,
        const MdOnePassBoundaryPlan *boundary_plan,
        uint64_t *local_records,
        int *range_has_records,
        int *range_first_tid, int *range_first_pos,
        int *range_last_tid, int *range_last_pos,
        MpiMarkdupStats *stats) {
    MdInputPass source(&reader);
    return MdMarkdupOnePassToSink(
        source, global_block_begin, local_block_count,
        rank, comm_size, include_fails, max_read_length,
        memory_limit, compress_level, sink, boundary_plan, local_records,
        range_has_records, range_first_tid, range_first_pos,
        range_last_tid, range_last_pos, stats);
}

static int MdMarkdupOnePassToSink(
        const swbam::BamInputBackend &backend,
        const swbam::BgzfBlockSpan *spans, size_t span_count,
        long long global_block_begin, int rank, int comm_size,
        int include_fails, int max_read_length, size_t memory_limit,
        int compress_level, swbam::RankBodySink *sink,
        const MdOnePassBoundaryPlan *boundary_plan,
        uint64_t *local_records,
        int *range_has_records,
        int *range_first_tid, int *range_first_pos,
        int *range_last_tid, int *range_last_pos,
        MpiMarkdupStats *stats) {
    MdInputPass source(&backend, spans, span_count);
    return MdMarkdupOnePassToSink(
        source, global_block_begin, (long long)span_count,
        rank, comm_size, include_fails, max_read_length,
        memory_limit, compress_level, sink, boundary_plan, local_records,
        range_has_records, range_first_tid, range_first_pos,
        range_last_tid, range_last_pos, stats);
}

struct MdOutputLayout {
    size_t size;
    char *header_data;
    size_t header_size;
    long long total_body_size;
    std::vector<long long> body_sizes;
    std::vector<long long> body_offsets;
};

struct MdOutputStageCosts {
    double stage44;
    double stage45_header;
    double stage46;
};

static int MdPrepareOutputLayout(
        const swbam::RankBodySink &local_body, sam_hdr_t *header,
        int compress_level, int rank, int comm_size,
        MpiMarkdupStats *stats, MdOutputLayout *output,
        MdOutputStageCosts *costs) {
    memset(costs, 0, sizeof(*costs));
    output->size = 0;
    output->header_data = nullptr;
    output->header_size = 0;
    output->total_body_size = 0;
    output->body_sizes.clear();
    output->body_offsets.clear();

    // 4.4 统计各 rank 输出大小 / gather layout
    double stage44_t0 = GetTime();
    if (local_body.size() > (uint64_t)LLONG_MAX) return -1;
    long long local_size = (long long)local_body.size();
    std::vector<long long> sizes((size_t)comm_size, 0);
    MPI_Allgather(&local_size, 1, MPI_LONG_LONG,
                  sizes.data(), 1, MPI_LONG_LONG,
                  MPI_COMM_WORLD);
    std::vector<long long> offsets((size_t)comm_size, 0);
    long long total_body = 0;
    for (int i = 0; i < comm_size; ++i) {
        offsets[(size_t)i] = total_body;
        if (sizes[(size_t)i] < 0 ||
            total_body > LLONG_MAX - sizes[(size_t)i]) {
            return -1;
        }
        total_body += sizes[(size_t)i];
    }
    output->body_sizes = sizes;
    output->body_offsets = offsets;
    output->total_body_size = total_body;
    costs->stage44 = MdReduceMax(GetTime() - stage44_t0);
    if (rank == 0) {
        printf("Complete the 4.4 gather body sizes cost %lf\n",
               costs->stage44);
    }

    char *header_memory = nullptr;
    size_t header_size = 0;
    int local_ok = 1;
    double stage45_header_t0 = 0.0;

    // 4.5a header/layout
    stage45_header_t0 = GetTime();
    if (rank == 0) {
        if (sam_hdr_add_pg(header, "SWBAM",
                           "PN", "SWBAM",
                           "CL", "markdup",
                           NULL) != 0) {
            fprintf(stderr,
                    "ERROR: failed to add markdup @PG header line.\n");
            local_ok = 0;
        }
        if (local_ok &&
            MpiCommonBuildBamHeaderMemory(
                header, compress_level,
                &header_memory, &header_size) != 0) {
            fprintf(stderr,
                    "ERROR: failed to build markdup BAM header.\n");
            local_ok = 0;
        }
        unsigned long long final_size =
            (unsigned long long)header_size +
            (unsigned long long)total_body +
            sizeof(kMarkdupBgzfEofBlock);
        if (local_ok &&
            (final_size > (unsigned long long)SIZE_MAX)) {
            fprintf(stderr,
                    "ERROR: markdup output exceeds addressable memory.\n");
            local_ok = 0;
        }
        if (local_ok) {
            output->size = (size_t)final_size;
            output->header_data = header_memory;
            output->header_size = header_size;
            header_memory = nullptr;
        }
    }
    local_ok = MdAllRanksOk(local_ok);
    costs->stage45_header =
        MdReduceMax(GetTime() - stage45_header_t0);
    if (rank == 0 && local_ok) {
        printf("Complete the 4.5a header/layout cost %lf\n",
               costs->stage45_header);
    }
    if (!local_ok) goto fail;


    if (header_memory) free(header_memory);
    return MdAllRanksOk(local_ok) ? 0 : -1;

fail:
    if (header_memory) free(header_memory);
    if (output->header_data) {
        free(output->header_data);
        output->header_data = nullptr;
        output->header_size = 0;
    }
    return -1;
}

static int MdWriteOutputMpiio(
        const swbam::RankBodySink &local_body,
        const MdOutputLayout &layout,
        const std::string &output_path,
        int rank) {

    unsigned long long sizes[2] = {0, 0};
    if (rank == 0) {
        sizes[0] = (unsigned long long)layout.header_size;
        sizes[1] = (unsigned long long)layout.size;
    }
    MPI_Bcast(sizes, 2, MPI_UNSIGNED_LONG_LONG, 0,
              MPI_COMM_WORLD);

    int local_ok = 1;
    swbam::mpi::MpiFileOutput output;
    if (output.Open(output_path, sizes[1]) != 0) local_ok = 0;
    if (MdAllRanksOk(local_ok)) {
        if (rank == 0 && layout.header_size > 0 &&
            output.WriteAt(0, layout.header_data,
                           layout.header_size) != 0) {
            local_ok = 0;
        }
        const uint64_t chunk_size = 8ull * 1024ull * 1024ull;
        std::vector<unsigned char> buffer;
        if (local_body.size() > 0) {
            try {
                buffer.resize((size_t)std::min<uint64_t>(
                    local_body.size(), chunk_size));
            } catch (...) {
                local_ok = 0;
            }
        }
        uint64_t copied = 0;
        while (local_ok && copied < local_body.size()) {
            const size_t count = (size_t)std::min<uint64_t>(
                local_body.size() - copied, buffer.size());
            if (local_body.ReadAt(
                    copied, buffer.data(), count) != 0 ||
                output.WriteAt(
                    sizes[0] + (unsigned long long)
                        layout.body_offsets[(size_t)rank] + copied,
                    buffer.data(), count) != 0) {
                local_ok = 0;
            }
            copied += count;
            if (local_ok && local_body.ReleaseBefore(copied) != 0) local_ok = 0;
        }
        if (rank == 0 &&
            output.WriteAt(
                sizes[0] +
                    (unsigned long long)layout.total_body_size,
                kMarkdupBgzfEofBlock,
                sizeof(kMarkdupBgzfEofBlock)) != 0) {
            local_ok = 0;
        }
    }
    const int writes_ok = MdAllRanksOk(local_ok);
    if (writes_ok && output.Sync() != 0) local_ok = 0;
    if (output.Close() != 0) local_ok = 0;
    return MdAllRanksOk(local_ok) ? 0 : -1;
}

static double MdAccountedOptimizedTime(const MpiMarkdupStats &stats) {
    return stats.t_candidate_decomp +
           stats.t_candidate_extract +
           stats.t_extract_status +
           stats.t_extract_merge +
           stats.t_extract_memcheck +
           stats.t_candidate_exchange +
           stats.t_boundary_check +
           stats.t_candidate_cleanup +
           stats.t_group_prepare +
           stats.t_group +
           stats.t_owner_cleanup +
           stats.t_result_exchange +
           stats.t_result_cleanup +
           stats.t_stream_buffer +
           stats.t_stream_finalize +
           stats.t_rewrite_decomp +
           stats.t_rewrite +
           stats.t_pack +
           stats.t_compress +
           stats.t_rank_sync;
}

static void MdPrintStats(const MpiMarkdupStats &stats,
                         int rank, int comm_size) {
    const int count = 16;
    long long local_counts[count] = {
        stats.input_blocks,
        stats.group_count,
        stats.total_records,
        stats.examined_records,
        stats.excluded_records,
        stats.pair_candidates,
        stats.single_candidates,
        stats.owner_candidates,
        stats.pair_duplicates,
        stats.single_duplicates,
        stats.marked_records,
        stats.cleared_records,
        stats.removed_records,
        stats.qname_bytes,
        stats.mpi_candidate_bytes,
        stats.mpi_result_bytes
    };
    const double accounted = MdAccountedOptimizedTime(stats);
    const double unaccounted = stats.t_optimized_total - accounted;

    const int kLineBytes = 8192;
    char local_lines[kLineBytes];
    memset(local_lines, 0, sizeof(local_lines));
    int used = snprintf(
        local_lines, sizeof(local_lines),
        "[rank %d] markdup_stats blocks=%lld groups=%lld "
        "records=%lld examined=%lld excluded=%lld "
        "pair_candidates=%lld single_candidates=%lld "
        "owner_candidates=%lld qname_bytes=%lld "
        "pair_duplicates=%lld single_duplicates=%lld "
        "marked=%lld cleared=%lld removed=%lld "
        "mpi_candidate_bytes=%lld result_bytes=%lld "
        "bgzf=%lld tracked_peak=%lld stream_pending_peak=%lld\n",
        rank, stats.input_blocks, stats.group_count,
        stats.total_records, stats.examined_records,
        stats.excluded_records, stats.pair_candidates,
        stats.single_candidates, stats.owner_candidates,
        stats.qname_bytes, stats.pair_duplicates,
        stats.single_duplicates, stats.marked_records,
        stats.cleared_records, stats.removed_records,
        stats.mpi_candidate_bytes, stats.mpi_result_bytes,
        stats.bgzf_blocks, stats.tracked_peak_bytes,
        stats.stream_pending_peak_bytes);
    if (used < 0) used = 0;
    if (used >= kLineBytes) used = kLineBytes - 1;
    snprintf(
        local_lines + used, (size_t)(kLineBytes - used),
        "[rank %d] markdup_timing candidate_decomp=%.6f "
        "candidate_extract=%.6f extract_status=%.6f "
        "extract_merge=%.6f extract_memcheck=%.6f "
        "candidate_pack=%.6f candidate_mpi=%.6f "
        "candidate_exchange=%.6f boundary=%.6f "
        "group_sort=%.6f group_scan=%.6f group=%.6f "
        "flat_init=%.6f flat_probe=%.6f result=%.6f "
        "stream_buffer=%.6f stream_finalize=%.6f "
        "rewrite_decomp=%.6f "
        "rewrite=%.6f pack=%.6f compress=%.6f "
        "read=%.6f write=%.6f rank_sync=%.6f "
        "cpe_sync=%.6f accounted=%.6f "
        "unaccounted=%.6f optimized=%.6f\n",
        rank, stats.t_candidate_decomp,
        stats.t_candidate_extract, stats.t_extract_status,
        stats.t_extract_merge, stats.t_extract_memcheck,
        stats.t_candidate_pack,
        stats.t_candidate_exchange_mpi,
        stats.t_candidate_exchange,
        stats.t_boundary_check, stats.t_group_sort,
        stats.t_group_scan, stats.t_group,
        stats.t_flat_init, stats.t_flat_probe,
        stats.t_result_exchange,
        stats.t_stream_buffer, stats.t_stream_finalize,
        stats.t_rewrite_decomp, stats.t_rewrite, stats.t_pack,
        stats.t_compress, stats.t_read, stats.t_write,
        stats.t_rank_sync, stats.t_cpe_sync, accounted,
        unaccounted, stats.t_optimized_total);
    std::vector<char> gathered_lines;
    if (rank == 0) {
        gathered_lines.resize((size_t)comm_size * kLineBytes);
    }
    MPI_Gather(local_lines, kLineBytes, MPI_CHAR,
               rank == 0 ? gathered_lines.data() : nullptr,
               kLineBytes, MPI_CHAR, 0, MPI_COMM_WORLD);

    long long sums[count] = {};
    MPI_Reduce(local_counts, sums, count, MPI_LONG_LONG,
               MPI_SUM, 0, MPI_COMM_WORLD);
    long long peak = stats.tracked_peak_bytes;
    long long peak_max = 0;
    MPI_Reduce(&peak, &peak_max, 1, MPI_LONG_LONG,
               MPI_MAX, 0, MPI_COMM_WORLD);
    long long pending_peak = stats.stream_pending_peak_bytes;
    long long pending_peak_max = 0;
    MPI_Reduce(&pending_peak, &pending_peak_max, 1, MPI_LONG_LONG,
               MPI_MAX, 0, MPI_COMM_WORLD);

    const int time_count = 28;
    double local_times[time_count] = {
        stats.t_candidate_decomp,
        stats.t_candidate_extract,
        stats.t_extract_status,
        stats.t_extract_merge,
        stats.t_extract_memcheck,
        stats.t_candidate_pack,
        stats.t_candidate_exchange_mpi,
        stats.t_candidate_exchange,
        stats.t_boundary_check,
        stats.t_group_sort,
        stats.t_group_scan,
        stats.t_group,
        stats.t_flat_init,
        stats.t_flat_probe,
        stats.t_result_exchange,
        stats.t_stream_buffer,
        stats.t_stream_finalize,
        stats.t_rewrite_decomp,
        stats.t_rewrite,
        stats.t_pack,
        stats.t_compress,
        stats.t_read,
        stats.t_write,
        stats.t_rank_sync,
        stats.t_cpe_sync,
        accounted,
        unaccounted,
        stats.t_optimized_total
    };
    double time_sums[time_count] = {};
    double time_max[time_count] = {};
    MPI_Reduce(local_times, time_sums, time_count, MPI_DOUBLE,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_times, time_max, time_count, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int r = 0; r < comm_size; ++r) {
            fputs(gathered_lines.data() + (size_t)r * kLineBytes,
                  stdout);
        }
        printf("OptimizedMarkdupMPI finished. ranks=%d blocks=%lld "
               "records=%lld examined=%lld excluded=%lld\n",
               comm_size, sums[0], sums[2], sums[3], sums[4]);
        printf("  candidates pair=%lld single=%lld owner=%lld "
               "qname_bytes=%lld\n",
               sums[5], sums[6], sums[7], sums[13]);
        printf("  duplicates pair=%lld single=%lld marked=%lld "
               "cleared=%lld removed=%lld\n",
               sums[8], sums[9], sums[10], sums[11], sums[12]);
        printf("  mpi candidate_bytes=%lld result_bytes=%lld "
               "tracked_peak_max=%lld stream_pending_peak_max=%lld\n",
               sums[14], sums[15], peak_max, pending_peak_max);
        printf("  timing_sum candidate_decomp=%.3f "
               "candidate_extract=%.3f extract_status=%.3f "
               "extract_merge=%.3f extract_memcheck=%.3f "
               "candidate_pack=%.3f candidate_mpi=%.3f "
               "candidate_exchange=%.3f boundary=%.3f "
               "group_sort=%.3f group_scan=%.3f group=%.3f "
               "flat_init=%.3f flat_probe=%.3f result=%.3f "
               "stream_buffer=%.3f stream_finalize=%.3f "
               "rewrite_decomp=%.3f "
               "rewrite=%.3f pack=%.3f compress=%.3f "
               "read=%.3f write=%.3f rank_sync=%.3f "
               "cpe_sync=%.3f accounted=%.3f "
               "unaccounted=%.3f optimized=%.3f\n",
               time_sums[0], time_sums[1], time_sums[2],
               time_sums[3], time_sums[4], time_sums[5],
               time_sums[6], time_sums[7], time_sums[8],
               time_sums[9], time_sums[10], time_sums[11],
               time_sums[12], time_sums[13], time_sums[14],
               time_sums[15], time_sums[16], time_sums[17],
               time_sums[18], time_sums[19], time_sums[20],
               time_sums[21], time_sums[22], time_sums[23],
               time_sums[24], time_sums[25], time_sums[26],
               time_sums[27]);
        printf("  timing_max candidate_decomp=%.3f "
               "candidate_extract=%.3f extract_status=%.3f "
               "extract_merge=%.3f extract_memcheck=%.3f "
               "candidate_pack=%.3f candidate_mpi=%.3f "
               "candidate_exchange=%.3f boundary=%.3f "
               "group_sort=%.3f group_scan=%.3f group=%.3f "
               "flat_init=%.3f flat_probe=%.3f result=%.3f "
               "stream_buffer=%.3f stream_finalize=%.3f "
               "rewrite_decomp=%.3f "
               "rewrite=%.3f pack=%.3f compress=%.3f "
               "read=%.3f write=%.3f rank_sync=%.3f "
               "cpe_sync=%.3f accounted=%.3f "
               "unaccounted=%.3f optimized=%.3f\n",
               time_max[0], time_max[1], time_max[2],
               time_max[3], time_max[4], time_max[5],
               time_max[6], time_max[7], time_max[8],
               time_max[9], time_max[10], time_max[11],
               time_max[12], time_max[13], time_max[14],
               time_max[15], time_max[16], time_max[17],
               time_max[18], time_max[19], time_max[20],
               time_max[21], time_max[22], time_max[23],
               time_max[24], time_max[25], time_max[26],
               time_max[27]);
        fflush(stdout);
    }
}

} // namespace


int ProcessMarkdupMPI(CmdInfo *cmd_info) {

    // 1. 初始化阶段：拿 rank、解析内存参数
    double total_t0 = GetTime();
    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int exit_code = 1;
    int local_ok = 1;
    swbam::mpi::MpiBamInput input_handle;
    swbam::BamInputBackend *input_backend = nullptr;
    swbam::mpi::MpiBamInputPlan input_plan;
    sam_hdr_t *header = nullptr;
    long long body_start = 0;
    long long n_blocks = 0;
    long long local_block_begin = 0;
    size_t rank_input_size = 0;
    swbam_usage::BodySink rank_body_sink;
    swbam_usage::OrderedBody ordered_body(
        &rank_body_sink, swbam_usage::CacheMemoryLimit() / 4);
    swbam::RankBodySink *body_sink = &rank_body_sink;
    uint64_t rank_body_memory_limit = 0;
    size_t memory_limit = 0;
    MpiMarkdupStats stats = {};
    MdOutputLayout output_layout = {};
    MdOutputStageCosts output_stage_costs = {};
    double stage41_cost = 0.0;
    double stage42_cost = 0.0;
    double stage43_cost = 0.0;
    double body_total_t0 = 0.0;
    body_sink = &ordered_body;

    uint64_t local_records = 0;
    int range_has_records = 0;
    int range_first_tid = -1;
    int range_first_pos = -1;
    int range_last_tid = -1;
    int range_last_pos = -1;

    if (MdParseMemory(cmd_info->markdup_memory_,
                      &memory_limit) != 0) {
        if (rank == 0) {
            fprintf(stderr,
                    "ERROR: invalid markdup memory limit '%s'.\n",
                    cmd_info->markdup_memory_.c_str());
        }
        local_ok = 0;
    }
    if (swbam::ParseByteSize(
            cmd_info->rank_body_memory_limit_,
            &rank_body_memory_limit) != 0 ||
        rank_body_memory_limit == 0) {
        if (rank == 0) {
            fprintf(stderr,
                    "ERROR: invalid rank body memory limit '%s'.\n",
                    cmd_info->rank_body_memory_limit_.c_str());
        }
        local_ok = 0;
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;
    if (rank == 0) {
        printf("111Complete the initialization cost %lf-----\n",
               MdReduceMax(GetTime() - total_t0));
    } else {
        MdReduceMax(GetTime() - total_t0);
    }


    // 2. 通过公共 backend 打开 BAM；按批次从磁盘读取。
    {
        if (swbam_usage::OpenInput(input_handle,
                cmd_info->in_file_name_, cmd_info->io_backend_) != 0) {
            fprintf(stderr,
                    "[rank %d] ERROR: cannot open markdup input %s.\n",
                    rank, cmd_info->in_file_name_.c_str());
            local_ok = 0;
        }
        input_backend = input_handle.backend();
        if (local_ok && (!input_backend ||
                         input_backend->format() != bam)) {
            local_ok = 0;
        }

        double cost = MdReduceMax(input_handle.data_open_cost());
        if (rank == 0 && local_ok) {
            printf("MPI BAM input backend=%s\n", input_handle.selected_backend().c_str());
            printf("222Complete the input open cost %lf--\n", cost);
        }
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;


    // 3. 读取 backend header 并检查排序方式。
    {
        double t0 = GetTime();
        if (local_ok && input_backend->header()) {
            header = sam_hdr_dup(input_backend->header());
        }
        if (!header) local_ok = 0;
        if (local_ok) {
            kstring_t sort_order = {0, 0, nullptr};
            if (sam_hdr_find_tag_hd(header, "SO",
                                    &sort_order) == 0 &&
                sort_order.s &&
                strcmp(sort_order.s, "queryname") == 0) {
                if (rank == 0) {
                    fprintf(stderr,
                            "ERROR: markdup input is queryname sorted; "
                            "coordinate sort it first.\n");
                }
                local_ok = 0;
            }
            free(sort_order.s);
        }
        if (local_ok) {
            body_start = (long long)input_backend->body_offset();
            if (body_start < 0 || (unsigned long long)body_start >
                    (unsigned long long)input_backend->size()) {
                local_ok = 0;
            }
        }
        double cost = MdReduceMax(
            input_handle.header_open_cost() + GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("333Complete the head cost %lf---\n", cost);
        }
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;

    body_total_t0 = GetTime();
    {
        // 4.1：扫描 BGZF block，并按 block 分给各 rank
        double t0 = GetTime();
        if (rank == 0) {
            printf("Enable MPI BAM MARKDUP mode (%d MPE + %d CPEs)!!!\n",
                   comm_size, comm_size * 64);
            printf("MPI BAM output compression level=%d\n",
                   cmd_info->compress_level_);
        }
        if (input_handle.PreparePlan(&input_plan) != 0) {
            local_ok = 0;
        }
        n_blocks = (long long)input_plan.blocks.size();
        local_block_begin = (long long)input_plan.rank_begin;
        rank_input_size = 0;
        for (size_t i = input_plan.rank_begin;
             local_ok && i < input_plan.rank_end; ++i) {
            const size_t block_size =
                input_plan.blocks[i].compressed_size;
            if (block_size > SIZE_MAX - rank_input_size) {
                local_ok = 0;
            } else {
                rank_input_size += block_size;
            }
        }

        stage41_cost = MdReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("MPI BAM scan complete. data_blocks=%lld "
                   "body_start=%lld\n",
                   n_blocks, body_start);
            printf("Complete the 4.1 scan/split/broadcast/select cost %lf\n",
                   stage41_cost);
        }
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;

    // 4.2: initialize rank-local output staging.
    {
        double t0 = GetTime();
        char prefix[96];
        snprintf(prefix, sizeof(prefix),
                 "swbam-markdup-rank-%d.body", rank);
        if (rank_body_sink.Open(
                cmd_info->rank_body_backend_,
                rank_body_memory_limit,
                rank_input_size ? rank_input_size :
                                  64 * 1024 * 1024,
                prefix, cmd_info->rank_body_temp_dir_) != 0) {
            local_ok = 0;
        }
        stage42_cost = MdReduceMax(GetTime() - t0);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.2 init reader/writer cost %lf\n",
                   stage42_cost);
        }
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;

    // 4.3 核心：提取候选、交换、找 duplicate、重写 BAM
    {
        if (rank == 0) {
            printf("MPI BAM markdup algorithm=%s max_read_length=%d\n",
                   "coordinate-stream-one-pass",
                   cmd_info->markdup_max_read_length_);
        }
        double optimized_t0 = GetTime();

        MdOnePassBoundaryPlan boundary_plan;
        if (MdBuildOnePassBoundaryPlan(
                *input_backend, input_plan.rank_spans(),
                input_plan.rank_block_count(), local_block_begin,
                rank, comm_size,
                cmd_info->markdup_include_fails_,
                cmd_info->markdup_max_read_length_, memory_limit,
                &boundary_plan, &stats) != 0) {
            local_ok = 0;
        }
        if (!MdAllRanksOkTimed(local_ok, &stats)) goto cleanup;
        const int stream_ret = MdMarkdupOnePassToSink(
                *input_backend, input_plan.rank_spans(),
                input_plan.rank_block_count(), local_block_begin,
                rank, comm_size,
                cmd_info->markdup_include_fails_,
                cmd_info->markdup_max_read_length_, memory_limit,
                cmd_info->compress_level_, body_sink,
                &boundary_plan,
                &local_records, &range_has_records,
                &range_first_tid, &range_first_pos,
                &range_last_tid, &range_last_pos, &stats);
        if (stream_ret != 0) local_ok = 0;
        if (!MdAllRanksOkTimed(local_ok, &stats)) goto cleanup;

        stats.t_optimized_total = GetTime() - optimized_t0;
        stage43_cost = MdReduceMax(stats.t_optimized_total);
        if (rank == 0) {
            printf("Complete the 4.3 OptimizedMarkdupMPI cost %lf\n",
                   stage43_cost);
        }
    }

    {
        const int local_spool = body_sink->is_memory() ? 0 : 1;
        int spool_ranks = 0;
        MPI_Reduce(&local_spool, &spool_ranks, 1, MPI_INT,
                   MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            printf("MPI rank body backend requested=%s "
                   "memory_ranks=%d spool_ranks=%d\n",
                   cmd_info->rank_body_backend_.c_str(),
                   comm_size - spool_ranks, spool_ranks);
        }
    }

    // Build the header and rank-local file offsets before distributed output.
    if (!body_sink || body_sink->Flush() != 0) local_ok = 0;
    if (!MdAllRanksOk(local_ok)) goto cleanup;
    if (MdPrepareOutputLayout(*body_sink, header,
                              cmd_info->compress_level_,
                              rank, comm_size, &stats,
                              &output_layout,
                              &output_stage_costs) != 0) {
        local_ok = 0;
    }
    if (!MdAllRanksOk(local_ok)) goto cleanup;

    if (rank == 0) {
        double core_total = stage41_cost + stage42_cost +
                            stage43_cost + output_stage_costs.stage44 +
                            output_stage_costs.stage45_header +
                            output_stage_costs.stage46;
        printf("Complete the total (4.1~4.6) cost %lf-----\n",
               core_total);
    }

    // 4.7 MdPrintStats打印统计信息
    {
        double stage47_t0 = GetTime();
        MdPrintStats(stats, rank, comm_size);
        double stage47_cost =
            MdReduceMax(GetTime() - stage47_t0);
        if (rank == 0 && local_ok) {
            printf("Complete the 4.7 rank stats reduce/print cost %lf\n",
                   stage47_cost);
        }
    }

    {
        double body_total_cost =
            MdReduceMax(GetTime() - body_total_t0);
        if (rank == 0 && local_ok) {
            printf("444Complete the total body cost %lf\n",
                   body_total_cost);
        }
    }


    {
        double write_t0 = GetTime();
        if (MdWriteOutputMpiio(
                *body_sink, output_layout, cmd_info->out_file_name_,
                rank) != 0) {
            local_ok = 0;
        }
        const double write_cost =
            MdReduceMax(GetTime() - write_t0);
        if (rank == 0 && local_ok) {
            printf("555MPI-IO distributed output cost %lf--\n",
                   write_cost);
        }
        if (!MdAllRanksOk(local_ok)) goto cleanup;
        exit_code = 0;
        goto cleanup;
    }

cleanup:
    if (output_layout.header_data) free(output_layout.header_data);
    if (header) sam_hdr_destroy(header);
    if (rank_body_sink.Close() != 0) exit_code = 1;
    input_handle.Close();
    double close_cost = MdReduceMax(GetTime() - total_t0);
    if (rank == 0) {
        printf("666markdup total process cost %lf-----\n",
               close_cost);
    }
    return exit_code;
}
