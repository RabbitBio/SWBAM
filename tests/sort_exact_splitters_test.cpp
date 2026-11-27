// Host MPI regression: the exact selector and temporary-run reader do not use CPEs.
#include <cstdlib>
#include <random>

#ifndef PLATFORM_SUNWAY
static void __real_athread_spawn(void *, void *, int) { std::abort(); }
static void athread_join() { std::abort(); }
#endif

#include "../swbam-app/algorithms/sort/swbam_sort_optimized_mpi.cpp"

static void Check(bool condition, const char *message) {
    if (condition) return;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    fprintf(stderr, "[rank %d] exact sort splitter test failed: %s\n", rank, message);
    MPI_Abort(MPI_COMM_WORLD, 1);
    std::abort();
}

static void CheckSplitters(const std::vector<MpiSortKey> &splitters,
                           const std::vector<MpiSortRecordMeta> &reference,
                           int ranks) {
    Check(splitters.size() == size_t(ranks - 1), "splitter count");
    if (reference.empty()) return;
    const size_t base = reference.size() / ranks, extra = reference.size() % ranks;
    for (size_t b = 0; b < splitters.size(); ++b) {
        const size_t boundary = (b + 1) * base + std::min(b + 1, extra);
        Check(MpiSortCompareMetaKey(reference[boundary - 1], splitters[b]) == 0,
              "splitter differs from full-sort oracle");
    }
}

static void TestInputs(const std::vector<std::vector<MpiSortRecordMeta>> &local_runs,
                       const std::vector<MpiSortRecordMeta> &reference,
                       int rank, int ranks) {
    std::vector<MpiSortKey> splitters;
    std::vector<MpiSortRecordMeta> combined;
    for (const auto &run : local_runs)
        combined.insert(combined.end(), run.begin(), run.end());
    std::sort(combined.begin(), combined.end(), MpiSortMetaLess());
    const auto memory_size = [&](size_t) { return combined.size(); };
    const auto memory_read = [&](size_t, size_t i, MpiSortKey *key) {
        *key = MpiSortKeyFromMeta(combined[i]);
        return true;
    };
    Check(swbam::sort_detail::ChooseExactSplitters(
              1, memory_size, memory_read, MpiSortKeyLess(), rank, ranks,
              MPI_COMM_WORLD, SIZE_MAX, &splitters) == 0, "in-memory selection");
    CheckSplitters(splitters, reference, ranks);

    FILE *file = tmpfile();
    Check(file != nullptr, "open temp file");
    MpiSortStats stats = {};
    MpiSortStrictRunArena arena;
    std::vector<MpiSortStrictRun> runs;
    off_t end = 0;
    for (size_t r = 0; r < local_runs.size(); ++r) {
        const auto &records = local_runs[r];
        MpiSortStrictRun run;
        if (r == 0 && (rank & 1) == 0 && !records.empty()) {
            run.memory_backed = true;
            run.extent.record_count = records.size();
            arena.data = reinterpret_cast<unsigned char *>(
                const_cast<MpiSortRecordMeta *>(records.data()));
        } else {
            run.extent.fd = fileno(file);
            run.extent.record_count = records.size();
            run.extent.meta_offset = end;
            const size_t bytes = records.size() * sizeof(MpiSortRecordMeta);
            Check(MpiSortStrictPwriteAll(run.extent.fd, records.data(), bytes, end, &stats) == 0,
                  "write external metadata");
            end += bytes;
        }
        runs.push_back(run);
    }
    Check(MpiSortChooseExactRunSplitters(runs, arena, rank, ranks, 64u << 10,
                                         &stats, &splitters) == 0, "external selection");
    CheckSplitters(splitters, reference, ranks);
    if (ranks > 1 && !reference.empty()) {
        Check(swbam::sort_detail::ChooseExactSplitters(
                  1, memory_size, memory_read, MpiSortKeyLess(), rank, ranks,
                  MPI_COMM_WORLD, 0, &splitters) != 0, "workspace cap");
        const auto fail_read = [&](size_t r, size_t i, MpiSortKey *key) {
            return rank != 0 && memory_read(r, i, key);
        };
        Check(swbam::sort_detail::ChooseExactSplitters(
                  1, memory_size, fail_read, MpiSortKeyLess(), rank, ranks,
                  MPI_COMM_WORLD, SIZE_MAX, &splitters) != 0, "collective read failure");
    }
    fclose(file);
}

static void TestLargeVirtualInput(int rank, int ranks) {
    const size_t local = 10000019;
    size_t probes = 0;
    const auto run_size = [&](size_t) { return local; };
    const auto read_key = [&](size_t, size_t i, MpiSortKey *key) {
        ++probes;
        *key = {};
        key->global_order = i * ranks + rank;
        return true;
    };
    std::vector<MpiSortKey> splitters;
    Check(swbam::sort_detail::ChooseExactSplitters(
              1, run_size, read_key, MpiSortKeyLess(), rank, ranks,
              MPI_COMM_WORLD, SIZE_MAX, &splitters) == 0, "large virtual selection");
    for (size_t i = 0; i < splitters.size(); ++i)
        Check(splitters[i].global_order == (i + 1) * local - 1, "large virtual boundary");
    Check(probes < local / 100, "selection reads too many keys");
    const unsigned long long local_probes = probes;
    unsigned long long max_probes = 0;
    MPI_Reduce(&local_probes, &max_probes, 1, MPI_UNSIGNED_LONG_LONG,
               MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
        printf("PASS %zu virtual keys: max %llu probes per rank\n", local * ranks, max_probes);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    std::mt19937 random(15);
    for (int test = 0; test < 36; ++test) {
        std::vector<std::vector<MpiSortRecordMeta>> runs(1 + rank % 4);
        std::vector<MpiSortRecordMeta> reference;
        const unsigned n = test < 2 ? unsigned(test) : 1001 + random() % 6000;
        for (unsigned i = 0; i < n; ++i) {
            MpiSortRecordMeta record = {};
            record.global_order = i;
            record.tid = test % 3 == 0 ? int(random() % 5) - 1 : 0;
            record.pos = test % 3 == 0 ? int(random() % 11) - 1 : 23;
            record.flag = (random() % 2) ? BAM_FREVERSE : 0;
            reference.push_back(record);
            const int owner = i == 0 || (test % 2 == 0 && i % 19 != 0)
                ? 0 : random() % ranks;
            const unsigned run = random() % unsigned(1 + owner % 4);
            if (rank == owner) runs[run].push_back(record);
        }
        std::sort(reference.begin(), reference.end(), MpiSortMetaLess());
        for (auto &run : runs) std::sort(run.begin(), run.end(), MpiSortMetaLess());
        TestInputs(runs, reference, rank, ranks);
    }
    TestLargeVirtualInput(rank, ranks);
    if (rank == 0) printf("PASS exact sort splitters (memory and external, ranks=%d)\n", ranks);
    MPI_Finalize();
}
