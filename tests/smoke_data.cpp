#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <mpi.h>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
const size_t kRecords = 20000;

void Require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

void MakeFixture(const std::string &path) {
    const char *header_text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000000\n";
    sam_hdr_t *header = sam_hdr_parse(strlen(header_text), header_text);
    samFile *out = sam_open(path.c_str(), "wb1");
    bam1_t *record = bam_init1();
    Require(header && out && record, "fixture allocation/open");
    Require(sam_hdr_write(out, header) == 0 && bgzf_flush(out->fp.bgzf) == 0, "fixture header");
    size_t records = 0;
    for (int group = 0; group < 5000; ++group) {
        for (int end = 0; end < 2; ++end) {
            for (int copy = 0; copy < 2; ++copy) {
                const int position = 100 + 300 * group;
                std::ostringstream line;
                line << "pair" << group << (copy ? ".b" : ".a") << '\t'
                     << (end ? 147 : 99) << "\tchr1\t" << position + 100 * end
                     << '\t' << (copy ? 20 : 60) << "\t50M\t=\t"
                     << position + 100 * (1-end) << '\t' << (end ? -150 : 150)
                     << '\t' << std::string(50, 'A') << '\t'
                     << std::string(50, copy ? '#' : 'I') << "\tNM:i:0";
                std::string text = line.str();
                kstring_t sam = {text.size(), text.size()+1, &text[0]};
                Require(sam_parse1(&sam, header, record) >= 0 &&
                        sam_write1(out, header, record) >= 0, "fixture record");
                // Exercise rank edges without allowing a record to cross a BGZF block.
                if (++records % 31 == 0) Require(bgzf_flush(out->fp.bgzf) == 0, "fixture flush");
            }
        }
    }
    bam_destroy1(record); sam_hdr_destroy(header);
    Require(sam_close(out) == 0, "fixture close");
}

std::vector<std::string> Read(const std::string &path, bool sorted = false,
                              bool collated = false, bool marked = false,
                              bool mates = false, int min_mapq = -1) {
    samFile *in = sam_open(path.c_str(), "r");
    Require(in != nullptr, "cannot open result");
    sam_hdr_t *header = sam_hdr_read(in);
    bam1_t *b = bam_init1();
    Require(header && b, "result header");
    kstring_t text = {0, 0, nullptr};
    std::vector<std::string> records;
    std::set<std::string> completed;
    std::string previous_name;
    int previous_tid = -1;
    hts_pos_t previous_pos = -1;
    int status = 0;
    while ((status = sam_read1(in, header, b)) >= 0) {
        const std::string name = bam_get_qname(b);
        if (sorted) {
            Require(b->core.tid > previous_tid ||
                    (b->core.tid == previous_tid && b->core.pos >= previous_pos), "coordinate order");
            previous_tid = b->core.tid; previous_pos = b->core.pos;
        }
        if (collated && name != previous_name) {
            if (!previous_name.empty()) completed.insert(previous_name);
            Require(!completed.count(name), "QNAME continuity");
            previous_name = name;
        }
        if (mates) {
            Require(bam_aux_get(b, "MC") && bam_aux_get(b, "MQ") && bam_aux_get(b, "ms"), "missing mate tags");
            Require(std::string(bam_aux2Z(bam_aux_get(b, "MC"))) == "50M", "mate CIGAR");
            Require(bam_aux2i(bam_aux_get(b, "MQ")) == (name.back() == 'b' ? 20 : 60), "mate MAPQ");
            Require(bam_aux2i(bam_aux_get(b, "ms")) == (name.back() == 'b' ? 0 : 2000), "mate score");
        }
        if (marked) {
            Require(bool(b->core.flag & BAM_FDUP) == (name.back() == 'b'), "duplicate membership");
            b->core.flag &= ~BAM_FDUP;
        }
        if (min_mapq >= 0 && b->core.qual < min_mapq) continue;
        Require(sam_format1(header, b, &text) >= 0, "SAM formatting");
        records.push_back(std::string(text.s, text.l));
    }
    Require(status == -1, "truncated or malformed output");
    free(text.s); bam_destroy1(b); sam_hdr_destroy(header);
    Require(sam_close(in) == 0, "result close");
    std::sort(records.begin(), records.end());
    return records;
}

void Verify(const std::string &directory) {
    const auto path = [&](const std::string &name) { return directory + '/' + name; };
    const auto original = Read(path("input.bam"), true);
    Require(original.size() == kRecords, "fixture count");
    Require(Read(path("rewrite.bam")) == original, "rewrite fields/AUX");
    Require(Read(path("from-sam.bam")) == original, "SAM roundtrip fields/AUX");
    Require(Read(path("sdk-rewrite.bam")) == original, "SDK rewrite fields/AUX");
    Require(Read(path("filtered.bam")) == Read(path("input.bam"), false, false, false, false, 30), "filter membership");
    Require(Read(path("names.bam"), false, true) == original, "collate multiset");
    Require(Read(path("names-external.bam"), false, true) == original, "external collate multiset");
    const auto mates = Read(path("mates.bam"), false, true, false, true);
    Require(mates.size() == kRecords, "fixmate count");
    Require(Read(path("sorted.bam"), true) == mates, "sort multiset");
    Require(Read(path("sorted-external.bam"), true) == mates, "external sort multiset");
    Require(Read(path("marked.bam"), true, false, true) == mates, "markdup fields/AUX");
    Require(Read(path("workflow.bam"), true, false, true) == mates, "workflow fields/AUX");
    puts("SWBAM_SMOKE PASS records=20000 duplicates=10000");
}
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int status = 0;
    if (rank == 0) {
        try {
            Require(argc == 3, "Usage: swbam-smoke-data make|verify DIRECTORY");
            if (std::string(argv[1]) == "make") MakeFixture(std::string(argv[2]) + "/input.bam");
            else if (std::string(argv[1]) == "verify") Verify(argv[2]);
            else throw std::runtime_error("unknown smoke-data mode");
        } catch (const std::exception &e) {
            fprintf(stderr, "SWBAM_SMOKE FAIL: %s\n", e.what()); status = 1;
        }
    }
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    return status;
}
