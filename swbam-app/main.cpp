#include "CLI11.hpp"
#include "swbam_mpi.h"
#include <athread.h>
#include <cstdio>
#include <mpi.h>

static int RunCommand(int argc, char **argv, int rank) {
    CmdInfo cmd;
    CLI::App app("Distributed BAM/SAM processing on Sunway", "SWBAM");
    app.require_subcommand(1, 1);
    const auto input = [&](CLI::App *sub) {
        sub->add_option("-i,--inFile", cmd.in_file_name_, "Input BAM/SAM file")
            ->required()->check(CLI::ExistingFile);
        sub->add_option("--io-backend", cmd.io_backend_, "Input backend: posix or mpiio")
            ->check(CLI::IsMember(std::vector<std::string>{"posix", "mpiio"}));
        sub->add_flag("--verbose", cmd.verbose_, "Print detailed diagnostics");
    };
    const auto output = [&](CLI::App *sub) {
        sub->add_option("-o,--outFile", cmd.out_file_name_, "Output BAM/SAM file")->required();
        sub->add_option("--compress-level", cmd.compress_level_, "BAM compression level (paper implementation: 1)")
            ->check(CLI::IsMember(std::vector<int>{1}));
        sub->add_option("--rank-body-backend", cmd.rank_body_backend_, "Rank-local staging: memory, spool, or auto")
            ->check(CLI::IsMember(std::vector<std::string>{"memory", "spool", "auto"}));
        sub->add_option("--rank-body-memory-limit", cmd.rank_body_memory_limit_, "Per-rank staging memory budget");
        sub->add_option("--rank-body-temp-dir", cmd.rank_body_temp_dir_, "Directory for rank-local staging files");
    };
    CLI::App *view = app.add_subcommand("view", "Filter, convert, or rewrite BAM/SAM files");
    input(view); output(view);
    view->add_option("--min-mapq", cmd.min_mapq_, "Keep records with MAPQ at least this value");
    view->add_option("--max-mapq", cmd.max_mapq_, "Keep records with MAPQ at most this value");
    view->add_option("--require-flag", cmd.require_flag_, "Require all specified FLAG bits");
    view->add_option("--exclude-flag", cmd.exclude_flag_, "Exclude any specified FLAG bit");
    view->add_option("--ref-name", cmd.ref_name_, "Keep the specified reference sequence");
    view->add_option("--min-read-len", cmd.min_read_len_, "Minimum read length");
    view->add_option("--max-read-len", cmd.max_read_len_, "Maximum read length");
    CLI::App *flagstat = app.add_subcommand("flagstat", "Alignment flag statistics");
    input(flagstat);
    CLI::App *sort = app.add_subcommand("sort", "Distributed coordinate sorting");
    input(sort); output(sort);
    sort->add_option("-m,--memory", cmd.sort_memory_, "Algorithm workspace per rank, e.g. 8G");
    sort->add_option("-T,--temp-prefix", cmd.sort_temp_prefix_, "External-run file prefix");
    const auto sort_splitters = [&](CLI::App *command) {
        auto *exact = command->add_flag("--exact-splitters", cmd.sort_exact_splitters_,
                                         "Use exact global splitters (default)");
        auto *sample = command->add_flag_callback("--sample-splitters",
            [&]() { cmd.sort_exact_splitters_ = false; }, "Use sampled splitters instead");
        exact->excludes(sample);
    };
    sort_splitters(sort);
    CLI::App *collate = app.add_subcommand("collate", "Distributed QNAME collation");
    input(collate); output(collate);
    collate->add_option("-m,--memory", cmd.collate_memory_, "Algorithm workspace per rank, e.g. 8G");
    collate->add_option("-T,--temp-prefix", cmd.collate_temp_prefix_, "External-run file prefix");
    auto *collate_bins = collate->add_option("-n,--bins", cmd.collate_bins_, "Logical QNAME hash bins");
    CLI::App *fixmate = app.add_subcommand("fixmate", "Update mate fields and tags in name-collated BAM");
    input(fixmate); output(fixmate);
    fixmate->add_flag("-m", cmd.fixmate_mate_score_, "Calculate mate-score tags")->required();
    CLI::App *markdup = app.add_subcommand("markdup", "One-pass coordinate-stream duplicate marking");
    input(markdup); output(markdup);
    markdup->add_flag("--stream", "Single-pass streaming (always enabled)");
    markdup->add_option("-l,--max-read-length", cmd.markdup_max_read_length_, "Read-length bound for window retirement")
        ->check(CLI::PositiveNumber);
    markdup->add_option("-m,--memory", cmd.markdup_memory_, "Algorithm workspace per rank, e.g. 8G");
    CLI::App *workflow = app.add_subcommand("dedup-workflow", "Collate, fixmate, sort, and markdup using intermediate BAM files");
    input(workflow); output(workflow);
    workflow->add_option("-T,--temp-prefix", cmd.workflow_temp_prefix_, "Unique prefix for intermediate files")->required();
    workflow->add_option("-m,--memory", cmd.pipeline_memory_, "Algorithm workspace per rank")->required();
    sort_splitters(workflow);
    auto *workflow_bins = workflow->add_option("-n,--bins", cmd.collate_bins_, "Logical QNAME hash bins");
    workflow->add_option("-l,--max-read-length", cmd.markdup_max_read_length_, "Read-length bound for duplicate marking")
        ->check(CLI::PositiveNumber);
    workflow->add_flag("--keep-intermediates", cmd.workflow_keep_intermediates_, "Keep intermediate BAM files");

    int parse_exit = 0, parse_status = 0;
    try { app.parse(argc, argv); }
    catch (const CLI::ParseError &error) {
        parse_exit = 1;
        parse_status = error.get_exit_code();
        if (rank == 0) app.exit(error);
    }
    int any_exit = 0, result = 0;
    MPI_Allreduce(&parse_exit, &any_exit, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (any_exit) {
        MPI_Allreduce(&parse_status, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        return result;
    }
    cmd.collate_bins_explicit_ = collate_bins->count() || workflow_bins->count();
    if (!*workflow) {
        const size_t work = swbam_usage::WorkMemoryLimit();
        const size_t reserve = std::max<size_t>(256u << 20, work / 8);
        std::string *limit = *sort ? &cmd.sort_memory_ :
            *collate ? &cmd.collate_memory_ : *markdup ? &cmd.markdup_memory_ : nullptr;
        uint64_t requested = 0;
        int ok = work > reserve;
        if (limit) {
            if (limit->empty() && ok) *limit = std::to_string(work - reserve);
            ok = ok && swbam::ParseByteSize(*limit, &requested) == 0 && requested <= work - reserve;
        }
        if (!swbam::mpi::AllRanksOk(ok)) {
            if (rank == 0) fprintf(stderr, "ERROR: algorithm workspace exceeds --work-memory or node allowance\n");
            return 2;
        }
    }
    return *workflow ? ProcessDedupWorkflowMPI(&cmd) :
           *sort ? ProcessSortMPI(&cmd) : *collate ? ProcessCollateMPI(&cmd) :
           *fixmate ? ProcessFixmateMPI(&cmd) : *markdup ? ProcessMarkdupMPI(&cmd) :
           *flagstat ? ProcessFlagstatMPI(&cmd) : ProcessSwBamMPI(&cmd);
}

int main(int argc, char **argv) {
    swbam_usage::Session runtime;
    if (runtime.Prepare(&argc, &argv)) return 2;
    const double start = swbam_usage::WallTime();
    MPI_Init(&argc, &argv);
    const double processing_start = swbam_usage::WallTime();
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int status = 2;
    if (swbam::mpi::AllRanksOk(runtime.EnterMpi() == 0)) {
        athread_init();
        status = RunCommand(argc, argv, rank);
        athread_halt();
    }
    const double processing = swbam_usage::WallTime() - processing_start;
    MPI_Finalize();
    return runtime.Finish(status, swbam_usage::WallTime() - start, processing);
}
