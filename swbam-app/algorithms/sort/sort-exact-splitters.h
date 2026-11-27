#ifndef SWBAM_SORT_EXACT_SPLITTERS_H
#define SWBAM_SORT_EXACT_SPLITTERS_H

#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

#include <mpi.h>

namespace swbam {
namespace sort_detail {

using std::size_t;

template <class Key>
struct SplitterPivot {
    Key key;
    unsigned long long weight;
    size_t run;
    size_t index;
    int rank;
};

template <class Pivot, class Less>
Pivot WeightedMedian(std::vector<Pivot> *pivots, Less less) {
    std::sort(pivots->begin(), pivots->end(),
              [&](const Pivot &a, const Pivot &b) { return less(a.key, b.key); });
    unsigned long long total = 0;
    for (const auto &pivot : *pivots) total += pivot.weight;
    unsigned long long left = total / 2;
    for (auto pivot : *pivots) {
        if (left < pivot.weight) {
            pivot.weight = total;
            return pivot;
        }
        left -= pivot.weight;
    }
    return Pivot{};
}

// Each run must be sorted, and complete keys must be globally unique.
// Return inclusive upper bounds; the first N % ranks buckets get one extra record.
template <class Key, class Less, class RunSize, class ReadKey>
int ChooseExactSplitters(size_t run_count, RunSize run_size, ReadKey read_key,
                         Less less, int rank, int ranks, MPI_Comm comm,
                         size_t workspace_limit, std::vector<Key> *splitters) {
    using Pivot = SplitterPivot<Key>;
    struct Window { size_t lo, hi, upper; };
    using Count = unsigned long long;

    if (ranks < 2) {
        splitters->clear();
        return 0;
    }
    Count local_total = 0, total = 0;
    for (size_t r = 0; r < run_count; ++r) local_total += run_size(r);
    if (MPI_Allreduce(&local_total, &total, 1, MPI_UNSIGNED_LONG_LONG,
                      MPI_SUM, comm) != MPI_SUCCESS) return -1;
    const int goals = static_cast<int>(std::min<Count>(ranks - 1, total));

    // Bound scratch space even when there are many ranks or external runs.
    const size_t capacity = static_cast<size_t>(ranks - 1);
    const size_t left_after_result = capacity <= workspace_limit / sizeof(Key)
        ? workspace_limit - capacity * sizeof(Key) : 0;
    const size_t scratch_size = std::max(run_count, rank == 0 ? size_t(ranks) : 0u);
    int width = 0;
    if (left_after_result >= sizeof(Count) &&
        scratch_size <= (left_after_result - sizeof(Count)) / sizeof(Pivot)) {
        const size_t left = left_after_result - sizeof(Count) - scratch_size * sizeof(Pivot);
        const size_t pivot_slots = 2 + (rank == 0 ? size_t(ranks) : 0);
        const size_t fixed = pivot_slots * sizeof(Pivot) + 2 * sizeof(Count);
        if (left >= fixed && run_count <= (left - fixed) / sizeof(Window)) {
            width = static_cast<int>(std::min<size_t>(
                std::max(1, std::min(8, goals)), left / (fixed + run_count * sizeof(Window))));
        }
    }
    if (!goals && capacity <= workspace_limit / sizeof(Key)) width = 1;
    int batch_width = 0;
    if (MPI_Allreduce(&width, &batch_width, 1, MPI_INT, MPI_MIN, comm) != MPI_SUCCESS ||
        batch_width == 0) return -1;

    std::vector<Window> windows;
    std::vector<Pivot> proposals, pivots, gathered, scratch;
    std::vector<Count> targets, counts;
    int ok = 1, all_ok = 0;
    try {
        splitters->assign(capacity, Key{});
        if (goals != 0) {
            windows.resize(run_count * batch_width);
            proposals.resize(batch_width);
            pivots.resize(batch_width);
            if (rank == 0) gathered.resize(size_t(ranks) * batch_width);
            scratch.reserve(scratch_size);
            targets.resize(batch_width);
            counts.resize(batch_width + 1);
        }
    } catch (...) {
        ok = 0;
    }
    if (MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, comm) != MPI_SUCCESS ||
        !all_ok) return -1;

    for (int first = 0; first < goals; ) {
        const int active = std::min(batch_width, goals - first);
        for (int s = 0; s < active; ++s) {
            const Count j = static_cast<Count>(first) + s + 1;
            targets[s] = j * (total / ranks) + std::min(j, total % ranks);
            for (size_t r = 0; r < run_count; ++r)
                windows[size_t(s) * run_count + r] = {0, run_size(r), 0};
        }
        for (int unresolved = active; unresolved != 0; ) {
            // Weight midpoints by their remaining window lengths. The root
            // sees only one such representative per rank and target.
            for (int s = 0; s < active; ++s) {
                proposals[s] = Pivot{};
                if (!targets[s]) continue;
                scratch.clear();
                for (size_t r = 0; r < run_count; ++r) {
                    const Window &w = windows[size_t(s) * run_count + r];
                    if (w.lo == w.hi) continue;
                    Pivot p = {};
                    p.weight = w.hi - w.lo;
                    p.run = r;
                    p.index = w.lo + (w.hi - w.lo) / 2;
                    p.rank = rank;
                    if (!read_key(r, p.index, &p.key)) {
                        proposals[s].rank = -1;
                        break;
                    }
                    scratch.push_back(p);
                }
                if (proposals[s].rank >= 0) proposals[s] = WeightedMedian(&scratch, less);
            }
            const int pivot_bytes = static_cast<int>(active * sizeof(Pivot));
            if (MPI_Gather(proposals.data(), pivot_bytes, MPI_BYTE, gathered.data(),
                           pivot_bytes, MPI_BYTE, 0, comm) != MPI_SUCCESS) return -1;
            if (rank == 0) {
                for (int s = 0; s < active; ++s) {
                    pivots[s] = Pivot{};
                    if (!targets[s]) continue;
                    scratch.clear();
                    for (int r = 0; r < ranks; ++r) {
                        const Pivot &p = gathered[size_t(r) * active + s];
                        if (p.rank < 0) {
                            pivots[s].rank = -1;
                            break;
                        }
                        if (p.weight) scratch.push_back(p);
                    }
                    if (pivots[s].rank >= 0) pivots[s] = WeightedMedian(&scratch, less);
                }
            }
            if (MPI_Bcast(pivots.data(), pivot_bytes, MPI_BYTE, 0, comm) != MPI_SUCCESS)
                return -1;
            std::fill(counts.begin(), counts.begin() + active + 1, 0);
            for (int s = 0; s < active; ++s) {
                if (!targets[s]) continue;
                const Pivot &p = pivots[s];
                if (!p.weight || p.rank < 0) return -1;
                for (size_t r = 0; r < run_count; ++r) {
                    Window &w = windows[size_t(s) * run_count + r];
                    if (p.rank == rank && p.run == r) {
                        w.upper = p.index + 1;
                    } else {
                        size_t lo = w.lo, hi = w.hi;
                        while (lo < hi) {
                            const size_t mid = lo + (hi - lo) / 2;
                            Key key;
                            if (!read_key(r, mid, &key)) {
                                counts[active] = 1;
                                break;
                            }
                            if (less(p.key, key)) hi = mid;
                            else lo = mid + 1;
                        }
                        w.upper = lo;
                    }
                    counts[s] += w.upper;
                }
            }
            if (MPI_Allreduce(MPI_IN_PLACE, counts.data(), active + 1,
                              MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm) != MPI_SUCCESS ||
                counts[active]) return -1;
            for (int s = 0; s < active; ++s) {
                if (!targets[s]) continue;
                const Pivot &p = pivots[s];
                if (counts[s] == targets[s]) {
                    (*splitters)[first + s] = p.key;
                    targets[s] = 0;
                    --unresolved;
                } else {
                    // The absolute upper bounds count keys outside the active window too.
                    for (size_t r = 0; r < run_count; ++r) {
                        Window &w = windows[size_t(s) * run_count + r];
                        if (counts[s] < targets[s]) w.lo = w.upper;
                        else w.hi = w.upper - (p.rank == rank && p.run == r ? 1 : 0);
                    }
                }
            }
        }
        first += active;
    }
    if (goals != 0)
        std::fill(splitters->begin() + goals, splitters->end(), (*splitters)[goals - 1]);
    return 0;
}

} // namespace sort_detail
} // namespace swbam

#endif
