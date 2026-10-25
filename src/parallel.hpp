#pragma once

#include <algorithm>
#include <thread>
#include <vector>

namespace isp {

// Splits rows [0, rows) into contiguous bands, one per thread, and runs fn(y0, y1) on each band.
// Every stage writes only to its own output rows, so bands never overlap and no locking is needed.
// The calling thread processes the first band itself instead of sitting idle.
template <typename Fn>
void parallel_rows(int rows, int threads, Fn&& fn) {
    threads = std::max(1, std::min(threads, rows));
    if (threads == 1) {
        fn(0, rows);
        return;
    }
    const int band = (rows + threads - 1) / threads;
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(threads - 1));
    for (int t = 1; t < threads; ++t) {
        const int y0 = t * band;
        const int y1 = std::min(rows, y0 + band);
        if (y0 >= y1) break;
        workers.emplace_back([&fn, y0, y1] { fn(y0, y1); });
    }
    fn(0, std::min(rows, band));
    for (std::thread& w : workers) w.join();
}

inline int hardware_threads() {
    const unsigned n = std::thread::hardware_concurrency();
    return n == 0 ? 1 : static_cast<int>(n);
}

}  // namespace isp
