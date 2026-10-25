#include "pipeline.hpp"

#include <algorithm>
#include <chrono>
#include <stdexcept>

namespace isp {

PipelineResult run_pipeline(const RawImage& raw, const RawMeta& meta, const PipelineConfig& cfg,
                            bool keep_intermediates) {
    if (raw.width != meta.width || raw.height != meta.height) {
        throw std::runtime_error("raw image is " + std::to_string(raw.width) + "x" + std::to_string(raw.height) +
                                 " but metadata says " + std::to_string(meta.width) + "x" + std::to_string(meta.height));
    }
    if (raw.width < 4 || raw.height < 4) throw std::runtime_error("raw image must be at least 4x4");

    using Clock = std::chrono::steady_clock;
    PipelineResult result;
    auto timed = [&](const char* name, auto&& stage) {
        const auto t0 = Clock::now();
        stage();
        const std::chrono::duration<double, std::milli> elapsed = Clock::now() - t0;
        result.timings.push_back({name, elapsed.count()});
    };

    Plane mosaic;
    timed("Black Level", [&] { mosaic = black_level_correct(raw, meta, cfg.threads); });
    if (cfg.defect_correction) {
        timed("Defect Correction", [&] {
            result.defects_corrected = correct_defective_pixels(mosaic, cfg.defect_threshold, cfg.threads);
        });
    }
    timed("White Balance", [&] {
        result.wb_gains = cfg.gray_world ? gray_world_gains(mosaic, meta.cfa) : meta.wb_gains;
        apply_white_balance(mosaic, meta.cfa, result.wb_gains, cfg.threads);
    });
    if (keep_intermediates) result.mosaic_after_wb = mosaic;

    RgbImage rgb;
    timed("Demosaic", [&] { rgb = demosaic(mosaic, meta.cfa, cfg.demosaic, cfg.threads); });
    if (keep_intermediates) result.camera_rgb = rgb;

    timed("Color Correction", [&] { apply_color_correction(rgb, meta.ccm, cfg.threads); });
    timed("Gamma", [&] { result.output = gamma_encode(rgb, cfg.threads); });
    return result;
}

std::vector<BenchRow> benchmark(const RawImage& raw, const RawMeta& meta, PipelineConfig cfg,
                                const std::vector<int>& thread_counts, int repeats) {
    std::vector<BenchRow> rows;
    for (int t : thread_counts) {
        cfg.threads = t;
        run_pipeline(raw, meta, cfg);  // warm-up: page in buffers and build the gamma table
        std::vector<std::vector<double>> samples;
        std::vector<std::string> names;
        for (int r = 0; r < repeats; ++r) {
            const PipelineResult res = run_pipeline(raw, meta, cfg);
            if (samples.empty()) {
                samples.resize(res.timings.size());
                for (const StageTiming& s : res.timings) names.push_back(s.name);
            }
            for (size_t i = 0; i < res.timings.size(); ++i) samples[i].push_back(res.timings[i].ms);
        }
        BenchRow row;
        row.threads = t;
        for (size_t i = 0; i < samples.size(); ++i) {
            std::vector<double>& s = samples[i];
            std::sort(s.begin(), s.end());
            const double median = s[s.size() / 2];
            row.stage_ms.push_back({names[i], median});
            row.total_ms += median;
        }
        rows.push_back(row);
    }
    return rows;
}

}  // namespace isp
