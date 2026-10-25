#pragma once

#include <string>
#include <vector>

#include "image.hpp"
#include "parallel.hpp"
#include "raw_meta.hpp"
#include "stages.hpp"

namespace isp {

struct PipelineConfig {
    DemosaicMethod demosaic = DemosaicMethod::Malvar;
    bool defect_correction = true;
    float defect_threshold = 0.2f;   // in normalised units, where 1.0 is the full sensor range
    bool gray_world = false;         // estimate white balance from the image instead of using the metadata gains
    int threads = hardware_threads();
};

struct StageTiming {
    std::string name;
    double ms = 0.0;
};

struct PipelineResult {
    Rgb8 output;
    std::vector<StageTiming> timings;
    int defects_corrected = 0;
    Gains wb_gains{1.0f, 1.0f, 1.0f};
    // Filled only when intermediates are requested, for the stage-by-stage images.
    Plane mosaic_after_wb;
    RgbImage camera_rgb;
};

// RAW -> black level -> defect correction -> white balance -> demosaic -> colour correction -> gamma -> sRGB.
PipelineResult run_pipeline(const RawImage& raw, const RawMeta& meta, const PipelineConfig& cfg,
                            bool keep_intermediates = false);

struct BenchRow {
    int threads = 1;
    std::vector<StageTiming> stage_ms;  // median per stage over the repeats
    double total_ms = 0.0;              // sum of the per-stage medians
};

// Runs the full pipeline `repeats` times at each thread count (after one warm-up run) and keeps the medians.
std::vector<BenchRow> benchmark(const RawImage& raw, const RawMeta& meta, PipelineConfig cfg,
                                const std::vector<int>& thread_counts, int repeats);

}  // namespace isp
