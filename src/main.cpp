// camera-isp-sim: command-line front end.
//
//   isp_sim demo <out_dir>                    generate a scene, simulate raws, run and score the pipeline, benchmark it
//   isp_sim scene <out.ppm>                   write the synthetic test chart
//   isp_sim make-raw <in.ppm> <out_prefix>    simulate a sensor raw (<prefix>.pgm + <prefix>.txt) from an sRGB image
//   isp_sim run <raw.pgm> <meta.txt> <out>    process a raw into an sRGB .bmp or .ppm
//   isp_sim bench <raw.pgm> <meta.txt>        time each stage across thread counts

#include <algorithm>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "analysis.hpp"
#include "image.hpp"
#include "pipeline.hpp"
#include "synth.hpp"

namespace fs = std::filesystem;
using namespace isp;

namespace {

struct Args {
    std::vector<std::string> positional;
    std::map<std::string, std::string> options;  // --key value, or --flag mapped to "1"

    bool has(const std::string& key) const { return options.count(key) != 0; }
    std::string get(const std::string& key, const std::string& fallback) const {
        auto it = options.find(key);
        return it == options.end() ? fallback : it->second;
    }
    int get_int(const std::string& key, int fallback) const { return has(key) ? std::stoi(options.at(key)) : fallback; }
    float get_float(const std::string& key, float fallback) const {
        return has(key) ? std::stof(options.at(key)) : fallback;
    }
};

const std::vector<std::string> kFlags = {"--no-dpc", "--gray-world"};

Args parse_args(int argc, char** argv, int first) {
    Args a;
    for (int i = first; i < argc; ++i) {
        const std::string s = argv[i];
        if (s.rfind("--", 0) == 0) {
            const bool is_flag = std::find(kFlags.begin(), kFlags.end(), s) != kFlags.end();
            if (is_flag) {
                a.options[s] = "1";
            } else {
                if (i + 1 >= argc) throw std::runtime_error("option " + s + " needs a value");
                a.options[s] = argv[++i];
            }
        } else {
            a.positional.push_back(s);
        }
    }
    return a;
}

void need_positional(const Args& a, size_t n, const char* usage) {
    if (a.positional.size() < n) throw std::runtime_error(std::string("usage: ") + usage);
}

DemosaicMethod parse_demosaic(const std::string& s) {
    if (s == "bilinear") return DemosaicMethod::Bilinear;
    if (s == "malvar" || s == "mhc") return DemosaicMethod::Malvar;
    throw std::runtime_error("unknown demosaic method '" + s + "' (use bilinear or malvar)");
}

std::string fmt(double v, int decimals) {
    if (v == std::numeric_limits<double>::infinity()) return "inf";
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(decimals) << v;
    return ss.str();
}

std::string compiler_name() {
#if defined(_MSC_VER)
    return "MSVC " + std::to_string(_MSC_VER);
#elif defined(__clang__)
    return std::string("Clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("GCC ") + __VERSION__;
#else
    return "unknown compiler";
#endif
}

std::vector<int> default_thread_counts() {
    const int max_t = hardware_threads();
    std::vector<int> counts;
    for (int t = 1; t < max_t; t *= 2) counts.push_back(t);
    counts.push_back(max_t);
    return counts;
}

void print_timings(const std::vector<StageTiming>& timings) {
    double total = 0.0;
    for (const StageTiming& s : timings) {
        std::cout << "  " << std::left << std::setw(18) << s.name << std::right << std::setw(9) << fmt(s.ms, 2) << " ms\n";
        total += s.ms;
    }
    std::cout << "  " << std::left << std::setw(18) << "Total" << std::right << std::setw(9) << fmt(total, 2) << " ms\n";
}

void write_raw_files(const std::string& prefix, const RawGenResult& r) {
    write_pgm16(prefix + ".pgm", r.raw, r.meta.white_level);
    write_meta(prefix + ".txt", r.meta);
}

std::string bench_table(const std::vector<BenchRow>& rows, double megapixels) {
    std::ostringstream md;
    md << "| Threads | Total (ms) | Speedup | Throughput (MP/s) |\n| --- | --- | --- | --- |\n";
    for (const BenchRow& r : rows) {
        md << "| " << r.threads << " | " << fmt(r.total_ms, 1) << " | " << fmt(rows.front().total_ms / r.total_ms, 2)
           << "x | " << fmt(megapixels / (r.total_ms / 1000.0), 1) << " |\n";
    }
    md << "\nPer-stage medians (ms):\n\n| Stage |";
    for (const BenchRow& r : rows) md << ' ' << r.threads << (r.threads == 1 ? " thread |" : " threads |");
    md << "\n| --- |";
    for (size_t i = 0; i < rows.size(); ++i) md << " --- |";
    md << '\n';
    for (size_t s = 0; s < rows.front().stage_ms.size(); ++s) {
        md << "| " << rows.front().stage_ms[s].name << " |";
        for (const BenchRow& r : rows) md << ' ' << fmt(r.stage_ms[s].ms, 2) << " |";
        md << '\n';
    }
    return md.str();
}

// ---------------------------------------------------------------------------------------------------------------

int cmd_scene(const Args& a) {
    need_positional(a, 1, "isp_sim scene <out.ppm|out.bmp> [--width 1280 --height 960]");
    const Rgb8 scene = make_test_scene(a.get_int("--width", 1280), a.get_int("--height", 960));
    write_image(a.positional[0], scene);
    std::cout << "Wrote " << a.positional[0] << " (" << scene.width << "x" << scene.height << ")\n";
    return 0;
}

int cmd_make_raw(const Args& a) {
    need_positional(a, 2, "isp_sim make-raw <in.ppm> <out_prefix> [--noise 1] [--hot-fraction 0.0005] [--seed 7] "
                          "[--pattern RGGB]");
    const Rgb8 scene = read_ppm(a.positional[0]);
    RawGenConfig cfg;
    cfg.noise_scale = a.get_float("--noise", 1.0f);
    cfg.hot_pixel_fraction = a.get_float("--hot-fraction", 0.0005f);
    cfg.seed = static_cast<uint64_t>(a.get_int("--seed", 7));
    cfg.meta.cfa = Cfa::parse(a.get("--pattern", "RGGB"));
    const RawGenResult r = make_raw(scene, cfg);
    write_raw_files(a.positional[1], r);
    std::cout << "Wrote " << a.positional[1] << ".pgm and " << a.positional[1] << ".txt (" << r.raw.width << "x"
              << r.raw.height << ", " << r.meta.cfa.name() << ", " << r.hot_pixels << " hot pixels)\n";
    return 0;
}

PipelineConfig config_from_args(const Args& a) {
    PipelineConfig cfg;
    cfg.demosaic = parse_demosaic(a.get("--demosaic", "malvar"));
    cfg.defect_correction = !a.has("--no-dpc");
    cfg.defect_threshold = a.get_float("--dpc-threshold", cfg.defect_threshold);
    cfg.gray_world = a.has("--gray-world");
    cfg.threads = a.get_int("--threads", hardware_threads());
    return cfg;
}

int cmd_run(const Args& a) {
    need_positional(a, 3, "isp_sim run <raw.pgm> <meta.txt> <out.bmp|out.ppm> [--demosaic bilinear|malvar] "
                          "[--threads N] [--no-dpc] [--dpc-threshold 0.25] [--gray-world] [--ref truth.ppm] "
                          "[--dump-dir DIR]");
    const RawImage raw = read_pgm16(a.positional[0]);
    const RawMeta meta = read_meta(a.positional[1]);
    const PipelineConfig cfg = config_from_args(a);
    const bool dump = a.has("--dump-dir");
    const PipelineResult r = run_pipeline(raw, meta, cfg, dump);
    write_image(a.positional[2], r.output);

    std::cout << "Processed " << raw.width << "x" << raw.height << " " << meta.cfa.name() << " raw with "
              << cfg.threads << " thread(s)\n";
    print_timings(r.timings);
    if (cfg.defect_correction) std::cout << "  Defective pixels corrected: " << r.defects_corrected << '\n';
    std::cout << "  White balance gains (R, G, B): " << fmt(r.wb_gains[0], 3) << ", " << fmt(r.wb_gains[1], 3) << ", "
              << fmt(r.wb_gains[2], 3) << '\n';
    if (a.has("--ref")) {
        const Rgb8 ref = read_ppm(a.get("--ref", ""));
        std::cout << "  PSNR vs reference (4 px border excluded): " << fmt(psnr(r.output, ref, 4), 2) << " dB\n";
    }
    if (dump) {
        const fs::path dir = a.get("--dump-dir", "");
        fs::create_directories(dir);
        write_bmp((dir / "1_mosaic.bmp").string(), mosaic_as_gray(r.mosaic_after_wb));
        write_bmp((dir / "2_mosaic_color.bmp").string(), mosaic_as_color(r.mosaic_after_wb, meta.cfa));
        write_bmp((dir / "3_demosaic_camera_rgb.bmp").string(), linear_to_display(r.camera_rgb));
        write_bmp((dir / "4_final.bmp").string(), r.output);
        std::cout << "  Stage images written to " << dir.string() << '\n';
    }
    std::cout << "Wrote " << a.positional[2] << '\n';
    return 0;
}

int cmd_bench(const Args& a) {
    need_positional(a, 2, "isp_sim bench <raw.pgm> <meta.txt> [--repeats 5] [--demosaic malvar]");
    const RawImage raw = read_pgm16(a.positional[0]);
    const RawMeta meta = read_meta(a.positional[1]);
    PipelineConfig cfg = config_from_args(a);
    const auto rows = benchmark(raw, meta, cfg, default_thread_counts(), a.get_int("--repeats", 5));
    const double mp = static_cast<double>(raw.width) * raw.height / 1e6;
    std::cout << bench_table(rows, mp);
    return 0;
}

int cmd_demo(const Args& a) {
    need_positional(a, 1, "isp_sim demo <out_dir> [--width 1280 --height 960] [--bench-width 4032 --bench-height 3024] "
                          "[--repeats 5]");
    const fs::path out = a.positional[0];
    fs::create_directories(out / "stages");
    const int width = a.get_int("--width", 1280);
    const int height = a.get_int("--height", 960);
    const int threads = hardware_threads();
    auto path = [&](const std::string& name) { return (out / name).string(); };

    std::cout << "[1/5] Generating " << width << "x" << height << " test scene and simulated raws\n";
    const Rgb8 scene = make_test_scene(width, height);
    write_ppm(path("scene_truth.ppm"), scene);
    write_bmp(path("scene_truth.bmp"), scene);

    RawGenConfig clean_cfg;
    clean_cfg.noise_scale = 0.0f;
    clean_cfg.hot_pixel_fraction = 0.0f;
    const RawGenResult clean = make_raw(scene, clean_cfg);
    write_raw_files(path("raw_clean"), clean);

    RawGenConfig noisy_cfg;
    noisy_cfg.noise_scale = 1.0f;
    noisy_cfg.hot_pixel_fraction = 0.0005f;
    const RawGenResult noisy = make_raw(scene, noisy_cfg);
    write_raw_files(path("raw_noisy"), noisy);

    const SceneLayout layout = scene_layout(width, height);
    const std::vector<Rect> fine_detail = {layout.zone_plate, layout.siemens_star};
    const std::vector<Rect> flat_color = {layout.color_patches};

    std::cout << "[2/5] Demosaic comparison on the noise-free raw\n";
    PipelineConfig cfg;
    cfg.threads = threads;
    cfg.defect_correction = false;  // only the demosaic method differs between these two runs
    cfg.demosaic = DemosaicMethod::Bilinear;
    const PipelineResult bil = run_pipeline(clean.raw, clean.meta, cfg);
    cfg.demosaic = DemosaicMethod::Malvar;
    const PipelineResult mhc = run_pipeline(clean.raw, clean.meta, cfg, true);
    write_bmp(path("out_bilinear.bmp"), bil.output);
    write_bmp(path("out_malvar.bmp"), mhc.output);

    const double bil_full = psnr(bil.output, scene, 4), mhc_full = psnr(mhc.output, scene, 4);
    const double bil_fine = psnr(bil.output, scene, fine_detail), mhc_fine = psnr(mhc.output, scene, fine_detail);
    const double bil_flat = psnr(bil.output, scene, flat_color), mhc_flat = psnr(mhc.output, scene, flat_color);

    // Colour accuracy check: away from edges, a noise-free raw must reproduce every patch exactly.
    int max_patch_error = 0;
    for (const Rect& r : color_patch_rects(width, height, 4)) {
        for (int y = r.y0; y < r.y1; ++y)
            for (int x = r.x0; x < r.x1; ++x)
                for (int c = 0; c < 3; ++c) {
                    const size_t i = (static_cast<size_t>(y) * width + x) * 3 + static_cast<size_t>(c);
                    max_patch_error = std::max(max_patch_error, std::abs(mhc.output.data[i] - scene.data[i]));
                }
    }

    std::cout << "[3/5] Full pipeline on the noisy raw (with " << noisy.hot_pixels << " injected hot pixels)\n";
    cfg.defect_correction = false;
    const PipelineResult noisy_no_dpc = run_pipeline(noisy.raw, noisy.meta, cfg);
    cfg.defect_correction = true;
    const PipelineResult noisy_dpc = run_pipeline(noisy.raw, noisy.meta, cfg);
    const PipelineResult clean_dpc = run_pipeline(clean.raw, clean.meta, cfg);  // no defects: anything fixed is a false positive
    write_bmp(path("out_noisy_no_dpc.bmp"), noisy_no_dpc.output);
    write_bmp(path("out_noisy.bmp"), noisy_dpc.output);
    // Score defect correction away from the zone plate, whose outer ring holds detail finer than the
    // same-colour sampling grid; there, real detail and single-pixel defects cannot be told apart.
    const std::vector<Rect> outside_zone = {layout.color_patches, layout.ramps, layout.siemens_star};
    const double noisy_no_dpc_psnr = psnr(noisy_no_dpc.output, scene, outside_zone);
    const double noisy_dpc_psnr = psnr(noisy_dpc.output, scene, outside_zone);
    // False positives: on the defect-free raw, any output pixel that defect correction changed.
    int changed_outside_zone = 0, changed_total = 0;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t i = (static_cast<size_t>(y) * width + x) * 3;
            const bool changed = clean_dpc.output.data[i] != mhc.output.data[i] ||
                                 clean_dpc.output.data[i + 1] != mhc.output.data[i + 1] ||
                                 clean_dpc.output.data[i + 2] != mhc.output.data[i + 2];
            if (!changed) continue;
            ++changed_total;
            const Rect& z = layout.zone_plate;
            if (x < z.x0 || x >= z.x1 || y < z.y0 || y >= z.y1) ++changed_outside_zone;
        }
    }

    cfg.gray_world = true;
    const PipelineResult gw = run_pipeline(clean.raw, clean.meta, cfg);
    cfg.gray_world = false;

    std::cout << "[4/5] Writing stage images and zoomed comparisons\n";
    write_bmp((out / "stages" / "1_mosaic.bmp").string(), mosaic_as_gray(mhc.mosaic_after_wb));
    write_bmp((out / "stages" / "2_mosaic_color.bmp").string(), mosaic_as_color(mhc.mosaic_after_wb, clean.meta.cfa));
    write_bmp((out / "stages" / "3_demosaic_camera_rgb.bmp").string(), linear_to_display(mhc.camera_rgb));
    write_bmp((out / "stages" / "4_final.bmp").string(), mhc.output);

    // Truth | bilinear | Malvar, zoomed 4x, at the star centre and the zone plate's high-frequency edge.
    const Rect& star = layout.siemens_star;
    const int scx = (star.x0 + star.x1) / 2, scy = (star.y0 + star.y1) / 2;
    const Rect star_crop{scx - 48, scy - 48, scx + 48, scy + 48};
    const Rect& zp = layout.zone_plate;
    const int zr = std::min(zp.x1 - zp.x0, zp.y1 - zp.y0) / 2;
    const int zcx = (zp.x0 + zp.x1) / 2, zcy = (zp.y0 + zp.y1) / 2;
    const Rect zone_crop{zcx + zr - 100, zcy - 48, zcx + zr - 4, zcy + 48};
    for (const auto& [name, rect] : {std::pair<const char*, Rect>{"compare_star.bmp", star_crop},
                                     std::pair<const char*, Rect>{"compare_zone_plate.bmp", zone_crop}}) {
        write_bmp(path(name), side_by_side({crop_zoom(scene, rect, 4), crop_zoom(bil.output, rect, 4),
                                            crop_zoom(mhc.output, rect, 4)}, 8));
    }

    const int bw = a.get_int("--bench-width", 4032), bh = a.get_int("--bench-height", 3024);
    const int repeats = a.get_int("--repeats", 5);
    std::cout << "[5/5] Benchmarking a " << bw << "x" << bh << " raw (" << repeats << " runs per thread count)\n";
    const Rgb8 big_scene = make_test_scene(bw, bh);
    const RawGenResult big = make_raw(big_scene, noisy_cfg);
    const auto rows = benchmark(big.raw, big.meta, PipelineConfig{}, default_thread_counts(), repeats);
    const double mp = static_cast<double>(bw) * bh / 1e6;

    // ---- report ----
    std::ostringstream md;
    md << "# Results\n\n";
    md << "Generated by `isp_sim demo` on a machine with " << threads << " hardware threads (" << compiler_name()
       << ").\n\n";
    md << "## Demosaic Quality\n\n";
    md << "PSNR of the final sRGB output against the ground-truth scene, from a noise-free " << width << "x" << height
       << " raw. Higher is better.\n\n";
    md << "| Region | Bilinear (dB) | Malvar-He-Cutler (dB) | Gain (dB) |\n| --- | --- | --- | --- |\n";
    md << "| Whole image (4 px border excluded) | " << fmt(bil_full, 2) << " | " << fmt(mhc_full, 2) << " | "
       << fmt(mhc_full - bil_full, 2) << " |\n";
    md << "| Fine detail (zone plate and Siemens star) | " << fmt(bil_fine, 2) << " | " << fmt(mhc_fine, 2) << " | "
       << fmt(mhc_fine - bil_fine, 2) << " |\n";
    md << "| Colour patch quadrant, including patch edges | " << fmt(bil_flat, 2) << " | " << fmt(mhc_flat, 2) << " | "
       << fmt(mhc_flat - bil_flat, 2) << " |\n\n";
    md << "The chart deliberately contains detail finer than a Bayer sensor can resolve (red and blue are sampled at "
          "half resolution, and the zone plate reaches the full-resolution Nyquist limit), so absolute PSNR is low for "
          "every method; the gap between methods is the meaningful number.\n\n";
    md << "Colour accuracy check: inside the 24 patches (4 px away from edges) the largest error is "
       << max_patch_error << " code values out of 255.\n\n";
    md << "## Noise and Defective Pixels\n\n";
    md << "Noisy raw (read noise " << noisy_cfg.read_noise_dn << " DN plus shot noise) with " << noisy.hot_pixels
       << " injected hot pixels, Malvar demosaic.\n\n";
    md << "PSNR is measured outside the zone plate (see below).\n\n";
    md << "| Run | PSNR (dB) | Pixels corrected |\n| --- | --- | --- |\n";
    md << "| Without defect correction | " << fmt(noisy_no_dpc_psnr, 2) << " | 0 |\n";
    md << "| With defect correction | " << fmt(noisy_dpc_psnr, 2) << " | " << noisy_dpc.defects_corrected << " |\n\n";
    md << "False positives: on the noise-free raw, which has no defects, the stage still flagged "
       << clean_dpc.defects_corrected << " pixels. Of the output pixels this changed, " << changed_outside_zone
       << " of " << changed_total << " lie outside the zone plate; the rest are in its outer ring, where the "
          "pattern is finer than the same-colour sampling grid and a single-channel check cannot separate detail "
          "from defects.\n\n";
    md << "## White Balance\n\n";
    md << "| Gains | R | G | B |\n| --- | --- | --- | --- |\n";
    md << "| True (from metadata) | " << fmt(clean.meta.wb_gains[0], 3) << " | " << fmt(clean.meta.wb_gains[1], 3)
       << " | " << fmt(clean.meta.wb_gains[2], 3) << " |\n";
    md << "| Gray-world estimate | " << fmt(gw.wb_gains[0], 3) << " | " << fmt(gw.wb_gains[1], 3) << " | "
       << fmt(gw.wb_gains[2], 3) << " |\n\n";
    md << "Gray-world output PSNR: " << fmt(psnr(gw.output, scene, 4), 2)
       << " dB. The chart is not grey on average, so gray-world is expected to miss here.\n\n";
    md << "## Performance\n\n";
    md << "Full pipeline on a " << bw << "x" << bh << " (" << fmt(mp, 1) << " MP) raw, median of " << repeats
       << " runs, CPU only.\n\n";
    md << bench_table(rows, mp);
    const std::string report = md.str();

    std::ofstream(path("results.md")) << report;
    std::cout << '\n' << report << "\nAll outputs are in " << out.string() << '\n';
    return 0;
}

void print_usage() {
    std::cout << "camera-isp-sim: a camera image signal processor pipeline in C++17\n\n"
                 "Commands:\n"
                 "  isp_sim demo <out_dir>                   end-to-end demo: scene, raws, quality scores, benchmark\n"
                 "  isp_sim scene <out.ppm|out.bmp>          write the synthetic test chart\n"
                 "  isp_sim make-raw <in.ppm> <out_prefix>   simulate a sensor raw from an sRGB image\n"
                 "  isp_sim run <raw.pgm> <meta.txt> <out>   process a raw into an sRGB .bmp or .ppm\n"
                 "  isp_sim bench <raw.pgm> <meta.txt>       time each stage across thread counts\n\n"
                 "Run a command without arguments to see its options.\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage();
        return 1;
    }
    const std::string cmd = argv[1];
    try {
        const Args args = parse_args(argc, argv, 2);
        if (cmd == "demo") return cmd_demo(args);
        if (cmd == "scene") return cmd_scene(args);
        if (cmd == "make-raw") return cmd_make_raw(args);
        if (cmd == "run") return cmd_run(args);
        if (cmd == "bench") return cmd_bench(args);
        print_usage();
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
