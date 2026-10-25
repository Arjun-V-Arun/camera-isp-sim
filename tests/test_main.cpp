// Unit tests for camera-isp-sim. No framework needed: each CHECK records a failure and the run exits non-zero.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "analysis.hpp"
#include "image.hpp"
#include "pipeline.hpp"
#include "stages.hpp"
#include "synth.hpp"

using namespace isp;

namespace {

int g_failures = 0;
int g_checks = 0;

#define CHECK(cond)                                                                  \
    do {                                                                             \
        ++g_checks;                                                                  \
        if (!(cond)) {                                                               \
            ++g_failures;                                                            \
            std::cerr << __FILE__ << ':' << __LINE__ << ": CHECK failed: " #cond "\n"; \
        }                                                                            \
    } while (0)

#define CHECK_NEAR(a, b, tol) CHECK(std::fabs(static_cast<double>(a) - static_cast<double>(b)) <= (tol))

const char* kPatterns[] = {"RGGB", "BGGR", "GRBG", "GBRG"};

// A mosaic sampled from an image whose R, G and B are given functions of (x, y).
template <typename F>
Plane mosaic_from(int w, int h, const Cfa& cfa, F value_at) {
    Plane p;
    p.width = w;
    p.height = h;
    p.data.resize(static_cast<size_t>(w) * h);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) p.data[static_cast<size_t>(y) * w + x] = value_at(x, y, cfa.at(x, y));
    return p;
}

void test_reflect101() {
    CHECK(reflect101(-1, 10) == 1);
    CHECK(reflect101(-2, 10) == 2);
    CHECK(reflect101(10, 10) == 8);
    CHECK(reflect101(11, 10) == 7);
    CHECK(reflect101(5, 10) == 5);
    // Parity is preserved, so a reflected Bayer sample keeps its colour.
    for (int i = -4; i < 14; ++i) CHECK(((reflect101(i, 10) - i) % 2) == 0);
}

void test_cfa_parse() {
    for (const char* name : kPatterns) CHECK(Cfa::parse(name).name() == name);
    bool threw = false;
    try {
        Cfa::parse("RGBG");  // greens not on a diagonal
    } catch (const std::exception&) {
        threw = true;
    }
    CHECK(threw);
    CHECK(Cfa::parse("RGGB").row_color(0) == kRed);
    CHECK(Cfa::parse("RGGB").row_color(1) == kBlue);
}

void test_black_level() {
    RawMeta m;
    m.width = 4;
    m.height = 1;
    RawImage raw{4, 1, {static_cast<uint16_t>(m.black_level), static_cast<uint16_t>(m.white_level), 0,
                        static_cast<uint16_t>((m.black_level + m.white_level) / 2)}};
    const Plane p = black_level_correct(raw, m, 1);
    CHECK_NEAR(p.data[0], 0.0, 1e-6);
    CHECK_NEAR(p.data[1], 1.0, 1e-6);
    CHECK_NEAR(p.data[2], 0.0, 1e-6);  // below black is clamped
    CHECK_NEAR(p.data[3], 0.5, 1e-3);
}

void test_white_balance_positions() {
    for (const char* name : kPatterns) {
        const Cfa cfa = Cfa::parse(name);
        Plane p = mosaic_from(6, 6, cfa, [](int, int, int) { return 0.25f; });
        apply_white_balance(p, cfa, Gains{2.0f, 1.0f, 3.0f}, 2);
        for (int y = 0; y < 6; ++y) {
            for (int x = 0; x < 6; ++x) {
                const float expect = 0.25f * (cfa.at(x, y) == kRed ? 2.0f : cfa.at(x, y) == kBlue ? 3.0f : 1.0f);
                CHECK_NEAR(p.data[static_cast<size_t>(y) * 6 + x], expect, 1e-6);
            }
        }
    }
}

// A uniform colour must come back exactly, everywhere, for both methods and all four patterns.
void test_demosaic_flat_field() {
    const float rgb[3] = {0.7f, 0.4f, 0.2f};
    for (const char* name : kPatterns) {
        const Cfa cfa = Cfa::parse(name);
        const Plane p = mosaic_from(16, 12, cfa, [&](int, int, int ch) { return rgb[ch]; });
        for (DemosaicMethod m : {DemosaicMethod::Bilinear, DemosaicMethod::Malvar}) {
            const RgbImage out = demosaic(p, cfa, m, 3);
            double worst = 0.0;
            for (size_t i = 0; i < out.data.size(); ++i) worst = std::max(worst, static_cast<double>(std::fabs(out.data[i] - rgb[i % 3])));
            CHECK(worst < 1e-5);
        }
    }
}

// Both methods are exact on a grey linear ramp away from the border (the kernels are symmetric
// and their gradient terms cancel on a plane).
void test_demosaic_linear_ramp() {
    for (const char* name : kPatterns) {
        const Cfa cfa = Cfa::parse(name);
        const Plane p = mosaic_from(20, 20, cfa, [](int x, int y, int) { return 0.1f + 0.02f * x + 0.01f * y; });
        for (DemosaicMethod m : {DemosaicMethod::Bilinear, DemosaicMethod::Malvar}) {
            const RgbImage out = demosaic(p, cfa, m, 1);
            double worst = 0.0;
            for (int y = 2; y < 18; ++y)
                for (int x = 2; x < 18; ++x)
                    for (int c = 0; c < 3; ++c) {
                        const double truth = 0.1 + 0.02 * x + 0.01 * y;
                        worst = std::max(worst, std::fabs(static_cast<double>(out.data[(static_cast<size_t>(y) * 20 + x) * 3 + c]) - truth));
                    }
            CHECK(worst < 1e-5);
        }
    }
}

void test_defect_correction() {
    const Cfa cfa = Cfa::parse("RGGB");
    Plane p = mosaic_from(10, 10, cfa, [](int, int, int ch) { return ch == kGreen ? 0.4f : 0.3f; });
    // Nothing to fix in a clean flat field.
    Plane clean = p;
    CHECK(correct_defective_pixels(clean, 0.25f, 2) == 0);
    // One hot and one dead pixel.
    p.data[5 * 10 + 4] = 1.0f;  // (4, 5) is green in RGGB
    p.data[2 * 10 + 2] = 0.0f;  // (2, 2) is red in RGGB
    const int fixed = correct_defective_pixels(p, 0.25f, 2);
    CHECK(fixed == 2);
    CHECK_NEAR(p.data[5 * 10 + 4], 0.4, 1e-6);
    CHECK_NEAR(p.data[2 * 10 + 2], 0.3, 1e-6);
    // Vertical stripes two pixels wide reach the Bayer sampling limit, yet nothing is flagged,
    // because along each stripe the same-colour neighbours match the centre.
    Plane stripes = mosaic_from(12, 12, cfa, [](int x, int, int) { return (x / 2) % 2 ? 0.9f : 0.1f; });
    CHECK(correct_defective_pixels(stripes, 0.2f, 1) == 0);
}

void test_srgb_transfer() {
    CHECK_NEAR(srgb_encode(0.0f), 0.0, 1e-7);
    CHECK_NEAR(srgb_encode(1.0f), 1.0, 1e-6);
    CHECK_NEAR(srgb_encode(0.0031308f), 0.04045, 1e-5);  // the two segments meet here
    CHECK_NEAR(srgb_encode(0.18f), 0.4614, 1e-3);        // 18% grey encodes to about 46%
    for (int i = 0; i <= 100; ++i) {
        const float v = static_cast<float>(i) / 100.0f;
        CHECK_NEAR(srgb_encode(srgb_decode(v)), v, 1e-5);
    }
}

void test_matrix() {
    const RawMeta m;
    const Mat3 prod = mat3_multiply(m.ccm, mat3_inverse(m.ccm));
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) CHECK_NEAR(prod[static_cast<size_t>(r * 3 + c)], r == c ? 1.0 : 0.0, 1e-5);
    // Rows sum to one, so grey stays grey through the colour matrix.
    RgbImage grey{1, 1, {0.5f, 0.5f, 0.5f}};
    apply_color_correction(grey, m.ccm, 1);
    for (float v : grey.data) CHECK_NEAR(v, 0.5, 1e-6);
}

void test_threads_are_deterministic() {
    const Rgb8 scene = make_test_scene(256, 192);
    RawGenConfig cfg;
    cfg.hot_pixel_fraction = 0.001f;
    const RawGenResult r = make_raw(scene, cfg);
    PipelineConfig pc;
    pc.threads = 1;
    const PipelineResult one = run_pipeline(r.raw, r.meta, pc);
    pc.threads = 7;  // an odd count, so the bands are uneven
    const PipelineResult many = run_pipeline(r.raw, r.meta, pc);
    CHECK(one.output.data == many.output.data);
    CHECK(one.defects_corrected == many.defects_corrected);
}

void test_end_to_end_quality() {
    const Rgb8 scene = make_test_scene(320, 240);
    RawGenConfig cfg;
    cfg.noise_scale = 0.0f;
    cfg.hot_pixel_fraction = 0.0f;
    const RawGenResult r = make_raw(scene, cfg);
    PipelineConfig pc;
    pc.defect_correction = false;
    const SceneLayout L = scene_layout(320, 240);
    // Away from edges, flat colour patches must come back exactly from a noise-free raw.
    pc.demosaic = DemosaicMethod::Malvar;
    const PipelineResult mhc = run_pipeline(r.raw, r.meta, pc);
    CHECK(psnr(mhc.output, scene, color_patch_rects(320, 240, 3)) > 50.0);
    // Malvar should beat bilinear on fine detail.
    pc.demosaic = DemosaicMethod::Bilinear;
    const PipelineResult bil = run_pipeline(r.raw, r.meta, pc);
    const std::vector<Rect> fine = {L.zone_plate, L.siemens_star};
    CHECK(psnr(mhc.output, scene, fine) > psnr(bil.output, scene, fine));
}

void test_file_roundtrip() {
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "isp_sim_tests";
    std::filesystem::create_directories(dir);

    RawImage raw{3, 2, {0, 1, 255, 256, 4095, 65535}};
    write_pgm16((dir / "t.pgm").string(), raw, 65535);
    CHECK(read_pgm16((dir / "t.pgm").string()).data == raw.data);

    Rgb8 img{2, 2, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}};
    write_ppm((dir / "t.ppm").string(), img);
    CHECK(read_ppm((dir / "t.ppm").string()).data == img.data);

    RawMeta m;
    m.width = 8;
    m.height = 6;
    m.cfa = Cfa::parse("GBRG");
    write_meta((dir / "t.txt").string(), m);
    const RawMeta back = read_meta((dir / "t.txt").string());
    CHECK(back.cfa.name() == "GBRG");
    CHECK(back.width == 8 && back.height == 6);
    CHECK_NEAR(back.ccm[7], m.ccm[7], 1e-6);
    std::filesystem::remove_all(dir);
}

}  // namespace

int main() {
    test_reflect101();
    test_cfa_parse();
    test_black_level();
    test_white_balance_positions();
    test_demosaic_flat_field();
    test_demosaic_linear_ramp();
    test_defect_correction();
    test_srgb_transfer();
    test_matrix();
    test_threads_are_deterministic();
    test_end_to_end_quality();
    test_file_roundtrip();
    std::cout << (g_checks - g_failures) << '/' << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
