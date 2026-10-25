#include "synth.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "stages.hpp"

namespace isp {
namespace {

constexpr double kPi = 3.14159265358979323846;

// SplitMix64 generator with Box-Muller normals. Written out by hand so the noise is identical on every compiler
// (std::normal_distribution is implementation-defined, so MSVC and GCC would generate different noise).
class Rng {
public:
    explicit Rng(uint64_t seed) : state_(seed) {}

    uint64_t next() {
        uint64_t z = (state_ += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }

    double uniform() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }  // [0, 1)

    double normal() {
        const double u1 = 1.0 - uniform();  // (0, 1], safe for log
        const double u2 = uniform();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * kPi * u2);
    }

private:
    uint64_t state_;
};

void fill(Rgb8& img, const Rect& r, uint8_t red, uint8_t green, uint8_t blue) {
    for (int y = r.y0; y < r.y1; ++y) {
        for (int x = r.x0; x < r.x1; ++x) {
            uint8_t* px = &img.data[(static_cast<size_t>(y) * img.width + x) * 3];
            px[0] = red;
            px[1] = green;
            px[2] = blue;
        }
    }
}

void set_gray(Rgb8& img, int x, int y, double v01) {
    const uint8_t v = static_cast<uint8_t>(std::lround(std::clamp(v01, 0.0, 1.0) * 255.0));
    uint8_t* px = &img.data[(static_cast<size_t>(y) * img.width + x) * 3];
    px[0] = px[1] = px[2] = v;
}

// HSV to RGB with h in [0, 360) and s, v in [0, 1].
std::array<double, 3> hsv_to_rgb(double h, double s, double v) {
    const double c = v * s;
    const double hp = h / 60.0;
    const double x = c * (1.0 - std::fabs(std::fmod(hp, 2.0) - 1.0));
    double r = 0, g = 0, b = 0;
    if (hp < 1) { r = c; g = x; }
    else if (hp < 2) { r = x; g = c; }
    else if (hp < 3) { g = c; b = x; }
    else if (hp < 4) { g = x; b = c; }
    else if (hp < 5) { r = x; b = c; }
    else { r = c; b = x; }
    const double m = v - c;
    return {r + m, g + m, b + m};
}

// Approximate sRGB values of the 24 X-Rite ColorChecker Classic patches, in chart order.
constexpr uint8_t kPatches[24][3] = {
    {115, 82, 68},   {194, 150, 130}, {98, 122, 157},  {87, 108, 67},   {133, 128, 177}, {103, 189, 170},
    {214, 126, 44},  {80, 91, 166},   {193, 90, 99},   {94, 60, 108},   {157, 188, 64},  {224, 163, 46},
    {56, 61, 150},   {70, 148, 73},   {175, 54, 60},   {231, 199, 31},  {187, 86, 149},  {8, 133, 161},
    {243, 243, 242}, {200, 200, 200}, {160, 160, 160}, {122, 122, 121}, {85, 85, 85},    {52, 52, 52}};

}  // namespace

SceneLayout scene_layout(int width, int height) {
    const int hw = width / 2;
    const int hh = height / 2;
    return SceneLayout{{0, 0, hw, hh}, {hw, 0, width, hh}, {0, hh, hw, height}, {hw, hh, width, height}};
}

std::vector<Rect> color_patch_rects(int width, int height, int inset) {
    const Rect q = scene_layout(width, height).color_patches;
    const int qw = q.x1 - q.x0, qh = q.y1 - q.y0;
    const int gap = std::max(2, qw / 80);
    const int patch = std::min((qw - 7 * gap) / 6, (qh - 5 * gap) / 4);
    const int ox = q.x0 + (qw - (6 * patch + 5 * gap)) / 2;
    const int oy = q.y0 + (qh - (4 * patch + 3 * gap)) / 2;
    std::vector<Rect> rects;
    for (int i = 0; i < 24; ++i) {
        const int x0 = ox + (i % 6) * (patch + gap);
        const int y0 = oy + (i / 6) * (patch + gap);
        rects.push_back({x0 + inset, y0 + inset, x0 + patch - inset, y0 + patch - inset});
    }
    return rects;
}

Rgb8 make_test_scene(int width, int height) {
    if (width < 64 || height < 64 || width % 2 != 0 || height % 2 != 0) {
        throw std::runtime_error("scene width and height must be even and at least 64");
    }
    Rgb8 img;
    img.width = width;
    img.height = height;
    img.data.assign(static_cast<size_t>(width) * height * 3, 0);
    const SceneLayout L = scene_layout(width, height);

    // Colour patches: a 6 x 4 grid on a dark grey surround.
    {
        fill(img, L.color_patches, 40, 40, 40);
        const std::vector<Rect> patches = color_patch_rects(width, height);
        for (size_t i = 0; i < patches.size(); ++i) fill(img, patches[i], kPatches[i][0], kPatches[i][1], kPatches[i][2]);
    }

    // Zone plate: cos(pi r^2 / (2 R)) has local frequency r / (2R) cycles per pixel, reaching 0.5 (Nyquist) at r = R.
    {
        const Rect& q = L.zone_plate;
        const double cx = (q.x0 + q.x1 - 1) / 2.0, cy = (q.y0 + q.y1 - 1) / 2.0;
        const double radius = std::min(q.x1 - q.x0, q.y1 - q.y0) / 2.0 - 2.0;
        for (int y = q.y0; y < q.y1; ++y) {
            for (int x = q.x0; x < q.x1; ++x) {
                const double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
                const double v = r2 <= radius * radius ? 0.5 + 0.5 * std::cos(kPi * r2 / (2.0 * radius)) : 0.5;
                set_gray(img, x, y, 0.08 + 0.84 * v);
            }
        }
    }

    // Ramps: a grey ramp on top, and below it every hue left to right with saturation falling top to bottom.
    {
        const Rect& q = L.ramps;
        const int qw = q.x1 - q.x0;
        const int mid = q.y0 + (q.y1 - q.y0) / 3;
        for (int y = q.y0; y < q.y1; ++y) {
            for (int x = q.x0; x < q.x1; ++x) {
                const double t = static_cast<double>(x - q.x0) / (qw - 1);
                if (y < mid) {
                    set_gray(img, x, y, t);
                } else {
                    const double s = 1.0 - 0.85 * (y - mid) / static_cast<double>(q.y1 - mid - 1);
                    const auto rgb = hsv_to_rgb(std::min(t * 360.0, 359.999), s, 0.9);
                    uint8_t* px = &img.data[(static_cast<size_t>(y) * width + x) * 3];
                    for (int c = 0; c < 3; ++c) px[c] = static_cast<uint8_t>(std::lround(rgb[static_cast<size_t>(c)] * 255.0));
                }
            }
        }
    }

    // Siemens star: 36 alternating sectors, so spatial frequency rises without limit towards the centre.
    {
        const Rect& q = L.siemens_star;
        const double cx = (q.x0 + q.x1 - 1) / 2.0, cy = (q.y0 + q.y1 - 1) / 2.0;
        const double radius = std::min(q.x1 - q.x0, q.y1 - q.y0) / 2.0 - 8.0;
        const double sector = 2.0 * kPi / 36.0;
        for (int y = q.y0; y < q.y1; ++y) {
            for (int x = q.x0; x < q.x1; ++x) {
                const double dx = x - cx, dy = y - cy;
                const double r = std::sqrt(dx * dx + dy * dy);
                double v = 0.5;
                if (r < radius && r > 4.0) {
                    const int k = static_cast<int>(std::floor((std::atan2(dy, dx) + kPi) / sector));
                    v = (k % 2 == 0) ? 0.1 : 0.9;
                }
                set_gray(img, x, y, v);
            }
        }
    }
    return img;
}

RawGenResult make_raw(const Rgb8& scene, const RawGenConfig& cfg) {
    RawGenResult out;
    out.meta = cfg.meta;
    out.meta.width = scene.width;
    out.meta.height = scene.height;
    const RawMeta& m = out.meta;

    const Mat3 inv_ccm = mat3_inverse(m.ccm);
    const float range = static_cast<float>(m.white_level - m.black_level);

    // Inverse gamma for every 8-bit code value.
    std::array<float, 256> to_linear{};
    for (int i = 0; i < 256; ++i) to_linear[static_cast<size_t>(i)] = srgb_decode(static_cast<float>(i) / 255.0f);

    Rng rng(cfg.seed);
    out.raw.width = scene.width;
    out.raw.height = scene.height;
    out.raw.data.resize(static_cast<size_t>(scene.width) * scene.height);

    for (int y = 0; y < scene.height; ++y) {
        for (int x = 0; x < scene.width; ++x) {
            const size_t i = static_cast<size_t>(y) * scene.width + x;
            const uint8_t* px = &scene.data[i * 3];
            const int ch = m.cfa.at(x, y);
            // Only the channel this pixel's filter passes is needed: one row of the inverse colour matrix.
            float cam = 0.0f;
            for (int k = 0; k < 3; ++k) cam += inv_ccm[static_cast<size_t>(ch * 3 + k)] * to_linear[px[k]];
            // Out-of-gamut colours can map to negative camera values; a real sensor cannot record those.
            cam = std::max(cam, 0.0f);
            const float sensor = std::min(cam / m.wb_gains[static_cast<size_t>(ch)], 1.0f);  // undo white balance
            double dn_signal = static_cast<double>(sensor) * range;
            if (cfg.noise_scale > 0.0f) {
                const double var = static_cast<double>(cfg.read_noise_dn) * cfg.read_noise_dn +
                                   static_cast<double>(cfg.shot_noise_dn_per_dn) * dn_signal;
                dn_signal += cfg.noise_scale * std::sqrt(var) * rng.normal();
            }
            const double dn = std::round(m.black_level + dn_signal);
            out.raw.data[i] = static_cast<uint16_t>(std::clamp(dn, 0.0, static_cast<double>(m.white_level)));
        }
    }

    // Hot pixels: stuck at the saturation value regardless of the scene.
    const size_t n = out.raw.data.size();
    const size_t wanted = static_cast<size_t>(std::llround(static_cast<double>(cfg.hot_pixel_fraction) * n));
    std::vector<bool> hot(n, false);
    for (size_t k = 0; k < wanted; ++k) {
        const size_t idx = static_cast<size_t>(rng.next() % n);
        if (hot[idx]) continue;
        hot[idx] = true;
        out.raw.data[idx] = static_cast<uint16_t>(m.white_level);
        ++out.hot_pixels;
    }
    return out;
}

}  // namespace isp
