#pragma once

#include <vector>

#include "analysis.hpp"
#include "image.hpp"
#include "raw_meta.hpp"

namespace isp {

// Where each part of the generated test scene sits, so results can be measured per region.
struct SceneLayout {
    Rect color_patches;  // top left: 24 ColorChecker-style patches
    Rect zone_plate;     // top right: circular zone plate, frequency rising to Nyquist at its edge
    Rect ramps;          // bottom left: grey ramp above a hue and saturation sweep
    Rect siemens_star;   // bottom right: 36-sector Siemens star
};

SceneLayout scene_layout(int width, int height);

// The 24 colour patches, in chart order. `inset` shrinks each rectangle on every side,
// which is useful for measuring patch interiors away from edges.
std::vector<Rect> color_patch_rects(int width, int height, int inset = 0);

// An 8-bit sRGB test chart built to stress an ISP: flat colours for colour accuracy, and a zone plate and
// Siemens star whose fine detail exposes demosaicing artifacts. Width and height must be even and at least 64.
Rgb8 make_test_scene(int width, int height);

struct RawGenConfig {
    RawMeta meta;                     // pattern, levels, gains and matrix to simulate; width and height come from the scene
    float noise_scale = 1.0f;         // 0 gives a noise-free raw; 1 is the default sensor noise
    float read_noise_dn = 2.0f;       // standard deviation of read noise, in DN
    float shot_noise_dn_per_dn = 0.25f;  // shot-noise variance per DN of signal (a Gaussian approximation of Poisson noise)
    float hot_pixel_fraction = 0.0f;  // fraction of pixels forced to saturation
    uint64_t seed = 7;
};

struct RawGenResult {
    RawImage raw;
    RawMeta meta;
    int hot_pixels = 0;  // number of distinct pixels made hot
};

// Runs the ISP backwards to turn a ground-truth sRGB image into a plausible sensor raw:
// undo gamma, undo the colour matrix, undo white balance, sample through the Bayer filter,
// add black level and noise, quantise to the sensor bit depth, then inject hot pixels.
// Because the original image is known, the forward pipeline's output can be scored against it.
RawGenResult make_raw(const Rgb8& scene, const RawGenConfig& cfg);

}  // namespace isp
