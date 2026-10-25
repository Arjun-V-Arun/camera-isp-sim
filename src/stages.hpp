#pragma once

#include "image.hpp"
#include "raw_meta.hpp"

namespace isp {

enum class DemosaicMethod { Bilinear, Malvar };

// Mirrors an out-of-range index back into [0, n) without repeating the edge pixel
// (the same rule as OpenCV's BORDER_REFLECT_101). It preserves the parity of the index,
// so a reflected Bayer sample always has the same colour as the one it stands in for.
inline int reflect101(int i, int n) {
    if (n == 1) return 0;
    while (i < 0 || i >= n) {
        if (i < 0) i = -i;
        if (i >= n) i = 2 * (n - 1) - i;
    }
    return i;
}

// Copies a plane into a larger one with `border` pixels of reflect-101 padding on every side,
// so neighbourhood stages can index without bounds checks.
Plane pad_reflect101(const Plane& p, int border, int threads);

// Stage 1. Subtracts the sensor black level and scales to [0, 1].
Plane black_level_correct(const RawImage& raw, const RawMeta& meta, int threads);

// Stage 2. Replaces hot or dead pixels. A pixel is defective when it sits more than `threshold`
// above the maximum or below the minimum of its eight same-colour neighbours (the ring two pixels away),
// and is replaced by the median of those neighbours. Returns the number of pixels corrected.
int correct_defective_pixels(Plane& bayer, float threshold, int threads);

// Stage 3. White balance on the mosaic, before demosaicing, as hardware ISPs do.
Gains gray_world_gains(const Plane& bayer, const Cfa& cfa);
void apply_white_balance(Plane& bayer, const Cfa& cfa, const Gains& gains, int threads);

// Stage 4. Reconstructs full RGB at every pixel from the one sample the sensor measured.
RgbImage demosaic(const Plane& bayer, const Cfa& cfa, DemosaicMethod method, int threads);

// Stage 5. Camera RGB to linear sRGB through a 3x3 matrix, then clips to [0, 1].
void apply_color_correction(RgbImage& img, const Mat3& ccm, int threads);

// Stage 6. Linear light to 8-bit sRGB through a precomputed lookup table.
Rgb8 gamma_encode(const RgbImage& img, int threads);

// The sRGB transfer function and its inverse, on values in [0, 1].
float srgb_encode(float linear);
float srgb_decode(float encoded);

}  // namespace isp
