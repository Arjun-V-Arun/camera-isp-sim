#pragma once

#include <vector>

#include "image.hpp"
#include "raw_meta.hpp"

namespace isp {

// Half-open pixel rectangle [x0, x1) x [y0, y1).
struct Rect {
    int x0 = 0;
    int y0 = 0;
    int x1 = 0;
    int y1 = 0;
};

// Peak signal-to-noise ratio in dB between two 8-bit images, over the given regions combined.
// Returns infinity when the images are identical there.
double psnr(const Rgb8& a, const Rgb8& b, const std::vector<Rect>& regions);

// PSNR over the whole image except a `border`-pixel frame, the usual convention for demosaicing results.
double psnr(const Rgb8& a, const Rgb8& b, int border);

// Images for looking at intermediate stages.
Rgb8 mosaic_as_gray(const Plane& mosaic);                   // raw mosaic as a grey image
Rgb8 mosaic_as_color(const Plane& mosaic, const Cfa& cfa);  // each pixel shown in its filter colour
Rgb8 linear_to_display(const RgbImage& linear);             // gamma only, no colour correction

// Cuts out a region and enlarges it with nearest-neighbour scaling so individual pixels stay visible.
Rgb8 crop_zoom(const Rgb8& img, const Rect& region, int scale);

// Places images side by side with a gap between them.
Rgb8 side_by_side(const std::vector<Rgb8>& images, int gap);

}  // namespace isp
