#include "analysis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "stages.hpp"

namespace isp {

double psnr(const Rgb8& a, const Rgb8& b, const std::vector<Rect>& regions) {
    if (a.width != b.width || a.height != b.height) throw std::runtime_error("PSNR needs images of the same size");
    double sq_err = 0.0;
    double samples = 0.0;
    for (const Rect& r : regions) {
        for (int y = std::max(0, r.y0); y < std::min(a.height, r.y1); ++y) {
            for (int x = std::max(0, r.x0); x < std::min(a.width, r.x1); ++x) {
                const size_t i = (static_cast<size_t>(y) * a.width + x) * 3;
                for (int c = 0; c < 3; ++c) {
                    const double d = static_cast<double>(a.data[i + c]) - b.data[i + c];
                    sq_err += d * d;
                }
                samples += 3.0;
            }
        }
    }
    if (samples == 0.0) throw std::runtime_error("PSNR region is empty");
    const double mse = sq_err / samples;
    if (mse == 0.0) return std::numeric_limits<double>::infinity();
    return 10.0 * std::log10(255.0 * 255.0 / mse);
}

double psnr(const Rgb8& a, const Rgb8& b, int border) {
    return psnr(a, b, std::vector<Rect>{{border, border, a.width - border, a.height - border}});
}

Rgb8 mosaic_as_gray(const Plane& mosaic) {
    Rgb8 out;
    out.width = mosaic.width;
    out.height = mosaic.height;
    out.data.resize(mosaic.data.size() * 3);
    for (size_t i = 0; i < mosaic.data.size(); ++i) {
        const uint8_t v = static_cast<uint8_t>(std::lround(srgb_encode(mosaic.data[i]) * 255.0f));
        out.data[3 * i] = out.data[3 * i + 1] = out.data[3 * i + 2] = v;
    }
    return out;
}

Rgb8 mosaic_as_color(const Plane& mosaic, const Cfa& cfa) {
    Rgb8 out;
    out.width = mosaic.width;
    out.height = mosaic.height;
    out.data.assign(mosaic.data.size() * 3, 0);
    for (int y = 0; y < mosaic.height; ++y) {
        for (int x = 0; x < mosaic.width; ++x) {
            const size_t i = static_cast<size_t>(y) * mosaic.width + x;
            const float v = srgb_encode(mosaic.data[i]);
            out.data[3 * i + static_cast<size_t>(cfa.at(x, y))] = static_cast<uint8_t>(std::lround(v * 255.0f));
        }
    }
    return out;
}

Rgb8 linear_to_display(const RgbImage& linear) {
    Rgb8 out;
    out.width = linear.width;
    out.height = linear.height;
    out.data.resize(linear.data.size());
    for (size_t i = 0; i < linear.data.size(); ++i) {
        out.data[i] = static_cast<uint8_t>(std::lround(srgb_encode(linear.data[i]) * 255.0f));
    }
    return out;
}

Rgb8 crop_zoom(const Rgb8& img, const Rect& region, int scale) {
    const Rect r{std::max(0, region.x0), std::max(0, region.y0), std::min(img.width, region.x1),
                 std::min(img.height, region.y1)};
    Rgb8 out;
    out.width = (r.x1 - r.x0) * scale;
    out.height = (r.y1 - r.y0) * scale;
    if (out.width <= 0 || out.height <= 0) throw std::runtime_error("crop region is empty");
    out.data.resize(static_cast<size_t>(out.width) * out.height * 3);
    for (int y = 0; y < out.height; ++y) {
        for (int x = 0; x < out.width; ++x) {
            const size_t src = (static_cast<size_t>(r.y0 + y / scale) * img.width + (r.x0 + x / scale)) * 3;
            const size_t dst = (static_cast<size_t>(y) * out.width + x) * 3;
            for (int c = 0; c < 3; ++c) out.data[dst + c] = img.data[src + c];
        }
    }
    return out;
}

Rgb8 side_by_side(const std::vector<Rgb8>& images, int gap) {
    Rgb8 out;
    for (const Rgb8& im : images) {
        out.width += im.width;
        out.height = std::max(out.height, im.height);
    }
    out.width += gap * static_cast<int>(images.size() > 0 ? images.size() - 1 : 0);
    out.data.assign(static_cast<size_t>(out.width) * out.height * 3, 255);
    int x_off = 0;
    for (const Rgb8& im : images) {
        for (int y = 0; y < im.height; ++y) {
            std::copy_n(&im.data[static_cast<size_t>(y) * im.width * 3], static_cast<size_t>(im.width) * 3,
                        &out.data[(static_cast<size_t>(y) * out.width + x_off) * 3]);
        }
        x_off += im.width + gap;
    }
    return out;
}

}  // namespace isp
