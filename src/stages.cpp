#include "stages.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "parallel.hpp"

namespace isp {

Plane pad_reflect101(const Plane& p, int border, int threads) {
    Plane q;
    q.width = p.width + 2 * border;
    q.height = p.height + 2 * border;
    q.data.resize(static_cast<size_t>(q.width) * q.height);

    // Column lookup is the same for every row, so compute it once.
    std::vector<int> src_x(static_cast<size_t>(q.width));
    for (int x = 0; x < q.width; ++x) src_x[static_cast<size_t>(x)] = reflect101(x - border, p.width);

    parallel_rows(q.height, threads, [&](int y0, int y1) {
        for (int y = y0; y < y1; ++y) {
            const float* src = &p.data[static_cast<size_t>(reflect101(y - border, p.height)) * p.width];
            float* dst = &q.data[static_cast<size_t>(y) * q.width];
            for (int x = 0; x < q.width; ++x) dst[x] = src[src_x[static_cast<size_t>(x)]];
        }
    });
    return q;
}

Plane black_level_correct(const RawImage& raw, const RawMeta& meta, int threads) {
    Plane p;
    p.width = raw.width;
    p.height = raw.height;
    p.data.resize(raw.data.size());

    const float black = static_cast<float>(meta.black_level);
    const float scale = 1.0f / static_cast<float>(meta.white_level - meta.black_level);

    parallel_rows(raw.height, threads, [&](int y0, int y1) {
        for (int y = y0; y < y1; ++y) {
            const uint16_t* src = &raw.data[static_cast<size_t>(y) * raw.width];
            float* dst = &p.data[static_cast<size_t>(y) * raw.width];
            for (int x = 0; x < raw.width; ++x) {
                // Read noise can push dark pixels below black; clamp so later stages never see negative light.
                dst[x] = std::clamp((static_cast<float>(src[x]) - black) * scale, 0.0f, 1.0f);
            }
        }
    });
    return p;
}

int correct_defective_pixels(Plane& bayer, float threshold, int threads) {
    // In a Bayer mosaic, the ring of pixels two steps away (horizontally, vertically and diagonally) all share
    // the centre pixel's colour, whatever the CFA layout, so this stage does not need to know the pattern.
    // Using all 8 ring neighbours rather than only the 4 in line matters: fine stripes or a zone plate make
    // the centre an extreme along one direction, but along the stripe there is always a neighbour with a
    // similar value, so real detail is left alone. A stuck pixel stands out in every direction.
    const Plane padded = pad_reflect101(bayer, 2, threads);
    const int pw = padded.width;
    const int offsets[8] = {-2, 2, -2 * pw, 2 * pw, -2 * pw - 2, -2 * pw + 2, 2 * pw - 2, 2 * pw + 2};
    std::atomic<int> corrected{0};

    parallel_rows(bayer.height, threads, [&](int y0, int y1) {
        int local = 0;
        for (int y = y0; y < y1; ++y) {
            const float* row = &padded.data[static_cast<size_t>(y + 2) * pw + 2];
            for (int x = 0; x < bayer.width; ++x) {
                const float* c = row + x;
                float lo = c[offsets[0]], hi = lo;
                for (int k = 1; k < 8; ++k) {
                    lo = std::min(lo, c[offsets[k]]);
                    hi = std::max(hi, c[offsets[k]]);
                }
                if (c[0] > hi + threshold || c[0] < lo - threshold) {
                    float n[8];
                    for (int k = 0; k < 8; ++k) n[k] = c[offsets[k]];
                    std::sort(n, n + 8);
                    bayer.data[static_cast<size_t>(y) * bayer.width + x] = 0.5f * (n[3] + n[4]);  // median of 8
                    ++local;
                }
            }
        }
        corrected += local;  // one atomic add per band, not per pixel
    });
    return corrected.load();
}

Gains gray_world_gains(const Plane& bayer, const Cfa& cfa) {
    // Gray-world assumption: the scene averages to neutral grey, so scale R and B until their means match G.
    double sum[3] = {0.0, 0.0, 0.0};
    double count[3] = {0.0, 0.0, 0.0};
    for (int y = 0; y < bayer.height; ++y) {
        for (int x = 0; x < bayer.width; ++x) {
            const int ch = cfa.at(x, y);
            sum[ch] += bayer.data[static_cast<size_t>(y) * bayer.width + x];
            count[ch] += 1.0;
        }
    }
    const double mean_r = sum[kRed] / count[kRed];
    const double mean_g = sum[kGreen] / count[kGreen];
    const double mean_b = sum[kBlue] / count[kBlue];
    if (mean_r <= 0.0 || mean_b <= 0.0) return Gains{1.0f, 1.0f, 1.0f};
    return Gains{static_cast<float>(mean_g / mean_r), 1.0f, static_cast<float>(mean_g / mean_b)};
}

void apply_white_balance(Plane& bayer, const Cfa& cfa, const Gains& gains, int threads) {
    parallel_rows(bayer.height, threads, [&](int y0, int y1) {
        for (int y = y0; y < y1; ++y) {
            // Each row alternates between just two channels, so look the gains up once per row.
            const float g_even = gains[static_cast<size_t>(cfa.at(0, y))];
            const float g_odd = gains[static_cast<size_t>(cfa.at(1, y))];
            float* row = &bayer.data[static_cast<size_t>(y) * bayer.width];
            for (int x = 0; x < bayer.width; ++x) row[x] *= (x & 1) ? g_odd : g_even;
        }
    });
}

namespace {

// Bilinear: keep the measured channel, and fill each missing channel with the mean of the
// same-colour samples in the 3x3 neighbourhood (4 cross or 4 diagonal neighbours, or 2 in line).
void demosaic_bilinear_rows(const Plane& padded, const Cfa& cfa, RgbImage& out, int y0, int y1) {
    const int pw = padded.width;
    for (int y = y0; y < y1; ++y) {
        for (int x = 0; x < out.width; ++x) {
            const float* c = &padded.data[static_cast<size_t>(y + 2) * pw + (x + 2)];
            const int own = cfa.at(x, y);
            float sum[3] = {0.0f, 0.0f, 0.0f};
            int count[3] = {0, 0, 0};
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    if (dx == 0 && dy == 0) continue;
                    // +2 keeps the coordinate non-negative without changing its parity.
                    const int ch = cfa.at(x + dx + 2, y + dy + 2);
                    sum[ch] += c[dy * pw + dx];
                    ++count[ch];
                }
            }
            float* px = &out.data[(static_cast<size_t>(y) * out.width + x) * 3];
            for (int ch = 0; ch < 3; ++ch) px[ch] = (ch == own) ? c[0] : sum[ch] / static_cast<float>(count[ch]);
        }
    }
}

// Malvar, He and Cutler (2004), "High-quality linear interpolation for demosaicing of Bayer-patterned
// color images". Bilinear interpolation plus a correction from the measured channel's local
// gradient (its Laplacian), which removes most of bilinear's blur and colour fringing at edges.
// Each 5x5 kernel's weights sum to 8, hence the final division by 8.
void demosaic_malvar_rows(const Plane& padded, const Cfa& cfa, RgbImage& out, int y0, int y1) {
    const int pw = padded.width;
    for (int y = y0; y < y1; ++y) {
        const int row_color = cfa.row_color(y);  // the R or B channel that shares this row with green
        for (int x = 0; x < out.width; ++x) {
            const float* c = &padded.data[static_cast<size_t>(y + 2) * pw + (x + 2)];
            auto v = [&](int dx, int dy) { return c[dy * pw + dx]; };
            const float center = c[0];
            const float cross1 = v(-1, 0) + v(1, 0) + v(0, -1) + v(0, 1);
            const float cross2 = v(-2, 0) + v(2, 0) + v(0, -2) + v(0, 2);
            const float diag1 = v(-1, -1) + v(1, -1) + v(-1, 1) + v(1, 1);

            float rgb[3];
            const int own = cfa.at(x, y);
            if (own == kGreen) {
                // Kernel for the colour found left and right of this pixel.
                const float horizontal = (5.0f * center + 4.0f * (v(-1, 0) + v(1, 0)) - (v(-2, 0) + v(2, 0)) - diag1 +
                                          0.5f * (v(0, -2) + v(0, 2))) / 8.0f;
                // Kernel for the colour found above and below (the transpose of the one above).
                const float vertical = (5.0f * center + 4.0f * (v(0, -1) + v(0, 1)) - (v(0, -2) + v(0, 2)) - diag1 +
                                        0.5f * (v(-2, 0) + v(2, 0))) / 8.0f;
                const int col_color = (row_color == kRed) ? kBlue : kRed;
                rgb[kGreen] = center;
                rgb[row_color] = horizontal;
                rgb[col_color] = vertical;
            } else {
                // At a red or blue pixel: green from the cross, the opposite colour from the diagonals.
                const int opposite = (own == kRed) ? kBlue : kRed;
                rgb[own] = center;
                rgb[kGreen] = (4.0f * center + 2.0f * cross1 - cross2) / 8.0f;
                rgb[opposite] = (6.0f * center + 2.0f * diag1 - 1.5f * cross2) / 8.0f;
            }
            float* px = &out.data[(static_cast<size_t>(y) * out.width + x) * 3];
            // The gradient correction can undershoot next to strong edges; light cannot be negative.
            for (int ch = 0; ch < 3; ++ch) px[ch] = std::max(rgb[ch], 0.0f);
        }
    }
}

}  // namespace

RgbImage demosaic(const Plane& bayer, const Cfa& cfa, DemosaicMethod method, int threads) {
    // Padding by 2 (an even number) keeps every sample's colour, so padded (x + 2, y + 2) has the CFA colour of (x, y).
    const Plane padded = pad_reflect101(bayer, 2, threads);
    RgbImage out;
    out.width = bayer.width;
    out.height = bayer.height;
    out.data.resize(static_cast<size_t>(bayer.width) * bayer.height * 3);
    parallel_rows(bayer.height, threads, [&](int y0, int y1) {
        if (method == DemosaicMethod::Bilinear) {
            demosaic_bilinear_rows(padded, cfa, out, y0, y1);
        } else {
            demosaic_malvar_rows(padded, cfa, out, y0, y1);
        }
    });
    return out;
}

void apply_color_correction(RgbImage& img, const Mat3& m, int threads) {
    parallel_rows(img.height, threads, [&](int y0, int y1) {
        for (int y = y0; y < y1; ++y) {
            float* px = &img.data[static_cast<size_t>(y) * img.width * 3];
            for (int x = 0; x < img.width; ++x, px += 3) {
                const float r = px[0], g = px[1], b = px[2];
                px[0] = std::clamp(m[0] * r + m[1] * g + m[2] * b, 0.0f, 1.0f);
                px[1] = std::clamp(m[3] * r + m[4] * g + m[5] * b, 0.0f, 1.0f);
                px[2] = std::clamp(m[6] * r + m[7] * g + m[8] * b, 0.0f, 1.0f);
            }
        }
    });
}

float srgb_encode(float linear) {
    if (linear <= 0.0f) return 0.0f;
    if (linear >= 1.0f) return 1.0f;
    return linear <= 0.0031308f ? 12.92f * linear : 1.055f * std::pow(linear, 1.0f / 2.4f) - 0.055f;
}

float srgb_decode(float encoded) {
    if (encoded <= 0.0f) return 0.0f;
    if (encoded >= 1.0f) return 1.0f;
    return encoded <= 0.04045f ? encoded / 12.92f : std::pow((encoded + 0.055f) / 1.055f, 2.4f);
}

namespace {

// 65536 entries: linear values are quantised to 16 bits before lookup, which is far finer than the 8-bit output.
// A function-local static is initialised once and thread-safely (C++11 guarantees this).
const std::vector<uint8_t>& gamma_lut() {
    static const std::vector<uint8_t> lut = [] {
        std::vector<uint8_t> t(65536);
        for (size_t i = 0; i < t.size(); ++i) {
            const float encoded = srgb_encode(static_cast<float>(i) / 65535.0f);
            t[i] = static_cast<uint8_t>(std::lround(encoded * 255.0f));
        }
        return t;
    }();
    return lut;
}

}  // namespace

Rgb8 gamma_encode(const RgbImage& img, int threads) {
    const std::vector<uint8_t>& lut = gamma_lut();
    Rgb8 out;
    out.width = img.width;
    out.height = img.height;
    out.data.resize(img.data.size());
    parallel_rows(img.height, threads, [&](int y0, int y1) {
        const size_t begin = static_cast<size_t>(y0) * img.width * 3;
        const size_t end = static_cast<size_t>(y1) * img.width * 3;
        for (size_t i = begin; i < end; ++i) {
            const float v = std::clamp(img.data[i], 0.0f, 1.0f);
            out.data[i] = lut[static_cast<size_t>(v * 65535.0f + 0.5f)];
        }
    });
    return out;
}

}  // namespace isp
