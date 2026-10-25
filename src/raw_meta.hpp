#pragma once

#include <array>
#include <cstddef>
#include <string>

namespace isp {

enum Channel : int { kRed = 0, kGreen = 1, kBlue = 2 };

using Mat3 = std::array<float, 9>;   // row-major 3x3 matrix
using Gains = std::array<float, 3>;  // per-channel multipliers, R, G, B

// Colour filter array: which channel each pixel of the repeating 2x2 tile sees.
// Index into the tile is (x & 1) + 2 * (y & 1). Callers pass non-negative coordinates.
struct Cfa {
    std::array<int, 4> tile{kRed, kGreen, kGreen, kBlue};  // RGGB by default

    int at(int x, int y) const { return tile[static_cast<size_t>(((y & 1) << 1) | (x & 1))]; }

    // The non-green channel that appears in row y (R or B).
    int row_color(int y) const {
        const int a = at(0, y);
        return a != kGreen ? a : at(1, y);
    }

    static Cfa parse(const std::string& name);  // "RGGB", "BGGR", "GRBG" or "GBRG"
    std::string name() const;
};

// Everything the pipeline needs to know about a raw frame. Stored next to the raw file as key=value text,
// playing the role that EXIF/DNG tags play for a real camera file.
struct RawMeta {
    int width = 0;
    int height = 0;
    Cfa cfa;
    int black_level = 256;   // DN the sensor reports for zero light
    int white_level = 4095;  // saturation DN (12-bit sensor)
    Gains wb_gains{2.0f, 1.0f, 1.6f};
    // Camera RGB to linear sRGB. Each row sums to 1 so neutral greys stay neutral.
    Mat3 ccm{1.6f, -0.4f, -0.2f,
             -0.3f, 1.5f, -0.2f,
             0.0f, -0.5f, 1.5f};
};

RawMeta read_meta(const std::string& path);
void write_meta(const std::string& path, const RawMeta& meta);

Mat3 mat3_inverse(const Mat3& m);
Mat3 mat3_multiply(const Mat3& a, const Mat3& b);

}  // namespace isp
