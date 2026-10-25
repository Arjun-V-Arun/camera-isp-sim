#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace isp {

// Sensor output: one unsigned digital number (DN) per pixel, stored row by row.
struct RawImage {
    int width = 0;
    int height = 0;
    std::vector<uint16_t> data;
};

// Single-channel float image. Holds the Bayer mosaic once black level is removed (values in [0, 1]).
struct Plane {
    int width = 0;
    int height = 0;
    std::vector<float> data;
};

// Three-channel float image, interleaved R, G, B, in linear light.
struct RgbImage {
    int width = 0;
    int height = 0;
    std::vector<float> data;
};

// Three-channel 8-bit image, interleaved R, G, B, display-ready sRGB.
struct Rgb8 {
    int width = 0;
    int height = 0;
    std::vector<uint8_t> data;
};

// Netpbm PGM (P5). 16-bit samples are big-endian, as the format requires.
RawImage read_pgm16(const std::string& path);
void write_pgm16(const std::string& path, const RawImage& img, int maxval);

// Netpbm PPM (P6), 8 bits per channel.
Rgb8 read_ppm(const std::string& path);
void write_ppm(const std::string& path, const Rgb8& img);

// 24-bit uncompressed BMP, which Windows can open without extra software.
void write_bmp(const std::string& path, const Rgb8& img);

// Chooses PPM or BMP from the file extension.
void write_image(const std::string& path, const Rgb8& img);

}  // namespace isp
