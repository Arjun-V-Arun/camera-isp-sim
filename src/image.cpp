#include "image.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>

namespace isp {
namespace {

std::ifstream open_input(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + path);
    return in;
}

std::ofstream open_output(const std::string& path) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + path);
    return out;
}

// Reads one whitespace-separated header token, skipping '#' comments.
// When it returns, exactly one whitespace byte after the token has been consumed,
// which is where Netpbm binary data starts after the final header field.
std::string next_token(std::istream& in) {
    std::string token;
    char c;
    while (in.get(c)) {
        if (token.empty() && c == '#') {
            std::string comment;
            std::getline(in, comment);
            continue;
        }
        if (std::isspace(static_cast<unsigned char>(c))) {
            if (!token.empty()) break;
            continue;
        }
        token.push_back(c);
    }
    return token;
}

struct NetpbmHeader {
    std::string magic;
    int width = 0;
    int height = 0;
    int maxval = 0;
};

NetpbmHeader read_header(std::istream& in, const std::string& path) {
    NetpbmHeader h;
    h.magic = next_token(in);
    try {
        h.width = std::stoi(next_token(in));
        h.height = std::stoi(next_token(in));
        h.maxval = std::stoi(next_token(in));
    } catch (const std::exception&) {
        throw std::runtime_error(path + ": malformed Netpbm header");
    }
    if (h.width <= 0 || h.height <= 0 || h.maxval <= 0 || h.maxval > 65535) {
        throw std::runtime_error(path + ": invalid Netpbm dimensions or maxval");
    }
    return h;
}

void put_u16(std::ostream& out, uint16_t v) {
    out.put(static_cast<char>(v & 0xFF));
    out.put(static_cast<char>(v >> 8));
}

void put_u32(std::ostream& out, uint32_t v) {
    for (int i = 0; i < 4; ++i) out.put(static_cast<char>((v >> (8 * i)) & 0xFF));
}

}  // namespace

RawImage read_pgm16(const std::string& path) {
    std::ifstream in = open_input(path);
    const NetpbmHeader h = read_header(in, path);
    if (h.magic != "P5") throw std::runtime_error(path + ": not a binary PGM (P5) file");

    RawImage img;
    img.width = h.width;
    img.height = h.height;
    const size_t n = static_cast<size_t>(h.width) * h.height;
    img.data.resize(n);

    if (h.maxval < 256) {
        std::vector<unsigned char> bytes(n);
        in.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(n));
        if (static_cast<size_t>(in.gcount()) != n) throw std::runtime_error(path + ": truncated pixel data");
        for (size_t i = 0; i < n; ++i) img.data[i] = bytes[i];
    } else {
        std::vector<unsigned char> bytes(2 * n);
        in.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(2 * n));
        if (static_cast<size_t>(in.gcount()) != 2 * n) throw std::runtime_error(path + ": truncated pixel data");
        for (size_t i = 0; i < n; ++i) {
            img.data[i] = static_cast<uint16_t>((bytes[2 * i] << 8) | bytes[2 * i + 1]);
        }
    }
    return img;
}

void write_pgm16(const std::string& path, const RawImage& img, int maxval) {
    std::ofstream out = open_output(path);
    out << "P5\n" << img.width << ' ' << img.height << '\n' << maxval << '\n';
    const size_t n = static_cast<size_t>(img.width) * img.height;
    if (maxval < 256) {
        for (size_t i = 0; i < n; ++i) {
            out.put(static_cast<char>(std::min<int>(img.data[i], maxval)));
        }
    } else {
        for (size_t i = 0; i < n; ++i) {
            const uint16_t v = static_cast<uint16_t>(std::min<int>(img.data[i], maxval));
            out.put(static_cast<char>(v >> 8));
            out.put(static_cast<char>(v & 0xFF));
        }
    }
}

Rgb8 read_ppm(const std::string& path) {
    std::ifstream in = open_input(path);
    const NetpbmHeader h = read_header(in, path);
    if (h.magic != "P6") throw std::runtime_error(path + ": not a binary PPM (P6) file");
    if (h.maxval != 255) throw std::runtime_error(path + ": only 8-bit PPM (maxval 255) is supported");

    Rgb8 img;
    img.width = h.width;
    img.height = h.height;
    const size_t n = static_cast<size_t>(h.width) * h.height * 3;
    img.data.resize(n);
    in.read(reinterpret_cast<char*>(img.data.data()), static_cast<std::streamsize>(n));
    if (static_cast<size_t>(in.gcount()) != n) throw std::runtime_error(path + ": truncated pixel data");
    return img;
}

void write_ppm(const std::string& path, const Rgb8& img) {
    std::ofstream out = open_output(path);
    out << "P6\n" << img.width << ' ' << img.height << "\n255\n";
    out.write(reinterpret_cast<const char*>(img.data.data()), static_cast<std::streamsize>(img.data.size()));
}

void write_bmp(const std::string& path, const Rgb8& img) {
    std::ofstream out = open_output(path);
    const uint32_t row_bytes = static_cast<uint32_t>(img.width) * 3;
    const uint32_t padding = (4 - row_bytes % 4) % 4;
    const uint32_t data_size = (row_bytes + padding) * static_cast<uint32_t>(img.height);

    // BITMAPFILEHEADER (14 bytes)
    out.put('B');
    out.put('M');
    put_u32(out, 54 + data_size);
    put_u32(out, 0);
    put_u32(out, 54);
    // BITMAPINFOHEADER (40 bytes); a positive height means rows are stored bottom-up.
    put_u32(out, 40);
    put_u32(out, static_cast<uint32_t>(img.width));
    put_u32(out, static_cast<uint32_t>(img.height));
    put_u16(out, 1);
    put_u16(out, 24);
    put_u32(out, 0);
    put_u32(out, data_size);
    put_u32(out, 2835);
    put_u32(out, 2835);
    put_u32(out, 0);
    put_u32(out, 0);

    const char zeros[3] = {0, 0, 0};
    for (int y = img.height - 1; y >= 0; --y) {
        const uint8_t* row = &img.data[static_cast<size_t>(y) * row_bytes];
        for (int x = 0; x < img.width; ++x) {
            out.put(static_cast<char>(row[3 * x + 2]));
            out.put(static_cast<char>(row[3 * x + 1]));
            out.put(static_cast<char>(row[3 * x + 0]));
        }
        out.write(zeros, padding);
    }
}

void write_image(const std::string& path, const Rgb8& img) {
    std::string ext = path.size() >= 4 ? path.substr(path.size() - 4) : "";
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (ext == ".bmp") {
        write_bmp(path, img);
    } else if (ext == ".ppm") {
        write_ppm(path, img);
    } else {
        throw std::runtime_error(path + ": output must end in .bmp or .ppm");
    }
}

}  // namespace isp
