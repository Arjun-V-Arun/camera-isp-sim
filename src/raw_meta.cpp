#include "raw_meta.hpp"

#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace isp {
namespace {

std::string trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

template <size_t N>
std::array<float, N> parse_floats(const std::string& key, const std::string& value) {
    std::array<float, N> out{};
    std::istringstream ss(value);
    for (size_t i = 0; i < N; ++i) {
        if (!(ss >> out[i])) throw std::runtime_error("metadata: '" + key + "' needs " + std::to_string(N) + " numbers");
    }
    return out;
}

}  // namespace

Cfa Cfa::parse(const std::string& name) {
    if (name.size() != 4) throw std::runtime_error("CFA pattern must have 4 letters, got '" + name + "'");
    Cfa cfa;
    int greens = 0;
    for (size_t i = 0; i < 4; ++i) {
        switch (name[i]) {
            case 'R': cfa.tile[i] = kRed; break;
            case 'G': cfa.tile[i] = kGreen; ++greens; break;
            case 'B': cfa.tile[i] = kBlue; break;
            default: throw std::runtime_error("CFA pattern may only contain R, G, B: '" + name + "'");
        }
    }
    // A Bayer tile has its two greens on a diagonal and one red and one blue.
    const bool diagonal_greens = (cfa.tile[0] == kGreen && cfa.tile[3] == kGreen) ||
                                 (cfa.tile[1] == kGreen && cfa.tile[2] == kGreen);
    if (greens != 2 || !diagonal_greens || name.find('R') == std::string::npos || name.find('B') == std::string::npos) {
        throw std::runtime_error("not a Bayer pattern: '" + name + "' (use RGGB, BGGR, GRBG or GBRG)");
    }
    return cfa;
}

std::string Cfa::name() const {
    static const char letters[3] = {'R', 'G', 'B'};
    std::string s;
    for (int c : tile) s.push_back(letters[c]);
    return s;
}

RawMeta read_meta(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        const size_t eq = line.find('=');
        if (eq == std::string::npos) throw std::runtime_error(path + ": expected key=value, got '" + line + "'");
        kv[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
    }
    auto need = [&](const std::string& key) -> const std::string& {
        auto it = kv.find(key);
        if (it == kv.end()) throw std::runtime_error(path + ": missing key '" + key + "'");
        return it->second;
    };

    RawMeta m;
    m.width = std::stoi(need("width"));
    m.height = std::stoi(need("height"));
    m.cfa = Cfa::parse(need("pattern"));
    m.black_level = std::stoi(need("black_level"));
    m.white_level = std::stoi(need("white_level"));
    m.wb_gains = parse_floats<3>("wb_gains", need("wb_gains"));
    m.ccm = parse_floats<9>("ccm", need("ccm"));
    if (m.white_level <= m.black_level) throw std::runtime_error(path + ": white_level must exceed black_level");
    return m;
}

void write_meta(const std::string& path, const RawMeta& m) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << "# Raw frame metadata for camera-isp-sim\n";
    out << "width=" << m.width << '\n';
    out << "height=" << m.height << '\n';
    out << "pattern=" << m.cfa.name() << '\n';
    out << "black_level=" << m.black_level << '\n';
    out << "white_level=" << m.white_level << '\n';
    out << "wb_gains=" << m.wb_gains[0] << ' ' << m.wb_gains[1] << ' ' << m.wb_gains[2] << '\n';
    out << "ccm=";
    for (size_t i = 0; i < 9; ++i) out << m.ccm[i] << (i + 1 < 9 ? " " : "\n");
}

Mat3 mat3_inverse(const Mat3& m) {
    // Adjugate divided by determinant, computed in double for accuracy.
    const double a = m[0], b = m[1], c = m[2];
    const double d = m[3], e = m[4], f = m[5];
    const double g = m[6], h = m[7], i = m[8];
    const double co00 = e * i - f * h, co01 = -(d * i - f * g), co02 = d * h - e * g;
    const double det = a * co00 + b * co01 + c * co02;
    if (std::fabs(det) < 1e-12) throw std::runtime_error("colour matrix is singular");
    const double inv = 1.0 / det;
    return Mat3{
        static_cast<float>(co00 * inv), static_cast<float>(-(b * i - c * h) * inv), static_cast<float>((b * f - c * e) * inv),
        static_cast<float>(co01 * inv), static_cast<float>((a * i - c * g) * inv), static_cast<float>(-(a * f - c * d) * inv),
        static_cast<float>(co02 * inv), static_cast<float>(-(a * h - b * g) * inv), static_cast<float>((a * e - b * d) * inv)};
}

Mat3 mat3_multiply(const Mat3& a, const Mat3& b) {
    Mat3 r{};
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k) sum += static_cast<double>(a[row * 3 + k]) * b[k * 3 + col];
            r[row * 3 + col] = static_cast<float>(sum);
        }
    }
    return r;
}

}  // namespace isp
