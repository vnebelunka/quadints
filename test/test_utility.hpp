#ifndef TEST_UTILITY_HPP
#define TEST_UTILITY_HPP

#include <array>
#include <cmath>

struct point2d {
    std::array<double, 2> coords;
    auto operator+(const point2d& other) const -> point2d {
        return {coords[0] + other.coords[0], coords[1] + other.coords[1]};
    }
    auto operator*(double scalar) const -> point2d { return {coords[0] * scalar, coords[1] * scalar}; }
    auto operator+=(const point2d& other) -> point2d& {
        coords[0] += other.coords[0];
        coords[1] += other.coords[1];
        return *this;
    }
    auto operator*=(double scalar) -> point2d& {
        coords[0] *= scalar;
        coords[1] *= scalar;
        return *this;
    }
    auto operator==(const point2d& other) const -> bool {
        return std::abs(coords[0] - other.coords[0]) < 1e-12 && std::abs(coords[1] - other.coords[1]) < 1e-12;
    }
    auto operator!=(const point2d& other) const -> bool { return !(*this == other); }
};

inline auto operator*(double scalar, const point2d& p) -> point2d {
    return {scalar * p.coords[0], scalar * p.coords[1]};
}

inline auto norm(const point2d& p) -> double {
    return std::sqrt((p.coords[0] * p.coords[0]) + (p.coords[1] * p.coords[1]));
}

struct Triangle {
    using point_type = point2d;
    std::array<point2d, 3> _vertices;
    [[nodiscard]] auto mes() const -> double {
        // Compute the area of the triangle using the determinant formula
        const auto& A = _vertices[0].coords;
        const auto& B = _vertices[1].coords;
        const auto& C = _vertices[2].coords;
        return 0.5 * std::abs((A[0] * (B[1] - C[1])) + (B[0] * (C[1] - A[1])) + (C[0] * (A[1] - B[1])));
    }
    [[nodiscard]] constexpr auto vertices() const -> const std::array<point2d, 3>& { return _vertices; }
};

#endif  // TEST_UTILITY_HPP
