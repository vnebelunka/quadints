#ifndef TEST_UTILITY_HPP
#define TEST_UTILITY_HPP

#include "quads_triangle.hpp"
#include <array>
#include <cmath>

struct point2d {
  std::array<double, 2> coords;
  point2d operator+(const point2d &other) const {
    return {coords[0] + other.coords[0], coords[1] + other.coords[1]};
  }
  point2d operator*(double scalar) const {
    return {coords[0] * scalar, coords[1] * scalar};
  }
  point2d &operator+=(const point2d &other) {
    coords[0] += other.coords[0];
    coords[1] += other.coords[1];
    return *this;
  }
  point2d &operator*=(double scalar) {
    coords[0] *= scalar;
    coords[1] *= scalar;
    return *this;
  }
  bool operator==(const point2d &other) const {
    return std::abs(coords[0] - other.coords[0]) < 1e-12 &&
           std::abs(coords[1] - other.coords[1]) < 1e-12;
  }
  bool operator!=(const point2d &other) const { return !(*this == other); }
};

inline point2d operator*(double scalar, const point2d &p) {
  return {scalar * p.coords[0], scalar * p.coords[1]};
}

inline double norm(const point2d &p) {
  return std::sqrt(p.coords[0] * p.coords[0] + p.coords[1] * p.coords[1]);
}

struct Triangle {
  using point_type = point2d;
  std::array<point2d, 3> _vertices;
  constexpr double mes() const {
    // Compute the area of the triangle using the determinant formula
    const auto &A = _vertices[0].coords;
    const auto &B = _vertices[1].coords;
    const auto &C = _vertices[2].coords;
    return 0.5 * std::abs(A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) +
                          C[0] * (A[1] - B[1]));
  }
  constexpr const std::array<point2d, 3> &vertices() const { return _vertices; }
};

#endif // TEST_UTILITY_HPP
