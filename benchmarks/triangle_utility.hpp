#ifndef QUADINTS_BENCHMARKS_TRIANGLE_UTILITY_HPP
#define QUADINTS_BENCHMARKS_TRIANGLE_UTILITY_HPP

struct point2d {
    std::array<double, 2> coords;
    point2d operator+(const point2d& other) const { return {coords[0] + other.coords[0], coords[1] + other.coords[1]}; }
    point2d operator*(double scalar) const { return {coords[0] * scalar, coords[1] * scalar}; }
    point2d& operator+=(const point2d& other) {
        coords[0] += other.coords[0];
        coords[1] += other.coords[1];
        return *this;
    }
    point2d& operator*=(double scalar) {
        coords[0] *= scalar;
        coords[1] *= scalar;
        return *this;
    }
    point2d& operator-=(const point2d& other) {
        coords[0] -= other.coords[0];
        coords[1] -= other.coords[1];
        return *this;
    }
    bool operator==(const point2d& other) const {
        return std::abs(coords[0] - other.coords[0]) < 1e-12 && std::abs(coords[1] - other.coords[1]) < 1e-12;
    }
    bool operator!=(const point2d& other) const { return !(*this == other); }
};

inline point2d operator*(double scalar, const point2d& p) { return {scalar * p.coords[0], scalar * p.coords[1]}; }

inline point2d operator-(const point2d& p1, const point2d& p2) {
    return {p1.coords[0] - p2.coords[0], p1.coords[1] - p2.coords[1]};
}

inline double norm(const point2d& p) { return std::sqrt(p.coords[0] * p.coords[0] + p.coords[1] * p.coords[1]); }

struct Triangle {
    using point_type = point2d;
    std::array<point2d, 3> _vertices;
    constexpr double mes() const {
        // Compute the area of the triangle using the determinant formula
        const auto& A = _vertices[0].coords;
        const auto& B = _vertices[1].coords;
        const auto& C = _vertices[2].coords;
        return 0.5 * std::abs(A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]));
    }
    constexpr const std::array<point2d, 3>& vertices() const { return _vertices; }
};

static std::vector<Triangle> generate_random_triangles(size_t count, double min_coord = 0.0, double max_coord = 10.0) {
    std::mt19937 rng(12345);  // Fixed seed for reproducibility
    std::uniform_real_distribution<double> dist(min_coord, max_coord);

    std::vector<Triangle> triangles;
    triangles.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        Triangle T;
        T._vertices[0] = point2d{dist(rng), dist(rng)};
        T._vertices[1] = point2d{dist(rng), dist(rng)};
        T._vertices[2] = point2d{dist(rng), dist(rng)};
        triangles.push_back(T);
    }
    return triangles;
}

#endif  // QUADINTS_BENCHMARKS_TRIANGLE_UTILITY_HPP
