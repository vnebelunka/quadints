#ifndef TRIANGLE_QUADS_HPP
#define TRIANGLE_QUADS_HPP

#include <array>
#include <cassert>
#include <cmath>

#include "interface.hpp"

namespace quadints {
template <typename T, typename Scalar>
concept triangle = requires(T t) {
    typename T::point_type;
    requires banach_vec<Scalar, typename T::point_type>;
    t.vertices();
    t.vertices()[0], t.vertices()[1], t.vertices()[2];
    requires t.vertices().size() == 3;
};

template <typename Scalar>
struct barycentric_triangle {
    std::array<Scalar, 2> coords;
    constexpr barycentric_triangle(Scalar x, Scalar y) : coords{x, y} {}
    constexpr barycentric_triangle(Scalar x, Scalar y, Scalar z) : coords{x, y} {
        assert(std::abs(x + y + z - 1) < 1e-10);
    }
    constexpr barycentric_triangle() : coords{0, 0} {}
    constexpr Scalar x() const { return coords[0]; }
    constexpr Scalar y() const { return coords[1]; }
    constexpr Scalar z() const { return 1 - coords[0] - coords[1]; }

    template <typename triangle_t>
        requires triangle<triangle_t, Scalar>
    constexpr auto to_domain(const triangle_t& t) const {
        const auto& vertices = t.vertices();
        const auto& A = vertices[0];
        const auto& B = vertices[1];
        const auto& C = vertices[2];
        typename triangle_t::point_type D = A * coords[0] + B * coords[1] + C * (1 - coords[0] - coords[1]);
        return D;
    }
};
template <size_t N_POINTS, typename triangle_t, typename Scalar = triangle_t::point_type>
constexpr auto to_domain(const triangle_t& t, const std::array<barycentric_triangle<Scalar>, N_POINTS>& coords) {
    const auto& vertices = t.vertices();
    const auto& A = vertices[0];
    const auto& B = vertices[1];
    const auto& C = vertices[2];
    const auto CA = A + C * -1;
    const auto CB = B + C * -1;
    std::array<typename triangle_t::point_type, N_POINTS> result;
    for (size_t i = 0; i < N_POINTS; ++i) {
        result[i] = C + CA * coords[i].x() + CB * coords[i].y();
    }
    return result;
}

template <typename Scalar>
struct barycentric_direction {
    std::array<Scalar, 2> dirs;
    constexpr barycentric_direction(Scalar dx, Scalar dy) : dirs{dx, dy} {}
    constexpr barycentric_direction(Scalar dx, Scalar dy, Scalar dz) : dirs{dx, dy} {
        assert(std::abs(dx + dy + dz) < 1e-10);
    }
};

template <typename Scalar>
constexpr barycentric_triangle<Scalar> operator+(const barycentric_triangle<Scalar>& c,
                                                 const barycentric_direction<Scalar>& d) {
    return {c.coords[0] + d.dirs[0], c.coords[1] + d.dirs[1]};
}
template <typename Scalar>
constexpr barycentric_triangle<Scalar> operator+(const barycentric_direction<Scalar>& d,
                                                 const barycentric_triangle<Scalar>& c) {
    return c + d;
}

template <typename Scalar>
constexpr barycentric_triangle<Scalar>& operator+=(barycentric_triangle<Scalar>& c,
                                                   const barycentric_direction<Scalar>& d) {
    c.coords[0] += d.dirs[0];
    c.coords[1] += d.dirs[1];
    return c;
}

template <typename Scalar, unsigned int n_points>
struct TriangleQuadrature {
    static_assert("We have no quadrature of such dimension and order");
};

template <typename Scalar>
struct TriangleQuadrature<Scalar, 1> {
    static constexpr size_t n_points = 1;
    using point_type = barycentric_triangle<Scalar>;
    static constexpr std::array<point_type, 1> points{{{1. / 3, 1. / 3}}};
    static constexpr std::array<Scalar, 1> weights{1.};
};

template <typename Scalar>
struct TriangleQuadrature<Scalar, 3> {
    static constexpr size_t n_points = 3;
    using point_type = barycentric_triangle<Scalar>;
    static constexpr std::array<point_type, 3> points{point_type{0.5, 0.}, point_type{0.5, 0.5}, point_type{0., 0.5}};
    static constexpr std::array<Scalar, 3> weights{1. / 3, 1. / 3, 1. / 3};
};

template <typename Scalar>
struct TriangleQuadrature<Scalar, 7> {
    static constexpr size_t n_points = 7;
    using point_type = barycentric_triangle<Scalar>;

    // Precompute constants using std::sqrt (requires C++26 or a constexpr sqrt
    // implementation)
    static constexpr Scalar sqrt15 =
        static_cast<Scalar>(3.872983346207416885179265399782399610832921705291590826587573766113483091936996);
    static constexpr Scalar l1 = (6 + sqrt15) / 21;
    static constexpr Scalar l2 = (6 - sqrt15) / 21;
    static constexpr Scalar w1 = (155 + sqrt15) / 1200;
    static constexpr Scalar w2 = (155 - sqrt15) / 1200;

    static constexpr std::array<point_type, 7> points{point_type{Scalar(1.0 / 3.0), Scalar(1.0 / 3.0)},  // centroid
                                                      point_type{l1, l1},
                                                      point_type{l1, 1 - 2 * l1},
                                                      point_type{1 - 2 * l1, l1},
                                                      point_type{l2, l2},
                                                      point_type{l2, 1 - 2 * l2},
                                                      point_type{1 - 2 * l2, l2}};

    static constexpr std::array<Scalar, 7> weights{Scalar(0.225),  // = 9/40
                                                   w1,
                                                   w1,
                                                   w1,
                                                   w2,
                                                   w2,
                                                   w2};
};
}  // namespace quadints
#endif  // TRIANGLE_QUADS_HPP
