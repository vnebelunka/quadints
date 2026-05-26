#ifndef SEGMENT_QUADS_HPP
#define SEGMENT_QUADS_HPP
#include <array>
#include <concepts>

#include "interface.hpp"

namespace quadints {

template <typename T, typename Scalar>
concept segment = requires(T s) {
    typename T::point_type;
    requires banach_vec<typename T::point_type, Scalar>;
    { s.start() } -> std::convertible_to<typename T::point_type>;
    { s.end() } -> std::convertible_to<typename T::point_type>;
};

template <typename Scalar, unsigned int n_points>
struct SegmentQuadrature {
    static_assert("We have no quadrature of such segment quadrature");
};

template <typename Scalar>
struct barycentric_segment {
    using point_type = std::array<Scalar, 1>;
    point_type coords;
    template <typename segment_t>
        requires segment<segment_t, Scalar>
    constexpr auto to_domain(const segment_t& s) const {
        return s.start() + (s.end() - s.start()) * coords[0];
    }
};

template <size_t N_points, typename segment_t, typename Scalar = segment_t::point_type>
constexpr auto to_domain(const segment_t& s, const std::array<barycentric_segment<Scalar>, N_points>& coords) {
    std::array<Scalar, N_points> result;
    auto len = s.end() - s.start();
    for (size_t i = 0; i < N_points; ++i) {
        result[i] = s.start() + len * coords[i].coords[0];
    }
    return result;
}

template <typename Scalar>
struct SegmentQuadrature<Scalar, 1> {
    using point_type = barycentric_segment<Scalar>;
    static constexpr std::array<point_type, 1> points{{0.5}};
    static constexpr std::array<Scalar, 1> weights{1.};
    static constexpr size_t n_points = 1;
};

template <typename Scalar>
struct SegmentQuadrature<Scalar, 2> {
   private:
    static constexpr Scalar sqrt_3 =
        static_cast<Scalar>(1.732050807568877293527446341505872366942805253810380628055806979451933016908798);

   public:
    using point_type = barycentric_segment<Scalar>;
    static constexpr std::size_t n_points = 2;
    static constexpr std::array<point_type, 2> points{point_type{{1. / 2 * (1 - 1. / sqrt_3)}},
                                                      point_type{{1. / 2 * (1 + 1. / sqrt_3)}}};
    static constexpr std::array<Scalar, 2> weights{1. / 2, 1. / 2};
};

template <typename Scalar>
struct SegmentQuadrature<Scalar, 3> {
   private:
    static constexpr Scalar sqrt_3_over_5 =
        static_cast<Scalar>(0.7745966692414833770358530799564799221665843410583181653175147532226966183873991);

   public:
    using point_type = barycentric_segment<Scalar>;
    static constexpr std::size_t n_points = 3;
    static constexpr std::array<point_type, n_points> points{point_type{{1. / 2 - (1. / 2) * sqrt_3_over_5}},
                                                             point_type{{1. / 2 + (1. / 2) * sqrt_3_over_5}},
                                                             point_type{{1. / 2}}};
    static constexpr std::array<Scalar, n_points> weights{5. / 18, 5. / 18, 4. / 9};
};

}  // namespace quadints
#endif  // SEGMENT_QUADS_HPP
