#ifndef SEGMENT_QUADS_HPP
#define SEGMENT_QUADS_HPP
#include "interface.hpp"
#include <array>
#include <concepts>

template <typename T, typename Scalar>
concept segment = requires(T s) {
  typename T::point_type;
  requires banach_vec<typename T::point_type, Scalar>;
  { s.start() } -> std::convertible_to<typename T::point_type>;
  { s.end() } -> std::convertible_to<typename T::point_type>;
};

template <typename Scalar, unsigned int n_points> struct SegmentQuadrature {
  static_assert("We have no quadrature of such segment quadrature");
};

template <typename Scalar> struct barycentric_segment {
  std::array<Scalar, 1> coords;
};

template <typename Scalar, typename segment_t>
  requires segment<segment_t, Scalar>
constexpr auto
from_reference_domain(const barycentric_segment<Scalar> &reference_point,
                      const segment_t &s) -> segment_t::point_type {
  return s.start() + (s.end() - s.start()) * reference_point.coords[0];
}

template <typename Scalar> struct SegmentQuadrature<Scalar, 1> {
  static constexpr std::array<barycentric_segment<Scalar>, 1> points{{0.5}};
  static constexpr std::array<Scalar, 1> weights{1.};
  static constexpr int n_points = 1;
};

template <typename Scalar> struct SegmentQuadrature<Scalar, 2> {
private:
  static constexpr Scalar sqrt_3 = static_cast<Scalar>(
      1.732050807568877293527446341505872366942805253810380628055806979451933016908798);

public:
  static constexpr std::array<barycentric_segment<Scalar>, 2> points{
      barycentric_segment<Scalar>{{1. / 2 * (1 - 1. / sqrt_3)}},
      barycentric_segment<Scalar>{{1. / 2 * (1 + 1. / sqrt_3)}}};
  static constexpr std::array<Scalar, 2> weights{1. / 2, 1. / 2};
  static constexpr std::size_t n_points = 2;
};

template <typename Scalar> struct SegmentQuadrature<Scalar, 3> {
private:
  static constexpr Scalar sqrt_3_over_5 = static_cast<Scalar>(
      0.7745966692414833770358530799564799221665843410583181653175147532226966183873991);

public:
  static constexpr std::array<barycentric_segment<Scalar>, 3> points{
      barycentric_segment<Scalar>{{1. / 2 - (1. / 2) * sqrt_3_over_5}},
      barycentric_segment<Scalar>{{1. / 2 + (1. / 2) * sqrt_3_over_5}},
      barycentric_segment<Scalar>{{1. / 2}}};
  static constexpr std::array<Scalar, 3> weights{5. / 18, 5. / 18, 4. / 9};
  static constexpr std::size_t n_points = 3;
};

#endif // SEGMENT_QUADS_HPP
