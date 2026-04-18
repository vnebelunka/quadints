#ifndef BARYCENTRIC_PRODUCT_HPP
#define BARYCENTRIC_PRODUCT_HPP

#include "interface.hpp"
#include <array>
#include <cmath>
#include <type_traits>

template <typename QuadratureX, typename QuadratureY,
          typename Scalar =
              std::remove_cvref_t<decltype(QuadratureX::weights[0])>>
struct Quadrature_product {
  using BarycentricX = std::remove_cvref_t<decltype(QuadratureX::points[0])>;
  using BarycentricY = std::remove_cvref_t<decltype(QuadratureY::points[0])>;
  static constexpr size_t dim_x = QuadratureX::points[0].coords.size();
  static constexpr size_t dim_y = QuadratureY::points[0].coords.size();

  static constexpr std::size_t n_points =
      QuadratureX::n_points * QuadratureY::n_points;
  struct barycentric_coord {
    std::array<Scalar, dim_x + dim_y> coords;
    template <typename Domain>
      requires requires(barycentric_coord c, Domain d) {
        d.from_reference_domain(c);
      }
    constexpr auto to_domain(const Domain &domain) const {
      return domain.from_reference_domain(*this);
    }
  };

  constexpr static auto points = []() {
    std::array<barycentric_coord, n_points> arr{};
    size_t idx = 0;
    for (size_t i = 0; i < QuadratureX::n_points; ++i) {
      for (size_t j = 0; j < QuadratureY::n_points; ++j) {
        std::copy(QuadratureX::points[i].coords.begin(),
                  QuadratureX::points[i].coords.end(), arr[idx].coords.begin());
        std::copy(QuadratureY::points[j].coords.begin(),
                  QuadratureY::points[j].coords.end(),
                  arr[idx].coords.begin() + dim_x);
        ++idx;
      }
    }
    return arr;
  }();

  constexpr static auto weights = []() {
    std::array<Scalar, n_points> arr{};
    size_t idx = 0;
    for (size_t i = 0; i < QuadratureX::n_points; ++i) {
      for (size_t j = 0; j < QuadratureY::n_points; ++j) {
        arr[idx++] = QuadratureX::weights[i] * QuadratureY::weights[j];
      }
    }
    return arr;
  }();
};

#endif // BARYCENTRIC_PRODUCT_HPP
