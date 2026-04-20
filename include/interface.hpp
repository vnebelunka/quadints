#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include <concepts>
#include <functional>
#include <iterator>
#include <ranges>

namespace quadints {

template <typename T, typename scalar>
concept has_norm = requires(T t) {
  { norm(t) } -> std::convertible_to<scalar>;
};

template <typename T, typename scalar>
concept has_abs = requires(T t) {
  { std::abs(t) } -> std::convertible_to<scalar>;
};

template <typename T, typename scalar>
concept has_magnitude = has_norm<T, scalar> || has_abs<T, scalar>;

template <typename scalar, typename T>
concept banach_vec = requires(scalar alpha, T t, T u) {
  t + u;
  t *alpha;
  alpha *t;
  t += u;
  t *= alpha;
  requires has_magnitude<T, scalar>;
};

template <typename scalar, typename T>
concept domain = requires(T t) {
  typename T::point_type;
  { t.mes() } -> std::convertible_to<scalar>;
};

template <typename T, typename D, typename scalar>
concept quadrature_rule = requires() {
  requires domain<scalar, D>;
  { T::n_points } -> std::convertible_to<size_t>;
  T::points;
  requires std::ranges::forward_range<decltype(T::points)>;
  T::weights;
  requires std::ranges::forward_range<decltype(T::weights)>;
  { *T::weights.begin() } -> std::convertible_to<scalar>;
  requires T::n_points == T::points.size() && T::n_points == T::weights.size();
};

template <typename F, typename arg, typename Scalar>
concept integrable =
    requires(F f, arg x) { requires banach_vec<Scalar, decltype(f(x))>; };

template <typename F, typename arg1, typename arg2, typename Scalar>
concept integrable2 = requires(F f, arg1 x, arg2 y) {
  requires banach_vec<Scalar, decltype(f(x, y))>;
};

namespace detail_ {
template <typename QuadRule, typename Domain, typename Func, typename Scalar,
          size_t... Idx>
requires quadrature_rule<QuadRule, Domain, Scalar> &&
         integrable<Func, typename Domain::point_type, Scalar>
constexpr decltype(auto) impl_integrate_(Func &&f, const Domain &cell,
                               std::index_sequence<Idx...> /**/) {
  using return_type =
      decltype(std::declval<
                   std::invoke_result_t<Func, typename Domain::point_type>>() *
               std::declval<Scalar>());

  return ((std::invoke(std::forward<Func>(f),
                       QuadRule::points[Idx].to_domain(cell)) *
           QuadRule::weights[Idx]) +
          ... + return_type{}) *
         cell.mes();
}
} // namespace detail_

template <typename QuadRule, typename Domain, typename Func,
          typename Scalar = decltype(*QuadRule::weights.begin())>
constexpr decltype(auto) integrate(Func &&f, const Domain &cell) {
  return detail_::impl_integrate_<QuadRule, Domain, Func, Scalar>(
      std::forward<Func>(f), cell,
      std::make_index_sequence<QuadRule::n_points>{});
}

template <typename QuadRule1, typename QuadRule2 = QuadRule1, typename Domain1,
          typename Domain2, typename Func,
          typename Scalar = decltype(*QuadRule1::weights.begin())>
requires quadrature_rule<QuadRule1, Domain1, Scalar> &&
         quadrature_rule<QuadRule2, Domain2, Scalar> &&
         integrable2<Func, typename Domain1::point_type,
                     typename Domain2::point_type, Scalar>
constexpr auto integrate2(Func &&f, const Domain1 &cell1,
                          const Domain2 &cell2) {
  using return_type =
      decltype(std::declval<
                   std::invoke_result_t<Func, typename Domain1::point_type,
                                        typename Domain2::point_type>>() *
               std::declval<Scalar>());

  return_type res{};
  auto points1_it = std::ranges::begin(QuadRule1::points);
  auto weights1_it = std::ranges::begin(QuadRule1::weights);
  auto points1_end = std::ranges::end(QuadRule1::points);
  for (; points1_it != points1_end; ++points1_it, ++weights1_it) {
    const auto &ref_point1 = *points1_it;
    const auto weight1 = static_cast<Scalar>(*weights1_it);
    const auto point1 = ref_point1.to_domain(cell1);
    for (auto weights2_it = std::begin(QuadRule2::weights),
              points2_it = std::begin(QuadRule2::points);
         points2_it != std::end(QuadRule2::points);
         ++points2_it, ++weights2_it) {
      const auto &ref_point2 = *points2_it;
      const auto weight2 = static_cast<Scalar>(*weights2_it);
      const auto point2 = ref_point2.to_domain(cell2);
      res += std::invoke(f, point1, point2) * weight1 * weight2;
    }
  }
  return res * cell1.mes() * cell2.mes();
};

} // namespace quadints
#endif // INTERFACE_HPP
