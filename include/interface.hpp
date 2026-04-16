#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include <concepts>
#include <functional>
#include <iterator>
#include <ranges>

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
  t * alpha;
  alpha * t;
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

template <typename QuadRule, typename Domain, typename Func,
          typename Scalar = decltype(*(QuadRule::weights.begin()))>
  requires quadrature_rule<QuadRule, Domain, Scalar> &&
           integrable<Func, typename Domain::point_type, Scalar>
constexpr auto integrate(Func &&f, const Domain &cell)
    -> std::invoke_result_t<Func, typename Domain::point_type> {
  using return_type = std::invoke_result_t<Func, typename Domain::point_type>;
  return_type res{};
  auto points_it = std::ranges::begin(QuadRule::points);
  auto weights_it = std::ranges::begin(QuadRule::weights);
  const auto points_end = std::ranges::end(QuadRule::points);

  for (; points_it != points_end; ++points_it, ++weights_it) {
    const auto &ref_point = *points_it;
    const auto weight = static_cast<Scalar>(*weights_it); // convert to Scalar
    const auto point = ref_point.to_domain(cell);
    res += std::invoke(f, point) * weight;
  }
  return res * cell.mes();
}

#endif // INTERFACE_HPP
