#ifndef QUADINTS_CONCEPTS_HPP
#define QUADINTS_CONCEPTS_HPP

#include <cmath>
#include <concepts>
#include <ranges>

namespace quadints {
template <typename T, typename scalar>
concept has_norm = requires(T t) {
    { norm(t) } -> std::convertible_to<scalar>;
};

template <typename T, typename scalar>
concept has_abs = requires(T t) {
    { abs(t) } -> std::convertible_to<scalar>;
};

template <typename T, typename scalar>
concept has_magnitude = has_norm<T, scalar> || has_abs<T, scalar>;

template <typename scalar, typename T>
concept banach_vec = requires(scalar alpha, T t, T u) {
    t + u;
    t* alpha;
    alpha* t;
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
    requires std::ranges::contiguous_range<decltype(T::points)>;
    T::points;
    T::weights;
    requires std::ranges::contiguous_range<decltype(T::weights)>;
    { T::weights[0] } -> std::convertible_to<scalar>;
    requires T::n_points == T::points.size() && T::n_points == T::weights.size();
};

template <typename F, typename arg, typename Scalar>
concept integrable = requires(F f, arg x) { requires banach_vec<Scalar, decltype(f(x))>; };

template <typename F, typename arg1, typename arg2, typename Scalar>
concept integrable2 = requires(F f, arg1 x, arg2 y) {
    { f(x, y) };
    requires banach_vec<Scalar, decltype(f(x, y))>;
};

template <typename F, typename DomainT, typename Scalar>
concept integrable_over_domain = domain<Scalar, DomainT> && integrable<F, typename DomainT::point_type, Scalar>;

template <typename F, typename DomainXT, typename DomainYT, typename Scalar>
concept integrable2_over_domain = domain<Scalar, DomainXT> && domain<Scalar, DomainYT> &&
                                  integrable2<F, typename DomainXT::point_type, typename DomainYT::point_type, Scalar>;

}  // namespace quadints

#endif  // QUADINTS_CONCEPTS_HPP
