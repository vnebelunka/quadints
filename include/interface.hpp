#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include <concepts>
#include <functional>
#include <iostream>
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
    requires std::ranges::contiguous_range<decltype(T::points)>;
    T::weights;
    requires std::ranges::contiguous_range<decltype(T::weights)>;
    { T::weights[0] } -> std::convertible_to<scalar>;
    requires T::n_points == T::points.size() && T::n_points == T::weights.size();
};

template <typename F, typename arg, typename Scalar>
concept integrable = requires(F f, arg x) { requires banach_vec<Scalar, decltype(f(x))>; };

template <typename F, typename arg1, typename arg2, typename Scalar>
concept integrable2 = requires(F f, arg1 x, arg2 y) { requires banach_vec<Scalar, decltype(f(x, y))>; };

namespace detail {

template <typename QuadRule, typename Domain, typename Func, typename Scalar = decltype(*(QuadRule::weights.begin()))>
    requires quadrature_rule<QuadRule, Domain, Scalar> && integrable<Func, typename Domain::point_type, Scalar>
constexpr auto integrate_iter(Func&& f, const Domain& cell) -> std::invoke_result_t<Func, typename Domain::point_type> {
    using return_type = std::invoke_result_t<Func, typename Domain::point_type>;
    return_type res{};
    auto points_it = std::ranges::begin(QuadRule::points);
    auto weights_it = std::ranges::begin(QuadRule::weights);
    const auto points_end = std::ranges::end(QuadRule::points);

    for (; points_it != points_end; ++points_it, ++weights_it) {
        const auto& ref_point = *points_it;
        const auto weight = static_cast<Scalar>(*weights_it);  // convert to Scalar
        const auto point = ref_point.to_domain(cell);
        res += std::invoke(f, point) * weight;
    }
    return res * cell.mes();
}

template <typename Domain, typename RefPoint, typename DomainPoint, std::size_t N>
concept has_batch_to_domain = requires(const Domain& cell, const std::array<RefPoint, N>& ref_pts) {
    { to_domain<N>(cell, ref_pts) } -> std::convertible_to<std::array<DomainPoint, N>>;
};

template <typename QuadRule, typename Domain>
auto get_domain_points(const Domain& cell) -> std::array<typename Domain::point_type, QuadRule::n_points> {
    std::array<typename Domain::point_type, QuadRule::n_points> res;
    if constexpr (has_batch_to_domain<Domain, typename QuadRule::point_type, typename Domain::point_type,
                                      QuadRule::n_points>) {
        res = to_domain<QuadRule::n_points>(cell, QuadRule::points);
    } else {
        for (std::size_t i = 0; i < QuadRule::n_points; ++i) {
            res[i] = QuadRule::points[i].to_domain(cell);
        }
    }
    return res;
}

template <typename Func, typename DomainPoint, std::size_t N>
concept has_batch_func = requires(Func&& f, const std::array<DomainPoint, N>& pts) {
    { f(pts) } -> std::convertible_to<std::array<typename DomainPoint::scalar_type, N>>;
};

template <typename Func, typename DomainPointX, typename DomainPointY, std::size_t N>
concept has_batch_func2 =
    requires(Func&& f, const std::array<DomainPointX, N>& ptsX, const std::array<DomainPointY, N>& ptsY) {
        { f(ptsX, ptsY) } -> std::convertible_to<std::array<typename DomainPointX::scalar_type, N>>;
    };

template <typename QuadRule, typename Domain, typename Func,
          typename Scalar = std::decay_t<decltype(*(QuadRule::weights.begin()))>>
    requires quadrature_rule<QuadRule, Domain, Scalar> && integrable<Func, typename Domain::point_type, Scalar>
constexpr auto integrate_collect(Func&& f, const Domain& cell)
    -> std::invoke_result_t<Func, typename Domain::point_type> {
    using return_type = std::invoke_result_t<Func, typename Domain::point_type>;
    return_type res{};

    std::array<Scalar, QuadRule::n_points> weights_arr;
    std::array<return_type, QuadRule::n_points> func_arr;
    auto domain_points = get_domain_points<QuadRule>(cell);
    if constexpr (has_batch_func<Func, typename Domain::point_type, QuadRule::n_points>) {
        func_arr = std::invoke(f, domain_points);
    } else {
        for (size_t i = 0; i < QuadRule::n_points; ++i) {
            func_arr[i] = std::invoke(f, domain_points[i]);
        }
    }
    for (size_t i = 0; i < QuadRule::n_points; ++i) {
        res += func_arr[i] * QuadRule::weights[i];
    }
    return res * cell.mes();
}

template <typename QuadRule1, typename QuadRule2 = QuadRule1, typename Domain1, typename Domain2, typename Func,
          typename Scalar = decltype(*(QuadRule1::weights.begin()))>
    requires quadrature_rule<QuadRule1, Domain1, Scalar> && quadrature_rule<QuadRule2, Domain2, Scalar> &&
             integrable2<Func, typename Domain1::point_type, typename Domain2::point_type, Scalar>
constexpr auto integrate2_iter(Func&& f, const Domain1& cell1, const Domain2& cell2) {
    using return_type = std::invoke_result_t<Func, typename Domain1::point_type, typename Domain2::point_type>;
    return_type res{};
    auto points1_it = std::ranges::begin(QuadRule1::points);
    auto weights1_it = std::ranges::begin(QuadRule1::weights);
    auto points1_end = std::ranges::end(QuadRule1::points);
    for (; points1_it != points1_end; ++points1_it, ++weights1_it) {
        const auto& ref_point1 = *points1_it;
        const auto weight1 = static_cast<Scalar>(*weights1_it);
        const auto point1 = ref_point1.to_domain(cell1);
        for (auto weights2_it = std::begin(QuadRule2::weights), points2_it = std::begin(QuadRule2::points);
             points2_it != std::end(QuadRule2::points); ++points2_it, ++weights2_it) {
            const auto& ref_point2 = *points2_it;
            const auto weight2 = static_cast<Scalar>(*weights2_it);
            const auto point2 = ref_point2.to_domain(cell2);
            res += std::invoke(f, point1, point2) * weight1 * weight2;
        }
    }
    return res * cell1.mes() * cell2.mes();
};

template <typename QuadRule1, typename QuadRule2 = QuadRule1, typename Domain1, typename Domain2, typename Func,
          typename Scalar = decltype(*(QuadRule1::weights.begin()))>
    requires quadrature_rule<QuadRule1, Domain1, Scalar> && quadrature_rule<QuadRule2, Domain2, Scalar> &&
             integrable2<Func, typename Domain1::point_type, typename Domain2::point_type, Scalar>
constexpr auto integrate2_collect(Func&& f, const Domain1& cell1, const Domain2& cell2) {
    using return_type = std::invoke_result_t<Func, typename Domain1::point_type, typename Domain2::point_type>;
    return_type res{};
    constexpr size_t n_points = QuadRule1::n_points * QuadRule2::n_points;

    std::array<return_type, n_points> func_arr{};
    std::array<typename Domain1::point_type, QuadRule1::n_points> domain_points1 = get_domain_points<QuadRule1>(cell1);
    std::array<typename Domain2::point_type, QuadRule2::n_points> domain_points2 = get_domain_points<QuadRule2>(cell2);
    if constexpr (has_batch_func2<Func, typename Domain1::point_type, typename Domain2::point_type, n_points>) {
        func_arr = std::invoke(f, domain_points1, domain_points2);
    } else {
        for (size_t j = 0; j < QuadRule2::n_points; ++j) {
            for (size_t i = 0; i < QuadRule1::n_points; ++i) {
                func_arr[j * QuadRule1::n_points + i] = std::invoke(f, domain_points1[i], domain_points2[j]);
            }
        }
    }
    for (size_t index_y = 0; index_y < QuadRule2::n_points; ++index_y) {
        return_type resj = 0;
        Scalar wj = QuadRule2::weights[index_y];
        for (size_t index_x = 0; index_x < QuadRule1::n_points; ++index_x) {
            resj += func_arr[index_y * QuadRule1::n_points + index_x] * QuadRule1::weights[index_x];
        }
        res += resj * wj;
    }
    return res;
};

}  // namespace detail

template <typename QuadRule1, typename QuadRule2 = QuadRule1, typename Domain1, typename Domain2, typename Func,
          typename Scalar = decltype(*(QuadRule1::weights.begin()))>
    requires quadrature_rule<QuadRule1, Domain1, Scalar> && quadrature_rule<QuadRule2, Domain2, Scalar> &&
             integrable2<Func, typename Domain1::point_type, typename Domain2::point_type, Scalar>
constexpr auto integrate2(Func&& f, const Domain1& cell1, const Domain2& cell2) {
    return detail::integrate2_collect<QuadRule1, QuadRule2>(f, cell1, cell2);
}

template <typename QuadRule, typename Domain, typename Func, typename Scalar = decltype(*(QuadRule::weights.begin()))>
    requires quadrature_rule<QuadRule, Domain, Scalar> && integrable<Func, typename Domain::point_type, Scalar>
constexpr auto integrate(Func&& f, const Domain& cell) -> std::invoke_result_t<Func, typename Domain::point_type> {
    return detail::integrate_collect<QuadRule>(f, cell);
}

}  // namespace quadints
#endif  // INTERFACE_HPP
