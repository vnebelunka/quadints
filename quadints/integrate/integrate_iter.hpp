#ifndef QUADINTS_INTEGRATE_ITER_HPP
#define QUADINTS_INTEGRATE_ITER_HPP

#include <functional>
#include <ranges>

#include "quadints/common/common.hpp"
#include "quadints/concepts.hpp"
namespace quadints::detail {

template <typename QuadRule, typename Domain, typename Func, typename Scalar = decltype(*(QuadRule::weights.begin()))>
    requires quadrature_rule<QuadRule, Domain, Scalar> && integrable_over_domain<Func, Domain, Scalar>
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

template <typename Domain, typename Func, typename QuadRulePointType, typename Scalar>
constexpr auto integrate_iter(Func&& f, const Domain& cell, const std::vector<QuadRulePointType>& quad_rule_points,
                              const std::vector<Scalar>& quad_rule_weights) -> return_type<Func, Domain> {
    using return_type = return_type<Func, Domain>;
    return_type res{};
    auto points_it = std::ranges::begin(quad_rule_points);
    auto weights_it = std::ranges::begin(quad_rule_weights);
    const auto points_end = std::ranges::end(quad_rule_points);
    for (; points_it != points_end; ++points_it, ++weights_it) {
        const auto& ref_point = *points_it;
        const auto weight = static_cast<Scalar>(*weights_it);  // convert to Scalar
        const auto point = ref_point.to_domain(cell);
        res += std::invoke(f, point) * weight;
    }
    return res * cell.mes();
}
template <typename Domain, typename Func, typename RuntimeQuadRule, typename Scalar = decltype(Domain().mes())>
constexpr auto integrate_iter(Func&& f, const Domain& cell, const RuntimeQuadRule& quad) {
    return integrate_iter(std::forward<Func>(f), cell, quad.points(), quad.weights());
}

template <typename Domain1, typename Domain2, typename Func, typename QuadRule1PointType1, typename QuadRule2PointType2,
          typename Scalar = decltype(Domain2().mes())>
    requires integrable2<Func, typename Domain1::point_type, typename Domain2::point_type, Scalar>
constexpr auto integrate2_iter(Func&& f, const Domain1& cell_x, const Domain2& cell_y,
                               const std::vector<QuadRule1PointType1>& quad_rule_points_x,
                               const std::vector<Scalar>& quad_rule_weights_x,
                               const std::vector<QuadRule2PointType2>& quad_rule_points_y,
                               const std::vector<Scalar>& quad_rule_weights_y) {
    using return_type = std::invoke_result_t<Func, typename Domain1::point_type, typename Domain2::point_type>;
    return_type res{};
    auto points1_it = std::ranges::begin(quad_rule_points_x);
    auto weights1_it = std::ranges::begin(quad_rule_weights_x);
    auto points1_end = std::ranges::end(quad_rule_points_x);
    for (; points1_it != points1_end; ++points1_it, ++weights1_it) {
        const auto& ref_point1 = *points1_it;
        const auto weight1 = static_cast<Scalar>(*weights1_it);
        const auto point1 = ref_point1.to_domain(cell_x);
        for (auto weights2_it = std::begin(quad_rule_weights_y), points2_it = std::begin(quad_rule_points_y);
             points2_it != std::end(quad_rule_points_y); ++points2_it, ++weights2_it) {
            const auto& ref_point2 = *points2_it;
            const auto weight2 = static_cast<Scalar>(*weights2_it);
            const auto point2 = ref_point2.to_domain(cell_y);
            res += std::invoke(f, point1, point2) * weight1 * weight2;
        }
    }
    return res * cell_x.mes() * cell_y.mes();
}

template <typename Domain1, typename Domain2, typename Func, typename DynamicQuadRuleX, typename DynamicWaudRuleY,
          typename Scalar = decltype(Domain2().mes())>
constexpr auto integrate2_iter(Func&& f, const Domain1& cellx, const Domain2& celly, const DynamicQuadRuleX& quadx,
                               const DynamicWaudRuleY& quady) {
    return integrate2_iter(std::forward<Func>(f), cellx, celly, quadx.points(), quadx.weights(), quady.points(),
                           quady.weights());
}
}  // namespace quadints::detail

#endif  // QUADINTS_INTEGRATE_ITER_HPP
