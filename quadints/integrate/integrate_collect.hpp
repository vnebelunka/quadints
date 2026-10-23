#ifndef QUADINTS_INTEGRATE_COLLECT_HPP
#define QUADINTS_INTEGRATE_COLLECT_HPP
#include <functional>
#include <utility>

#include "quadints/common/common.hpp"
#include "quadints/common/to_domain.hpp"

namespace quadints::detail {
template <typename QuadRule, typename Domain, typename Func, typename Scalar = std::decay_t<decltype(Domain().mes())>>
    requires quadrature_rule<QuadRule, Domain, Scalar> && integrable_over_domain<Func, Domain, Scalar>
constexpr auto integrate_collect(Func&& f, const Domain& cell) -> return_type<Func, Domain> {
    using return_type = return_type<Func, Domain>;
    return_type res{};
    std::array<return_type, QuadRule::n_points> func_arr{};
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

template <typename Domain, typename Func, typename QuadRulePointType, typename Scalar>
constexpr auto integrate_collect(Func&& f, const Domain& cell, const std::vector<QuadRulePointType>& quad_rule_points,
                                 const std::vector<Scalar>& quad_rule_weights) -> return_type<Func, Domain> {
    using return_type = std::invoke_result_t<Func, typename Domain::point_type>;
    return_type res{};
    auto domain_points = get_domain_points(cell, quad_rule_points);
    std::vector<return_type> func_arr(quad_rule_points.size());
    for (std::size_t i = 0; i < quad_rule_points.size(); ++i) {
        func_arr[i] = std::invoke(f, domain_points[i]);
    }
    for (std::size_t i = 0; i < quad_rule_weights.size(); ++i) {
        res += func_arr[i] * quad_rule_weights[i];
    }
    return res * cell.mes();
}

template <typename Domain, typename Func, typename QuadRulePointType, typename BufferT,
          typename Scalar = decltype(Domain().mes())>
constexpr auto integrate_collect(Func&& f, const Domain& cell, const std::vector<QuadRulePointType>& quad_rule_points,
                                 const std::vector<Scalar>& quad_rule_weights,
                                 BufferT& buff) -> return_type<Func, Domain> {
    using return_type = return_type<Func, Domain>;
    return_type res{};
    auto domain_points = get_domain_points(cell, quad_rule_points);
    for (size_t i = 0; i < quad_rule_points.size(); ++i) {
        buff[i] = std::invoke(f, domain_points[i]);
    }
    for (size_t i = 0; i < quad_rule_weights.size(); ++i) {
        res += buff[i] * quad_rule_weights[i];
    }
    return res * cell.mes();
}

template <typename Domain, typename Func, typename RuntimeQuadRule, typename Scalar = decltype(Domain().mes())>
constexpr auto integrate_collect(Func&& f, const Domain& cell, const RuntimeQuadRule& quad) {
    return integrate_collect(std::forward<Func>(f), cell, quad.points(), quad.weights());
}

template <typename Domain, typename Func, typename RuntimeQuadRule, typename BufferT,
          typename Scalar = decltype(Domain().mes())>
constexpr auto integrate_collect(Func&& f, const Domain& cell, const RuntimeQuadRule& quad, BufferT& buff) {
    return integrate_collect(std::forward<Func>(f), cell, quad.points(), quad.weights(), buff);
}

template <typename QuadRule1, typename QuadRule2 = QuadRule1, typename Domain1, typename Domain2, typename Func,
          typename Scalar = decltype(*(QuadRule1::weights.begin()))>
constexpr auto integrate2_collect(Func&& f, const Domain1& cell1, const Domain2& cell2) {
    using return_type = std::invoke_result_t<Func, typename Domain1::point_type, typename Domain2::point_type>;
    return_type res{};
    constexpr size_t n_points = QuadRule1::n_points * QuadRule2::n_points;
    constexpr std::size_t DEFAULT_STACK_BYTES = static_cast<size_t>(8) * 1024 * 1024;
    std::array<return_type, n_points> func_arr{};
    static_assert(n_points * sizeof(return_type) < DEFAULT_STACK_BYTES, "func_arr too large");
    std::array<typename Domain1::point_type, QuadRule1::n_points> domain_points1 = get_domain_points<QuadRule1>(cell1);
    std::array<typename Domain2::point_type, QuadRule2::n_points> domain_points2 = get_domain_points<QuadRule2>(cell2);
    if constexpr (has_batch_func2<Func, typename Domain1::point_type, typename Domain2::point_type, n_points>) {
        func_arr = std::invoke(f, domain_points1, domain_points2);
    } else {
        for (size_t j = 0; j < QuadRule2::n_points; ++j) {
            for (size_t i = 0; i < QuadRule1::n_points; ++i) {
                func_arr[(j * QuadRule1::n_points) + i] = std::invoke(f, domain_points1[i], domain_points2[j]);
            }
        }
    }
    for (size_t index_y = 0; index_y < QuadRule2::n_points; ++index_y) {
        return_type resj = 0;
        Scalar wj = QuadRule2::weights[index_y];
        for (size_t index_x = 0; index_x < QuadRule1::n_points; ++index_x) {
            resj += func_arr[(index_y * QuadRule1::n_points) + index_x] * QuadRule1::weights[index_x];
        }
        res += resj * wj;
    }
    return res * cell1.mes() * cell2.mes();
};

template <typename Domain1, typename Domain2, typename Func, typename QuadRule1PointType1, typename QuadRule2PointType2,
          typename Scalar = decltype(Domain2().mes())>
    requires integrable2<Func, typename Domain1::point_type, typename Domain2::point_type, Scalar>
constexpr auto integrate2_collect(Func&& f, const Domain1& cell_x, const Domain2& cell_y,
                                  const std::vector<QuadRule1PointType1>& quad_rule_points_x,
                                  const std::vector<Scalar>& quad_rule_weights_x,
                                  const std::vector<QuadRule2PointType2>& quad_rule_points_y,
                                  const std::vector<Scalar>& quad_rule_weights_y) {
    using rt = std::invoke_result_t<Func, typename Domain1::point_type, typename Domain2::point_type>;
    auto domain_points_x = get_domain_points(cell_x, quad_rule_points_x);
    auto domain_points_y = get_domain_points(cell_y, quad_rule_points_y);
    std::vector<rt> func_arr(quad_rule_points_x.size() * quad_rule_points_y.size());
    size_t const Nx = quad_rule_points_x.size();
    size_t const Ny = quad_rule_points_y.size();

    for (size_t j = 0; j < Ny; ++j) {
        for (size_t i = 0; i < Nx; ++i) {
            func_arr[(j * Ny) + i] = std::invoke(f, domain_points_x[i], domain_points_y[j]);
        }
    }

    rt res{};

    for (size_t index_y = 0; index_y < Ny; ++index_y) {
        rt resj{};
        Scalar wj = quad_rule_weights_y[index_y];
        for (size_t index_x = 0; index_x < Nx; ++index_x) {
            resj += func_arr[(index_y * Nx) + index_x] * quad_rule_weights_x[index_x];
        }
        res += resj * wj;
    }
    return res * cell_x.mes() * cell_y.mes();
}

template <typename Domain1, typename Domain2, typename Func, typename DynamicQuadRuleX, typename DynamicWaudRuleY,
          typename Scalar = decltype(Domain2().mes())>
constexpr auto integrate2_collect(Func&& f, const Domain1& cellx, const Domain2& celly, const DynamicQuadRuleX& quadx,
                                  const DynamicWaudRuleY& quady) {
    return integrate2_collect(std::forward<Func>(f), cellx, celly, quadx.points(), quadx.weights(), quady.points(),
                              quady.weights());
}

}  // namespace quadints::detail

#endif  // QUADINTS_INTEGRATE_COLLECT_HPP
