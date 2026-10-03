#ifndef QUADINTS_ADAPTIVE_HPP
#define QUADINTS_ADAPTIVE_HPP
#include <memory>

#include "StopCriterion.hpp"
#include "interface.hpp"
#include "splitter_triangle.hpp"
namespace quadints {

static constexpr size_t comptime_max_level = 6;

template <typename Func, typename Domain>
using return_type = std::invoke_result_t<Func, typename Domain::point_type>;
template <typename Func, typename DomainX, typename DomainY>
using return_type_2d = std::invoke_result_t<Func, typename DomainX::point_type, typename DomainY::point_type>;

template <size_t Level>
struct IntegrateDispatcher {
    template <typename Quadrule, typename Func, typename Domain, typename Scalar>
    static auto call(size_t level, Func&& f, const Domain& cell) -> return_type<Func, Domain> {
        if (level == Level) {
            return quadints::integrate<CollectedQuadratureStatic<Quadrule, Level, Scalar>>(std::forward<Func>(f), cell);
        } else {
            return IntegrateDispatcher<Level - 1>::template call<Quadrule, Func, Domain, Scalar>(
                level, std::forward<Func>(f), cell);
        }
    }
};

template <>
struct IntegrateDispatcher<0> {
    template <typename Quadrule, typename Func, typename Domain, typename Scalar>
    static auto call(size_t level, Func&& f, const Domain& cell) -> return_type<Func, Domain> {
        if (level == 0) {
            return quadints::integrate<CollectedQuadratureStatic<Quadrule, 0, Scalar>>(std::forward<Func>(f), cell);
        } else {
            throw std::runtime_error("reached comptime max level (" + std::to_string(comptime_max_level) + ")");
        }
    }
};

template <size_t Level>
struct IntegrateDispatcher2d {
    template <typename QuadruleX, typename QuadruleY, typename DomainX, typename DomainY, typename Func,
              typename Scalar>
    static auto call(size_t level, Func&& f, const DomainX& cellx, const DomainY& celly)
        -> return_type_2d<Func, DomainX, DomainY> {
        if (level == Level) {
            using q1 = CollectedQuadratureStatic<QuadruleX, Level, Scalar>;
            using q2 = CollectedQuadratureStatic<QuadruleY, Level, Scalar>;
            return quadints::integrate2<q1, q2, DomainX, DomainY, Func, Scalar>(std::forward<Func>(f), cellx, celly);
        } else {
            return IntegrateDispatcher2d<Level - 1>::template call<QuadruleX, QuadruleY, DomainX, DomainY, Func,
                                                                   Scalar>(level, std::forward<Func>(f), cellx, celly);
        }
    }
};

template <>
struct IntegrateDispatcher2d<0> {
    template <typename QuadruleX, typename QuadruleY, typename DomainX, typename DomainY, typename Func,
              typename Scalar>
    static auto call(size_t level, Func&& f, const DomainX& cellx, const DomainY& celly)
        -> return_type_2d<Func, DomainX, DomainY> {
        if (level == 0) {
            using q1 = CollectedQuadratureStatic<QuadruleX, 0, Scalar>;
            using q2 = CollectedQuadratureStatic<QuadruleY, 0, Scalar>;
            return quadints::integrate2<q1, q2, DomainX, DomainY, Func, Scalar>(std::forward<Func>(f), cellx, celly);
        } else {
            throw std::runtime_error("Invalid level");
        }
    }
};

static constexpr size_t comptime_max_level_2d = 3;

/**
 * Concept for an integrator callback function.
 *
 * @tparam Callback The callback function type.
 * @tparam return_type The return type of the .
 * @note signature of callback function is:
 *
 * @param current_integral The current integral value.
 * @param previous_integral The previous integral value.
 * @param level The current level of the integration.
 * @param converged Whether the integration has converged.
 */
template <typename Callback, typename Integrand_return_type>
concept IntegratorCallback = std::invocable<Callback&, CurrentIntegral<Integrand_return_type>,
                                            PreviousIntegral<Integrand_return_type>, size_t, bool>;

template <typename Quadrule, typename Domain, typename Criterion>
class AdaptiveIntegrator {
    using Scalar = decltype(std::declval<const Domain&>().mes());
    static_assert(quadrature_rule<Quadrule, Domain, Scalar>,
                  "Quadrule must be a quadrature rule for the given domain and scalar type");
    static_assert(StopCriterion<Criterion, Scalar>, "Criterion must be a stop criterion");
    template <typename Func>
    constexpr auto integrate_over_level(Func&& f, const Domain& cell, size_t level) -> return_type<Func, Domain> {
        if (level <= comptime_max_level) {
            return IntegrateDispatcher<comptime_max_level>::call<Quadrule, Func, Domain, Scalar>(
                level, std::forward<Func>(f), cell);
        } else {
            auto q = CollectedQuadratureDynamic<Quadrule, Scalar>(level);
            return quadints::integrate<q>(std::forward<Func>(f), cell);
        }
    }
    Criterion criterion;
    size_t max_level;
    size_t start_level;

   public:
    AdaptiveIntegrator(Criterion criterion, size_t max_level, size_t start_level)
        : criterion(criterion), max_level(max_level), start_level(start_level) {}

    template <typename Func>
        requires integrable<Func, typename Domain::point_type, Scalar>
    constexpr auto integrate(Func&& f, const Domain& cell) -> return_type<Func, Domain> {
        using rt = return_type<Func, Domain>;
        rt cur_integral = integrate_over_level(f, cell, start_level);
        rt prev_integral{};

        for (size_t curlevel = start_level + 1; curlevel < max_level; ++curlevel) {
            prev_integral = cur_integral;
            cur_integral = integrate_over_level(f, cell, curlevel);
            if (criterion(CurrentIntegral(cur_integral), PreviousIntegral(prev_integral)) == StopCriterionType::STOP) {
                break;
            }
        }
        return cur_integral;
    }

    template <typename Func, typename Callback>
        requires integrable<Func, typename Domain::point_type, Scalar> &&
                 IntegratorCallback<Callback, return_type<Func, Domain>>
    constexpr auto integrate(Func&& f, const Domain& cell, Callback&& callback) -> return_type<Func, Domain> {
        using rt = return_type<Func, Domain>;
        rt cur_integral = integrate_over_level(f, cell, start_level);
        rt prev_integral{};
        bool stoped_early = false;
        size_t cur_level;
        for (cur_level = start_level + 1; cur_level < max_level; ++cur_level) {
            prev_integral = cur_integral;
            cur_integral = integrate_over_level(f, cell, cur_level);
            if (criterion(CurrentIntegral(cur_integral), PreviousIntegral(prev_integral)) == StopCriterionType::STOP) {
                stoped_early = true;
                break;
            }
        }
        callback(CurrentIntegral(cur_integral), PreviousIntegral(prev_integral), cur_level, stoped_early);
        return cur_integral;
    }
};

template <typename Quadrule, typename Domain, typename Criterion>
auto make_adaptive_integrator(Criterion&& c, std::size_t max_lvl, std::size_t start_lvl) {
    return AdaptiveIntegrator<Quadrule, Domain, Criterion>(std::move(c), max_lvl, start_lvl);
}

template <typename QuadruleX, typename QuadruleY, typename DomainX, typename DomainY, typename Criterion>
class AdaptiveIntegrator2d {
    using Scalar = decltype(std::declval<const DomainX&>().mes());
    static_assert(std::is_same_v<Scalar, decltype(std::declval<const DomainY&>().mes())>,
                  "DomainX and DomainY must have the same scalar type");
    static_assert(quadrature_rule<QuadruleX, DomainX, Scalar>, "QuadruleX must be a quadrature rule for DomainX");
    static_assert(quadrature_rule<QuadruleY, DomainY, Scalar>, "QuadruleY must be a quadrature rule for DomainY");
    static_assert(StopCriterion<Criterion, Scalar>, "Criterion must be a stop criterion for Scalar");
    template <typename Func>
    constexpr auto integrate_over_level(Func&& f, const DomainX& cellx, const DomainY& celly, size_t level)
        -> return_type_2d<Func, DomainX, DomainY> {
        if (level <= comptime_max_level_2d) {
            return IntegrateDispatcher2d<comptime_max_level_2d>::call<QuadruleX, QuadruleY, DomainX, DomainY, Func,
                                                                      Scalar>(level, std::forward<Func>(f), cellx,
                                                                              celly);
        } else {
            auto qx = CollectedQuadratureDynamic<QuadruleX, Scalar>(level);
            auto qy = CollectedQuadratureDynamic<QuadruleY, Scalar>(level);
            return quadints::integrate2<qx, qy, DomainX, DomainY, Func, Scalar>(std::forward<Func>(f), cellx, celly);
        }
    }
    Criterion criterion;
    size_t max_level;
    size_t start_level;

   public:
    AdaptiveIntegrator2d(Criterion criterion, size_t max_level, size_t start_level)
        : criterion(criterion), max_level(max_level), start_level(start_level) {}

    template <typename Func>
        requires integrable2<Func, typename DomainX::point_type, typename DomainY::point_type, Scalar>
    constexpr auto integrate(Func&& f, const DomainX& cellx, const DomainY& celly)
        -> return_type_2d<Func, DomainX, DomainY> {
        using rt = return_type_2d<Func, DomainX, DomainY>;
        rt cur_integral = integrate_over_level(f, cellx, celly, start_level);
        rt prev_integral{};
        for (size_t curlevel = start_level + 1; curlevel <= max_level; ++curlevel) {
            prev_integral = cur_integral;
            cur_integral = integrate_over_level(f, cellx, celly, curlevel);
            if (criterion(CurrentIntegral(cur_integral), PreviousIntegral(prev_integral)) == StopCriterionType::STOP) {
                break;
            }
        }
        return cur_integral;
    }
    template <typename Func, typename Callback>
        requires integrable2<Func, typename DomainX::point_type, typename DomainY::point_type, Scalar> &&
                 IntegratorCallback<Callback, return_type_2d<Func, DomainX, DomainY>>
    constexpr auto integrate(Func&& f, const DomainX& cellx, const DomainY& celly, Callback&& callback)
        -> return_type_2d<Func, DomainX, DomainY> {
        using rt = return_type_2d<Func, DomainX, DomainY>;
        rt cur_integral = integrate_over_level(f, cellx, celly, start_level);
        rt prev_integral{};
        bool stopped_early = false;
        size_t curlevel;
        for (curlevel = start_level + 1; curlevel <= max_level; ++curlevel) {
            prev_integral = cur_integral;
            cur_integral = integrate_over_level(f, cellx, celly, curlevel);
            if (criterion(CurrentIntegral(cur_integral), PreviousIntegral(prev_integral)) == StopCriterionType::STOP) {
                stopped_early = true;
                break;
            }
        }
        callback(CurrentIntegral(cur_integral), PreviousIntegral(prev_integral), curlevel, stopped_early);
        return cur_integral;
    }
};

template <typename QuadRuleX, typename QuadRuleY, typename DomainX, typename DomainY, typename Criterion>
auto make_adaptive_integrator2d(Criterion c, size_t max_level, size_t start_level) {
    return AdaptiveIntegrator2d<QuadRuleX, QuadRuleY, DomainX, DomainY, Criterion>(std::move(c), max_level,
                                                                                   start_level);
}
}  // namespace quadints

#endif  // QUADINTS_ADAPTIVE_HPP
