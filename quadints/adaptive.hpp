#ifndef QUADINTS_ADAPTIVE_HPP
#define QUADINTS_ADAPTIVE_HPP
#include <memory>

#include "interface.hpp"
#include "splitter_triangle.hpp"
namespace quadints {
template <typename Scalar>
struct IntegrationParams {
    Scalar rtol;
    Scalar atol;
    size_t start_level;
    size_t max_level;
};

template <size_t Level>
struct IntegrateDispatcher {
    template <typename Quadrule, typename Func, typename Domain, typename Scalar>
    static auto call(size_t level, Func&& f, const Domain& cell) {
        if (level == Level) {
            return quadints::integrate<CollectedQuadrature<Quadrule, Level, Scalar>>(std::forward<Func>(f), cell);
        } else {
            return IntegrateDispatcher<Level - 1>::template call<Quadrule, Func, Domain, Scalar>(
                level, std::forward<Func>(f), cell);
        }
    }
};
template <>
struct IntegrateDispatcher<0> {
    template <typename Quadrule, typename Func, typename Domain, typename Scalar>
    static auto call(size_t level, Func&& f, const Domain& cell) {
        if (level == 0) {
            return quadints::integrate<CollectedQuadrature<Quadrule, 0, Scalar>>(std::forward<Func>(f), cell);
        } else {
            throw std::runtime_error("Invalid level");
        }
    }
};

template <size_t Level>
struct IntegrateDispatcher2d {
    template <typename QuadruleX, typename QuadruleY, typename DomainX, typename DomainY, typename Func,
              typename Scalar>
    static Scalar call(size_t level, Func&& f, const DomainX& cellx, const DomainY& celly) {
        if (level == Level) {
            using q1 = CollectedQuadrature<QuadruleX, Level, Scalar>;
            using q2 = CollectedQuadrature<QuadruleY, Level, Scalar>;
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
    static Scalar call(size_t level, Func&& f, const DomainX& cellx, const DomainY& celly) {
        if (level == 0) {
            using q1 = CollectedQuadrature<QuadruleX, 0, Scalar>;
            using q2 = CollectedQuadrature<QuadruleY, 0, Scalar>;
            return quadints::integrate2<q1, q2, DomainX, DomainY, Func, Scalar>(std::forward<Func>(f), cellx, celly);
        } else {
            throw std::runtime_error("Invalid level");
        }
    }
};

static constexpr size_t comptime_max_level = 7;
static constexpr size_t comptime_max_level_2d = 4;

template <typename Func, typename Domain>
using return_type = std::invoke_result_t<Func, typename Domain::point_type>;
template <typename Func, typename DomainX, typename DomainY>
using return_type_2d = std::invoke_result_t<Func, typename DomainX::point_type, typename DomainY::point_type>;

template <typename Quadrule, typename Domain, typename Scalar>
class AdaptiveIntegrator {
    template <typename Func>
    constexpr auto integrate_over_level(Func&& f, const Domain& cell, size_t level) -> return_type<Func, Domain> {
        return IntegrateDispatcher<comptime_max_level>::call<Quadrule, Func, Domain, Scalar>(
            level, std::forward<Func>(f), cell);
    }

   public:
    template <typename Func>
        requires quadints::quadrature_rule<Quadrule, Domain, Scalar> &&
                 integrable<Func, typename Domain::point_type, Scalar>
    constexpr auto integrate(Func&& f, const Domain& cell, const IntegrationParams<Scalar>& params)
        -> return_type<Func, Domain> {
        using rt = return_type<Func, Domain>;
        rt cur_integral = integrate_over_level(f, cell, params.start_level);
        rt prev_integral{};

        for (size_t curlevel = params.start_level + 1; curlevel < params.max_level; ++curlevel) {
            prev_integral = cur_integral;
            cur_integral = integrate_over_level(f, cell, curlevel);
            std::cerr << "Level " << curlevel << " integral: " << cur_integral << std::endl;
            if (std::abs(cur_integral - prev_integral) < params.atol + params.rtol * std::abs(cur_integral)) {
                std::cerr << "Converged at level " << curlevel << std::endl;
                break;
            }
        }
        return cur_integral;
    }
};

template <typename QuadruleX, typename QuadruleY, typename DomainX, typename DomainY, typename Scalar>
class AdaptiveIntegrator2d {
    template <typename Func>
    constexpr auto integrate_over_level(Func&& f, const DomainX& cellx, const DomainY& celly, size_t level)
        -> return_type_2d<Func, DomainX, DomainY> {
        return IntegrateDispatcher2d<comptime_max_level_2d>::call<QuadruleX, QuadruleY, DomainX, DomainY, Func, Scalar>(
            level, std::forward<Func>(f), cellx, celly);
    }

   public:
    template <typename Func>
    constexpr auto integrate(Func&& f, const DomainX& cellx, const DomainY& celly,
                             const IntegrationParams<Scalar>& params) -> return_type_2d<Func, DomainX, DomainY> {
        using rt = return_type_2d<Func, DomainX, DomainY>;
        rt cur_integral = integrate_over_level(f, cellx, celly, params.start_level);
        rt prev_integral{};
        for (size_t curlevel = params.start_level + 1; curlevel <= comptime_max_level; ++curlevel) {
            prev_integral = cur_integral;
            cur_integral = integrate_over_level(f, cellx, celly, curlevel);
            std::cerr << "level " << curlevel << " integral " << cur_integral << std::endl;
            if (std::abs(cur_integral - prev_integral) < params.atol + params.rtol * std::abs(cur_integral)) {
                std::cerr << "Converged at level " << curlevel << std::endl;
                std::cerr << "diff: " << std::abs(cur_integral - prev_integral) << std::endl;
                std::cerr << "integral: " << cur_integral << std::endl;
                std::cerr << "prev_integral: " << prev_integral << std::endl;
                break;
            }
        }
        return cur_integral;
    }
};

}  // namespace quadints

#endif  // QUADINTS_ADAPTIVE_HPP
