#ifndef QUADINTS_ADAPTIVE_HPP
#define QUADINTS_ADAPTIVE_HPP
#include "splitter_triangle.hpp"

#include "interface.hpp"
namespace quadints {
template <typename Scalar>
struct IntegrationParams {
    Scalar rtol;
    Scalar atol;
    size_t start_level;
    size_t max_level;
};

template <typename Quadrule, typename Domain, typename Scalar>
class AdaptiveIntegrator {
    template <typename Func>
    constexpr auto integrate_over_level_runtime(Func&& f, const Domain& cell, size_t level)
        -> std::invoke_result_t<Func, typename Domain::point_type> {
        using return_type = std::invoke_result_t<Func, typename Domain::point_type>;
        return_type cur_res{};
        // 4^6 = 4096 is limit for constexpr evaluation of quadrature.
        switch (level) {
            case 0:
                cur_res = quadints::integrate<CollectedQuadrature<Quadrule, 0, Scalar>>(f, cell);
                break;
            case 1:
                cur_res = quadints::integrate<CollectedQuadrature<Quadrule, 1, Scalar>>(f, cell);
                break;
            case 2:
                cur_res = quadints::integrate<CollectedQuadrature<Quadrule, 2, Scalar>>(f, cell);
                break;
            case 3:
                cur_res = quadints::integrate<CollectedQuadrature<Quadrule, 3, Scalar>>(f, cell);
                break;
            case 4:
                cur_res = quadints::integrate<CollectedQuadrature<Quadrule, 4, Scalar>>(f, cell);
                break;
            case 5:
                cur_res = quadints::integrate<CollectedQuadrature<Quadrule, 5, Scalar>>(f, cell);
                break;
            case 6:
                cur_res = quadints::integrate<CollectedQuadrature<Quadrule, 6, Scalar>>(f, cell);
                break;
            case 7:
                cur_res = quadints::integrate<CollectedQuadrature<Quadrule, 7, Scalar>>(f, cell);
                break;
            default:
                throw std::runtime_error("Invalid level");
        }
        return cur_res;
    }

   public:
    template <typename Func>
        requires quadints::quadrature_rule<Quadrule, Domain, Scalar> &&
                 integrable<Func, typename Domain::point_type, Scalar>
    constexpr auto integrate(Func&& f, const Domain& cell, const IntegrationParams<Scalar>& params)
        -> std::invoke_result_t<Func, typename Domain::point_type> {
        using return_type = std::invoke_result_t<Func, typename Domain::point_type>;
        return_type cur_integral = integrate_over_level_runtime(f, cell, params.start_level);
        return_type prev_integral{};

        for (size_t curlevel = params.start_level + 1; curlevel < params.max_level; ++curlevel) {
            prev_integral = cur_integral;
            cur_integral = integrate_over_level_runtime(f, cell, curlevel);
            std::cerr << "Level " << curlevel << " integral: " << cur_integral << std::endl;
            if (std::abs(cur_integral - prev_integral) < params.atol + params.rtol * std::abs(cur_integral)) {
                std::cerr << "Converged at level " << curlevel << std::endl;
                break;
            }
        }
        return cur_integral;
    }
};

}  // namespace quadints

#endif  // QUADINTS_ADAPTIVE_HPP
