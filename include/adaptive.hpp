#ifndef QUADINTS_ADAPTIVE_HPP
#define QUADINTS_ADAPTIVE_HPP
#include <splitter_triangle.hpp>

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
    template <size_t L, typename Func>
        requires quadrature_rule<Quadrule, Domain, Scalar> && integrable<Func, typename Domain::point_type, Scalar>
    constexpr auto integrate_over_level(Func&& f, const Domain& cell)
        -> std::invoke_result_t<Func, typename Domain::point_type> {
        using return_type = std::invoke_result_t<Func, typename Domain::point_type>;
        using cur_quad = CollectedQuadrature<Quadrule, L, Scalar>;
        auto domain_points = quadints::detail::get_domain_points<cur_quad>(cell);
        std::array<return_type, cur_quad::n_points> func_arr;
        if constexpr (quadints::detail::has_batch_func<Func, typename Domain::point_type, cur_quad::n_points>) {
            func_arr = std::invoke(f, domain_points);
        } else {
            for (size_t i = 0; i < cur_quad::n_points; ++i) {
                func_arr[i] = std::invoke(f, domain_points[i]);
            }
        }
        return_type cur_res{};
        for (size_t i = 0; i < cur_quad::n_points; ++i) {
            cur_res += func_arr[i] * cur_quad::weights[i];
        }
        return cur_res * cell.mes();
    }

    template <typename Func>
    constexpr auto integrate_over_level_runtime(Func&& f, const Domain& cell, size_t level)
        -> std::invoke_result_t<Func, typename Domain::point_type> {
        using return_type = std::invoke_result_t<Func, typename Domain::point_type>;
        return_type cur_res{};
        switch (level) {
            case 0:
                cur_res = integrate_over_level<0>(f, cell);
                break;
            case 1:
                cur_res = integrate_over_level<1>(f, cell);
                break;
            case 2:
                cur_res = integrate_over_level<2>(f, cell);
                break;
            case 3:
                cur_res = integrate_over_level<3>(f, cell);
                break;
            case 4:
                cur_res = integrate_over_level<4>(f, cell);
                break;
            case 5:
                cur_res = integrate_over_level<5>(f, cell);
                break;
            case 6:
                cur_res = integrate_over_level<6>(f, cell);
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
