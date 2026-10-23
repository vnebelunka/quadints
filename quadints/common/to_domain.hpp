#ifndef QUADINTS_TO_DOMAIN_HPP
#define QUADINTS_TO_DOMAIN_HPP

#include <array>
#include <concepts>
#include <cstdlib>
#include <vector>

namespace quadints {
template <std::size_t N, typename Domain, typename RefPoint>
auto to_domain(const Domain& cell, const std::array<RefPoint, N>& ref_pts)
    -> std::enable_if_t<false, std::array<typename Domain::point_type, N>>;

template <typename Domain, typename RefPoint, typename DomainPoint, std::size_t N>
concept has_batch_to_domain = requires(const Domain& cell, const std::array<RefPoint, N>& ref_pts) {
    { to_domain<N>(cell, ref_pts) } -> std::convertible_to<std::array<DomainPoint, N>>;
};

template <typename Domain, typename RefPoint, typename DomainPoint>
concept has_batch_to_domain_dynamic = requires(const Domain& cell, const std::vector<RefPoint>& ref_pts) {
    { to_domain(cell, ref_pts) } -> std::convertible_to<std::vector<DomainPoint>>;
};

template <typename QuadRule, typename Domain>
auto get_domain_points(const Domain& cell) -> std::array<typename Domain::point_type, QuadRule::n_points> {
    std::array<typename Domain::point_type, QuadRule::n_points> res{};
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

template <typename BarPointType, typename Domain>
auto get_domain_points(const Domain& cell, const std::vector<BarPointType>& points) {
    if constexpr (has_batch_to_domain_dynamic<Domain, BarPointType, typename Domain::point_type>) {
        return to_domain(cell, points);
    } else {
        std::vector<typename Domain::point_type> res(points.size());
        for (std::size_t i = 0; i < points.size(); ++i) {
            res[i] = points[i].to_domain(cell);
        }
        return res;
    }
}
}  // namespace quadints

#endif  // QUADINTS_TO_DOMAIN_HPP
