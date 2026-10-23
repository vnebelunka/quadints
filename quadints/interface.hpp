#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include <concepts>
#include <functional>
#include <iostream>
#include <iterator>
#include <ranges>
#include <utility>

#include "common/common.hpp"
#include "common/to_domain.hpp"
#include "concepts.hpp"
#include "integrate/integrate_collect.hpp"
#include "integrate/integrate_iter.hpp"

namespace quadints {
template <typename QuadRule1, typename QuadRule2 = QuadRule1, typename Domain1, typename Domain2, typename Func,
          typename Scalar = decltype(*(QuadRule1::weights.begin()))>
    requires quadrature_rule<QuadRule1, Domain1, Scalar> && quadrature_rule<QuadRule2, Domain2, Scalar> &&
             integrable2<Func, typename Domain1::point_type, typename Domain2::point_type, Scalar>
constexpr auto integrate2(Func&& f, const Domain1& cell1, const Domain2& cell2) {
    return detail::integrate2_collect<QuadRule1, QuadRule2>(std::forward<Func>(f), cell1, cell2);
}

template <typename QuadRule, typename Domain, typename Func, typename Scalar = decltype(*(QuadRule::weights.begin()))>
    requires quadrature_rule<QuadRule, Domain, Scalar> && integrable<Func, typename Domain::point_type, Scalar>
constexpr auto integrate(Func&& f, const Domain& cell) -> std::invoke_result_t<Func, typename Domain::point_type> {
    return detail::integrate_collect<QuadRule>(std::forward<Func>(f), cell);
}

}  // namespace quadints
#endif  // INTERFACE_HPP
