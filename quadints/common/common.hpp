#ifndef QUADINTS_COMMON_HPP
#define QUADINTS_COMMON_HPP

#include "quadints/concepts.hpp"

namespace quadints {
template <typename T, typename scalar>
    requires has_magnitude<T, scalar>
auto magnitude(T t) -> scalar {
    if constexpr (has_norm<T, scalar>) {
        return norm(t);
    } else {
        return std::abs(t);
    }
}

template <typename Func, typename Domain>
using return_type = std::invoke_result_t<Func, typename Domain::point_type>;
template <typename Func, typename DomainX, typename DomainY>
using return_type_2d = std::invoke_result_t<Func, typename DomainX::point_type, typename DomainY::point_type>;

template <typename Func, typename DomainPoint, std::size_t N>
concept has_batch_func = requires(Func&& f, const std::array<DomainPoint, N>& pts) {
    { f(pts) } -> std::convertible_to<std::array<std::invoke_result_t<Func, DomainPoint>, N>>;
};

template <typename Func, typename DomainPointX, typename DomainPointY, std::size_t N>
concept has_batch_func2 =
    requires(Func&& f, const std::array<DomainPointX, N>& ptsX, const std::array<DomainPointY, N>& ptsY) {
        { f(ptsX, ptsY) } -> std::convertible_to<std::array<std::invoke_result_t<Func, DomainPointX, DomainPointY>, N>>;
    };

}  // namespace quadints

#endif  // QUADINTS_COMMON_HPP
