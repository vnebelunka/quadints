#ifndef BARYCENTRIC_PRODUCT_HPP
#define BARYCENTRIC_PRODUCT_HPP

#include <array>
#include <cmath>

namespace quadints {

namespace utils {
// Реализация тензорного произведение
template<typename t1, typename t2, std::size_t ... Idx1, std::size_t ... Idx2>
constexpr decltype(auto) impl_tensor_product(std::array<t1, sizeof...(Idx1)> a1, std::array<t2, sizeof...(Idx2)> a2,
                                 std::index_sequence<Idx1...>, std::index_sequence<Idx2...>) {
  using ret_type = decltype(std::declval<t1>() * std::declval<t2>());
  std::array<std::array<ret_type, sizeof...(Idx1)>, sizeof...(Idx2)> result{
    [&a1](auto&& value) { return std::array{a1[Idx1] * value ...}; }(a2[Idx2]) ...
};
  return result;
}

// Интерфейс тензорного произведения
template<typename t1, typename t2, std::size_t N, std::size_t M>
constexpr decltype(auto) tensor_product(std::array<t1, N> a1, std::array<t2, M> a2) {
  return impl_tensor_product(a1, a2, std::make_index_sequence<N>{}, std::make_index_sequence<M>{});
}

// Реализация декартового произведение
template<typename t1, typename t2, std::size_t ... Idx1, std::size_t ... Idx2>
constexpr decltype(auto) impl_cartesian_product(std::array<t1, sizeof...(Idx1)> a1, std::array<t2, sizeof...(Idx2)> a2,
                                 std::index_sequence<Idx1...>, std::index_sequence<Idx2...>) {
  std::array<std::array<std::pair<t1, t2>, sizeof...(Idx1)>, sizeof...(Idx2)> result{
    [&a1](auto&& value) { return std::array{std::make_pair(a1[Idx1], value) ...}; }(a2[Idx2]) ...
};
  return result;
}

// Интерфейс декартового произведения
template<typename t1, typename t2, std::size_t N, std::size_t M>
constexpr decltype(auto) carteian_product(std::array<t1, N> a1, std::array<t2, M> a2) {
  return impl_cartesian_product(a1, a2, std::make_index_sequence<N>{}, std::make_index_sequence<M>{});
}

// Приведение к одномерному массиву
// p.s. тут не общо написано, потому что не обобщается на произвольное
// количество уровней вложенности массива
//
// Обобщение будет делаться по шаблонной рекурсии на основе этих функций для 2D
template<typename T, std::size_t N, std::size_t M>
constexpr decltype(auto) flatten(std::array<std::array<T, N>, M> arr) {
  std::array<T, N * M> result{};
  for (std::size_t i = 0; i < M; i++)
    for (std::size_t j = 0; j < N; j++) {
      result[i * N + j] = arr[i][j];
    }
  return result;
}

} // namespace utils

template <typename QuadratureX, typename QuadratureY,
          typename Scalar =
              std::remove_cvref_t<decltype(QuadratureX::weights[0])>>
struct QuadratureCartesianProduct {
  using BarycentricPointX = std::remove_cvref_t<decltype(QuadratureX::points[0])>;
  using BarycentricPointY = std::remove_cvref_t<decltype(QuadratureY::points[0])>;
  static constexpr size_t dim_x = QuadratureX::points[0].coords.size();
  static constexpr size_t dim_y = QuadratureY::points[0].coords.size();
};

}


#endif // BARYCENTRIC_PRODUCT_HPP
