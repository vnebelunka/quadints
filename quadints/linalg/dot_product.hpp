#ifndef QUADINTS_LINALG_DOT_HPP
#define QUADINTS_LINALG_DOT_HPP
#if defined(__AVX512F__) || defined(HAS_AVX512)
#include <immintrin.h>
#endif  // defined(__AVX512F__) || defined(HAS_AVX512)
#include <array>
#include <complex>
#include <cstdlib>

#include "quadints/common/aligned_array.hpp"

namespace quadints::linalg::detail {

template <typename T>
struct static_extent;

template <typename T, std::size_t N>
struct static_extent<std::array<T, N>> : std::integral_constant<std::size_t, N> {};

template <typename T, std::size_t N, std::size_t Alignment>
struct static_extent<common::aligned_array<T, N, Alignment>> : std::integral_constant<std::size_t, N> {};

template <typename T>
inline constexpr std::size_t static_extent_v = static_extent<std::remove_cvref_t<T>>::value;

template <typename T>
concept fixed_contiguous_array =
    requires {
        typename std::remove_cvref_t<T>::value_type;
        static_extent_v<T>;
    } && std::ranges::contiguous_range<const std::remove_cvref_t<T>> &&
    std::ranges::sized_range<const std::remove_cvref_t<T>>;

template <typename A, typename B>
concept compatible_arrays =
    fixed_contiguous_array<A> && fixed_contiguous_array<B> && (static_extent_v<A> == static_extent_v<B>);

}  // namespace quadints::linalg::detail

namespace quadints::linalg {

template <typename FuncArray, typename WeightArray>
    requires detail::compatible_arrays<FuncArray, WeightArray>
inline auto dot_product_fallback(const FuncArray& func_arr, const WeightArray& weights) {
    using ReturnType = std::remove_cvref_t<decltype(func_arr[0] * weights[0])>;

    constexpr std::size_t N = detail::static_extent_v<FuncArray>;

    ReturnType res{};

    for (std::size_t i = 0; i < N; ++i) {
        res += func_arr[i] * weights[i];
    }

    return res;
}

}  // namespace quadints::linalg

namespace quadints::linalg::detail {

template <typename Array>
inline constexpr std::size_t declared_alignment = [] {
    using T = std::remove_cvref_t<Array>;

    if constexpr (requires { T::alignment; }) {
        return T::alignment;
    } else {
        return std::size_t{0};
    }
}();

#if defined(__AVX512F__) || defined(HAS_AVX512)

template <typename Array>
inline __m512d load_pd(const double* ptr) noexcept {
    if constexpr (declared_alignment<Array> >= 64) {
        return _mm512_load_pd(ptr);
    } else {
        return _mm512_loadu_pd(ptr);
    }
}

template <typename Array>
inline __m512d maskz_load_pd(__mmask8 mask, const double* ptr) noexcept {
    if constexpr (declared_alignment<Array> >= 64) {
        return _mm512_maskz_load_pd(mask, ptr);
    } else {
        return _mm512_maskz_loadu_pd(mask, ptr);
    }
}

#endif

}  // namespace quadints::linalg::detail

#if defined(__AVX512F__) || defined(HAS_AVX512)

namespace quadints::linalg::detail {

template <typename FuncArray, typename WeightArray>
    requires compatible_arrays<FuncArray, WeightArray> && std::same_as<typename FuncArray::value_type, double> &&
             std::same_as<typename WeightArray::value_type, double>
inline double dot_product_double_avx512(const FuncArray& func_arr, const WeightArray& weights) {
    constexpr std::size_t N = static_extent_v<FuncArray>;

    if constexpr (N <= 8) {
        return dot_product_fallback(func_arr, weights);
    } else {
        constexpr std::size_t VEC_SIZE = 8;
        constexpr std::size_t UNROLL = 4;
        constexpr std::size_t BLOCK_SIZE = VEC_SIZE * UNROLL;

        const double* ptr_a = func_arr.data();
        const double* ptr_b = weights.data();

        __m512d sum0 = _mm512_setzero_pd();
        __m512d sum1 = _mm512_setzero_pd();
        __m512d sum2 = _mm512_setzero_pd();
        __m512d sum3 = _mm512_setzero_pd();

        std::size_t i = 0;

        for (; i + BLOCK_SIZE <= N; i += BLOCK_SIZE) {
            sum0 = _mm512_fmadd_pd(load_pd<FuncArray>(ptr_a + i), load_pd<WeightArray>(ptr_b + i), sum0);

            sum1 = _mm512_fmadd_pd(load_pd<FuncArray>(ptr_a + i + 8), load_pd<WeightArray>(ptr_b + i + 8), sum1);

            sum2 = _mm512_fmadd_pd(load_pd<FuncArray>(ptr_a + i + 16), load_pd<WeightArray>(ptr_b + i + 16), sum2);

            sum3 = _mm512_fmadd_pd(load_pd<FuncArray>(ptr_a + i + 24), load_pd<WeightArray>(ptr_b + i + 24), sum3);
        }

        sum0 = _mm512_add_pd(sum0, sum1);
        sum2 = _mm512_add_pd(sum2, sum3);
        sum0 = _mm512_add_pd(sum0, sum2);

        for (; i + VEC_SIZE <= N; i += VEC_SIZE) {
            sum0 = _mm512_fmadd_pd(load_pd<FuncArray>(ptr_a + i), load_pd<WeightArray>(ptr_b + i), sum0);
        }

        if constexpr (N % VEC_SIZE != 0) {
            if (i < N) {
                const auto remainder = N - i;
                const __mmask8 mask = static_cast<__mmask8>((1u << remainder) - 1u);

                const __m512d va = maskz_load_pd<FuncArray>(mask, ptr_a + i);

                const __m512d vb = maskz_load_pd<WeightArray>(mask, ptr_b + i);

                sum0 = _mm512_fmadd_pd(va, vb, sum0);
            }
        }

        return _mm512_reduce_add_pd(sum0);
    }
}

template <typename FuncArray, typename WeightArray>
    requires compatible_arrays<FuncArray, WeightArray>
inline std::complex<double> dot_product_complex_avx512(const FuncArray& a, const WeightArray& weights) {
    constexpr std::size_t N = static_extent_v<FuncArray>;
    static_assert(sizeof(std::complex<double>) == 2 * sizeof(double));
    if constexpr (N < 8) {
        return dot_product_fallback(a, weights);
    } else {
        constexpr std::size_t VEC_SIZE = 8;
        constexpr std::size_t UNROLL = 2;
        constexpr std::size_t BLOCK_SIZE = VEC_SIZE * UNROLL;

        const double* ptr_a = reinterpret_cast<const double*>(a.data());
        const double* ptr_b = weights.data();
        // [w0, w0, w1, w1, w2, w2, w3, w3]
        const __m512i idx_lo = _mm512_set_epi64(3, 3, 2, 2, 1, 1, 0, 0);
        // [w4, w4, w5, w5, w6, w6, w7, w7]
        const __m512i idx_hi = _mm512_set_epi64(7, 7, 6, 6, 5, 5, 4, 4);
        __m512d sum0 = _mm512_setzero_pd();
        __m512d sum1 = _mm512_setzero_pd();
        __m512d sum2 = _mm512_setzero_pd();
        __m512d sum3 = _mm512_setzero_pd();
        std::size_t i = 0;
        for (; i + BLOCK_SIZE <= N; i += BLOCK_SIZE) {
            const __m512d w0 = load_pd<WeightArray>(ptr_b + i);
            const __m512d w0_lo = _mm512_permutexvar_pd(idx_lo, w0);
            const __m512d w0_hi = _mm512_permutexvar_pd(idx_hi, w0);
            const __m512d a0 = load_pd<FuncArray>(ptr_a + 2 * i);
            const __m512d a1 = load_pd<FuncArray>(ptr_a + 2 * i + 8);
            sum0 = _mm512_fmadd_pd(a0, w0_lo, sum0);
            sum1 = _mm512_fmadd_pd(a1, w0_hi, sum1);
            const __m512d w1 = load_pd<WeightArray>(ptr_b + i + 8);
            const __m512d w1_lo = _mm512_permutexvar_pd(idx_lo, w1);
            const __m512d w1_hi = _mm512_permutexvar_pd(idx_hi, w1);
            const __m512d a2 = load_pd<FuncArray>(ptr_a + 2 * i + 16);
            const __m512d a3 = load_pd<FuncArray>(ptr_a + 2 * i + 24);
            sum2 = _mm512_fmadd_pd(a2, w1_lo, sum2);
            sum3 = _mm512_fmadd_pd(a3, w1_hi, sum3);
        }
        for (; i + VEC_SIZE <= N; i += VEC_SIZE) {
            const __m512d w = load_pd<WeightArray>(ptr_b + i);
            const __m512d w_lo = _mm512_permutexvar_pd(idx_lo, w);
            const __m512d w_hi = _mm512_permutexvar_pd(idx_hi, w);
            const __m512d a0 = load_pd<FuncArray>(ptr_a + 2 * i);
            const __m512d a1 = load_pd<FuncArray>(ptr_a + 2 * i + 8);
            sum0 = _mm512_fmadd_pd(a0, w_lo, sum0);
            sum1 = _mm512_fmadd_pd(a1, w_hi, sum1);
        }
        sum0 = _mm512_add_pd(sum0, sum1);
        sum2 = _mm512_add_pd(sum2, sum3);
        sum0 = _mm512_add_pd(sum0, sum2);
        const double real = _mm512_reduce_add_pd(_mm512_maskz_mov_pd(0x55, sum0));
        const double imag = _mm512_reduce_add_pd(_mm512_maskz_mov_pd(0xAA, sum0));
        std::complex<double> result{real, imag};
        for (; i < N; ++i) {
            result += a[i] * weights[i];
        }
        return result;
    }
}
}  // namespace quadints::linalg::detail

#endif

namespace quadints::linalg {

template <typename FuncArray, typename WeightArray>
    requires detail::compatible_arrays<FuncArray, WeightArray>
inline auto dot_product(const FuncArray& func_arr, const WeightArray& weights) {
#if defined(__AVX512F__) || defined(HAS_AVX512)
    using FuncType = typename std::remove_cvref_t<FuncArray>::value_type;
    using WeightType = typename std::remove_cvref_t<WeightArray>::value_type;
    if constexpr (std::same_as<FuncType, double> && std::same_as<WeightType, double>) {
        return detail::dot_product_double_avx512(func_arr, weights);
    } else if constexpr (std::same_as<FuncType, std::complex<double>> && std::same_as<WeightType, double>) {
        return detail::dot_product_complex_avx512(func_arr, weights);
    } else {
        return dot_product_fallback(func_arr, weights);
    }
#else
    return dot_product_fallback(func_arr, weights);
#endif
}

}  // namespace quadints::linalg

#endif  // QUADINTS_LINALG_DOT_HPP
