#ifndef QUADINTS_LINALG_DOT_HPP
#define QUADINTS_LINALG_DOT_HPP
#if defined(__AVX512F__) || defined(HAS_AVX512)
#include <immintrin.h>
#endif  // defined(__AVX512F__) || defined(HAS_AVX512)
#include <array>
#include <complex>
#include <cstdlib>
template <size_t N, typename ReturnType, typename ScalarType>
inline ReturnType dot_product_fallback(const std::array<ReturnType, N>& func_arr,
                                       const std::array<ScalarType, N>& weights) {
    ReturnType res{};
    for (size_t i = 0; i < N; ++i) {
        res += func_arr[i] * weights[i];
    }
    return res;
}

namespace quadints::linalg {
template <size_t N, typename ReturnType, typename ScalarType>
inline ReturnType dot_product(const std::array<ReturnType, N>& func_arr, const std::array<ScalarType, N>& weights) {
    return dot_product_fallback(func_arr, weights);
}

#if defined(__AVX512F__) || defined(HAS_AVX512)
template <size_t N>
inline double dot_product(const std::array<double, N>& func_arr, const std::array<double, N>& weights) {
    if constexpr (N <= 8) {
        return dot_product_fallback(func_arr, weights);
    }
    constexpr size_t VEC_SIZE = 8;
    constexpr size_t UNROLL = 4;
    constexpr size_t BLOCK_SIZE = VEC_SIZE * UNROLL;  // 32 элемента (256 байт)

    const double* ptr_a = func_arr.data();
    const double* ptr_b = weights.data();

    // 4 независимых аккумулятора для параллельной загрузки FMA-конвейеров
    __m512d sum0 = _mm512_setzero_pd();
    __m512d sum1 = _mm512_setzero_pd();
    __m512d sum2 = _mm512_setzero_pd();
    __m512d sum3 = _mm512_setzero_pd();

    size_t i = 0;

    // 1. Основной развернутый цикл: 32 элемента за итерацию
    for (; i + BLOCK_SIZE <= N; i += BLOCK_SIZE) {
        __m512d va0 = _mm512_loadu_pd(ptr_a + i + 0);
        __m512d vb0 = _mm512_loadu_pd(ptr_b + i + 0);
        sum0 = _mm512_fmadd_pd(va0, vb0, sum0);

        __m512d va1 = _mm512_loadu_pd(ptr_a + i + 8);
        __m512d vb1 = _mm512_loadu_pd(ptr_b + i + 8);
        sum1 = _mm512_fmadd_pd(va1, vb1, sum1);

        __m512d va2 = _mm512_loadu_pd(ptr_a + i + 16);
        __m512d vb2 = _mm512_loadu_pd(ptr_b + i + 16);
        sum2 = _mm512_fmadd_pd(va2, vb2, sum2);

        __m512d va3 = _mm512_loadu_pd(ptr_a + i + 24);
        __m512d vb3 = _mm512_loadu_pd(ptr_b + i + 24);
        sum3 = _mm512_fmadd_pd(va3, vb3, sum3);
    }

    // Сворачиваем 4 аккумулятора в один
    sum0 = _mm512_add_pd(sum0, sum1);
    sum2 = _mm512_add_pd(sum2, sum3);
    sum0 = _mm512_add_pd(sum0, sum2);

    // 2. Вторичный цикл по 8 элементов (для промежуточных остатков)
    for (; i + VEC_SIZE <= N; i += VEC_SIZE) {
        __m512d va = _mm512_loadu_pd(ptr_a + i);
        __m512d vb = _mm512_loadu_pd(ptr_b + i);
        sum0 = _mm512_fmadd_pd(va, vb, sum0);
    }

    // 3. Обработка хвоста (< 8 элементов)
    // if constexpr гарантирует, что для N кратном 8 этот блок вообще не попадет в бинарник
    if constexpr (N % VEC_SIZE != 0) {
        if (i < N) {
            size_t remainder = N - i;
            __mmask8 mask = (1u << remainder) - 1;

            __m512d va = _mm512_maskz_loadu_pd(mask, ptr_a + i);
            __m512d vb = _mm512_maskz_loadu_pd(mask, ptr_b + i);

            sum0 = _mm512_fmadd_pd(va, vb, sum0);
        }
    }

    return _mm512_reduce_add_pd(sum0);
}
template <size_t N>
inline std::complex<double> dot_product(const std::array<std::complex<double>, N>& a, const std::array<double, N>& b) {
    // 1. Для малых N — наивная скалярная реализация (без оверхеда на маски и permute)
    if constexpr (N < 8) {
        return dot_product_fallback(a, b);
    } else {
        constexpr size_t VEC_SIZE = 8;  // 8 комплексных чисел = 16 double (2 x AVX-512)
        constexpr size_t UNROLL = 2;    // 16 комплексных чисел за итерацию
        constexpr size_t BLOCK_SIZE = VEC_SIZE * UNROLL;

        const double* ptr_a = reinterpret_cast<const double*>(a.data());
        const double* ptr_b = b.data();

        // Маски индексов для дублирования элементов b_i:
        // idx0 превращает [b0..b7] в [b0, b0, b1, b1, b2, b2, b3, b3]
        // idx1 превращает [b0..b7] в [b4, b4, b5, b5, b6, b6, b7, b7]
        const __m512i idx0 = _mm512_set_epi64(3, 3, 2, 2, 1, 1, 0, 0);
        const __m512i idx1 = _mm512_set_epi64(7, 7, 6, 6, 5, 5, 4, 4);

        __m512d sum0 = _mm512_setzero_pd();
        __m512d sum1 = _mm512_setzero_pd();
        __m512d sum2 = _mm512_setzero_pd();
        __m512d sum3 = _mm512_setzero_pd();

        size_t i = 0;

        // 2. Развернутый основной цикл (по 16 комплексных элементов)
        for (; i + BLOCK_SIZE <= N; i += BLOCK_SIZE) {
            // Пачка 0 (элементы i ... i+7)
            __m512d va0 = _mm512_loadu_pd(ptr_a + 2 * i + 0);  // Re0, Im0, Re1, Im1, Re2, Im2, Re3, Im3
            __m512d va1 = _mm512_loadu_pd(ptr_a + 2 * i + 8);  // Re4, Im4, Re5, Im5, Re6, Im6, Re7, Im7
            __m512d vb_raw0 = _mm512_loadu_pd(ptr_b + i);      // b0..b7

            __m512d vb0 = _mm512_permutexvar_pd(idx0, vb_raw0);
            __m512d vb1 = _mm512_permutexvar_pd(idx1, vb_raw0);

            sum0 = _mm512_fmadd_pd(va0, vb0, sum0);
            sum1 = _mm512_fmadd_pd(va1, vb1, sum1);

            // Пачка 1 (элементы i+8 ... i+15)
            __m512d va2 = _mm512_loadu_pd(ptr_a + 2 * i + 16);
            __m512d va3 = _mm512_loadu_pd(ptr_a + 2 * i + 24);
            __m512d vb_raw1 = _mm512_loadu_pd(ptr_b + i + 8);

            __m512d vb2 = _mm512_permutexvar_pd(idx0, vb_raw1);
            __m512d vb3 = _mm512_permutexvar_pd(idx1, vb_raw1);

            sum2 = _mm512_fmadd_pd(va2, vb2, sum2);
            sum3 = _mm512_fmadd_pd(va3, vb3, sum3);
        }

        // 3. Вторичный цикл (по 8 комплексных элементов)
        for (; i + VEC_SIZE <= N; i += VEC_SIZE) {
            __m512d va0 = _mm512_loadu_pd(ptr_a + 2 * i + 0);
            __m512d va1 = _mm512_loadu_pd(ptr_a + 2 * i + 8);
            __m512d vb_raw = _mm512_loadu_pd(ptr_b + i);

            __m512d vb0 = _mm512_permutexvar_pd(idx0, vb_raw);
            __m512d vb1 = _mm512_permutexvar_pd(idx1, vb_raw);

            sum0 = _mm512_fmadd_pd(va0, vb0, sum0);
            sum1 = _mm512_fmadd_pd(va1, vb1, sum1);
        }

        // Сворачиваем векторные аккумуляторы
        sum0 = _mm512_add_pd(sum0, sum1);
        sum2 = _mm512_add_pd(sum2, sum3);
        sum0 = _mm512_add_pd(sum0, sum2);

        // 4. Покомпонентная векторная редукция (сворачивание 512b -> 256b -> 128b)
        __m256d lo = _mm512_extractf64x4_pd(sum0, 0);
        __m256d hi = _mm512_extractf64x4_pd(sum0, 1);
        __m256d sum256 = _mm256_add_pd(lo, hi);

        __m128d lo128 = _mm256_extractf128_pd(sum256, 0);
        __m128d hi128 = _mm256_extractf128_pd(sum256, 1);
        __m128d res128 = _mm_add_pd(lo128, hi128);
        alignas(16) double res_arr[2];
        _mm_store_pd(res_arr, res128);
        std::complex<double> total_res{res_arr[0], res_arr[1]};

        // 5. Обработка оставшегося хвоста (< 8 элементов)
        if constexpr (N % VEC_SIZE != 0) {
            for (; i < N; ++i) {
                total_res += a[i] * b[i];
            }
        }

        return total_res;
    }
}

#endif  // defined(__AVX512F__) || defined(HAS_AVX512)

}  // namespace quadints::linalg

#endif  // QUADINTS_LINALG_DOT_HPP
