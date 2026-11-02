#ifndef QUADINTS_LINALG_QUADRATIC_FORM
#define QUADINTS_LINALG_QUADRATIC_FORM
#include <type_traits>
#if defined(__AVX512F__) || defined(HAS_AVX512)
#include <immintrin.h>
#endif  // defined(__AVX512F__) || defined(HAS_AVX512)
#include <array>
#include <complex>
#include <cstdlib>
namespace quadints::linalg {
template <typename QuadRuleX, typename QuadRuleY, typename ReturnType>
inline ReturnType quadratic_form_fallback(const std::array<ReturnType, QuadRuleX::n_points * QuadRuleY::n_points>& A) {
    constexpr size_t N = QuadRuleX::n_points;
    constexpr size_t M = QuadRuleY::n_points;
    ReturnType res{};
    for (size_t col = 0; col < M; ++col) {
        ReturnType resj = 0;
        for (size_t row = 0; row < N; ++row) {
            resj += A[col * N + row] * QuadRuleX::weights[row];
        }
        res += resj * QuadRuleY::weights[col];
    }
    return res;
}

#if defined(__AVX512F__) || defined(HAS_AVX512)
template <typename QuadRuleX, typename QuadRuleY>
inline double quadratic_form_double(const std::array<double, QuadRuleX::n_points * QuadRuleY::n_points>& A) {
    constexpr size_t N = QuadRuleX::n_points;
    constexpr size_t M = QuadRuleY::n_points;
    if constexpr (N < 8 || N * M < 32) {
        return quadratic_form_fallback<QuadRuleX, QuadRuleY>(A);
    } else {
        constexpr size_t VEC_SIZE = 8;
        constexpr size_t COL_UNROLL = 4;

        const double* ptr_x = QuadRuleX::weights.data();
        const double* ptr_A = A.data();
        const double* ptr_y = QuadRuleY::weights.data();

        // Единый векторный аккумулятор для всей квадратичной формы
        __m512d acc[COL_UNROLL];
        for (int k = 0; k < COL_UNROLL; ++k) {
            acc[k] = _mm512_setzero_pd();
        }
        constexpr size_t n_remainder = N % VEC_SIZE;
        constexpr size_t n_vec_end = (N / VEC_SIZE) * VEC_SIZE;
        constexpr __mmask8 tail_mask = (n_remainder > 0) ? ((1u << n_remainder) - 1) : 0;

        size_t j = 0;

        // 1. Основной цикл: обработка пачек по 4 столбца
        for (; j + COL_UNROLL <= M; j += COL_UNROLL) {
            const double* ptr_Ak[COL_UNROLL];
            for (int i = 0; i < COL_UNROLL; ++i) {
                ptr_Ak[i] = ptr_A + (j + i) * N;
            }
            __m512d sum[COL_UNROLL];
            for (int i = 0; i < COL_UNROLL; ++i) {
                sum[i] = _mm512_setzero_pd();
            }
            size_t i = 0;
            for (; i < n_vec_end; i += VEC_SIZE) {
                __m512d vx = _mm512_loadu_pd(ptr_x + i);
                for (int k = 0; k < COL_UNROLL; ++k) {
                    sum[k] = _mm512_fmadd_pd(_mm512_loadu_pd(ptr_Ak[k] + i), vx, sum[k]);
                }
            }

            if constexpr (n_remainder > 0) {
                __m512d vx = _mm512_maskz_loadu_pd(tail_mask, ptr_x + n_vec_end);

                for (int k = 0; k < COL_UNROLL; ++k) {
                    sum[k] = _mm512_fmadd_pd(_mm512_maskz_loadu_pd(tail_mask, ptr_Ak[k] + n_vec_end), vx, sum[k]);
                }
            }
            for (int k = 0; k < COL_UNROLL; ++k) {
                acc[k] = _mm512_fmadd_pd(sum[k], _mm512_set1_pd(ptr_y[j + k]), acc[k]);
            }
        }

        if constexpr (M % 4 != 0) {
            for (; j < M; ++j) {
                const double* ptr_A0 = ptr_A + j * N;
                __m512d sum0 = _mm512_setzero_pd();

                size_t i = 0;
                for (; i < n_vec_end; i += VEC_SIZE) {
                    __m512d vx = _mm512_loadu_pd(ptr_x + i);
                    __m512d va0 = _mm512_loadu_pd(ptr_A0 + i);
                    sum0 = _mm512_fmadd_pd(va0, vx, sum0);
                }
                if constexpr (n_remainder > 0) {
                    __m512d vx = _mm512_maskz_loadu_pd(tail_mask, ptr_x + n_vec_end);
                    __m512d va0 = _mm512_maskz_loadu_pd(tail_mask, ptr_A0 + n_vec_end);
                    sum0 = _mm512_fmadd_pd(va0, vx, sum0);
                }

                __m512d vy0 = _mm512_set1_pd(ptr_y[j]);
                acc[0] = _mm512_fmadd_pd(sum0, vy0, acc[0]);
            }
        }
        __m512d acc_vec = _mm512_add_pd(_mm512_add_pd(acc[0], acc[1]), _mm512_add_pd(acc[2], acc[3]));
        return _mm512_reduce_add_pd(acc_vec);
    }
}

template <typename QuadRule>
inline double quadratic_form_double(const std::array<double, QuadRule::n_points * QuadRule::n_points>& A) {
    constexpr size_t N = QuadRule::n_points;
    if constexpr (N < 8) {
        return quadratic_form_fallback<QuadRule, QuadRule>(A);
    } else {
        constexpr size_t VEC_SIZE = 8;
        constexpr size_t COL_UNROLL = 4;

        const double* ptr_x = QuadRule::weights.data();
        const double* ptr_A = A.data();

        __m512d acc[COL_UNROLL];
        for (size_t k = 0; k < COL_UNROLL; ++k) {
            acc[k] = _mm512_setzero_pd();
        }

        constexpr size_t n_remainder = N % VEC_SIZE;
        constexpr size_t n_vec_end = (N / VEC_SIZE) * VEC_SIZE;
        constexpr __mmask8 tail_mask = (n_remainder > 0) ? ((1u << n_remainder) - 1) : 0;

        size_t j = 0;

        for (; j + COL_UNROLL <= N; j += COL_UNROLL) {
            const double* ptr_Ak[COL_UNROLL];
            __m512d sum[COL_UNROLL];

            for (size_t k = 0; k < COL_UNROLL; ++k) {
                ptr_Ak[k] = ptr_A + (j + k) * N;
                sum[k] = _mm512_setzero_pd();
            }

            size_t i = 0;
            for (; i < n_vec_end; i += VEC_SIZE) {
                __m512d vx = _mm512_loadu_pd(ptr_x + i);
                for (size_t k = 0; k < COL_UNROLL; ++k) {
                    sum[k] = _mm512_fmadd_pd(_mm512_loadu_pd(ptr_Ak[k] + i), vx, sum[k]);
                }
            }

            if constexpr (n_remainder > 0) {
                __m512d vx = _mm512_maskz_loadu_pd(tail_mask, ptr_x + n_vec_end);
                for (size_t k = 0; k < COL_UNROLL; ++k) {
                    sum[k] = _mm512_fmadd_pd(_mm512_maskz_loadu_pd(tail_mask, ptr_Ak[k] + n_vec_end), vx, sum[k]);
                }
            }

            for (size_t k = 0; k < COL_UNROLL; ++k) {
                acc[k] = _mm512_fmadd_pd(sum[k], _mm512_set1_pd(ptr_x[j + k]), acc[k]);
            }
        }

        if constexpr (N % 4 != 0) {
            for (; j < N; ++j) {
                const double* ptr_A0 = ptr_A + j * N;
                __m512d sum0 = _mm512_setzero_pd();

                size_t i = 0;
                for (; i < n_vec_end; i += VEC_SIZE) {
                    __m512d vx = _mm512_loadu_pd(ptr_x + i);
                    __m512d va0 = _mm512_loadu_pd(ptr_A0 + i);
                    sum0 = _mm512_fmadd_pd(va0, vx, sum0);
                }

                if constexpr (n_remainder > 0) {
                    __m512d vx = _mm512_maskz_loadu_pd(tail_mask, ptr_x + n_vec_end);
                    __m512d va0 = _mm512_maskz_loadu_pd(tail_mask, ptr_A0 + n_vec_end);
                    sum0 = _mm512_fmadd_pd(va0, vx, sum0);
                }

                acc[0] = _mm512_fmadd_pd(sum0, _mm512_set1_pd(ptr_x[j]), acc[0]);
            }
        }

        __m512d acc_vec = _mm512_add_pd(_mm512_add_pd(acc[0], acc[1]), _mm512_add_pd(acc[2], acc[3]));
        return _mm512_reduce_add_pd(acc_vec);
    }
}

#endif  // defined(__AVX512F__) || defined(HAS_AVX512)

template <typename QuadRuleX, typename QuadRuleY, typename ReturnType>
inline ReturnType quadratic_form(const std::array<ReturnType, QuadRuleX::n_points * QuadRuleY::n_points>& A) {
#if defined(__AVX512F__) || defined(HAS_AVX512)
    if constexpr (std::is_same_v<ReturnType, double>) {
        if constexpr (std::is_same_v<QuadRuleX, QuadRuleY>) {
            return quadratic_form_double<QuadRuleX>(A);
        }
        return quadratic_form_double<QuadRuleX, QuadRuleY>(A);
    }
#endif
    return quadratic_form_fallback<QuadRuleX, QuadRuleY>(A);
}

}  // namespace quadints::linalg
#endif  // QUADINTS_LINALG_QUADRATIC_FORM
