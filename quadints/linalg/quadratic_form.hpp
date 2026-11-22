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
template <typename QuadRuleX, typename QuadRuleY, typename MatrixArray>
inline auto quadratic_form_fallback(const MatrixArray& A) {
    using ReturnType = std::remove_cvref_t<decltype(A[0])>;
    constexpr std::size_t N = QuadRuleX::n_points;
    constexpr std::size_t M = QuadRuleY::n_points;
    ReturnType result{};
    for (std::size_t col = 0; col < M; ++col) {
        ReturnType sum{};
        for (std::size_t row = 0; row < N; ++row) {
            sum += A[col * N + row] * QuadRuleX::weights[row];
        }
        result += sum * QuadRuleY::weights[col];
    }
    return result;
}

#if defined(__AVX512F__) || defined(HAS_AVX512)

template <typename Array>
inline constexpr std::size_t declared_alignment = [] {
    using T = std::remove_cvref_t<Array>;

    if constexpr (requires { T::alignment; }) {
        return T::alignment;
    } else {
        return std::size_t{0};
    }
}();

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

#if defined(__AVX512F__) || defined(HAS_AVX512)

template <typename QuadRuleX, typename QuadRuleY, typename MatrixArray>
inline double quadratic_form_double(const MatrixArray& A) {
    constexpr size_t N = QuadRuleX::n_points;
    constexpr size_t M = QuadRuleY::n_points;
    if constexpr (N % 8 != 0 || N * M < 32) {
        return quadratic_form_fallback<QuadRuleX, QuadRuleY>(A);
    } else {
        constexpr size_t VEC_SIZE = 8;
        constexpr size_t COL_UNROLL = 4;
        const double* ptr_x = QuadRuleX::weights.data();
        const double* ptr_y = QuadRuleY::weights.data();
        const double* ptr_A = A.data();
        __m512d acc[COL_UNROLL];
        for (size_t k = 0; k < COL_UNROLL; ++k) {
            acc[k] = _mm512_setzero_pd();
        }
        size_t j = 0;
        // Main loop: 4 columns.
        for (; j + COL_UNROLL <= M; j += COL_UNROLL) {
            const double* ptr_Ak[COL_UNROLL];
            __m512d sum[COL_UNROLL];
            for (size_t k = 0; k < COL_UNROLL; ++k) {
                ptr_Ak[k] = ptr_A + (j + k) * N;
                sum[k] = _mm512_setzero_pd();
            }
            for (size_t i = 0; i < N; i += VEC_SIZE) {
                const __m512d vx = load_pd<decltype(QuadRuleX::weights)>(ptr_x + i);
                for (size_t k = 0; k < COL_UNROLL; ++k) {
                    const __m512d va = load_pd<MatrixArray>(ptr_Ak[k] + i);
                    sum[k] = _mm512_fmadd_pd(va, vx, sum[k]);
                }
            }
            for (size_t k = 0; k < COL_UNROLL; ++k) {
                acc[k] = _mm512_fmadd_pd(sum[k], _mm512_set1_pd(ptr_y[j + k]), acc[k]);
            }
        }
        for (; j < M; ++j) {
            const double* ptr_A0 = ptr_A + j * N;
            __m512d sum = _mm512_setzero_pd();
            for (size_t i = 0; i < N; i += VEC_SIZE) {
                const __m512d vx = load_pd<decltype(QuadRuleX::weights)>(ptr_x + i);
                const __m512d va = load_pd<MatrixArray>(ptr_A0 + i);
                sum = _mm512_fmadd_pd(va, vx, sum);
            }
            acc[0] = _mm512_fmadd_pd(sum, _mm512_set1_pd(ptr_y[j]), acc[0]);
        }
        const __m512d acc_vec = _mm512_add_pd(_mm512_add_pd(acc[0], acc[1]), _mm512_add_pd(acc[2], acc[3]));
        return _mm512_reduce_add_pd(acc_vec);
    }
}

template <typename QuadRuleX, typename QuadRuleY, typename MatrixArray>
inline std::complex<double> quadratic_form_complex(const MatrixArray& A) {
    constexpr size_t N = QuadRuleX::n_points;
    constexpr size_t M = QuadRuleY::n_points;
    if constexpr (N % 8 != 0 || N * M < 32) {
        return quadratic_form_fallback<QuadRuleX, QuadRuleY>(A);
    } else {
        constexpr size_t VEC_SIZE = 8;
        constexpr size_t COL_UNROLL = 4;

        const double* ptr_x = QuadRuleX::weights.data();
        const double* ptr_y = QuadRuleY::weights.data();
        const double* ptr_A = reinterpret_cast<const double*>(A.data());
        // Indices for extracting real and imaginary components.
        const __m512i idx_re = _mm512_setr_epi64(0, 2, 4, 6, 8, 10, 12, 14);
        const __m512i idx_im = _mm512_setr_epi64(1, 3, 5, 7, 9, 11, 13, 15);
        __m512d acc_re[COL_UNROLL];
        __m512d acc_im[COL_UNROLL];
        for (size_t k = 0; k < COL_UNROLL; ++k) {
            acc_re[k] = _mm512_setzero_pd();
            acc_im[k] = _mm512_setzero_pd();
        }
        size_t j = 0;
        // Main loop: 4 columns.
        for (; j + COL_UNROLL <= M; j += COL_UNROLL) {
            const double* ptr_Ak[COL_UNROLL];
            __m512d sum_re[COL_UNROLL];
            __m512d sum_im[COL_UNROLL];
            for (size_t k = 0; k < COL_UNROLL; ++k) {
                ptr_Ak[k] = ptr_A + (j + k) * N * 2;
                sum_re[k] = _mm512_setzero_pd();
                sum_im[k] = _mm512_setzero_pd();
            }
            for (size_t i = 0; i < N; i += VEC_SIZE) {
                const __m512d vx = load_pd<decltype(QuadRuleX::weights)>(ptr_x + i);
                for (size_t k = 0; k < COL_UNROLL; ++k) {
                    // 8 complex numbers = 16 doubles.
                    const __m512d v0 = load_pd<MatrixArray>(ptr_Ak[k] + 2 * i);
                    const __m512d v1 = load_pd<MatrixArray>(ptr_Ak[k] + 2 * i + 8);
                    // Extract 8 real and 8 imaginary components.
                    const __m512d va_re = _mm512_permutex2var_pd(v0, idx_re, v1);
                    const __m512d va_im = _mm512_permutex2var_pd(v0, idx_im, v1);
                    // Multiply by real quadrature weights.
                    sum_re[k] = _mm512_fmadd_pd(va_re, vx, sum_re[k]);
                    sum_im[k] = _mm512_fmadd_pd(va_im, vx, sum_im[k]);
                }
            }
            // Multiply column sums by Y weights.
            for (size_t k = 0; k < COL_UNROLL; ++k) {
                const __m512d vy = _mm512_set1_pd(ptr_y[j + k]);
                acc_re[k] = _mm512_fmadd_pd(sum_re[k], vy, acc_re[k]);
                acc_im[k] = _mm512_fmadd_pd(sum_im[k], vy, acc_im[k]);
            }
        }
        for (; j < M; ++j) {
            const double* ptr_A0 = ptr_A + j * N * 2;

            __m512d sum_re = _mm512_setzero_pd();
            __m512d sum_im = _mm512_setzero_pd();

            for (size_t i = 0; i < N; i += VEC_SIZE) {
                const __m512d vx = load_pd<decltype(QuadRuleX::weights)>(ptr_x + i);
                const __m512d v0 = load_pd<MatrixArray>(ptr_A0 + 2 * i);
                const __m512d v1 = load_pd<MatrixArray>(ptr_A0 + 2 * i + 8);
                const __m512d va_re = _mm512_permutex2var_pd(v0, idx_re, v1);
                const __m512d va_im = _mm512_permutex2var_pd(v0, idx_im, v1);
                sum_re = _mm512_fmadd_pd(va_re, vx, sum_re);
                sum_im = _mm512_fmadd_pd(va_im, vx, sum_im);
            }
            const __m512d vy = _mm512_set1_pd(ptr_y[j]);
            acc_re[0] = _mm512_fmadd_pd(sum_re, vy, acc_re[0]);
            acc_im[0] = _mm512_fmadd_pd(sum_im, vy, acc_im[0]);
        }
        const __m512d acc_re_vec =
            _mm512_add_pd(_mm512_add_pd(acc_re[0], acc_re[1]), _mm512_add_pd(acc_re[2], acc_re[3]));
        const __m512d acc_im_vec =
            _mm512_add_pd(_mm512_add_pd(acc_im[0], acc_im[1]), _mm512_add_pd(acc_im[2], acc_im[3]));
        return {_mm512_reduce_add_pd(acc_re_vec), _mm512_reduce_add_pd(acc_im_vec)};
    }
}

#endif

template <typename QuadRuleX, typename QuadRuleY>
inline std::complex<double> quadratic_form_complex(
    const std::array<std::complex<double>, QuadRuleX::n_points * QuadRuleY::n_points>& A) {
    constexpr size_t N = QuadRuleX::n_points;
    constexpr size_t M = QuadRuleY::n_points;

    const double* ptr_x = QuadRuleX::weights.data();
    const double* ptr_y = QuadRuleY::weights.data();

    // Скалярный вариант для небольших размеров.
    if constexpr (N < 8 || N * M < 32) {
        std::complex<double> result{};

        for (size_t j = 0; j < M; ++j) {
            std::complex<double> sum{};

            for (size_t i = 0; i < N; ++i) {
                sum += A[j * N + i] * ptr_x[i];
            }

            result += sum * ptr_y[j];
        }

        return result;
    } else {
        constexpr size_t VEC_SIZE = 8;  // 8 комплексных чисел
        constexpr size_t COL_UNROLL = 4;

        constexpr size_t n_vec_end = (N / VEC_SIZE) * VEC_SIZE;

        // Индексы для разделения real и imag:
        //
        // v0 = [re0, im0, re1, im1, re2, im2, re3, im3]
        // v1 = [re4, im4, re5, im5, re6, im6, re7, im7]
        //
        // re = [re0, re1, re2, re3, re4, re5, re6, re7]
        // im = [im0, im1, im2, im3, im4, im5, im6, im7]

        const __m512i idx_re = _mm512_setr_epi64(0, 2, 4, 6, 8, 10, 12, 14);

        const __m512i idx_im = _mm512_setr_epi64(1, 3, 5, 7, 9, 11, 13, 15);

        // Отдельные аккумуляторы для действительной и мнимой частей.
        __m512d acc_re[COL_UNROLL];
        __m512d acc_im[COL_UNROLL];

        for (size_t k = 0; k < COL_UNROLL; ++k) {
            acc_re[k] = _mm512_setzero_pd();
            acc_im[k] = _mm512_setzero_pd();
        }

        // Накопление скалярных хвостов N % 8.
        double tail_re = 0.0;
        double tail_im = 0.0;

        auto accumulate_tail = [&](const double* ptr_Ak, double wy) {
            double sum_re = 0.0;
            double sum_im = 0.0;

            for (size_t i = n_vec_end; i < N; ++i) {
                sum_re += ptr_Ak[2 * i] * ptr_x[i];
                sum_im += ptr_Ak[2 * i + 1] * ptr_x[i];
            }

            tail_re += sum_re * wy;
            tail_im += sum_im * wy;
        };

        size_t j = 0;

        // Основной цикл: обработка по 4 столбца.
        for (; j + COL_UNROLL <= M; j += COL_UNROLL) {
            const double* ptr_Ak[COL_UNROLL];

            for (size_t k = 0; k < COL_UNROLL; ++k) {
                // std::complex<double> хранится как [real, imag].
                ptr_Ak[k] = reinterpret_cast<const double*>(A.data() + (j + k) * N);
            }

            __m512d sum_re[COL_UNROLL];
            __m512d sum_im[COL_UNROLL];

            for (size_t k = 0; k < COL_UNROLL; ++k) {
                sum_re[k] = _mm512_setzero_pd();
                sum_im[k] = _mm512_setzero_pd();
            }

            for (size_t i = 0; i < n_vec_end; i += VEC_SIZE) {
                // Восемь вещественных весов.
                const __m512d vx = _mm512_loadu_pd(ptr_x + i);

                for (size_t k = 0; k < COL_UNROLL; ++k) {
                    // Загружаем 8 комплексных чисел (16 double).
                    const __m512d v0 = _mm512_loadu_pd(ptr_Ak[k] + 2 * i);

                    const __m512d v1 = _mm512_loadu_pd(ptr_Ak[k] + 2 * i + 8);

                    // Разделяем действительную и мнимую части.
                    const __m512d va_re = _mm512_permutex2var_pd(v0, idx_re, v1);

                    const __m512d va_im = _mm512_permutex2var_pd(v0, idx_im, v1);

                    // Умножение на вещественные веса X.
                    sum_re[k] = _mm512_fmadd_pd(va_re, vx, sum_re[k]);

                    sum_im[k] = _mm512_fmadd_pd(va_im, vx, sum_im[k]);
                }
            }

            // Умножение результата на веса Y.
            for (size_t k = 0; k < COL_UNROLL; ++k) {
                const __m512d vy = _mm512_set1_pd(ptr_y[j + k]);

                acc_re[k] = _mm512_fmadd_pd(sum_re[k], vy, acc_re[k]);

                acc_im[k] = _mm512_fmadd_pd(sum_im[k], vy, acc_im[k]);

                // Обработка N % 8 оставшихся элементов.
                accumulate_tail(ptr_Ak[k], ptr_y[j + k]);
            }
        }

        // Оставшиеся столбцы, если M не кратно 4.
        for (; j < M; ++j) {
            const double* ptr_A0 = reinterpret_cast<const double*>(A.data() + j * N);

            __m512d sum_re = _mm512_setzero_pd();
            __m512d sum_im = _mm512_setzero_pd();

            for (size_t i = 0; i < n_vec_end; i += VEC_SIZE) {
                const __m512d vx = _mm512_loadu_pd(ptr_x + i);

                const __m512d v0 = _mm512_loadu_pd(ptr_A0 + 2 * i);

                const __m512d v1 = _mm512_loadu_pd(ptr_A0 + 2 * i + 8);

                const __m512d va_re = _mm512_permutex2var_pd(v0, idx_re, v1);

                const __m512d va_im = _mm512_permutex2var_pd(v0, idx_im, v1);

                sum_re = _mm512_fmadd_pd(va_re, vx, sum_re);

                sum_im = _mm512_fmadd_pd(va_im, vx, sum_im);
            }

            const __m512d vy = _mm512_set1_pd(ptr_y[j]);

            acc_re[0] = _mm512_fmadd_pd(sum_re, vy, acc_re[0]);

            acc_im[0] = _mm512_fmadd_pd(sum_im, vy, acc_im[0]);

            accumulate_tail(ptr_A0, ptr_y[j]);
        }

        // Объединяем четыре векторных аккумулятора.
        const __m512d acc_re_vec =
            _mm512_add_pd(_mm512_add_pd(acc_re[0], acc_re[1]), _mm512_add_pd(acc_re[2], acc_re[3]));

        const __m512d acc_im_vec =
            _mm512_add_pd(_mm512_add_pd(acc_im[0], acc_im[1]), _mm512_add_pd(acc_im[2], acc_im[3]));

        // Горизонтальное суммирование и добавление хвостов.
        const double result_re = _mm512_reduce_add_pd(acc_re_vec) + tail_re;

        const double result_im = _mm512_reduce_add_pd(acc_im_vec) + tail_im;

        return {result_re, result_im};
    }
}

#endif  // defined(__AVX512F__) || defined(HAS_AVX512)

template <typename QuadRuleX, typename QuadRuleY, typename MatrixArray>
inline auto quadratic_form(const MatrixArray& A) {
    using ReturnType = std::remove_cvref_t<decltype(A[0])>;
#if defined(__AVX512F__) || defined(HAS_AVX512)
    if constexpr (std::same_as<ReturnType, double>) {
        return quadratic_form_double<QuadRuleX, QuadRuleY>(A);
    } else if constexpr (std::same_as<ReturnType, std::complex<double>>) {
        return quadratic_form_complex<QuadRuleX, QuadRuleY>(A);
    } else {
        return quadratic_form_fallback<QuadRuleX, QuadRuleY>(A);
    }

#else

    return detail::quadratic_form_fallback<QuadRuleX, QuadRuleY>(A);

#endif
}

}  // namespace quadints::linalg
#endif  // QUADINTS_LINALG_QUADRATIC_FORM
