#include <gtest/gtest.h>
#include <immintrin.h>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

#include "quadints/adaptive/adaptive.hpp"
#include "quadints/quadratures/quads_triangle.hpp"
#include "sleef.h"
#include "test_utility.hpp"
struct cos10 {
    double w = 10;
    double operator()(const point2d& point) const { return cos(w * (point.coords[0] + point.coords[1])); }
    double operator()(const point2d& pointx, const point2d& pointy) const {
        return cos(w * (pointx.coords[0] + pointy.coords[1]));
    }

    template <std::size_t N>
    quadints::common::vectorized_array<double, N> operator()(const std::array<point2d, N>& pts) const {
        quadints::common::vectorized_array<double, N> result;
        std::size_t i = 0;

#if defined(__AVX512F__) || defined(HAS_AVX512)
        const __m512d vec_w = _mm512_set1_pd(w);

        // Индексы для деинтерливинга x (четные) и y (нечетные) из двух регистров
        const __m512i idx_x = _mm512_setr_epi64(0, 2, 4, 6, 8, 10, 12, 14);
        const __m512i idx_y = _mm512_setr_epi64(1, 3, 5, 7, 9, 11, 13, 15);

        // Основной цикл: 8 точек (16 double) за итерацию
        for (; i + 7 < N; i += 8) {
            // Загрузка двух блоков по 4 точки: [x0,y0..x3,y3] и [x4,y4..x7,y7]
            __m512d v0 = _mm512_loadu_pd(&pts[i].coords[0]);
            __m512d v1 = _mm512_loadu_pd(&pts[i + 4].coords[0]);

            // Разделение данных на вектор X [x0..x7] и вектор Y [y0..y7]
            __m512d vec_x = _mm512_permutex2var_pd(v0, idx_x, v1);
            __m512d vec_y = _mm512_permutex2var_pd(v0, idx_y, v1);

            // Вычисление w * (x + y)
            __m512d sum = _mm512_add_pd(vec_x, vec_y);
            __m512d arg = _mm512_mul_pd(vec_w, sum);

            // Вычисление векторного косинуса SLEEF (8 x double, 1.0 ULP)
            __m512d res_vec = Sleef_cosd8_u35avx512f(arg);

            // Сохранение 8 результатов
            _mm512_store_pd(&result[i], res_vec);
        }
#endif
        for (; i < N; ++i) {
            result[i] = operator()(pts[i]);
        }
        return result;
    }
    template <std::size_t N>
    quadints::common::vectorized_array<double, N * N> operator()(const std::array<point2d, N>& ptsX,
                                                                 const std::array<point2d, N>& ptsY) const {
        quadints::common::vectorized_array<double, N * N> result;

#if defined(__AVX512F__) || defined(HAS_AVX512)
        const __m512d vec_w = _mm512_set1_pd(w);
        // Маска нечетных индексов для извлечения Y-координат (coords[1]) из ptsY
        const __m512i idx_y = _mm512_setr_epi64(1, 3, 5, 7, 9, 11, 13, 15);

        // Внешний цикл перебирает точки X
        for (std::size_t i = 0; i < N; ++i) {
            // Broadcast X-координаты текущей точки ptsX[i] на весь 512-битный регистр
            const __m512d vec_x = _mm512_set1_pd(ptsX[i].coords[0]);

            // Внутренний векторизованный цикл по координатам Y (по 8 элементов за итерацию)
            std::size_t j = 0;
            for (; j + 7 < N; j += 8) {
                __m512d vy0 = _mm512_loadu_pd(&ptsY[j].coords[0]);
                __m512d vy1 = _mm512_loadu_pd(&ptsY[j + 4].coords[0]);

                __m512d vec_y = _mm512_permutex2var_pd(vy0, idx_y, vy1);

                __m512d sum = _mm512_add_pd(vec_x, vec_y);
                __m512d arg = _mm512_mul_pd(vec_w, sum);

                __m512d res_vec = Sleef_cosd8_u35avx512f(arg);

                // Запись в строку i матрицы результатов: result[i * N + j]
                if constexpr (N % 8 == 0) {
                    _mm512_store_pd(&result[i * N + j], res_vec);
                } else {
                    _mm512_storeu_pd(&result[i * N + j], res_vec);
                }
            }

            // Скалярная обработка остатка по оси Y (если N % 8 != 0)
            for (; j < N; ++j) {
                result[i * N + j] = operator()(ptsX[i], ptsY[j]);
            }
        }
#else
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                result[i * N + j] = operator()(ptsX[i], ptsY[j]);
            }
        }
#endif
        return result;
    }
};

namespace {

template <std::size_t N>
void VerifyCos10Batch() {
    cos10 func;
    std::array<point2d, N> ptsX{};
    std::array<point2d, N> ptsY{};

    // Заполняем тестовыми точками с разными знаками и масштабами
    for (std::size_t i = 0; i < N; ++i) {
        const double t = static_cast<double>(i);
        ptsX[i] = point2d{0.13 * t - 0.7, -0.21 * t + 0.4};
        ptsY[i] = point2d{-0.17 * t + 0.9, 0.29 * t - 0.5};
    }

    // 1. Проверка 1D batch-версии: operator()(ptsX)
    const auto batch_1d = func(ptsX);
    for (std::size_t i = 0; i < N; ++i) {
        const double expected = func(ptsX[i]);
        EXPECT_NEAR(batch_1d[i], expected, 1e-14) << "1D mismatch at N=" << N << ", i=" << i;
    }

    // 2. Проверка 2D batch-версии: operator()(ptsX, ptsY)
    const auto batch_2d = func(ptsX, ptsY);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            const double expected = func(ptsX[i], ptsY[j]);
            EXPECT_NEAR(batch_2d[i * N + j], expected, 1e-14) << "2D mismatch at N=" << N << ", i=" << i << ", j=" << j;
        }
    }
}

}  // namespace

TEST(Cos10Test, BatchMatchesScalarAcrossSizes) {
    VerifyCos10Batch<1>();   // граничный случай: 1 элемент
    VerifyCos10Batch<5>();   // N < 8: только скалярная ветка
    VerifyCos10Batch<8>();   // N == 8: ровно 1 регистр AVX-512, ветка (N % 8 == 0)
    VerifyCos10Batch<13>();  // N > 8, N % 8 != 0: AVX-512 + _mm512_storeu_pd + скалярный хвост
    VerifyCos10Batch<16>();  // N == 16: 2 регистра AVX-512 без остатка
}

TEST(Cos10Test, AdaptiveIntegralConvergesToAnalytical) {
    using namespace quadints;
    cos10 func;
    const double w = func.w;

    const double a = (1.0 - std::cos(w)) / (w * w);
    const double b = (w - std::sin(w)) / (w * w);
    const std::complex<double> expected{a * a - b * b, 0.0};  // ~ -0.01077941970978211

    {
        using quad = TriangleQuadrature<double, 3>;
        auto integrator =
            make_adaptive_integrator2d<quad, quad, Triangle, Triangle>(DefaultCriterion(Atol(0.), Rtol(1e-4)), 10, 0);
        auto tri_unit1 = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
        auto result3 = integrator.integrate(func, tri_unit1, tri_unit1);

        // std::real / std::imag работают как для std::complex<double>, так и для double
        EXPECT_NEAR(std::real(result3), expected.real(), 2e-4);
        EXPECT_NEAR(std::imag(result3), expected.imag(), 2e-4);
    }
}
