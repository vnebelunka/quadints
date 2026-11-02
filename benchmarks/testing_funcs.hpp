#ifndef QUADINTS_BENCHES_TESTING_FUNCS
#define QUADINTS_BENCHES_TESTING_FUNCS
#include <benchmark/benchmark.h>
#include <gtest/gtest.h>
#include <immintrin.h>
#include <sleef.h>

#include <complex>

#include "quadints/common/common.hpp"
#include "triangle_utility.hpp"

template <size_t p>
struct xpyp {
    double operator()(const point2d& point) const { return std::pow(point.coords[0] + point.coords[1], p); }
    auto operator()(const point2d& pointx, const point2d& pointy) const -> double {
        return std::pow(pointx.coords[0] + pointy.coords[1], p);
    }
};
struct constant {
    double operator()(const point2d& point) const { return 1.0; }
    auto operator()(const point2d& pointx, const point2d& pointy) const -> double { return 1.0; }
};

struct cos10 {
    double w = 10;
    double operator()(const point2d& point) const { return cos(w * (point.coords[0] + point.coords[1])); }
    double operator()(const point2d& pointx, const point2d& pointy) const {
        return cos(w * (pointx.coords[0] + pointy.coords[1]));
    }
    template <std::size_t N>
    std::array<double, N> operator()(const std::array<point2d, N>& pts) const {
        std::array<double, N> result;
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
            _mm512_storeu_pd(&result[i], res_vec);
        }
#endif
        for (; i < N; ++i) {
            result[i] = operator()(pts[i]);
        }
        return result;
    }
    template <std::size_t N>
    std::array<double, N * N> operator()(const std::array<point2d, N>& ptsX, const std::array<point2d, N>& ptsY) const {
        std::array<double, N * N> result;

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
                _mm512_storeu_pd(&result[i * N + j], res_vec);
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

static_assert(quadints::has_batch_func2<cos10, point2d, point2d, 16>);

struct expir {
    auto operator()(const point2d& point) -> std::complex<double> {
        return std::exp(norm(point) * std::complex<double>{0., 1.});
    }
    auto operator()(const point2d& pointx, const point2d& pointy) -> std::complex<double> {
        auto r = norm(pointx - pointy);
        return std::exp(r * std::complex<double>{0., 1.});
    }
};

struct norm_func {
    double operator()(const point2d& pointx, const point2d& pointy) const {
        return std::sqrt(std::pow(pointx.coords[0] - pointy.coords[0], 2) +
                         std::pow(pointx.coords[1] - pointy.coords[1], 2));
    }
};

struct constant_func {
    double operator()(const point2d& pointx, const point2d& pointy) const { return 1.0; }
};

struct expir_func {
    auto operator()(const point2d& pointx, const point2d& pointy) -> std::complex<double> {
        return std::exp(norm(pointx - pointy) * std::complex<double>{0., 1.});
    }
};

#endif  // QUADINTS_BENCHES_TESTING_FUNCS
