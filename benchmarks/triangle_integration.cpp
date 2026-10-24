#include <benchmark/benchmark.h>

#include <cmath>
#include <complex>
#include <vector>

#include "quadints/interface.hpp"
#include "quadints/quadratures/collected_quadrature.hpp"
#include "quadints/quadratures/quads_triangle.hpp"
#include "testing_funcs.hpp"
#include "triangle_utility.hpp"

using namespace quadints;

// --- Алиасы квадратурных правил ---
using Quad1 = TriangleQuadrature<double, 1>;
using Quad3 = TriangleQuadrature<double, 3>;
using Quad7 = TriangleQuadrature<double, 7>;

// ============================================================================
// 1D Бенчмарк: Одномерное интегрирование (Static Depth)
// ============================================================================

template <typename QuadRule, typename Integrand, size_t depth>
static void BM_triangle_1d(benchmark::State& state) {
    constexpr size_t n = 10'000;
    static const std::vector<Triangle> triangles = generate_random_triangles(n);

    using quad = CollectedQuadratureStatic<QuadRule, depth, double>;
    constexpr size_t evals_per_triangle = quad::n_points;

    size_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(quadints::detail::integrate_collect<quad>(Integrand{}, triangles[i]));
        i = (i + 1) % n;
    }

    // Полное количество вычислений функции в секунду
    state.SetItemsProcessed(state.iterations() * evals_per_triangle);
}

// ============================================================================
// 2D Бенчмарк: Двумерное интегрирование по парам треугольников (Static Depth)
// ============================================================================

template <typename QuadRule, typename Integrand, size_t depth>
static void BM_static_quadrature_xy(benchmark::State& state) {
    constexpr size_t n = 50;
    static const std::vector<Triangle> triangles = generate_random_triangles(n);

    using quad = CollectedQuadratureStatic<QuadRule, depth, double>;

    for (auto _ : state) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                benchmark::DoNotOptimize(
                    detail::integrate2_collect<quad, quad>(Integrand{}, triangles[i], triangles[j]));
            }
        }
    }

    // Время на обработку одной пары элементов
    state.counters["time_per_elem"] =
        benchmark::Counter(n * n, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);

    // Полное число вычислений подинтегральной функции: N_pairs * (N_pts * N_pts)
    constexpr size_t total_evals_per_iter = (n * n) * (quad::n_points * quad::n_points);
    state.SetItemsProcessed(state.iterations() * total_evals_per_iter);
}

// ============================================================================
// Автоматическая регистрация статически типизированных бенчмарков
// ============================================================================

// --- 1D Макросы (depth: 1..4) ---
#define REGISTER_BENCH_1D_DEPTH(QuadRule, IntegrandName, IntegrandType, Depth) \
    BENCHMARK_TEMPLATE(BM_triangle_1d, QuadRule, IntegrandType, Depth)         \
        ->Name("BM_1D/" #IntegrandName "/" #QuadRule "/depth:" #Depth)

#define REGISTER_BENCH_1D(QuadRule, IntegrandName, IntegrandType)       \
    REGISTER_BENCH_1D_DEPTH(QuadRule, IntegrandName, IntegrandType, 1); \
    REGISTER_BENCH_1D_DEPTH(QuadRule, IntegrandName, IntegrandType, 2); \
    REGISTER_BENCH_1D_DEPTH(QuadRule, IntegrandName, IntegrandType, 3); \
    REGISTER_BENCH_1D_DEPTH(QuadRule, IntegrandName, IntegrandType, 4)

#define REGISTER_ALL_1D(Name, Type)       \
    REGISTER_BENCH_1D(Quad1, Name, Type); \
    REGISTER_BENCH_1D(Quad3, Name, Type); \
    REGISTER_BENCH_1D(Quad7, Name, Type);

// --- 2D Макросы (depth: 1..3) ---
#define REGISTER_BENCH_2D_DEPTH(QuadRule, IntegrandName, IntegrandType, Depth)  \
    BENCHMARK_TEMPLATE(BM_static_quadrature_xy, QuadRule, IntegrandType, Depth) \
        ->Name("BM_2D/" #IntegrandName "/" #QuadRule "/depth:" #Depth)

#define REGISTER_BENCH_2D(QuadRule, IntegrandName, IntegrandType)       \
    REGISTER_BENCH_2D_DEPTH(QuadRule, IntegrandName, IntegrandType, 1); \
    REGISTER_BENCH_2D_DEPTH(QuadRule, IntegrandName, IntegrandType, 2); \
    REGISTER_BENCH_2D_DEPTH(QuadRule, IntegrandName, IntegrandType, 3)

#define REGISTER_ALL_2D(Name, Type)       \
    REGISTER_BENCH_2D(Quad1, Name, Type); \
    REGISTER_BENCH_2D(Quad3, Name, Type); \
    REGISTER_BENCH_2D(Quad7, Name, Type);

// ============================================================================
// Запуск всех комбинаций
// ============================================================================

// 1D Бенчмарки
REGISTER_ALL_1D(xpyp_1, xpyp<1>);
REGISTER_ALL_1D(xpyp_2, xpyp<2>);
REGISTER_ALL_1D(xpyp_4, xpyp<4>);
REGISTER_ALL_1D(cos10, cos10);
REGISTER_ALL_1D(expir, expir);

// 2D Бенчмарки
REGISTER_ALL_2D(xpyp_1, xpyp<1>);
REGISTER_ALL_2D(xpyp_2, xpyp<2>);
REGISTER_ALL_2D(xpyp_4, xpyp<4>);
REGISTER_ALL_2D(cos10, cos10);
REGISTER_ALL_2D(expir, expir);

BENCHMARK_MAIN();
