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
// 1D Бенчмарк: Одномерное интегрирование
// ============================================================================

template <typename QuadRule, typename Integrand>
static void BM_triangle_1d(benchmark::State& state) {
    constexpr size_t n = 10'000;
    static const std::vector<Triangle> triangles = generate_random_triangles(n);
    size_t level = state.range(0);
    auto quad = CollectedQuadratureDynamic<QuadRule>(level);
    size_t evals_per_triangle = quad.points().size();
    size_t i = 0;
    const auto& points = quad.points();
    const auto& weights = quad.weights();

    for (auto _ : state) {
        benchmark::DoNotOptimize(quadints::detail::integrate_iter(Integrand{}, triangles[i], points, weights));
        i = (i + 1) % n;
    }

    // Отчет о полном числе вычислений подинтегральной функции в секунду
    state.SetItemsProcessed(state.iterations() * evals_per_triangle);
}

// ============================================================================
// 2D Бенчмарк: Двумерное интегрирование (по парам треугольников)
// ============================================================================

template <typename QuadRule, typename Integrand>
static void BM_dynamic_quadrature_xy(benchmark::State& state) {
    constexpr size_t n = 50;
    static const std::vector<Triangle> triangles = generate_random_triangles(n);
    size_t level = state.range(0);
    auto integrator = CollectedQuadratureDynamic<QuadRule>(level);

    const auto& points = integrator.points();
    const auto& weights = integrator.weights();

    for (auto _ : state) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                benchmark::DoNotOptimize(
                    detail::integrate2_iter(Integrand{}, triangles[i], triangles[j], points, weights, points, weights));
            }
        }
    }

    // Время в секундах/наносекундах на одну пару элементов (50 * 50 = 2500 элементов)
    state.counters["time_per_elem"] =
        benchmark::Counter(n * n, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);

    // Полное число вычислений функции за одну итерацию: N_pairs * (N_pts * N_pts)
    size_t total_evals = (n * n) * (points.size() * points.size());
    state.SetItemsProcessed(state.iterations() * total_evals);
}

// ============================================================================
// Автоматическая регистрация бенчмарков
// ============================================================================

// --- 1D регистрация (depth: 1..4) ---
#define REGISTER_BENCH_1D(QuadRule, IntegrandName, IntegrandType) \
    BENCHMARK_TEMPLATE(BM_triangle_1d, QuadRule, IntegrandType)   \
        ->DenseRange(1, 4, 1)                                     \
        ->ArgName("depth")                                        \
        ->Name("BM_1D/" #IntegrandName "/" #QuadRule)

#define REGISTER_ALL_1D(Name, Type)       \
    REGISTER_BENCH_1D(Quad1, Name, Type); \
    REGISTER_BENCH_1D(Quad3, Name, Type); \
    REGISTER_BENCH_1D(Quad7, Name, Type);

// --- 2D регистрация (depth: 1..3) ---
#define REGISTER_BENCH_2D(QuadRule, IntegrandName, IntegrandType)         \
    BENCHMARK_TEMPLATE(BM_dynamic_quadrature_xy, QuadRule, IntegrandType) \
        ->DenseRange(1, 3, 1)                                             \
        ->ArgName("depth")                                                \
        ->Name("BM_2D/" #IntegrandName "/" #QuadRule)

#define REGISTER_ALL_2D(Name, Type)       \
    REGISTER_BENCH_2D(Quad1, Name, Type); \
    REGISTER_BENCH_2D(Quad3, Name, Type); \
    REGISTER_BENCH_2D(Quad7, Name, Type);

// ============================================================================
// Запуск всех комбинаций
// ============================================================================

// 1D Бенчмарки
REGISTER_ALL_1D(constant, constant);
REGISTER_ALL_1D(xpyp_1, xpyp<1>);
REGISTER_ALL_1D(xpyp_2, xpyp<2>);
REGISTER_ALL_1D(xpyp_4, xpyp<4>);
REGISTER_ALL_1D(cos10, cos10);
REGISTER_ALL_1D(expir, expir);

// 2D Бенчмарки
REGISTER_ALL_2D(constant, constant);
REGISTER_ALL_2D(xpyp_1, xpyp<1>);
REGISTER_ALL_2D(xpyp_2, xpyp<2>);
REGISTER_ALL_2D(xpyp_4, xpyp<4>);
REGISTER_ALL_2D(cos10, cos10);
REGISTER_ALL_2D(expir, expir);

BENCHMARK_MAIN();
