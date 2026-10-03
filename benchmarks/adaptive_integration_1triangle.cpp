#include <benchmark/benchmark.h>
#include <gtest/gtest.h>

#include <complex>
#include <random>
#include <vector>

#include "quadints/adaptive.hpp"
#include "quadints/interface.hpp"
#include "quadints/quads_triangle.hpp"
#include "triangle_utility.hpp"
using namespace quadints;

template <size_t p>
struct xpyp {
    double operator()(const point2d& point) const { return std::pow(point.coords[0] * point.coords[1], p); }
};
struct constant {
    double operator()(const point2d& point) const { return 1.0; }
};

struct expir {
    auto operator()(const point2d& point) -> std::complex<double> {
        return std::exp(norm(point) * std::complex<double>{0., 1.});
    }
};

namespace {
template <typename QuadRule, typename Integrand>
void BM_adaptive_integration(benchmark::State& state) {
    size_t n = 10'000;
    std::vector<Triangle> triangles = generate_random_triangles(n);
    size_t depth;
    size_t total_depth = 0;
    auto get_depth = [&total_depth](CurrentIntegral<return_type<Integrand, Triangle>>,
                                    PreviousIntegral<return_type<Integrand, Triangle>>, size_t current_depth,
                                    bool) { total_depth += current_depth; };
    auto integrator = make_adaptive_integrator<QuadRule, Triangle>(DefaultCriterion{Atol(0.0), Rtol(1e-5)}, 6, 0);
    size_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(integrator.integrate(Integrand{}, triangles[i], get_depth));
        i = (i + 1) % n;
    }
    state.counters["AvgDepth"] = benchmark::Counter(total_depth, benchmark::Counter::kAvgIterations);
    auto sum_squares = [](size_t n) { return n * (n + 1) * (2 * n + 1) / 6; };
    state.counters["AvgIntegrandCalls"] = benchmark::Counter(
        QuadRule::n_points * QuadRule::n_points * sum_squares(total_depth), benchmark::Counter::kAvgIterations);
}
}  // namespace

BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 1>, constant>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 3>, constant>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 7>, constant>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 1>, xpyp<1>>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 3>, xpyp<1>>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 7>, xpyp<1>>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 1>, xpyp<3>>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 3>, xpyp<3>>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 7>, xpyp<3>>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 1>, expir>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 3>, expir>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 7>, expir>);

BENCHMARK_MAIN();
