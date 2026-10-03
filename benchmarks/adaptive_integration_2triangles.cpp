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

namespace {
template <typename QuadRule, typename Integrand>
void BM_adaptive_integration(benchmark::State& state) {
    size_t n = 10'000;
    std::vector<Triangle> trianglesx = generate_random_triangles(n);
    std::vector<Triangle> trianglesy = generate_random_triangles(n, 15, 25);
    size_t depth;
    size_t total_depth = 0;
    auto get_depth = [&total_depth](CurrentIntegral<return_type_2d<Integrand, Triangle, Triangle>>,
                                    PreviousIntegral<return_type_2d<Integrand, Triangle, Triangle>>,
                                    size_t current_depth, bool) { total_depth += current_depth; };
    auto integrator = make_adaptive_integrator2d<QuadRule, QuadRule, Triangle, Triangle>(
        DefaultCriterion{Atol(0.0), Rtol(1e-5)}, 6, 0);
    size_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(integrator.integrate(Integrand{}, trianglesx[i], trianglesy[i], get_depth));
        i = (i + 1) % n;
    }
    state.counters["AvgDepth"] = benchmark::Counter(total_depth, benchmark::Counter::kAvgIterations);
}
}  // namespace

BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 1>, constant_func>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 3>, constant_func>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 7>, constant_func>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 1>, norm_func>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 3>, norm_func>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 7>, norm_func>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 1>, expir_func>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 3>, expir_func>);
BENCHMARK(BM_adaptive_integration<TriangleQuadrature<double, 7>, expir_func>);

BENCHMARK_MAIN();
