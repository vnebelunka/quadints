#include <benchmark/benchmark.h>
#include <gtest/gtest.h>

#include <complex>
#include <random>
#include <vector>

#include "quadints/adaptive.hpp"
#include "quadints/interface.hpp"
#include "quadints/quads_triangle.hpp"
#include "quadints/splitter_triangle.hpp"
#include "triangle_utility.hpp"
using namespace quadints;

template <size_t p>
struct xpyp {
    double operator()(const point2d& pointx, const point2d& pointy) const {
        return std::pow(pointx.coords[0] + pointy.coords[1], p);
    }
};
struct constant {
    double operator()(const point2d& pointx, const point2d& pointy) const { return 1.0; }
};

struct expir {
    auto operator()(const point2d& pointx, const point2d& pointy) -> std::complex<double> {
        auto r = norm(pointx - pointy);
        return std::exp(r * std::complex<double>{0., 1.});
    }
};

namespace {
template <typename QuadRule, typename Integrand>
void BM_dynamic_quadrature(benchmark::State& state) {
    size_t n = 50;
    std::vector<Triangle> triangles = generate_random_triangles(n);
    size_t level = state.range(0);
    auto integrator = CollectedQuadratureDynamic<QuadRule>(level);
    for (auto _ : state) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                benchmark::DoNotOptimize(detail::integrate2_collect(Integrand{}, triangles[i], triangles[j],
                                                                    integrator.points(), integrator.weights(),
                                                                    integrator.points(), integrator.weights()));
            }
        }
    }
    state.counters["time_per_elem"] =
        benchmark::Counter(n * n, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
}

template <typename QuadRule, typename Integrand, size_t depth>
void BM_static_quadrature(benchmark::State& state) {
    size_t n = 50;
    std::vector<Triangle> triangles = generate_random_triangles(n);
    using quad = CollectedQuadratureStatic<QuadRule, depth, double>;
    for (auto _ : state) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                benchmark::DoNotOptimize(
                    detail::integrate2_collect<quad, quad>(Integrand{}, triangles[i], triangles[j]));
            }
        }
    }
    state.counters["time_per_elem"] =
        benchmark::Counter(n * n, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
}

}  // namespace

BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 1>, constant>)->Arg(1)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 3>, constant>)->Arg(1)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 7>, constant>)->Arg(1)->ArgName("depth");
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 1>, constant, 1>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 3>, constant, 1>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 7>, constant, 1>);

BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 1>, constant>)->Arg(2)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 3>, constant>)->Arg(2)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 7>, constant>)->Arg(2)->ArgName("depth");
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 1>, constant, 2>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 3>, constant, 2>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 7>, constant, 2>);

BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 1>, constant>)->Arg(3)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 3>, constant>)->Arg(3)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 7>, constant>)->Arg(3)->ArgName("depth");
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 1>, constant, 3>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 3>, constant, 3>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 7>, constant, 3>);

BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 1>, expir>)->Arg(1)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 3>, expir>)->Arg(1)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 7>, expir>)->Arg(1)->ArgName("depth");
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 1>, expir, 1>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 3>, expir, 1>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 7>, expir, 1>);

BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 1>, expir>)->Arg(2)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 3>, expir>)->Arg(2)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 7>, expir>)->Arg(2)->ArgName("depth");
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 1>, expir, 2>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 3>, expir, 2>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 7>, expir, 2>);

BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 1>, expir>)->Arg(3)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 3>, expir>)->Arg(3)->ArgName("depth");
BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 7>, expir>)->Arg(3)->ArgName("depth");
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 1>, expir, 3>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 3>, expir, 3>);
BENCHMARK(BM_static_quadrature<TriangleQuadrature<double, 7>, expir, 3>);

BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 1>, expir>)->Arg(5)->ArgName("depth");

BENCHMARK_MAIN();
