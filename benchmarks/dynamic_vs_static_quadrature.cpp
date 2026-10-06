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
void BM_dynamic_quadrature(benchmark::State& state) {
    size_t n = 10'000;
    std::vector<Triangle> triangles = generate_random_triangles(n);
    size_t level = state.range(0);
    auto integrator = CollectedQuadratureDynamic<QuadRule>(level);
    auto points = integrator.points();
    auto weights = integrator.weights();
    size_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(detail::integrate_collect(Integrand{}, triangles[i], points, weights));
        i = (i + 1) % n;
    }
}

template <typename QuadRule, typename Integrand, size_t depth>
void BM_static_quadrature(benchmark::State& state) {
    size_t n = 10'000;
    std::vector<Triangle> triangles = generate_random_triangles(n);
    size_t i = 0;
    using quad = CollectedQuadratureStatic<QuadRule, depth, double>;
    for (auto _ : state) {
        benchmark::DoNotOptimize(detail::integrate_collect<quad>(Integrand{}, triangles[i]));
        i = (i + 1) % n;
    }
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

BENCHMARK(BM_dynamic_quadrature<TriangleQuadrature<double, 1>, expir>)->Arg(10)->ArgName("depth");

BENCHMARK_MAIN();
