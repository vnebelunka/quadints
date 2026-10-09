#include <benchmark/benchmark.h>
#include <gtest/gtest.h>

#include <random>
#include <vector>

#include "quadints/interface.hpp"
#include "quadints/quads_triangle.hpp"
#include "triangle_utility.hpp"

using namespace quadints;

static const Triangle tri_unit{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
// The integrand used in the test
static auto quad_xy = [](point2d p) { return std::exp(p.coords[0] * p.coords[1]); };

template <typename QuadRule>
static void BM_triangle_vec_iter(benchmark::State& state) {
    size_t n = 10000;
    auto tvec = generate_random_triangles(n);
    size_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(quadints::detail::integrate_iter<QuadRule>(quad_xy, tvec[i]));
        i += 1;
        i %= n;
    }
}

template <typename QuadRule>
static void BM_triangle_vec_collect(benchmark::State& state) {
    size_t n = 10000;
    auto tvec = generate_random_triangles(n);
    size_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(quadints::detail::integrate_collect<QuadRule>(quad_xy, tvec[i]));
        i += 1;
        i %= n;
    }
}

template <typename QuadRule>
static void BM_triangle_uni_collect(benchmark::State& state) {
    for (auto _ : state) {
        benchmark::DoNotOptimize(quadints::detail::integrate_collect<QuadRule>(quad_xy, tri_unit));
    }
}
template <typename QuadRule>
static void BM_triangle_uni_iter(benchmark::State& state) {
    for (auto _ : state) {
        benchmark::DoNotOptimize(quadints::detail::integrate_iter<QuadRule>(quad_xy, tri_unit));
    }
}

BENCHMARK(BM_triangle_vec_iter<TriangleQuadrature<double, 1>>);
BENCHMARK(BM_triangle_vec_iter<TriangleQuadrature<double, 3>>);
BENCHMARK(BM_triangle_vec_iter<TriangleQuadrature<double, 7>>);
BENCHMARK(BM_triangle_vec_collect<TriangleQuadrature<double, 1>>);
BENCHMARK(BM_triangle_vec_collect<TriangleQuadrature<double, 3>>);
BENCHMARK(BM_triangle_vec_collect<TriangleQuadrature<double, 7>>);
BENCHMARK(BM_triangle_uni_collect<TriangleQuadrature<double, 1>>);
BENCHMARK(BM_triangle_uni_collect<TriangleQuadrature<double, 3>>);
BENCHMARK(BM_triangle_uni_collect<TriangleQuadrature<double, 7>>);
BENCHMARK(BM_triangle_vec_iter<TriangleQuadrature<double, 1>>);
BENCHMARK(BM_triangle_vec_iter<TriangleQuadrature<double, 3>>);
BENCHMARK(BM_triangle_vec_iter<TriangleQuadrature<double, 7>>);

constexpr size_t N = 128;

template <typename Scalar>
consteval std::array<Scalar, N> fill_array(Scalar val) {
    std::array<Scalar, N> res;
    res.fill(val);
    return res;
}

template <typename Scalar>
struct LongQuadratureEmul {
    static constexpr size_t n_points = N;
    using point_type = barycentric_triangle<Scalar>;
    static constexpr std::array<point_type, n_points> points = fill_array(point_type{0.5, 0.5});
    static constexpr std::array<Scalar, n_points> weights = fill_array(1.);
};

BENCHMARK(BM_triangle_vec_iter<LongQuadratureEmul<double>>);
BENCHMARK(BM_triangle_vec_iter<LongQuadratureEmul<double>>);
BENCHMARK(BM_triangle_vec_iter<LongQuadratureEmul<double>>);
BENCHMARK(BM_triangle_vec_collect<LongQuadratureEmul<double>>);
BENCHMARK(BM_triangle_vec_collect<LongQuadratureEmul<double>>);
BENCHMARK(BM_triangle_vec_collect<LongQuadratureEmul<double>>);
BENCHMARK_MAIN();
