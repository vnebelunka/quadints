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

constexpr size_t kNumTriangles = 10'000;

struct TriangleData {
    std::vector<Triangle> x;
    std::vector<Triangle> y;
    TriangleData() : x(generate_random_triangles(kNumTriangles)), y(generate_random_triangles(kNumTriangles, 15, 25)) {}
};

const TriangleData& shared_triangles() {
    static const TriangleData data;
    return data;
}

template <typename QuadRule, typename Integrand>
void adaptive_integration_triangle_pair(benchmark::State& state) {
    const auto& tri = shared_triangles();
    const size_t n = tri.x.size();

    const auto rtol = std::pow(10, -static_cast<double>(state.range(0)));

    size_t total_depth = 0;
    size_t total_calls = 0;
    size_t total_converged = 0;

    auto get_info = [&total_depth, &total_calls, &total_converged](
                        CurrentIntegral<return_type_2d<Integrand, Triangle, Triangle>>,
                        PreviousIntegral<return_type_2d<Integrand, Triangle, Triangle>>, size_t current_depth,
                        bool is_converged) {
        total_depth += current_depth;
        constexpr auto num_triangles = [](size_t depth) { return (std::pow(4, depth + 1) - 1) / 3; };
        total_calls += QuadRule::n_points * QuadRule::n_points * num_triangles(current_depth);
        total_converged += is_converged;
    };

    auto integrator = make_adaptive_integrator2d<QuadRule, QuadRule, Triangle, Triangle>(
        DefaultCriterion{Atol(0.0), Rtol(rtol)}, comptime_max_level_2d, 0);

    size_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(integrator.integrate(Integrand{}, tri.x[i], tri.y[i], get_info));
        i = (i + 1) % n;
    }

    state.counters["AvgDepth"] = benchmark::Counter(total_depth, benchmark::Counter::kAvgIterations);
    state.counters["AvgIntegrandCalls"] = benchmark::Counter(total_calls, benchmark::Counter::kAvgIterations);
    state.counters["AvgConverged"] = benchmark::Counter(total_converged, benchmark::Counter::kAvgIterations);
}

}  // namespace

BENCHMARK(adaptive_integration_triangle_pair<TriangleQuadrature<double, 1>, constant_func>)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->ArgName("-log10(rtol)")
    ->ArgName("-log10(rtol)");
BENCHMARK(adaptive_integration_triangle_pair<TriangleQuadrature<double, 3>, constant_func>)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->ArgName("-log10(rtol)");
BENCHMARK(adaptive_integration_triangle_pair<TriangleQuadrature<double, 7>, constant_func>)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->ArgName("-log10(rtol)");
BENCHMARK(adaptive_integration_triangle_pair<TriangleQuadrature<double, 1>, norm_func>)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->ArgName("-log10(rtol)");
BENCHMARK(adaptive_integration_triangle_pair<TriangleQuadrature<double, 3>, norm_func>)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->ArgName("-log10(rtol)");
BENCHMARK(adaptive_integration_triangle_pair<TriangleQuadrature<double, 7>, norm_func>)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->ArgName("-log10(rtol)");
BENCHMARK(adaptive_integration_triangle_pair<TriangleQuadrature<double, 1>, expir_func>)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->ArgName("-log10(rtol)");
BENCHMARK(adaptive_integration_triangle_pair<TriangleQuadrature<double, 3>, expir_func>)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->ArgName("-log10(rtol)");
BENCHMARK(adaptive_integration_triangle_pair<TriangleQuadrature<double, 7>, expir_func>)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->ArgName("-log10(rtol)");

BENCHMARK_MAIN();
