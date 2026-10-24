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

// --- Type Aliases to prevent macro comma errors ---
using Quad1 = TriangleQuadrature<double, 1>;
using Quad3 = TriangleQuadrature<double, 3>;
using Quad7 = TriangleQuadrature<double, 7>;

template <typename QuadRule, typename Integrand>
static void BM_triangle(benchmark::State& state) {
    constexpr size_t n = 10'000;
    static const std::vector<Triangle> triangles = generate_random_triangles(n);

    size_t level = state.range(0);
    auto quad = CollectedQuadratureDynamic<QuadRule>(level);
    size_t evals_per_triangle = quad.points().size();
    size_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(quadints::detail::integrate_iter(Integrand{}, triangles[i], quad));
        i = (i + 1) % n;
    }

    // Reports total integrand evaluations per second in benchmark output
    state.SetItemsProcessed(state.iterations() * evals_per_triangle);
}

// --- Helper Macro for Clean Global Registration ---
#define REGISTER_BENCH(QuadRule, IntegrandName, IntegrandType) \
    BENCHMARK_TEMPLATE(BM_triangle, QuadRule, IntegrandType)   \
        ->DenseRange(1, 4, 1)                                  \
        ->ArgName("depth")                                     \
        ->Name("BM_triangle/" #IntegrandName "/" #QuadRule)

// Register all combinations across Order 1, 3, 7 and Depth 1..4
#define REGISTER_ALL_ORDERS(Name, Type) \
    REGISTER_BENCH(Quad1, Name, Type);  \
    REGISTER_BENCH(Quad3, Name, Type);  \
    REGISTER_BENCH(Quad7, Name, Type);

REGISTER_ALL_ORDERS(xpyp_1, xpyp<1>);
REGISTER_ALL_ORDERS(xpyp_2, xpyp<2>);
REGISTER_ALL_ORDERS(xpyp_2, xpyp<4>);
REGISTER_ALL_ORDERS(cos10, cos10);
REGISTER_ALL_ORDERS(expir, expir);

BENCHMARK_MAIN();
