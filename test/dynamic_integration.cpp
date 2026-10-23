#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "quadints/adaptive/StopCriterion.hpp"
#include "quadints/adaptive/adaptive.hpp"
#include "quadints/integrate/integrate_collect.hpp"
#include "quadints/quadratures/collected_quadrature.hpp"
#include "quadints/quadratures/quads_triangle.hpp"
#include "test_utility.hpp"
using namespace quadints;

using bar_type = barycentric_triangle<double>;
namespace {
struct MidPointDynamic {
    static constexpr size_t n_points() { return 1; }
    using point_type = bar_type;
    static std::vector<point_type> points() { return {point_type(1. / 3., 1. / 3.)}; }
    static constexpr std::vector<double> weights() { return {1.0}; }
};
}  // namespace

TEST(TestDynamicIntegration, MidPointRule_constant) {
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    MidPointDynamic const integrator;
    auto constant = [](const point2d&) -> double { return 1.0; };
    auto result = detail::integrate_collect(constant, tri_unit, integrator);
    ASSERT_NEAR(result, 0.5, 1e-12);
}

TEST(TestDynamicIntegration, CollectedQuadrature_constant) {
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator = CollectedQuadratureDynamic<TriangleQuadrature<double, 1>>(1);
    auto constant = [](const point2d&) -> double { return 1.0; };
    auto result = detail::integrate_collect(constant, tri_unit, integrator);
    ASSERT_NEAR(result, 0.5, 1e-12);
}

TEST(TestDynamicIntegration, AdaptiveIntegration) {
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator = make_adaptive_integrator<TriangleQuadrature<double, 1>, Triangle>(NonAdaptiveCriterion{}, 8, 8);
    auto constant = [](const point2d&) -> double { return 1.0; };
    auto result = integrator.integrate(constant, tri_unit);
    ASSERT_NEAR(result, 0.5, 1e-12);
}

TEST(TestDynamicIntegratopn, MidPointRule_constant_xy) {
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    MidPointDynamic const integrator;
    auto constant = [](const point2d&, const point2d&) -> double { return 1.0; };
    auto result = detail::integrate2_collect(constant, tri_unit, tri_unit, integrator.points(), integrator.weights(),
                                             integrator.points(), integrator.weights());
    ASSERT_NEAR(result, 0.25, 1e-12);
}

TEST(TestDynamicIntegration, CollectedQuadrature_constant_xy) {
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator = CollectedQuadratureDynamic<TriangleQuadrature<double, 1>>(1);
    auto constant = [](const point2d&, const point2d&) -> double { return 1.0; };
    auto result = detail::integrate2_collect(constant, tri_unit, tri_unit, integrator.points(), integrator.weights(),
                                             integrator.points(), integrator.weights());
    ASSERT_NEAR(result, 0.25, 1e-12);
}

TEST(TestDynamicIntegration, AdaptiveIntegration_xy) {
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    using TriangleMidpointRule = TriangleQuadrature<double, 1>;
    auto integrator = make_adaptive_integrator2d<TriangleMidpointRule, TriangleMidpointRule, Triangle, Triangle>(
        DefaultCriterion(Atol(0.), Rtol(1e-5)), 10, 0);
    auto constant = [](const point2d&, const point2d&) -> double { return 1.0; };
    auto result = integrator.integrate(constant, tri_unit, tri_unit);
    ASSERT_NEAR(result, 0.25, 1e-12);
}
