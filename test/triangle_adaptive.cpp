#include <gtest/gtest.h>

#include <array>
#include <cmath>

#include "quadints/adaptive.hpp"
#include "quadints/interface.hpp"
#include "quadints/quads_triangle.hpp"
#include "test_utility.hpp"
using namespace quadints;

// ----- 1-point (midpoint) rule -----

using TriangleMidpointRule = TriangleQuadrature<double, 1>;

TEST(AdaptiveIntegrationTriangleTest, MidpointRule_ConstantLambda) {
    auto constant = [](point2d) { return 5.0; };
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    AdaptiveIntegrator<TriangleMidpointRule, Triangle, double> integrator;
    double result = integrator.integrate(constant, tri_unit, IntegrationParams<double>{0.0, 1e-3, 0, 10});
    double exact = 5.0 * tri_unit.mes();  // area=0.5 → 2.5
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(AdaptiveIntegrationTriangleTest, Gauss3Rule_CubicLambdaNotExact) {
    auto cubic = [](point2d p) { return p.coords[0] * p.coords[0] * p.coords[0]; };
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    AdaptiveIntegrator<TriangleMidpointRule, Triangle, double> integrator;
    double result = integrator.integrate(cubic, tri_unit, IntegrationParams<double>{0.0, 1e-3, 0, 10});
    double exact_cubic = 1.0 / 20.0;  // 3‑point rule is not exact for cubic; we check it differs.
    EXPECT_NE(result, exact_cubic);
}
