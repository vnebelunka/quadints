#include <gtest/gtest.h>

#include <array>
#include <cmath>

#include "gtest/gtest.h"
#include "quadints/adaptive.hpp"
#include "quadints/interface.hpp"
#include "quadints/quads_triangle.hpp"
#include "test_utility.hpp"
using namespace quadints;

// ----- 1-point (midpoint) rule -----

using TriangleMidpointRule = TriangleQuadrature<double, 1>;

TEST(AdaptiveIntegrationTriangleTest, MidPoint0degree) {
    auto constant = [](point2d) { return 5.0; };
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator =
        make_adaptive_integrator<TriangleMidpointRule, Triangle>(DefaultCriterion{Atol(0.0), Rtol(1e-3)}, 10, 0);
    double result = integrator.integrate(constant, tri_unit);
    double exact = 5.0 * tri_unit.mes();  // area=0.5 → 2.5
    EXPECT_NEAR(result, exact, 1e-12);
}

TEST(AdaptiveIntegrationTriangleTest, MidPoint1degree) {
    auto linear = [](point2d p) { return p.coords[0]; };
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator =
        make_adaptive_integrator<TriangleMidpointRule, Triangle>(DefaultCriterion{Atol(0.0), Rtol(1e-3)}, 10, 0);
    double result = integrator.integrate(linear, tri_unit);
    double exact = 1. / 6;
    EXPECT_NEAR(result, exact, std::abs(exact) * 3e-3);
}

TEST(AdaptiveIntegrationTriangleTest, MidPoint2degree) {
    auto square = [](point2d p) { return p.coords[0] * p.coords[1]; };
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator =
        make_adaptive_integrator<TriangleMidpointRule, Triangle>(DefaultCriterion{Atol(0.0), Rtol(1e-4)}, 10, 0);
    double result = integrator.integrate(square, tri_unit);
    double exact_square = 1. / 24;
    EXPECT_NEAR(result, exact_square, std::abs(exact_square) * 1e-3);
}

TEST(AdaptiveIntegrationTriangleTest, MidPoint3degree) {
    auto cubic = [](point2d p) { return p.coords[0] * p.coords[0] * p.coords[0]; };
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator =
        make_adaptive_integrator<TriangleMidpointRule, Triangle>(DefaultCriterion{Atol(0.0), Rtol(1e-3)}, 10, 0);
    double result = integrator.integrate(cubic, tri_unit);
    double exact_cubic = 1.0 / 20.0;  // 3‑point rule is not exact for cubic; we check it differs.
    EXPECT_NEAR(result, exact_cubic, std::abs(exact_cubic) * 3e-3);
}

TEST(AdaptiveIntegrationTriangleTest, MidPoint3degree1e4) {
    auto cubic = [](point2d p) { return p.coords[0] * p.coords[0] * p.coords[0]; };
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator =
        make_adaptive_integrator<TriangleMidpointRule, Triangle>(DefaultCriterion{Atol(0.0), Rtol(2e-4)}, 10, 0);
    double result = integrator.integrate(cubic, tri_unit);
    double exact_cubic = 1.0 / 20.0;  // 3‑point rule is not exact for cubic; we check it differs.
    EXPECT_NEAR(result, exact_cubic, std::abs(exact_cubic) * 4e-4);
}

TEST(AdaptiveIntegrationTriangleTest, MidPoint3degree1e6) {
    auto cubic = [](point2d p) { return p.coords[0] * p.coords[0] * p.coords[0]; };
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator =
        make_adaptive_integrator<TriangleMidpointRule, Triangle>(DefaultCriterion{Atol(0.0), Rtol(1e-6)}, 10, 0);
    EXPECT_THROW(integrator.integrate(cubic, tri_unit);, std::runtime_error);
}

TEST(AdaptiveIntegrationTriangle2Test, MidPoint0degree) {
    auto constant = [](point2d px, point2d py) { return 1.0; };
    auto tri_unit1 = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto tri_unit2 = Triangle{point2d{1.0, 0.0}, point2d{0.0, 1.0}, point2d{1.0, 1.0}};
    // AdaptiveIntegrator2d<TriangleMidpointRule, TriangleMidpointRule, Triangle, Triangle> integrator;
    auto integrator = make_adaptive_integrator2d<TriangleMidpointRule, TriangleMidpointRule, Triangle, Triangle>(
        DefaultCriterion{Atol(0.0), Rtol(1e-3)}, 10, 0);
    auto result = integrator.integrate(constant, tri_unit1, tri_unit2);
    EXPECT_NEAR(result, tri_unit1.mes() * tri_unit2.mes(), 1e-12);
}

TEST(AdaptiveIntegrationTriangle2Test, MidPoint1degree) {
    auto constant = [](point2d px, point2d py) { return px.coords[0] * py.coords[0] + px.coords[1] * py.coords[1]; };
    auto tri_unit1 = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto tri_unit2 = Triangle{point2d{1.0, 0.0}, point2d{0.0, 1.0}, point2d{1.0, 1.0}};
    // AdaptiveIntegrator2d<TriangleMidpointRule, TriangleMidpointRule, Triangle, Triangle> integrator;
    auto integrator = make_adaptive_integrator2d<TriangleMidpointRule, TriangleMidpointRule, Triangle, Triangle>(
        DefaultCriterion(Atol(0.), Rtol(1e-3)), 10, 0);
    auto result = integrator.integrate(constant, tri_unit1, tri_unit2);
    auto expected = 1. / 9;
    EXPECT_NEAR(result, expected, std::abs(expected) * 3e-3);
}

TEST(AdaptiveIntegrationTriangle2Test, MidPointDist) {
    auto dist = [](point2d px, point2d py) {
        return std::sqrt((px.coords[0] - py.coords[0]) * (px.coords[0] - py.coords[0]) +
                         (px.coords[1] - py.coords[1]) * (px.coords[1] - py.coords[1]));
    };
    auto tri_unit1 = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto tri_unit2 = Triangle{point2d{1.0, 0.0}, point2d{0.0, 1.0}, point2d{1.0, 1.0}};
    auto integrator = make_adaptive_integrator2d<TriangleMidpointRule, TriangleMidpointRule, Triangle, Triangle>(
        DefaultCriterion(Atol(0.), Rtol(5e-3)), 10, 0);
    auto result = integrator.integrate(dist, tri_unit1, tri_unit2);
    auto expected = 0.157129;
    EXPECT_NEAR(result, expected, 3e-3 * std::abs(expected));
}

TEST(AdaptiveIntegrationTriangle2Test, MidPointDist1triangle) {
    auto dist = [](point2d px, point2d py) {
        return std::sqrt((px.coords[0] - py.coords[0]) * (px.coords[0] - py.coords[0]) +
                         (px.coords[1] - py.coords[1]) * (px.coords[1] - py.coords[1]));
    };
    auto tri_unit1 = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator = make_adaptive_integrator2d<TriangleMidpointRule, TriangleMidpointRule, Triangle, Triangle>(
        DefaultCriterion(Atol(0.), Rtol(1e-2)), 10, 1);
    auto result = integrator.integrate(dist, tri_unit1, tri_unit1);
    auto expected = 0.103576;
    EXPECT_NEAR(result, expected, 3e-3 * std::abs(expected));
}

TEST(AdaptiveIntegrationTriangle2Test, MidPointDist2triangles1e5) {
    auto dist = [](point2d px, point2d py) {
        return std::sqrt((px.coords[0] - py.coords[0]) * (px.coords[0] - py.coords[0]) +
                         (px.coords[1] - py.coords[1]) * (px.coords[1] - py.coords[1]));
    };
    auto tri_unit1 = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto tri_unit2 = Triangle{point2d{1.0, 0.0}, point2d{0.0, 1.0}, point2d{1.0, 1.0}};
    auto integrator = make_adaptive_integrator2d<TriangleMidpointRule, TriangleMidpointRule, Triangle, Triangle>(
        DefaultCriterion(Atol(0.), Rtol(1e-5)), 10, 0);
    EXPECT_THROW(integrator.integrate(dist, tri_unit1, tri_unit2), std::runtime_error);
    // EXPECT_NEAR(result, expected, 3e-5 * std::abs(expected));
}

TEST(AdaptiveIntegrationTriangle2Test, SecondDegreeDist) {
    auto dist = [](point2d px, point2d py) {
        return std::sqrt((px.coords[0] - py.coords[0]) * (px.coords[0] - py.coords[0]) +
                         (px.coords[1] - py.coords[1]) * (px.coords[1] - py.coords[1]));
    };
    auto tri_unit1 = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto tri_unit2 = Triangle{point2d{1.0, 0.0}, point2d{0.0, 1.0}, point2d{1.0, 1.0}};
    using quad = TriangleQuadrature<double, 3>;
    auto integrator =
        make_adaptive_integrator2d<quad, quad, Triangle, Triangle>(DefaultCriterion(Atol(0.), Rtol(1e-4)), 10, 0);
    auto result = integrator.integrate(dist, tri_unit1, tri_unit2);
    auto expected = 0.157129;
    EXPECT_NEAR(result, expected, 3e-4 * std::abs(expected));
}
