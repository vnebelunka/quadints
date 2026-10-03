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
TEST(AdaptiveIntegrationTriangleTestCallback, MidPoint0degree) {
    auto constant = [](point2d) { return 5.0; };
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator =
        make_adaptive_integrator<TriangleMidpointRule, Triangle>(DefaultCriterion{Atol(0.0), Rtol(1e-3)}, 10, 0);
    bool converged;
    size_t depth;
    double result = integrator.integrate(constant, tri_unit, [&converged, &depth](CurrentIntegral<double>, PreviousIntegral<double>, size_t current_depth, bool is_converged) {
        converged = is_converged;
        depth = current_depth;
    });
    ASSERT_TRUE(converged);
    ASSERT_EQ(depth, 1);
    double exact = 5.0 * tri_unit.mes();  // area=0.5 → 2.5
    EXPECT_NEAR(result, exact, 1e-12);
}
