#include <gtest/gtest.h>

#include "quadints/StopCriterion.hpp"
#include "quadints/adaptive.hpp"
#include "quadints/interface.hpp"
#include "quadints/quads_triangle.hpp"
#include "quadints/splitter_triangle.hpp"
#include "test_utility.hpp"
using namespace quadints;

using bar_type = barycentric_triangle<double>;

TEST(TestDynamicIntegration, size_1_depth_1) {
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator = CollectedQuadratureDynamic<TriangleQuadrature<double, 1>>(1);
    auto quad_size = integrator.points().size();
    std::vector<double> buff(quad_size);
    auto constant = [](const point2d& p) { return 1.0; };
    auto result = detail::integrate_collect(constant, tri_unit, integrator, buff);
    ASSERT_NEAR(result, 0.5, 1e-12);
}

TEST(TestDynamicIntegration, size_7_depth_1) {
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator = CollectedQuadratureDynamic<TriangleQuadrature<double, 7>>(1);
    auto quad_size = integrator.points().size();
    std::vector<double> buff(quad_size);
    auto constant = [](const point2d& p) { return 1.0; };
    auto result = detail::integrate_collect(constant, tri_unit, integrator, buff);
    ASSERT_NEAR(result, 0.5, 1e-12);
}

TEST(TestDynamicIntegration, size_7_depth_2) {
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator = CollectedQuadratureDynamic<TriangleQuadrature<double, 7>>(2);
    auto quad_size = integrator.points().size();
    std::vector<double> buff(quad_size);
    auto constant = [](const point2d& p) { return 1.0; };
    auto result = detail::integrate_collect(constant, tri_unit, integrator, buff);
    ASSERT_NEAR(result, 0.5, 1e-12);
}

TEST(TestDynamicIntegration, size_7_depth_3) {
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator = CollectedQuadratureDynamic<TriangleQuadrature<double, 7>>(3);
    auto quad_size = integrator.points().size();
    std::vector<double> buff(quad_size);
    auto constant = [](const point2d& p) { return 1.0; };
    auto result = detail::integrate_collect(constant, tri_unit, integrator, buff);
    ASSERT_NEAR(result, 0.5, 1e-12);
}

TEST(TestDynamicIntegration, size_7_depth_0_array) {
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto integrator = CollectedQuadratureDynamic<TriangleQuadrature<double, 7>>(0);
    auto quad_size = integrator.points().size();
    ASSERT_EQ(quad_size, 7);
    std::array<double, 7> buff;
    auto constant = [](const point2d& p) { return 1.0; };
    auto result = detail::integrate_collect(constant, tri_unit, integrator, buff);
    ASSERT_NEAR(result, 0.5, 1e-12);
}
