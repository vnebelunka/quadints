#include <gtest/gtest.h>

#include <algorithm>
#include <quadints/splitter_triangle.hpp>

#include "gtest/gtest.h"
#include "quadints/quads_triangle.hpp"
#include "test_utility.hpp"

using namespace quadints;

TEST(split_triangle, depth0) {
    constexpr point2d a{0.0, 0.0};
    constexpr point2d b{1.0, 0.0};
    constexpr point2d c{0.0, 1.0};
    constexpr Triangle tri{a, b, c};
    auto qp = TriangleRange<double, 0>().to_vector();
    ASSERT_EQ(qp.size(), 1);
    for (auto tqp : qp) {
        std::cerr << "tqp:\n";
        std::cerr << tqp[0].x() << " " << tqp[0].y() << std::endl;
        std::cerr << tqp[1].x() << " " << tqp[1].y() << std::endl;
        std::cerr << tqp[2].x() << " " << tqp[2].y() << std::endl;
        std::cerr << "domain:\n";
        std::cerr << tqp[0].to_domain(tri).coords[0] << " " << tqp[0].to_domain(tri).coords[1] << std::endl;
        std::cerr << tqp[1].to_domain(tri).coords[0] << " " << tqp[1].to_domain(tri).coords[1] << std::endl;
        std::cerr << tqp[2].to_domain(tri).coords[0] << " " << tqp[2].to_domain(tri).coords[1] << std::endl;
        ASSERT_TRUE(tqp[0].to_domain(tri) == c);
        ASSERT_TRUE(tqp[1].to_domain(tri) == a);
        ASSERT_TRUE(tqp[2].to_domain(tri) == b);
    }
}

TEST(split_triangle, depth0_midpoint) {
    using TriangleMidpointRule = TriangleQuadrature<double, 1>;
    auto quad_points = TriangleRange<double, 0>().collect_quadrature_points<TriangleMidpointRule>();
    ASSERT_EQ(quad_points.size(), 1);
    for (auto qp : quad_points) {
        ASSERT_NEAR(qp.x(), 1. / 3, 1e-15);
        ASSERT_NEAR(qp.y(), 1. / 3, 1e-15);
        ASSERT_NEAR(qp.z(), 1. / 3, 1e-15);
    }
}

TEST(split_triangle, depth1_midpoint) {
    using TriangleMidpointRule = TriangleQuadrature<double, 1>;
    auto quad_points = TriangleRange<double, 1>().collect_quadrature_points<TriangleMidpointRule>();
    ASSERT_EQ(quad_points.size(), 4);
    std::sort(quad_points.begin(), quad_points.end(),
              [](const auto& a, const auto& b) { return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y()); });
    auto expected = std::vector<std::array<double, 3>>{
        {1. / 6, 1. / 6, 2. / 3},
        {1. / 6, 2. / 3, 1. / 6},
        {1. / 3, 1. / 3, 1. / 3},
        {2. / 3, 1. / 6, 1. / 6},
    };
    for (size_t i = 0; i < quad_points.size(); ++i) {
        std::cerr << "quad_points[" << i << "] = (" << quad_points[i].x() << ", " << quad_points[i].y() << ", "
                  << quad_points[i].z() << ")" << std::endl;
        std::cerr << "expected[" << i << "] = (" << expected[i][0] << ", " << expected[i][1] << ", " << expected[i][2]
                  << ")" << std::endl;

        ASSERT_NEAR(quad_points[i].x(), expected[i][0], 1e-15);
        ASSERT_NEAR(quad_points[i].y(), expected[i][1], 1e-15);
        ASSERT_NEAR(quad_points[i].z(), expected[i][2], 1e-15);
    }
}

TEST(split_triangle, depth0_midpoint_quadrature) {
    using quad = CollectedQuadrature<TriangleQuadrature<double, 1>, 0, double>;
    ASSERT_EQ(quad::n_points, 1);
    constexpr auto w = quad::weights[0];
    ASSERT_NEAR(w, 1., 1e-15);
    constexpr auto point = quad::points[0];
    ASSERT_NEAR(point.x(), 1. / 3, 1e-15);
    ASSERT_NEAR(point.y(), 1. / 3, 1e-15);
    ASSERT_NEAR(point.z(), 1. / 3, 1e-15);
}

TEST(split_triangle, depth1_midpoint_quadrature) {
    using quad = CollectedQuadrature<TriangleQuadrature<double, 1>, 1, double>;
    ASSERT_EQ(quad::n_points, 4);
    constexpr auto weights = quad::weights;
    constexpr auto check_if_constexpr_points = quad::points;
    auto points = quad::points;
    std::sort(points.begin(), points.end(),
              [](const barycentric_triangle<double>& a, const barycentric_triangle<double>& b) {
                  return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y());
              });
    std::array<double, 4> expected_weights = {1. / 4, 1. / 4, 1. / 4, 1. / 4};
    for (size_t i = 0; i < quad::n_points; ++i) {
        ASSERT_NEAR(weights[i], expected_weights[i], 1e-15);
    }
    auto expected = std::vector<std::array<double, 3>>{
        {1. / 6, 1. / 6, 2. / 3},
        {1. / 6, 2. / 3, 1. / 6},
        {1. / 3, 1. / 3, 1. / 3},
        {2. / 3, 1. / 6, 1. / 6},
    };
    for (size_t i = 0; i < quad::n_points; ++i) {
        ASSERT_NEAR(points[i].x(), expected[i][0], 1e-15);
        ASSERT_NEAR(points[i].y(), expected[i][1], 1e-15);
        ASSERT_NEAR(points[i].z(), expected[i][2], 1e-15);
    }
}

TEST(split_triangle, depth1_midpoint_quadrature_integration) {
    using quad = CollectedQuadrature<TriangleQuadrature<double, 1>, 1, double>;
    auto constant = [](point2d) { return 5.0; };
    auto tri_unit = Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    auto result = integrate<quad>(constant, tri_unit);
    double exact = 5.0 * tri_unit.mes();
    EXPECT_NEAR(result, exact, 1e-12);
}
