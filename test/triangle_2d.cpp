#include <array>
#include <cmath>

#include "interface.hpp"
#include "quads_triangle.hpp"
#include "test_utility.hpp"
#include <gtest/gtest.h>

using namespace quadints;

struct Triangle {
  using point_type = point2d;
  std::array<point2d, 3> vertices;
  constexpr double mes() const {
    // Compute the area of the triangle using the determinant formula
    const auto &A = vertices[0].coords;
    const auto &B = vertices[1].coords;
    const auto &C = vertices[2].coords;
    return 0.5 * std::abs(A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) +
                          C[0] * (A[1] - B[1]));
  }
};

// Test fixture
class TriangleIntegrationTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Unit right triangle: (0,0), (1,0), (0,1)
    tri_unit =
        Triangle{point2d{0.0, 0.0}, point2d{1.0, 0.0}, point2d{0.0, 1.0}};
    // Custom triangle: (1,1), (4,1), (1,3)
    tri_custom =
        Triangle{point2d{1.0, 1.0}, point2d{4.0, 1.0}, point2d{1.0, 3.0}};
  }
  Triangle tri_unit;
  Triangle tri_custom;
};

// ----- 1-point (midpoint) rule -----

using TriangleMidpointRule = TriangleQuadrature<double, 1>;

TEST_F(TriangleIntegrationTest, MidpointRule_ConstantLambda) {
  auto constant = [](point2d) { return 5.0; };
  double result = integrate<TriangleMidpointRule>(constant, tri_unit);
  double exact = 5.0 * tri_unit.mes(); // area=0.5 → 2.5
  EXPECT_NEAR(result, exact, 1e-12);
}

TEST_F(TriangleIntegrationTest, MidpointRule_LinearXLambda) {
  auto linear_x = [](point2d p) { return p.coords[0]; };
  double result = integrate<TriangleMidpointRule>(linear_x, tri_unit);
  EXPECT_NEAR(result, 1.0 / 6.0, 1e-12);
}

TEST_F(TriangleIntegrationTest, MidpointRule_LinearYLambda) {
  auto linear_y = [](point2d p) { return p.coords[1]; };
  double result = integrate<TriangleMidpointRule>(linear_y, tri_unit);
  EXPECT_NEAR(result, 1.0 / 6.0, 1e-12);
}

TEST_F(TriangleIntegrationTest, MidpointRule_QuadraticX2Lambda) {
  auto quad_x2 = [](point2d p) { return p.coords[0] * p.coords[0]; };
  double result = integrate<TriangleMidpointRule>(quad_x2, tri_unit);
  EXPECT_NEAR(result, 1.0 / 18.0, 1e-12);
}

// ----- 3‑point Gauss rule (exact for degree ≤2) -----
using TriangleGauss3Rule = TriangleQuadrature<double, 3>;

TEST_F(TriangleIntegrationTest, Gauss3Rule_QuadraticX2Lambda) {
  auto quad_x2 = [](point2d p) { return p.coords[0] * p.coords[0]; };
  double result = integrate<TriangleGauss3Rule>(quad_x2, tri_unit);
  EXPECT_NEAR(result, 1.0 / 12.0, 1e-12);
}

TEST_F(TriangleIntegrationTest, Gauss3Rule_QuadraticY2Lambda) {
  auto quad_y2 = [](point2d p) { return p.coords[1] * p.coords[1]; };
  double result = integrate<TriangleGauss3Rule>(quad_y2, tri_unit);
  EXPECT_NEAR(result, 1.0 / 12.0, 1e-12);
}

TEST_F(TriangleIntegrationTest, Gauss3Rule_QuadraticXYLambda) {
  auto quad_xy = [](point2d p) { return p.coords[0] * p.coords[1]; };
  double result = integrate<TriangleGauss3Rule>(quad_xy, tri_unit);
  EXPECT_NEAR(result, 1.0 / 24.0, 1e-12);
}

// ----- 3‑point Gauss rule (exact for degree ≤2) -----
TEST_F(TriangleIntegrationTest, Gauss3Rule_CubicLambdaNotExact) {
  auto cubic = [](point2d p) {
    return p.coords[0] * p.coords[0] * p.coords[0];
  };
  double result = integrate<TriangleGauss3Rule>(cubic, tri_unit);
  double exact_cubic =
      1.0 / 20.0; // 3‑point rule is not exact for cubic; we check it differs.
  EXPECT_NE(result, exact_cubic);
}

// ----- 4‑point Gauss rule (exact for degree ≤3) ----
using TriangleGauss4Rule = TriangleQuadrature<double, 7>;

TEST_F(TriangleIntegrationTest, Gauss4Rule_CubicLambda) {
  auto cubic = [](point2d p) {
    return p.coords[0] * p.coords[0] * p.coords[0];
  };
  double result = integrate<TriangleGauss4Rule>(cubic, tri_unit);
  EXPECT_NEAR(result, 1.0 / 20.0, 1e-12);
}

TEST_F(TriangleIntegrationTest, Gauss4Rule_QuarticLambdaNotExact) {
  auto quartic = [](point2d p) {
    return p.coords[0] * p.coords[0] * p.coords[0] * p.coords[0];
  };
  double result = integrate<TriangleGauss4Rule>(quartic, tri_unit);
  double exact_quartic = 1.0 / 30.0;
  EXPECT_NE(result, exact_quartic);
}
