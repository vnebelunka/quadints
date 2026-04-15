#include "interface.hpp"
#include "segment_quads.hpp"
#include <cmath>
#include <gtest/gtest.h>

struct Segment {
  using point_type = double;
  double _start;
  double _end;
  constexpr double mes() const { return _end - _start; }
  constexpr double start() const { return _start; }
  constexpr double end() const { return _end; }
};

constexpr double TOL = 1e-12;

// -------------------------------------------------------------------
// Test fixture for 1D integration
// -------------------------------------------------------------------
class Integration1DTest : public ::testing::Test {
protected:
  void SetUp() override {
    seg_unit = Segment{0.0, 1.0};
    seg_custom = Segment{2.0, 5.0};
  }

  Segment seg_unit;
  Segment seg_custom;
};

// -------------------------------------------------------------------
// Test cases using MidpointRule
// -------------------------------------------------------------------

using MidpointRule = SegmentQuadrature<double, 1>;

TEST_F(Integration1DTest, MidpointRule_ConstantLambda) {
  auto constant = [](double) { return 5.0; };
  double result = integrate<MidpointRule>(constant, seg_unit);
  EXPECT_NEAR(result, 5.0 * seg_unit.mes(), TOL);
}

TEST_F(Integration1DTest, MidpointRule_LinearLambda) {
  auto linear = [](double x) { return 2.0 * x + 1.0; };
  double result = integrate<MidpointRule>(linear, seg_unit);
  // ∫₀¹ (2x+1)dx = 2
  EXPECT_NEAR(result, 2.0, TOL);
}

TEST_F(Integration1DTest, MidpointRule_QuadraticLambda) {
  auto quadratic = [](double x) { return x * x; };
  double result = integrate<MidpointRule>(quadratic, seg_unit);
  // Exact = 1/3, but midpoint rule has error
  EXPECT_NEAR(result, 1.0 / 3.0, 0.1);
}

// -------------------------------------------------------------------
// Test cases using GaussLegendre2 (exact for polynomials up to degree 3)
// -------------------------------------------------------------------

using GaussLegendre2 = SegmentQuadrature<double, 2>;

TEST_F(Integration1DTest, GaussLegendre2_CubicLambda) {
  auto cubic = [](double x) { return x * x * x; };
  double result = integrate<GaussLegendre2>(cubic, seg_unit);
  EXPECT_NEAR(result, 0.25, TOL);
}

TEST_F(Integration1DTest, GaussLegendre2_SineLambda) {
  Segment seg_pi{0.0, M_PI};
  auto sine = [](double x) { return std::sin(x); };
  double result = integrate<GaussLegendre2>(sine, seg_pi);
  EXPECT_NEAR(result, 2.0, 1e-1);
}

// -------------------------------------------------------------------
// Test cases using GaussLegendre3 (exact for polynomials up to degree 5)
// -------------------------------------------------------------------

using GaussLegendre3 = SegmentQuadrature<double, 3>;
TEST_F(Integration1DTest, GaussLegendre3_QuarticLambda) {
  auto quartic = [](double x) { return x * x * x * x; };
  double result = integrate<GaussLegendre3>(quartic, seg_unit);
  EXPECT_NEAR(result, 1.0 / 5.0, TOL);
}

TEST_F(Integration1DTest, GaussLegendre3_QuinticLambda) {
  auto quintic = [](double x) { return x * x * x * x * x; };
  double result = integrate<GaussLegendre3>(quintic, seg_unit);
  EXPECT_NEAR(result, 1.0 / 6.0, TOL);
}

TEST_F(Integration1DTest, GaussLegendre3_SexticLambda_NotExact) {
  auto sextic = [](double x) { return x * x * x * x * x * x; };
  double result = integrate<GaussLegendre3>(sextic, seg_unit);
  double exact = 1.0 / 7.0;
  EXPECT_NE(result, exact);
  EXPECT_GT(std::abs(result - exact), 1e-10);
}

TEST_F(Integration1DTest, GaussLegendre3_SineOnPi_MoreAccurate) {
  Segment seg_pi{0.0, M_PI};
  auto sine = [](double x) { return std::sin(x); };
  double result = integrate<GaussLegendre3>(sine, seg_pi);
  double exact = 2.0;
  // 3‑point rule is much more accurate than 2‑point for sin.
  EXPECT_NEAR(result, exact, 1e-2);
}

TEST_F(Integration1DTest, GaussLegendre3_Exponential) {
  auto exp_func = [](double x) { return std::exp(x); };
  double result = integrate<GaussLegendre3>(exp_func, seg_unit);
  double exact = std::exp(1.0) -
                 1.0; // Error is small but not exact; use moderate tolerance.
  EXPECT_NEAR(result, exact, 1e-6);
}

// -------------------------------------------------------------------
// Test on a larger interval with a non‑polynomial function
// -------------------------------------------------------------------
TEST_F(Integration1DTest, GaussLegendre3_LogOnInterval) {
  Segment seg_log{1.0, 2.0};
  auto log_func = [](double x) { return std::log(x); };
  double result = integrate<GaussLegendre3>(log_func, seg_log);
  double exact = 2.0 * std::log(2.0) - 1.0;
  EXPECT_NEAR(result, exact, 1e-5);
}

// -------------------------------------------------------------------
// Test with vector‑valued integrand
// -------------------------------------------------------------------

struct vec3d {
  double x, y, z;
  vec3d operator+(const vec3d &other) const {
    return {x + other.x, y + other.y, z + other.z};
  }
  vec3d operator*(double scalar) const {
    return {x * scalar, y * scalar, z * scalar};
  }
  vec3d &operator+=(const vec3d &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }
  vec3d &operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
  }
};
vec3d operator*(double scalar, const vec3d &v) {
  return {scalar * v.x, scalar * v.y, scalar * v.z};
}

double norm(const vec3d &v) {
  return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

TEST_F(Integration1DTest, VectorValuedLambda) {
  auto vec_func = [](double x) {
    return vec3d{x, x * x, x * x * x}; // assuming vec3d(double,double,double)
  };
  auto result = integrate<GaussLegendre3>(vec_func, seg_unit);
  vec3d expected{0.5, 1.0 / 3.0, 0.25};
  EXPECT_NEAR(result.x, expected.x, TOL);
  EXPECT_NEAR(result.y, expected.y, TOL);
  EXPECT_NEAR(result.z, expected.z, TOL);
}
