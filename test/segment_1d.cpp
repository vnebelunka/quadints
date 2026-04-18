#include "interface.hpp"
#include "quads_segment.hpp"
#include <cmath>

#include <gtest/gtest.h>

using namespace quadints;

// -------------------------------------------------------------------
// Segment structure for 1D
// -------------------------------------------------------------------
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
// Constexpr helpers for analytic integration of polynomials
// -------------------------------------------------------------------

// constexpr integer power (works for non‑negative exponents)
constexpr double ipow(double base, int exp) {
  double result = 1.0;
  for (int i = 0; i < exp; ++i)
    result *= base;
  return result;
}

// ∫_a^b x^n dx = (b^{n+1} - a^{n+1}) / (n+1)
constexpr double integral_monomial(int n, double a = 0.0, double b = 1.0) {
  return (ipow(b, n + 1) - ipow(a, n + 1)) / static_cast<double>(n + 1);
}

// ∫_a^b (c_0 + c_1 x + ... + c_{N-1} x^{N-1}) dx
template <size_t N>
constexpr double integral_polynomial(const std::array<double, N> &coeffs,
                                     double a = 0.0, double b = 1.0) {
  double result = 0.0;
  for (size_t k = 0; k < N; ++k) {
    result += coeffs[k] * integral_monomial(static_cast<int>(k), a, b);
  }
  return result;
}

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
  // ∫₀¹ (2x+1) dx = [x² + x]₀¹ = 2
  constexpr double exact = 2.0; // can also be computed analytically
  EXPECT_NEAR(result, exact, TOL);
}

TEST_F(Integration1DTest, MidpointRule_QuadraticLambda) {
  auto quadratic = [](double x) { return x * x; };
  double result = integrate<MidpointRule>(quadratic, seg_unit);
  constexpr double exact = integral_monomial(2);
  EXPECT_NEAR(result, exact, 0.1); // midpoint rule has error ~0.083
}

// -------------------------------------------------------------------
// Test cases using GaussLegendre2 (exact for polynomials up to degree 3)
// -------------------------------------------------------------------

using GaussLegendre2 = SegmentQuadrature<double, 2>;

TEST_F(Integration1DTest, GaussLegendre2_CubicLambda) {
  auto cubic = [](double x) { return x * x * x; };
  double result = integrate<GaussLegendre2>(cubic, seg_unit);
  constexpr double exact = integral_monomial(3); // 1/4
  EXPECT_NEAR(result, exact, TOL);
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
  constexpr double exact = integral_monomial(4); // 1/5
  EXPECT_NEAR(result, exact, TOL);
}

TEST_F(Integration1DTest, GaussLegendre3_QuinticLambda) {
  auto quintic = [](double x) { return x * x * x * x * x; };
  double result = integrate<GaussLegendre3>(quintic, seg_unit);
  constexpr double exact = integral_monomial(5); // 1/6
  EXPECT_NEAR(result, exact, TOL);
}

TEST_F(Integration1DTest, GaussLegendre3_SexticLambda_NotExact) {
  auto sextic = [](double x) { return x * x * x * x * x * x; };
  double result = integrate<GaussLegendre3>(sextic, seg_unit);
  constexpr double exact = integral_monomial(6); // 1/7
  EXPECT_NE(result, exact);
  EXPECT_GT(std::abs(result - exact), 1e-10);
}

TEST_F(Integration1DTest, GaussLegendre3_SineOnPi_MoreAccurate) {
  Segment seg_pi{0.0, M_PI};
  auto sine = [](double x) { return std::sin(x); };
  double result = integrate<GaussLegendre3>(sine, seg_pi);
  double exact = 2.0;
  EXPECT_NEAR(result, exact, 1e-2);
}

TEST_F(Integration1DTest, GaussLegendre3_Exponential) {
  auto exp_func = [](double x) { return std::exp(x); };
  double result = integrate<GaussLegendre3>(exp_func, seg_unit);
  double exact = std::exp(1.0) - 1.0;
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
// Test with vector‑valued integrand (vec3d)
// -------------------------------------------------------------------

struct vec3d {
  double x, y, z;
  vec3d() : x(0), y(0), z(0) {}
  vec3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

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
  auto vec_func = [](double x) { return vec3d{x, x * x, x * x * x}; };
  auto result = integrate<GaussLegendre3>(vec_func, seg_unit);
  vec3d expected{integral_monomial(1),  // ∫₀¹ x dx = 1/2
                 integral_monomial(2),  // ∫₀¹ x² dx = 1/3
                 integral_monomial(3)}; // ∫₀¹ x³ dx = 1/4
  EXPECT_NEAR(result.x, expected.x, TOL);
  EXPECT_NEAR(result.y, expected.y, TOL);
  EXPECT_NEAR(result.z, expected.z, TOL);
}
