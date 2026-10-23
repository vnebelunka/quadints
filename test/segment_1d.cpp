#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>

#include "quadints/interface.hpp"
#include "quadints/quadratures/quads_segment.hpp"

using namespace quadints;

namespace {

// -------------------------------------------------------------------
// Segment structure for 1D
// -------------------------------------------------------------------
struct Segment {
    using point_type = double;
    double _start;
    double _end;
    [[nodiscard]] constexpr auto mes() const -> double { return _end - _start; }
    [[nodiscard]] constexpr auto start() const -> double { return _start; }
    [[nodiscard]] constexpr auto end() const -> double { return _end; }
};

constexpr double TOL = 1e-12;

// -------------------------------------------------------------------
// Constexpr helpers for analytic integration of polynomials
// -------------------------------------------------------------------

// constexpr integer power (works for non‑negative exponents)
constexpr auto ipow(double base, int exp) -> double {
    double result = 1.0;
    for (int i = 0; i < exp; ++i) {
        result *= base;
    }
    return result;
}

// ∫_a^b x^n dx = (b^{n+1} - a^{n+1}) / (n+1)
constexpr auto integral_monomial(int n, double a = 0.0, double b = 1.0) -> double {
    return (ipow(b, n + 1) - ipow(a, n + 1)) / static_cast<double>(n + 1);
}

// ∫_a^b (c_0 + c_1 x + ... + c_{N-1} x^{N-1}) dx
template <size_t N>
constexpr auto integral_polynomial(const std::array<double, N>& coeffs, double a = 0.0, double b = 1.0) -> double {
    double result = 0.0;
    for (size_t k = 0; k < N; ++k) {
        result += coeffs[k] * integral_monomial(static_cast<int>(k), a, b);
    }
    return result;
}
}  // namespace

// -------------------------------------------------------------------
// Test fixture for 1D integration
// -------------------------------------------------------------------
class Integration1DTest : public ::testing::Test {
   protected:
    void SetUp() override {
        seg_unit = Segment{._start = 0.0, ._end = 1.0};
        seg_custom = Segment{._start = 2.0, ._end = 5.0};
    }

    Segment seg_unit{};
    Segment seg_custom{};
};

// -------------------------------------------------------------------
// Test cases using MidpointRule
// -------------------------------------------------------------------

using MidpointRule = SegmentQuadrature<double, 1>;

TEST_F(Integration1DTest, MidpointRule_ConstantLambda) {
    auto constant = [](double) -> double { return 5.0; };
    double const result = integrate<MidpointRule>(constant, seg_unit);
    EXPECT_NEAR(result, 5.0 * seg_unit.mes(), TOL);
}

TEST_F(Integration1DTest, MidpointRule_LinearLambda) {
    auto linear = [](double x) -> double { return (2.0 * x) + 1.0; };
    double const result = integrate<MidpointRule>(linear, seg_unit);
    // ∫₀¹ (2x+1) dx = [x² + x]₀¹ = 2
    constexpr double exact = 2.0;  // can also be computed analytically
    EXPECT_NEAR(result, exact, TOL);
}

TEST_F(Integration1DTest, MidpointRule_QuadraticLambda) {
    auto quadratic = [](double x) -> double { return x * x; };
    double const result = integrate<MidpointRule>(quadratic, seg_unit);
    constexpr double exact = integral_monomial(2);
    EXPECT_NEAR(result, exact, 0.1);  // midpoint rule has error ~0.083
}

// -------------------------------------------------------------------
// Test cases using GaussLegendre2 (exact for polynomials up to degree 3)
// -------------------------------------------------------------------

using GaussLegendre2 = SegmentQuadrature<double, 2>;

TEST_F(Integration1DTest, GaussLegendre2_CubicLambda) {
    auto cubic = [](double x) -> double { return x * x * x; };
    double const result = integrate<GaussLegendre2>(cubic, seg_unit);
    constexpr double exact = integral_monomial(3);  // 1/4
    EXPECT_NEAR(result, exact, TOL);
}

TEST_F(Integration1DTest, GaussLegendre2_SineLambda) {
    Segment const seg_pi{._start = 0.0, ._end = std::numbers::pi};
    auto sine = [](double x) -> double { return std::sin(x); };
    double const result = integrate<GaussLegendre2>(sine, seg_pi);
    EXPECT_NEAR(result, 2.0, 1e-1);
}

// -------------------------------------------------------------------
// Test cases using GaussLegendre3 (exact for polynomials up to degree 5)
// -------------------------------------------------------------------

using GaussLegendre3 = SegmentQuadrature<double, 3>;

TEST_F(Integration1DTest, GaussLegendre3_QuarticLambda) {
    auto quartic = [](double x) -> double { return x * x * x * x; };
    double const result = integrate<GaussLegendre3>(quartic, seg_unit);
    constexpr double exact = integral_monomial(4);  // 1/5
    EXPECT_NEAR(result, exact, TOL);
}

TEST_F(Integration1DTest, GaussLegendre3_QuinticLambda) {
    auto quintic = [](double x) -> double { return x * x * x * x * x; };
    double const result = integrate<GaussLegendre3>(quintic, seg_unit);
    constexpr double exact = integral_monomial(5);  // 1/6
    EXPECT_NEAR(result, exact, TOL);
}

TEST_F(Integration1DTest, GaussLegendre3_SexticLambda_NotExact) {
    auto sextic = [](double x) -> double { return x * x * x * x * x * x; };
    double const result = integrate<GaussLegendre3>(sextic, seg_unit);
    constexpr double exact = integral_monomial(6);  // 1/7
    EXPECT_NE(result, exact);
    EXPECT_GT(std::abs(result - exact), 1e-10);
}

TEST_F(Integration1DTest, GaussLegendre3_SineOnPi_MoreAccurate) {
    Segment const seg_pi{._start = 0.0, ._end = std::numbers::pi};
    auto sine = [](double x) -> double { return std::sin(x); };
    double const result = integrate<GaussLegendre3>(sine, seg_pi);
    double const exact = 2.0;
    EXPECT_NEAR(result, exact, 1e-2);
}

TEST_F(Integration1DTest, GaussLegendre3_Exponential) {
    auto exp_func = [](double x) -> double { return std::exp(x); };
    double const result = integrate<GaussLegendre3>(exp_func, seg_unit);
    double const exact = std::numbers::e - 1.0;
    EXPECT_NEAR(result, exact, 1e-6);
}

// -------------------------------------------------------------------
// Test on a larger interval with a non‑polynomial function
// -------------------------------------------------------------------
TEST_F(Integration1DTest, GaussLegendre3_LogOnInterval) {
    Segment const seg_log{._start = 1.0, ._end = 2.0};
    auto log_func = [](double x) -> double { return std::log(x); };
    double const result = integrate<GaussLegendre3>(log_func, seg_log);
    double const exact = (2.0 * std::numbers::ln2) - 1.0;
    EXPECT_NEAR(result, exact, 1e-5);
}

// -------------------------------------------------------------------
// Test with vector‑valued integrand (vec3d)
// -------------------------------------------------------------------

struct vec3d {
    double x, y, z;
    vec3d() : x(0), y(0), z(0) {}
    vec3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    auto operator+(const vec3d& other) const -> vec3d { return {x + other.x, y + other.y, z + other.z}; }
    auto operator*(double scalar) const -> vec3d { return {x * scalar, y * scalar, z * scalar}; }
    auto operator+=(const vec3d& other) -> vec3d& {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
    auto operator*=(double scalar) -> vec3d& {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }
};

static auto operator*(double scalar, const vec3d& v) -> vec3d { return {scalar * v.x, scalar * v.y, scalar * v.z}; }

static auto norm(const vec3d& v) -> double { return std::sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z)); }

TEST_F(Integration1DTest, VectorValuedLambda) {
    auto vec_func = [](double x) -> vec3d { return vec3d{x, x * x, x * x * x}; };
    auto result = integrate<GaussLegendre3>(vec_func, seg_unit);
    vec3d const expected{integral_monomial(1),   // ∫₀¹ x dx = 1/2
                         integral_monomial(2),   // ∫₀¹ x² dx = 1/3
                         integral_monomial(3)};  // ∫₀¹ x³ dx = 1/4
    EXPECT_NEAR(result.x, expected.x, TOL);
    EXPECT_NEAR(result.y, expected.y, TOL);
    EXPECT_NEAR(result.z, expected.z, TOL);
}
