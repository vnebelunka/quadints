#include "quads_product.hpp"
#include "quads_segment.hpp"
#include "test_utility.hpp"
#include <gtest/gtest.h>

using namespace quadints;

using MidpointRule = SegmentQuadrature<double, 1>;
using SquareMidpointRule = Quadrature_product<MidpointRule, MidpointRule>;

using barycentric_square = SquareMidpointRule::barycentric_coord;

struct square {
  using point_type = point2d;
  point2d a;
  double dx = 1;
  auto from_reference_domain(const barycentric_square &c) const -> point2d {
    return {a.coords[0] + c.coords[0] * dx, a.coords[1] + c.coords[1] * dx};
  }
  double mes() const { return dx * dx; }
};

TEST(test_quadrature, square) {
  square sq{{0.0, 0.0}, 1.0};
  auto I = integrate<SquareMidpointRule>(
      [](point2d p) { return p.coords[0] * p.coords[1]; }, sq);
  ASSERT_NEAR(I, 1.0 / 4.0, 1e-12);
}
