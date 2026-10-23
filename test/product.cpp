#include <gtest/gtest.h>

#include "quadints/interface.hpp"
#include "quadints/quadratures/quads_product.hpp"
#include "quadints/quadratures/quads_segment.hpp"
#include "test_utility.hpp"

using namespace quadints;

using MidpointRule = SegmentQuadrature<double, 1>;
using SquareMidpointRule = Quadrature_product<MidpointRule, MidpointRule>;

using barycentric_square = SquareMidpointRule::barycentric_coord;

struct square {
    using point_type = point2d;
    point2d a;
    double dx = 1;
    [[nodiscard]] auto from_reference_domain(const barycentric_square& c) const -> point2d {
        return {a.coords[0] + (c.coords[0] * dx), a.coords[1] + (c.coords[1] * dx)};
    }
    [[nodiscard]] auto mes() const -> double { return dx * dx; }
};

TEST(test_quadrature, square) {
    square const sq{.a = {0.0, 0.0}, .dx = 1.0};
    auto I = integrate<SquareMidpointRule>([](point2d p) { return p.coords[0] * p.coords[1]; }, sq);
    ASSERT_NEAR(I, 1.0 / 4.0, 1e-12);
}
