#include "test_utility.hpp"
#include <gtest/gtest.h>
#include <splitter_triangle.hpp>

using namespace quadints;

TEST(split_triangle, depth1) {
  constexpr point2d a{0.0, 0.0};
  constexpr point2d b{1.0, 0.0};
  constexpr point2d c{0.0, 1.0};
  constexpr Triangle tri{a, b, c};
  auto qp = TriangleRange<double, 0>().to_vector();
  ASSERT_EQ(qp.size(), 1);
  for (auto tqp : qp) {
    ASSERT_TRUE(tqp[0].to_domain(tri) == a);
    ASSERT_TRUE(tqp[1].to_domain(tri) == b);
    ASSERT_TRUE(tqp[2].to_domain(tri) == c);
  }
}
