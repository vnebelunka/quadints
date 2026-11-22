
#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "quadints/common/aligned_array.hpp"
#include "quadints/integrate/integrate_collect.hpp"
#include "quadints/quadratures/collected_quadrature.hpp"
#include "quadints/quadratures/quads_triangle.hpp"
#include "test_utility.hpp"

using namespace quadints::common;

struct cos10 {
    double w = 10.0;

    double operator()(const point2d& point) const { return std::cos(w * (point.coords[0] + point.coords[1])); }

    template <std::size_t N>
    aligned_array<double, N, vectorization_alignment<double>> operator()(const std::array<point2d, N>& points) const {
        aligned_array<double, N, vectorization_alignment<double>> result{};

        for (std::size_t i = 0; i < N; ++i) result[i] = (*this)(points[i]);

        return result;
    }
};

struct cos10_scalar {
    double operator()(const point2d& point) const { return std::cos(10.0 * (point.coords[0] + point.coords[1])); }
};

struct cos10_unaligned {
    double operator()(const point2d& point) const { return std::cos(10.0 * (point.coords[0] + point.coords[1])); }

    template <std::size_t N>
    std::array<double, N> operator()(const std::array<point2d, N>& points) const {
        std::array<double, N> result{};

        for (std::size_t i = 0; i < N; ++i) result[i] = (*this)(points[i]);

        return result;
    }
};

using namespace quadints;

using bar_type = barycentric_triangle<double>;
using point_type = point2d;

constexpr auto alignment = vectorization_alignment<double>;

TEST(BatchAlignTest, DetectsAlignedBatchFunction) {
    EXPECT_TRUE((has_aligned_batch_func<cos10, point_type, 4, alignment>));
    EXPECT_TRUE((has_aligned_batch_func<cos10, point_type, 8, alignment>));
    EXPECT_TRUE((has_aligned_batch_func<cos10, point_type, 16, alignment>));
}

TEST(BatchAlignTest, RejectsScalarOnlyFunction) {
    EXPECT_FALSE((has_aligned_batch_func<cos10_scalar, point_type, 4, alignment>));
}

TEST(BatchAlignTest, RejectsUnalignedBatchFunction) {
    EXPECT_FALSE((has_aligned_batch_func<cos10_unaligned, point_type, 4, alignment>));
}

template <std::size_t N>
void check_batch_results() {
    cos10 f{10.0};

    std::array<point_type, N> points{};

    for (std::size_t i = 0; i < N; ++i) {
        points[i].coords[0] = 0.1 * i;
        points[i].coords[1] = 0.2 * i;
    }

    const auto batch_result = f(points);

    for (std::size_t i = 0; i < N; ++i) {
        const double expected = f(points[i]);

        EXPECT_NEAR(batch_result[i], expected, 1e-12) << "Mismatch at index " << i;
    }
}

TEST(BatchAlignTest, BatchResultsMatchScalarResults) {
    check_batch_results<1>();
    check_batch_results<4>();
    check_batch_results<5>();
    check_batch_results<8>();
    check_batch_results<17>();
    check_batch_results<128>();
}
