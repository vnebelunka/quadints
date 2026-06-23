#include <array>
#include <cstddef>
#include <iterator>
#include <type_traits>

#include "quads_triangle.hpp"

namespace quadints {

template <typename Scalar, size_t depth>
struct TriangleIterator {
    using point = barycentric_triangle<Scalar>;

    using dir = barycentric_direction<Scalar>;
    static constexpr const double step = (1. / static_cast<double>(1u << depth));
    static constexpr auto dir_right = dir{-step, step, 0};
    static constexpr auto dir_left = dir{step, -step, 0};
    static constexpr auto dir_up = dir{-step, 0, step};
    using bar_coords_t = std::array<point, 3>;
    bar_coords_t current_triangle = {};
    size_t pos = 0;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = bar_coords_t;
    using difference_type = std::ptrdiff_t;
    using pointer = const value_type*;
    using reference = const value_type&;

    constexpr TriangleIterator(const std::array<point, 3>& coords, size_t pos = 0)
        : current_triangle(coords), pos(pos) {}
    constexpr reference operator*() const { return current_triangle; }
    constexpr pointer operator->() const { return &current_triangle; }

    constexpr TriangleIterator& operator++() noexcept {
        ++pos;

        auto& cur_down_point = current_triangle[0];
        auto& cur_up_point = current_triangle[1];
        // if level != 1 and triangle was above edge
        if (current_triangle[2].y() > current_triangle[0].y() && current_triangle[0].y() != 0) {
            current_triangle[2] = cur_up_point + dir_left;
            return *this;
        }
        cur_down_point += dir_up;
        cur_up_point += dir_up;
        // can't go left -> go to next level
        if (cur_up_point.x() < 0) {
            cur_down_point = {1 - current_triangle[0].y() - step, current_triangle[0].y() + step, 0};
            cur_up_point = cur_down_point + dir_up;
        }
        current_triangle[2] = cur_down_point + dir_right;
        return *this;
    }
    constexpr bool operator==(const TriangleIterator& other) const { return pos == other.pos; }
    constexpr bool operator!=(const TriangleIterator& other) const { return !(*this == other); }
};

template <typename Scalar, size_t depth>
struct TriangleRange {
    constexpr static size_t size() { return 1u << (2 * depth); }

    constexpr TriangleIterator<Scalar, depth> begin() const {
        constexpr double step = 1. / static_cast<double>(1u << depth);
        return TriangleIterator<Scalar, depth>(
            {barycentric_triangle<Scalar>{1., 0., 0}, {1 - step, 0, step}, {1 - step, step, 0}}, 0);
    }
    constexpr TriangleIterator<Scalar, depth> end() const {
        return TriangleIterator<Scalar, depth>{{barycentric_triangle<Scalar>{0., 0.}, {0., 0.}, {0., 0.}}, size()};
    }
    constexpr auto to_vector() const -> std::array<std::array<barycentric_triangle<Scalar>, 3>, size()> {
        std::vector<std::array<barycentric_triangle<Scalar>, 3>> tmp_vec(this->begin(), this->end());
        std::array<std::array<barycentric_triangle<Scalar>, 3>, size()> result;
        std::copy(tmp_vec.begin(), tmp_vec.end(), result.begin());
        return result;
    }
    template <typename QuadRule>
    constexpr auto collect_quadrature_points() const {
        using bar_point = barycentric_triangle<Scalar>;
        using bar_dir = barycentric_direction<Scalar>;
        std::array<bar_point, QuadRule::n_points * size()> bar_points;
        size_t triangle_counter = 0;
        for (const auto& tri : *this) {
            const bar_point& A = tri[0];
            const bar_point& B = tri[1];
            const bar_point& C = tri[2];
            const bar_dir CA = A - C;
            const bar_dir CB = B - C;
            for (size_t i = 0; i < QuadRule::n_points; ++i) {
                bar_points[triangle_counter * QuadRule::n_points + i] =
                    C + CA * QuadRule::points[i].x() + CB * QuadRule::points[i].y();
            }
            ++triangle_counter;
        }
        return bar_points;
    }
};

template <typename BaseQuadRule, size_t depth, typename Scalar>
struct CollectedQuadrature {
    using range = TriangleRange<Scalar, depth>;
    static constexpr size_t n_points = BaseQuadRule::n_points * range::size();
    using point_type = barycentric_triangle<Scalar>;
    static constexpr std::array<point_type, n_points> make_points() {
        constexpr range r;
        return r.template collect_quadrature_points<BaseQuadRule>();
    }
    static constexpr std::array<Scalar, n_points> make_weights() {
        std::array<Scalar, n_points> result;
        for (size_t i = 0; i < n_points; ++i) {
            result[i] = BaseQuadRule::weights[i % BaseQuadRule::n_points] / range::size();
        }
        return result;
    }
    static constexpr std::array<point_type, n_points> points = make_points();
    static constexpr std::array<Scalar, n_points> weights = make_weights();
};
}  // namespace quadints
