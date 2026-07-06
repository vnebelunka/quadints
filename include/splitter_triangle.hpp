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
    static constexpr size_t levels = 1u << depth;
    static constexpr Scalar step = (1. / static_cast<Scalar>(levels));
    using bar_coords_t = std::array<point, 3>;
    bar_coords_t current_triangle = {};
    size_t cur_level = 0;
    size_t cur_pos = 0;
    bool upper = true;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = bar_coords_t;
    using difference_type = std::ptrdiff_t;
    using pointer = const value_type*;
    using reference = const value_type&;

    constexpr TriangleIterator(size_t cur_level, size_t cur_pos, bool upper)
        : cur_level(cur_level), cur_pos(cur_pos), upper(upper) {
        point A{cur_pos * step, cur_level * step};
        point B{(cur_pos + 1) * step, cur_level * step};
        point C{cur_pos * step, (cur_level + 1) * step};
        current_triangle = bar_coords_t{A, B, C};
    }
    constexpr reference operator*() const { return current_triangle; }
    constexpr pointer operator->() const { return &current_triangle; }

    constexpr TriangleIterator& operator++() noexcept {
        ++cur_pos;
        if (upper) {
            if (cur_pos + cur_level >= levels) {
                ++cur_level;
                cur_pos = 0;
                if (cur_level == levels) {
                    upper = false;
                    cur_level = 0;
                };
            }
        } else {
            if (cur_pos + cur_level >= levels - 1) {
                ++cur_level;
                cur_pos = 0;
            }
        }
        point A, B, C;
        if (upper) {
            A = {cur_pos * step, cur_level * step};
            B = {(cur_pos + 1) * step, cur_level * step};
            C = {cur_pos * step, (cur_level + 1) * step};
        } else {
            A = {(cur_pos + 1) * step, cur_level * step};
            B = {(cur_pos + 1) * step, (cur_level + 1) * step};
            C = {cur_pos * step, (cur_level + 1) * step};
        }
        current_triangle = bar_coords_t{A, B, C};
        return *this;
    }
    constexpr bool operator==(const TriangleIterator& other) const {
        return cur_level == other.cur_level && cur_pos == other.cur_pos && upper == other.upper;
    }
    constexpr bool operator!=(const TriangleIterator& other) const { return !(*this == other); }
};

template <typename Scalar, size_t depth>
struct TriangleRange {
    constexpr static size_t size() { return 1u << (2 * depth); }
    constexpr TriangleIterator<Scalar, depth> begin() const {
        constexpr double step = 1. / static_cast<double>(1u << depth);

        return TriangleIterator<Scalar, depth>(0, 0, true);
    }
    constexpr TriangleIterator<Scalar, depth> end() const {
        return TriangleIterator<Scalar, depth>{(1 << depth) - 1, 0, false};
    }
    constexpr auto to_vector() const -> std::array<std::array<barycentric_triangle<Scalar>, 3>, size()> {
        std::array<std::array<barycentric_triangle<Scalar>, 3>, size()> result;
        int i = 0;
        for (auto it = this->begin(); it != this->end(); ++it, ++i) {
            result[i] = *it;
        }
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
