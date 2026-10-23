#ifndef QUADINTS_SPLITTER_TRIANGLE_HPP
#define QUADINTS_SPLITTER_TRIANGLE_HPP
#include <array>
#include <cstddef>
#include <iterator>
#include <vector>

#include "quadints/quadratures/quads_triangle.hpp"

namespace {
constexpr void step_range(size_t& cur_pos, size_t& cur_level, bool& upper, size_t levels) {
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
}
template <typename Scalar>
using point = quadints::barycentric_triangle<Scalar>;

template <typename Scalar>
using bar_coords_t = std::array<point<Scalar>, 3>;

template <typename Scalar>
constexpr auto construct_triangle(size_t cur_pos, size_t cur_level, bool upper, Scalar step) -> bar_coords_t<Scalar> {
    point<Scalar> A;
    point<Scalar> B;
    point<Scalar> C;
    if (upper) {
        A = {cur_pos * step, cur_level * step};
        B = {(cur_pos + 1) * step, cur_level * step};
        C = {cur_pos * step, (cur_level + 1) * step};
    } else {
        A = {(cur_pos + 1) * step, cur_level * step};
        B = {(cur_pos + 1) * step, (cur_level + 1) * step};
        C = {cur_pos * step, (cur_level + 1) * step};
    }
    return {A, B, C};
}
}  // namespace

namespace quadints {

template <typename Scalar, size_t depth>
class TriangleIteratorStatic {
    using point = barycentric_triangle<Scalar>;

    using dir = barycentric_direction<Scalar>;
    static constexpr size_t levels = 1U << depth;
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

    constexpr TriangleIteratorStatic(size_t cur_level, size_t cur_pos, bool upper)
        : cur_level(cur_level),
          cur_pos(cur_pos),
          upper(upper),
          current_triangle(construct_triangle(cur_pos, cur_level, upper, step)) {}
    constexpr auto operator*() const -> reference { return current_triangle; }
    constexpr auto operator->() const -> pointer { return &current_triangle; }

    constexpr auto operator++() noexcept -> TriangleIteratorStatic& {
        step_range(cur_pos, cur_level, upper, levels);
        current_triangle = construct_triangle(cur_pos, cur_level, upper, step);
        return *this;
    }
    constexpr auto operator==(const TriangleIteratorStatic& other) const -> bool {
        return cur_level == other.cur_level && cur_pos == other.cur_pos && upper == other.upper;
    }
    constexpr auto operator!=(const TriangleIteratorStatic& other) const -> bool { return !(*this == other); }
};

template <typename Scalar>
class TriangleIteratorDynamic {
    using point = barycentric_triangle<Scalar>;
    using dir = barycentric_direction<Scalar>;
    size_t levels;
    Scalar step;
    using bar_coords_t = std::array<point, 3>;
    bar_coords_t current_triangle = {};
    size_t cur_level = 0;
    size_t cur_pos = 0;
    bool upper = true;
    size_t depth;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = bar_coords_t;
    using difference_type = std::ptrdiff_t;
    using reference = bar_coords_t&;
    using const_reference = const bar_coords_t&;
    using pointer = bar_coords_t*;
    using const_pointer = const bar_coords_t*;
    constexpr TriangleIteratorDynamic(size_t depth, size_t cur_level, size_t cur_pos, bool upper)
        : depth(depth),
          cur_level(cur_level),
          cur_pos(cur_pos),
          upper(upper),
          levels(1U << depth),
          step(1. / static_cast<Scalar>(levels)),
          current_triangle(construct_triangle(cur_pos, cur_level, upper, step)) {}
    constexpr auto operator*() const -> const_reference { return current_triangle; }
    constexpr auto operator->() const -> const_pointer { return &current_triangle; }

    constexpr auto operator++() noexcept -> TriangleIteratorDynamic& {
        step_range(cur_pos, cur_level, upper, levels);
        current_triangle = construct_triangle(cur_pos, cur_level, upper, step);
        return *this;
    }
    constexpr auto operator==(const TriangleIteratorDynamic& other) const -> bool {
        return cur_level == other.cur_level && cur_pos == other.cur_pos && upper == other.upper;
    }
    constexpr auto operator!=(const TriangleIteratorDynamic& other) const -> bool { return !(*this == other); }
};

template <typename Scalar, size_t depth>
struct TriangleRangeStatic {
    constexpr static auto size() -> size_t { return 1U << (2 * depth); }
    [[nodiscard]] constexpr auto begin() const -> TriangleIteratorStatic<Scalar, depth> {
        constexpr double step = 1. / static_cast<double>(1U << depth);

        return TriangleIteratorStatic<Scalar, depth>(0, 0, true);
    }
    [[nodiscard]] constexpr auto end() const -> TriangleIteratorStatic<Scalar, depth> {
        return TriangleIteratorStatic<Scalar, depth>{(1 << depth) - 1, 0, false};
    }
    [[nodiscard]] constexpr auto to_vector() const -> std::array<std::array<barycentric_triangle<Scalar>, 3>, size()> {
        std::array<std::array<barycentric_triangle<Scalar>, 3>, size()> result;
        int i = 0;
        for (auto it = this->begin(); it != this->end(); ++it, ++i) {
            result[i] = *it;
        }
        return result;
    }
    template <typename QuadRule>
    [[nodiscard]] constexpr auto collect_quadrature_points() const {
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

template <typename Scalar>
struct TriangleRangeDynamic {
    size_t depth;
    constexpr TriangleRangeDynamic(size_t depth) : depth(depth) {}
    [[nodiscard]] constexpr auto size() const -> size_t { return 1U << (2 * depth); }
    [[nodiscard]] constexpr auto begin() const -> TriangleIteratorDynamic<Scalar> {
        return TriangleIteratorDynamic<Scalar>(depth, 0, 0, true);
    }
    [[nodiscard]] constexpr auto end() const -> TriangleIteratorDynamic<Scalar> {
        return TriangleIteratorDynamic<Scalar>(depth, (1 << depth) - 1, 0, false);
    }
    [[nodiscard]] constexpr auto to_vector() const -> std::vector<std::array<barycentric_triangle<Scalar>, 3>> {
        std::vector<std::array<barycentric_triangle<Scalar>, 3>> result(size());
        int i = 0;
        for (auto it = this->begin(); it != this->end(); ++it, ++i) {
            result[i] = *it;
        }
        return result;
    }
    template <typename QuadRule>
    [[nodiscard]] constexpr auto collect_quadrature_points() const {
        using bar_point = barycentric_triangle<Scalar>;
        using bar_dir = barycentric_direction<Scalar>;
        std::vector<bar_point> bar_points(size() * QuadRule::n_points);
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
}  // namespace quadints
#endif  // QUADINTS_SPLITTER_TRIANGLE_HPP
