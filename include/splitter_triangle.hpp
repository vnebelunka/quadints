#include "quads_triangle.hpp"
#include <cstddef>
#include <iterator>
#include <type_traits>

namespace quadints {

template <typename Scalar, size_t depth> struct TriangleIterator {

  using point = barycentric_triangle<Scalar>;
  using dir = barycentric_direction<Scalar>;
  static constexpr const double step = (1. / static_cast<double>(1u << depth));
  static constexpr auto dir_left = dir{-step, step, 0};
  static constexpr auto dir_down = dir{0, +step, -step};
  static constexpr auto dir_up = dir{-step, 0, step};
  using bar_coords_t = std::array<point, 3>;
  bar_coords_t current_triangle = {};
  size_t pos = 0;

public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = bar_coords_t;
  using difference_type = std::ptrdiff_t;
  using pointer = const value_type *;
  using reference = const value_type &;

  constexpr TriangleIterator(const bar_coords_t &coords, size_t pos = 0)
      : current_triangle(coords), pos(pos) {}
  constexpr reference operator*() const { return current_triangle; }
  constexpr pointer operator->() const { return &current_triangle; }

  constexpr TriangleIterator &operator++() noexcept {
    ++pos;

    auto &cur_right_point = current_triangle[0];
    auto &cur_left_point = current_triangle[1];
    // if level != 1 and triangle was above edge
    if (current_triangle[2].z() > current_triangle[0].z() &&
        current_triangle[0].z() != 0) {
      current_triangle[2] = cur_right_point + dir_down;
      return *this;
    }
    cur_right_point += dir_left;
    cur_left_point += dir_left;
    // can't go left -> go to next level
    if (cur_left_point.x() < 0) {
      cur_right_point = {1 - current_triangle[0].z() - step, 0.,
                         current_triangle[0].z() + step};
      cur_left_point = cur_right_point + dir_left;
    }
    current_triangle[2] = cur_right_point + dir_up;
    return *this;
  }
  constexpr bool operator==(const TriangleIterator &other) const {
    return pos == other.pos;
  }
  constexpr bool operator!=(const TriangleIterator &other) const {
    return !(*this == other);
  }
};

template <typename Scalar, size_t depth> struct TriangleRange {
  constexpr size_t size() const { return 1u << (2 * depth); }

  constexpr TriangleIterator<Scalar, depth> begin() const {
    constexpr double step = 1. / static_cast<double>(1u << depth);
    return TriangleIterator<Scalar, depth>(
        {barycentric_triangle<Scalar>{1., 0., 0},
         {1 - step, step, 0.},
         {1 - step, 0., step}},
        0);
  }
  constexpr TriangleIterator<Scalar, depth> end() const {
    return TriangleIterator<Scalar, depth>{
        {barycentric_triangle<Scalar>{0., 0.}, {0., 0.}, {0., 0.}}, size()};
  }
  constexpr auto to_vector() const
      -> std::array<std::array<barycentric_triangle<Scalar>, 3>,
                    1u << (2 * depth)> {
    std::vector<std::array<barycentric_triangle<Scalar>, 3>> tmp_vec;
    for (const auto &tri : *this) {
      tmp_vec.push_back(tri);
    }
    std::array<std::array<barycentric_triangle<Scalar>, 3>, 1u << (2 * depth)>
        result;
    std::copy(tmp_vec.begin(), tmp_vec.end(), result.begin());
    return result;
  }
};
} // namespace quadints
