#include <array>
#include <iostream>

template <unsigned int dim, unsigned int n_points> struct GaussianQuadrature {
  static_assert("We have no quadrature of such dimension and order");
};

struct Segment {
  using point_type = double;
  double start;
  double end;
  constexpr double mes() const { return end - start; }
};

struct barycentric_coord {
  std::array<double, 1> coords;
};

constexpr auto from_reference_domain(const barycentric_coord &reference_point,
                                     const Segment &s) -> Segment::point_type {
  return s.start + (s.end - s.start) * reference_point.coords[0];
}

template <> struct GaussianQuadrature<1, 1> {
  static constexpr std::array<barycentric_coord, 1> points{{0.5}};
  static constexpr std::array<double, 1> weights{1.};
  static constexpr int n_points = 1;
};

template <> struct GaussianQuadrature<1, 3> {
  static constexpr std::array<barycentric_coord, 3> points{
      barycentric_coord{0.}, {0.5}, {1.}};
  static constexpr std::array<double, 3> weights{1. / 6, 4. / 6, 1. / 6};
  static constexpr size_t n_points = 3;
};

template <typename quadrature, typename domain, typename func>
constexpr auto integrate(func &&f, const domain &cell)
    -> decltype(f(typename domain::point_type())) {
  using return_type = decltype(f(typename domain::point_type()));
  return_type res{0};
  for (int i = 0; i < quadrature::n_points; ++i) {
    res += f(from_reference_domain(quadrature::points[i], cell)) *
           quadrature::weights[i];
  }
  return res * cell.mes();
}

int main() {
  auto f = [](double x) { return x * x; };
  constexpr Segment s{0., 1.};
  constexpr auto I1 = integrate<GaussianQuadrature<1, 1>>(f, s);
  std::cout << integrate<GaussianQuadrature<1, 1>>(f, s) << '\n';
  std::cout << integrate<GaussianQuadrature<1, 3>>(f, s) << '\n';
}
