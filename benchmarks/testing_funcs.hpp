#ifndef QUADINTS_BENCHES_TESTING_FUNCS
#define QUADINTS_BENCHES_TESTING_FUNCS
#include <complex>

#include "triangle_utility.hpp"

template <size_t p>
struct xpyp {
    double operator()(const point2d& point) const { return std::pow(point.coords[0] + point.coords[1], p); }
    auto operator()(const point2d& pointx, const point2d& pointy) const -> double {
        return std::pow(pointx.coords[0] + pointy.coords[1], p);
    }
};
struct constant {
    double operator()(const point2d& point) const { return 1.0; }
    auto operator()(const point2d& pointx, const point2d& pointy) const -> double { return 1.0; }
};

struct cos10 {
    double w = 10;
    double operator()(const point2d& point) const { return cos(w * (point.coords[0] + point.coords[1])); }
};

struct expir {
    auto operator()(const point2d& point) -> std::complex<double> {
        return std::exp(norm(point) * std::complex<double>{0., 1.});
    }
    auto operator()(const point2d& pointx, const point2d& pointy) -> std::complex<double> {
        auto r = norm(pointx - pointy);
        return std::exp(r * std::complex<double>{0., 1.});
    }
};

struct norm_func {
    double operator()(const point2d& pointx, const point2d& pointy) const {
        return std::sqrt(std::pow(pointx.coords[0] - pointy.coords[0], 2) +
                         std::pow(pointx.coords[1] - pointy.coords[1], 2));
    }
};

struct constant_func {
    double operator()(const point2d& pointx, const point2d& pointy) const { return 1.0; }
};

struct expir_func {
    auto operator()(const point2d& pointx, const point2d& pointy) -> std::complex<double> {
        return std::exp(norm(pointx - pointy) * std::complex<double>{0., 1.});
    }
};

#endif  // QUADINTS_BENCHES_TESTING_FUNCS
