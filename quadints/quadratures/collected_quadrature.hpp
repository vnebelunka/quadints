#ifndef QUADINTS_COLLECTED_QUADRATURE_HPP
#define QUADINTS_COLLECTED_QUADRATURE_HPP

#include "quadints/adaptive/splitter_triangle.hpp"
namespace quadints {
template <typename BaseQuadRule, size_t depth, typename Scalar>
struct CollectedQuadratureStatic {
    using range = TriangleRangeStatic<Scalar, depth>;
    static constexpr size_t n_points = BaseQuadRule::n_points * range::size();
    using point_type = barycentric_triangle<Scalar>;
    static constexpr auto make_points() -> std::array<point_type, n_points> {
        constexpr range r;
        return r.template collect_quadrature_points<BaseQuadRule>();
    }
    static constexpr auto make_weights() -> std::array<Scalar, n_points> {
        std::array<Scalar, n_points> result{};
        for (size_t i = 0; i < n_points; ++i) {
            result[i] = BaseQuadRule::weights[i % BaseQuadRule::n_points] / range::size();
        }
        return result;
    }
    static constexpr std::array<point_type, n_points> points = make_points();
    static constexpr std::array<Scalar, n_points> weights = make_weights();
};

template <typename BaseQuadRule, typename Scalar = std::decay_t<decltype(BaseQuadRule::weights[0])>>
class CollectedQuadratureDynamic {
    using range = TriangleRangeDynamic<Scalar>;
    size_t depth;
    size_t range_size;
    using point_type = barycentric_triangle<Scalar>;

    std::vector<point_type> _points;
    std::vector<Scalar> _weights;

   public:
    CollectedQuadratureDynamic(size_t depth) : depth(depth), range_size(range(depth).size()) {
        range const r(depth);
        _points = r.template collect_quadrature_points<BaseQuadRule>();
        _weights.resize(n_points());
        for (size_t i = 0; i < n_points(); ++i) {
            _weights[i] = BaseQuadRule::weights[i % BaseQuadRule::n_points] / range_size;
        }
    }
    [[nodiscard]] auto n_points() const -> size_t { return BaseQuadRule::n_points * range_size; }
    [[nodiscard]] auto points() const -> const std::vector<point_type>& { return _points; }
    [[nodiscard]] auto weights() const -> const std::vector<Scalar>& { return _weights; }
};

}  // namespace quadints

#endif  // QUADINTS_COLLECTED_QUADRATURE_HPP
