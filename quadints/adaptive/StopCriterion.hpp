#ifndef QUADINTS_STOP_CRITERION_HPP
#define QUADINTS_STOP_CRITERION_HPP
#include "quadints/interface.hpp"

namespace quadints {
template <typename Scalar>
class CurrentIntegral {
    Scalar value;

   public:
    explicit CurrentIntegral(Scalar value) : value(value) {}
    operator Scalar() const { return value; }
};

template <typename Scalar>
class PreviousIntegral {
    Scalar value;

   public:
    explicit PreviousIntegral(Scalar value) : value(value) {}
    operator Scalar() const { return value; }
};
template <typename Scalar>
class Atol {
    Scalar value;

   public:
    explicit Atol(Scalar value) : value(value) {}
    operator Scalar() const { return value; }
};
template <typename Scalar>
class Rtol {
    Scalar value;

   public:
    explicit Rtol(Scalar value) : value(value) {}
    operator Scalar() const { return value; }
};

enum class StopCriterionType : bool { STOP, CONTINUE };

template <typename T, typename Scalar>
concept StopCriterion = requires(T t, CurrentIntegral<Scalar> cur, PreviousIntegral<Scalar> prev) {
    { t(cur, prev) } -> std::convertible_to<StopCriterionType>;
};

struct NonAdaptiveCriterion {
    NonAdaptiveCriterion() = default;
    template <typename Scalar>
    auto operator()(CurrentIntegral<Scalar> /*unused*/,
                    PreviousIntegral<Scalar> /*unused*/) const -> StopCriterionType {
        return StopCriterionType::STOP;
    }
};

template <typename Scalar>
struct DefaultCriterion {
    Scalar atol;
    Scalar rtol;

    DefaultCriterion(Atol<Scalar> atol, Rtol<Scalar> rtol) : atol(atol), rtol(rtol) {}

    template <typename ReturnType>
    auto operator()(CurrentIntegral<ReturnType> current,
                    PreviousIntegral<ReturnType> previous) const -> StopCriterionType {
        auto error =
            magnitude<ReturnType, Scalar>(static_cast<ReturnType>(current) - static_cast<ReturnType>(previous));
        return error < atol + (rtol * magnitude<ReturnType, Scalar>(static_cast<ReturnType>(current)))
                   ? StopCriterionType::STOP
                   : StopCriterionType::CONTINUE;
    }
};
}  // namespace quadints

#endif  // QUADINTS_STOP_CRITERION_HPP
