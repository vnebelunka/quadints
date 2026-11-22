#ifndef QUADINTS_ALIGNED_ARRAY
#define QUADINTS_ALIGNED_ARRAY
#include <array>
#include <concepts>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>

namespace quadints::common {

#if defined(__AVX512F__) || defined(HAS_AVX512)
inline constexpr std::size_t simd_alignment = 64;
#else
inline constexpr std::size_t simd_alignment = 1;
#endif

template <typename T>
inline constexpr std::size_t vectorization_alignment = simd_alignment > alignof(T) ? simd_alignment : alignof(T);

template <typename T, std::size_t N, std::size_t Alignment>
class aligned_array {
    static_assert(N > 0);
    static_assert(Alignment >= alignof(T));

    alignas(Alignment) T elements[N];

   public:
    using value_type = T;
    using iterator = T*;
    using const_iterator = const T*;

    static constexpr std::size_t extent = N;
    static constexpr std::size_t alignment = Alignment;

    constexpr T* data() noexcept { return elements; }
    constexpr const T* data() const noexcept { return elements; }

    constexpr iterator begin() noexcept { return elements; }
    constexpr const_iterator begin() const noexcept { return elements; }

    constexpr iterator end() noexcept { return elements + N; }
    constexpr const_iterator end() const noexcept { return elements + N; }

    constexpr T& operator[](std::size_t i) noexcept { return elements[i]; }

    constexpr const T& operator[](std::size_t i) const noexcept { return elements[i]; }

    constexpr std::size_t size() const noexcept { return N; }
};

template <typename T, std::size_t N>
using vectorized_array = aligned_array<T, N, vectorization_alignment<T>>;

template <typename T>
struct aligned_array_traits {
    static constexpr bool is_aligned = false;
    static constexpr std::size_t extent = 0;
    static constexpr std::size_t alignment = 0;

    using value_type = void;
};

template <typename T, std::size_t N, std::size_t A>
struct aligned_array_traits<aligned_array<T, N, A>> {
    static constexpr bool is_aligned_array = true;
    static constexpr std::size_t extent = N;
    static constexpr std::size_t alignment = A;

    using value_type = T;
};

template <typename R, std::size_t N, std::size_t Alignment>
concept is_aligned_array = aligned_array_traits<std::remove_cvref_t<R>>::is_aligned_array &&
                           aligned_array_traits<std::remove_cvref_t<R>>::extent == N &&
                           aligned_array_traits<std::remove_cvref_t<R>>::alignment == Alignment &&
                           alignof(std::remove_cvref_t<R>) == Alignment;

}  // namespace quadints::common
#endif  // QUADINTS_ALIGNED_ARRAY;
