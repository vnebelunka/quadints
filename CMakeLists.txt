cmake_minimum_required(VERSION 3.24)
project(NumericalIntegration LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CheckCXXCompilerFlag)
include(CheckCXXSourceRuns)

include(CheckIPOSupported)
check_ipo_supported(RESULT lto_supported OUTPUT lto_error)

# ------------------------------
# Library (header‑only)
# ------------------------------

add_library(integration_lib INTERFACE)
target_include_directories(integration_lib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)
# 2. Определяем флаг компиляции под систему
if(MSVC)
    set(AVX512_FLAG "/arch:AVX512")
else()
    set(AVX512_FLAG "-mavx512f")
endif()

if(lto_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "LTO not supported, building without it: ${lto_error}")
endif()

# 3. Проверяем, поддерживается ли флаг компилятором
check_cxx_compiler_flag("${AVX512_FLAG}" COMPILER_SUPPORTS_AVX512)

# 4. Проверяем, может ли текущий процессор и ОС выполнить инструкцию AVX-512
if(COMPILER_SUPPORTS_AVX512)
    set(CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${AVX512_FLAG}")

    check_cxx_source_runs("
        #include <immintrin.h>
        int main() {
            #if defined(__GNUC__) || defined(__clang__)
                if (__builtin_cpu_supports(\"avx512f\")) {
                    __m512 v = _mm512_setzero_ps();
                    (void)v;
                    return 0;
                }
                return 1;
            #else
                __m512 v = _mm512_setzero_ps();
                (void)v;
                return 0;
            #endif
        }
    " HOST_SUPPORTS_AVX512)

    set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})
endif()

# 5. Если тест успешно пройден — пробрасываем флаги в INTERFACE
if(COMPILER_SUPPORTS_AVX512 AND HOST_SUPPORTS_AVX512)
    message(STATUS "AVX-512 успешно обнаружен и доступен для выполнения.")
    target_compile_options(integration_lib INTERFACE ${AVX512_FLAG})
    target_compile_definitions(integration_lib INTERFACE HAS_AVX512)
else()
    message(WARNING "AVX-512 не поддерживается компилятором или текущей системой.")
endif()

# ------------------------------
# Testing
# ------------------------------
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)

if(BUILD_TESTS)
  enable_testing()

  include(FetchContent)
  FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
    )
  FetchContent_MakeAvailable(googletest)

  add_executable(integration_tests
        test/segment_1d.cpp
        test/triangle_2d.cpp
        test/product.cpp
        test/split_triangle.cpp
        test/triangle_adaptive.cpp
        test/triangle_adaptive_callback.cpp
        test/dynamic_integration.cpp
        test/collect_with_buffer.cpp
    )
  target_link_libraries(integration_tests
        PRIVATE
            integration_lib
            gtest_main
    )

  # Optional: enable C++20 for the test executable explicitly
  target_compile_features(integration_tests PRIVATE cxx_std_20)

  # Discover tests
  include(GoogleTest)
  gtest_discover_tests(integration_tests)
endif()

if(BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()
