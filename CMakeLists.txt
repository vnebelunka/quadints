cmake_minimum_required(VERSION 3.20)
project(NumericalIntegration LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------------
# Library (header‑only)
# ------------------------------

add_library(integration_lib INTERFACE)
target_include_directories(integration_lib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# ------------------------------
# Testing
# ------------------------------
option(BUILD_TESTS "Build unit tests" ON)

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
