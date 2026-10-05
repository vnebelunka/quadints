cmake_minimum_required(VERSION 3.24) project(NumericalIntegration LANGUAGES CXX)

    set(CMAKE_CXX_STANDARD 20) set(CMAKE_CXX_STANDARD_REQUIRED ON) set(CMAKE_CXX_EXTENSIONS OFF)

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- --
#Library(header‑only)
#-- -- -- -- -- -- -- -- -- -- -- -- -- -- --

        add_library(integration_lib INTERFACE) target_include_directories(
            integration_lib INTERFACE $<BUILD_INTERFACE : ${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE : include>)

#-- -- -- -- -- -- -- -- -- -- -- -- -- -- --
#Testing
#-- -- -- -- -- -- -- -- -- -- -- -- -- -- --
            option(BUILD_TESTS "Build unit tests" ON) option(BUILD_BENCHMARKS "Build benchmarks" ON)

                if (BUILD_TESTS) enable_testing()

                    include(FetchContent) FetchContent_Declare(
                        googletest URL https:  // github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
                        ) FetchContent_MakeAvailable(googletest)

                        add_executable(integration_tests test / segment_1d.cpp test / triangle_2d.cpp test /
                                       product.cpp test / split_triangle.cpp test / triangle_adaptive.cpp test /
                                       triangle_adaptive_callback.cpp)
                            add_executable(dynamic_integration_test test / dynamic_integration.cpp)
                                target_link_libraries(integration_tests PRIVATE integration_lib gtest_main)
                                    target_link_libraries(dynamic_integration_test PRIVATE integration_lib gtest_main)

                                        target_compile_features(integration_tests PRIVATE cxx_std_20)

                                            include(GoogleTest) gtest_discover_tests(integration_tests)
                                                gtest_discover_tests(dynamic_integration_test) endif()

                                                    if (BUILD_BENCHMARKS) add_subdirectory(benchmarks) endif()
