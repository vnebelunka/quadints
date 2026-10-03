find_package(benchmark REQUIRED)

add_executable(triangle_integration_bench
    triangle_integration.cpp
)

target_link_libraries(triangle_integration_bench
    PRIVATE
        integration_lib
        benchmark::benchmark
)

add_executable(adaptive_integration_bench
    adaptive_integration.cpp
)

target_link_libraries(adaptive_integration_bench
    PRIVATE
        integration_lib
        benchmark::benchmark
)
