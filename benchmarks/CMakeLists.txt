find_package(benchmark REQUIRED)

add_executable(triangle_integration_bench
    triangle_integration.cpp
)

target_link_libraries(triangle_integration_bench
    PRIVATE
        integration_lib 
        benchmark::benchmark
)
