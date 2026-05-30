cmake_minimum_required(VERSION 3.14)
project(mapf_benchmark CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

add_executable(mapf_benchmark src/benchmark_main.cpp)
target_include_directories(mapf_benchmark PRIVATE include)
target_compile_options(mapf_benchmark PRIVATE -O2 -Wall)
