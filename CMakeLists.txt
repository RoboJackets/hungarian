cmake_minimum_required(VERSION 3.0.0)
project("Hungarian Algorithm")

set(CMAKE_CXX_STANDARD 11)

add_library(hungarian hungarian.cpp hungarian.hpp)

add_executable(demo test.cpp)
target_link_libraries(demo hungarian)
