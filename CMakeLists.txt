cmake_minimum_required(VERSION 3.21)
project(PollutantSimulation)

set(CMAKE_CXX_STANDARD 14)

add_executable(PollutantSimulation main.cpp)
add_library(Libraries STATIC Libraries.cpp)
target_link_libraries(${PROJECT_NAME} Libraries)