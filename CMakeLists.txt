cmake_minimum_required(VERSION 3.21)
project(PollutantSimulation)

set(CMAKE_CXX_STANDARD 14)

add_executable(PollutantSimulation main.cpp)
add_library(Libraries STATIC Libraries.cpp)
add_library(compute_FG STATIC compute_FG.cpp)
add_library(compute_poisson STATIC compute_poisson.cpp)
add_library(compute_uv STATIC compute_uv.cpp)
add_library(compute_phi STATIC compute_phi.cpp)
add_library(math_diff STATIC math_diff.cpp)

target_link_libraries(${PROJECT_NAME} compute_FG)
target_link_libraries(${PROJECT_NAME} compute_poisson)
target_link_libraries(${PROJECT_NAME} compute_uv)
target_link_libraries(${PROJECT_NAME} math_diff)
target_link_libraries(${PROJECT_NAME} compute_phi)
target_link_libraries(${PROJECT_NAME} Libraries)