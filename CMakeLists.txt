# cmake_minimum_required(VERSION <specify CMake version here>)
cmake_minimum_required(VERSION 3.13)

project(EvoSymbolicRegressionClion)

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/boost_1_71_0)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/boost_1_71_0)

add_executable(EvoSymbolicRegressionClion main.cpp linearRegModel.cpp linearRegModel.h)