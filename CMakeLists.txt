cmake_minimum_required(VERSION 3.8)
project(PMGBP_PDEVS)

set(CMAKE_CXX_STANDARD 14)

include_directories(../cadmium/include ../DESTimes atomic-models libs)

set(SOURCES libs/types.cpp libs/unit_definition.cpp main.cpp atomic-models)

add_executable(pmgbp ${SOURCES})