cmake_minimum_required(VERSION 3.8)
project(PMGBP_PDEVS)

set(CMAKE_CXX_STANDARD 14)

include_directories(
        ../cadmium/include
        ../DESTimes
        atomics
        libs
        structures
        vendor/tinyxml2
        vendor/DEVSDiagramer/model_json_exporter)

set(SOURCES
        structures/types.cpp
        structures/space.cpp
        libs/unit_definition.cpp
        vendor/tinyxml2/tinyxml2.cpp
        main.cpp
        atomics
        vendor/DEVSDiagramer/model_json_exporter)

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories(include ${Boost_INCLUDE_DIRS})

set(TEST_SOURCES libs structures test/cpp/libs)

enable_testing()
FILE(GLOB TestSources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/cpp/libs/*_test.cpp)
foreach(testSrc ${TestSources})
    get_filename_component(testName ${testSrc} NAME_WE)
    add_executable(${testName} test/cpp/libs/main-test.cpp ${testSrc} ${TEST_SOURCES})
    target_link_libraries(${testName} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    add_test(${testName} ${testName})
endforeach(testSrc)

add_executable(pmgbp ${SOURCES})
