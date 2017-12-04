cmake_minimum_required(VERSION 3.8)
project(PMGBP_PDEVS)

set(CMAKE_CXX_STANDARD 14)

include_directories(
        include
        vendor/tinyxml2
        vendor/cadmium/include
        vendor/DESTimes/include
        vendor/DEVSDiagrammer/model_json_exporter)

set(SOURCES
        src/structures/types.cpp
        src/structures/space.cpp
        main.cpp
        vendor/tinyxml2/tinyxml2.cpp
        vendor/DEVSDiagrammer/model_json_exporter)

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories(include ${Boost_INCLUDE_DIRS})

set(TEST_SOURCES include test/cpp/libs)

enable_testing()
FILE(GLOB TestSources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/cpp/libs/*_test.cpp)
foreach(testSrc ${TestSources})
    get_filename_component(testName ${testSrc} NAME_WE)
    add_executable(${testName} test/cpp/libs/main-test.cpp ${testSrc} ${TEST_SOURCES})
    target_link_libraries(${testName} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    add_test(${testName} ${testName})
endforeach(testSrc)

add_executable(pmgbp ${SOURCES})
