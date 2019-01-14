cmake_minimum_required(VERSION 3.8)
project(PMGBP_PDEVS)

set(CMAKE_CXX_STANDARD 17)

SET(SHOW_INFO "-D show_info")
SET(CMAKE_CXX_FLAGS  "${SHOW_INFO}")

include_directories(
        include/pmgbp
        vendor/tinyxml2
        vendor/cadmium/include
        vendor/DESTimes/include
        vendor/MeMoRe/include
        vendor/DEVSDiagrammer/model_json_exporter/include)

set(SOURCES
        src/pmgbp/structures/types.cpp
        src/pmgbp/structures/space.cpp
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

find_package(libmongocxx REQUIRED)
find_package(libbsoncxx REQUIRED)
include_directories(${LIBMONGOCXX_INCLUDE_DIR})
include_directories(${LIBBSONCXX_INCLUDE_DIRS})
include_directories("/usr/local/include/mongocxx/v_noabi")
include_directories("/usr/local/include/bsoncxx/v_noabi")
include_directories("/usr/local/include/libmongoc-1.0")
include_directories("/usr/local/include/libbson-1.0")
include_directories("/usr/local/lib")

add_executable(pmgbp ${SOURCES})

target_link_libraries(pmgbp ${LIBMONGOCXX_LIBRARIES})
target_link_libraries(pmgbp ${LIBBSONCXX_LIBRARIES})
