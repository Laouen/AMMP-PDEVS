cmake_minimum_required(VERSION 3.8)
project(PMGBP_PDEVS)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        atomic-models/biomass.hpp
        atomic-models/cdboost-reaction.hpp
        atomic-models/filter.hpp
        atomic-models/reaction.hpp
        atomic-models/space.hpp
        data-structures/message_types.hpp
        data-structures/randomNumbers.hpp
        data-structures/types.cpp
        data-structures/types.hpp
        data-structures/unit_definition.cpp
        data-structures/unit_definition.hpp
        model_generator/model_generator.hpp
        model_generator/model_generator.py
        model_generator/parser.py
        parser/parser.cpp
        parser/parser.hpp
        test-files/test_1.xml
        test-files/test_2.xml
        test-files/test_3.xml
        test-files/test_4.xml
        main-ammp.cpp
        main.cpp
        makefile
        model-engine.hpp
        msb201165-sup-0003.xml
        README.md)

add_executable(PMGBP_PDEVS ${SOURCE_FILES})