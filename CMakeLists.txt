cmake_minimum_required(VERSION 3.8)
project(PMGBP_PDEVS)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        atomic-models/biomass.hpp
        atomic-models/reaction.hpp
        atomic-models/filter.hpp
        atomic-models/reaction.hpp
        atomic-models/space.hpp
        libs/message_types.hpp
        libs/randomNumbers.hpp
        libs/types.cpp
        libs/types.hpp
        libs/unit_definition.cpp
        libs/unit_definition.hpp
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