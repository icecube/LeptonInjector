# CMakeLists.txt for pip install (scikit-build-core)
# This is an alternative CMake configuration optimized for pip installation
# that uses pybind11 instead of Boost.Python

cmake_minimum_required(VERSION 3.15...3.27)
project(LeptonInjector VERSION 1.1.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 99)

# Find required packages
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# HDF5
find_package(HDF5 REQUIRED)
if(HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_INCLUDE_DIRS}")
endif()

# Boost (headers only, not Boost.Python)
find_package(Boost REQUIRED)
if(Boost_FOUND)
    message(STATUS "Found Boost: ${Boost_INCLUDE_DIRS}")
endif()

# photospline
find_package(photospline REQUIRED)
if(photospline_FOUND)
    message(STATUS "Found photospline")
endif()

# Core library sources
set(CORE_SOURCES
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/BasicInjectionConfiguration.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/LeptonInjector.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/Controller.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/Coordinates.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/DataWriter.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/EventProps.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/Particle.cxx
    ${CMAKE_SOURCE_DIR}/private/LeptonInjector/Random.cxx
    ${CMAKE_SOURCE_DIR}/private/earthmodel-service/EarthModelCalculator.cxx
    ${CMAKE_SOURCE_DIR}/private/earthmodel-service/EarthModelService.cxx
    ${CMAKE_SOURCE_DIR}/private/phys-services/LICrossSection.cxx
)

# Build core library as static for embedding in the Python module
add_library(lepton_injector_static STATIC ${CORE_SOURCES})
target_include_directories(lepton_injector_static PUBLIC
    ${CMAKE_SOURCE_DIR}/public/
    ${HDF5_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)
target_link_libraries(lepton_injector_static PUBLIC
    ${HDF5_LIBRARIES}
    photospline
)
set_target_properties(lepton_injector_static PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

# Python module with pybind11
pybind11_add_module(LeptonInjector
    ${CMAKE_SOURCE_DIR}/private/pybindings_pybind11/LeptonInjector.cxx
)
target_link_libraries(LeptonInjector PRIVATE lepton_injector_static)
target_include_directories(LeptonInjector PRIVATE
    ${CMAKE_SOURCE_DIR}/public/
)

# Install the Python module
install(TARGETS LeptonInjector LIBRARY DESTINATION LeptonInjector)

# Install the Python package files
install(DIRECTORY ${CMAKE_SOURCE_DIR}/resources/python/LeptonInjector/
    DESTINATION LeptonInjector
    FILES_MATCHING PATTERN "*.py"
)

# Install resource files (earth model parameters, etc.)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/resources/earthparams
    DESTINATION LeptonInjector/resources
)
