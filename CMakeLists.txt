cmake_minimum_required(VERSION 3.18.4)

# set the project name
project(hetbuilder)
message("-- Project source directory: ${PROJECT_SOURCE_DIR}")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17 -fopenmp")

# specify sources
file(GLOB SOURCES    
    ${PROJECT_SOURCE_DIR}/backend/*.cpp
)

# target name
set(TARGET hetbuilder_backend)

if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE MATCHES RELEASE))
    set(CMAKE_BUILD_TYPE RELEASE)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

if((CMAKE_BUILD_TYPE MATCHES DEBUG) OR (CMAKE_BUILD_TYPE MATCHES Debug))
  set(CMAKE_BUILD_TYPE DEBUG)
  set(CMAKE_CXX_FLAGS_RELEASE "-O0")
endif()

message("-- CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

# looking for python
find_package(Python3 COMPONENTS Interpreter Development)
include_directories(${Python3_INCLUDE_DIRS})

# looking for git
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

# looking for pybind11
find_package(pybind11 CONFIG)
if(NOT pybind11_FOUND)
    message("-- Could not find an installation of pybind11.")
    message("-- Linking to submodule instead.")
    add_subdirectory(${PROJECT_SOURCE_DIR}/external/pybind11)
endif()

pybind11_add_module(${TARGET} ${PROJECT_SOURCE_DIR}/backend/pybindings.cpp ${SOURCES})

#looking for spglib
find_library(SPGLIB symspg libsymspg libsymspg.a ${Python3_LIBRARY})
if(NOT SPGLIB)
    message("-- Could not find an installation of the space group library (spglib).")
    message("-- Building from submodule instead. This might take a while.")
    set(SPGLIB_SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/spglib)    
    add_subdirectory(${SPGLIB_SOURCE_DIR})    
    target_include_directories(${TARGET} PRIVATE ${SPGLIB_SOURCE_DIR}/src)
    set(SPGLIB symspg)
else()
    message("-- Found spglib ${SPGLIB}")
endif()

# looking for OpenMP
option(USE_OPENMP "build OpenMP programs with OpenMP support enabled" ON)
if(USE_OPENMP)
  find_package(OpenMP 4.5 REQUIRED)
endif()

target_link_libraries(${TARGET} PRIVATE ${Python3_LIBRARY} OpenMP::OpenMP_CXX ${SPGLIB})


