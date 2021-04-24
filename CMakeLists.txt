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

# looking for spglib
find_library(SPGLIB symspg)
if(NOT SPGLIB)
  message(FATAL_ERROR "Could not find an installation of the space group library (spglib).")
else()
   message("-- Found spglib ${SPGLIB}")
endif()

# looking for OpenMP
option(USE_OPENMP "build OpenMP programs with OpenMP support enabled" ON)
if(USE_OPENMP)
  find_package(OpenMP 4.5 REQUIRED)
endif()

find_library(PYBIND symspg)
if(NOT SPGLIB)
  message(FATAL_ERROR "Could not find an installation of the space group library (spglib).")
else()
   message("-- Found spglib ${SPGLIB}")
endif()

find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(hetbuilder_backend ${PROJECT_SOURCE_DIR}/backend/pybindings.cpp ${SOURCES})

target_link_libraries(hetbuilder_backend PRIVATE ${Python3_LIBRARY} OpenMP::OpenMP_CXX ${SPGLIB})


