cmake_minimum_required(VERSION 3.17.2)

# set the project name
project(hetbuilder)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17 -fopenmp")

# specify sources
file(GLOB SOURCES    
    ${PROJECT_SOURCE_DIR}/hetbuilder/backend/*.cpp
)

if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE MATCHES RELEASE))
  set(CMAKE_BUILD_TYPE RELEASE)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  # add the dynamic library
  set(TARGET hetbuild)
  add_library(${TARGET} MODULE ${SOURCES})  
  message("-- TARGET LIBRARY: ${TARGET}")
endif()

if(CMAKE_BUILD_TYPE MATCHES TESTING)
  message("-- BUILD TYPE: ${CMAKE_BUILD_TYPE}")
  set(TARGET hetbuild.exe)
  add_executable(${TARGET} ${SOURCES})  
  message("-- TARGET EXECUTABLE: ${TARGET}")
endif()


find_package(Python3 COMPONENTS Interpreter Development)
#find_package(PythonLibs REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})

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


message("-- Project source directory: ${PROJECT_SOURCE_DIR}")

target_link_libraries(${TARGET} PUBLIC ${Python3_LIBRARY} OpenMP::OpenMP_CXX ${SPGLIB})

set(PYBINDDIR ${PROJECT_SOURCE_DIR}/external/pybind11)
message("-- Pybind11 source directory: ${PYBINDDIR}")
add_subdirectory(${PROJECT_SOURCE_DIR}/external/pybind11)
pybind11_add_module(coincidence_algorithm ${PROJECT_SOURCE_DIR}/hetbuilder/backend/pybindings.cpp)

if(CMAKE_BUILD_TYPE MATCHES RELEASE)
  set_target_properties(coincidence_algorithm
      PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/hetbuilder"
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/hetbuilder"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/hetbuilder"
  )
endif()