cmake_minimum_required(VERSION 3.17.2)

# set the project name
project(heterostructure-builder)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17 -fopenmp")

# looking for OpenMP
option(USE_OPENMP "build OpenMP programs with OpenMP support enabled" ON)
if(USE_OPENMP)
  find_package(OpenMP 4.5 REQUIRED)
endif()


# specify sources
file(GLOB SOURCES    
    backend/*.cpp
)

if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE MATCHES Release))
  set(CMAKE_BUILD_TYPE Release)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  # add the dynamic library
  set(TARGET coincidence_algorithm)
  add_library(${TARGET} MODULE ${SOURCES})  
  message("-- TARGET: ${TARGET}")
endif()

if(CMAKE_BUILD_TYPE MATCHES TESTING)
  message("-- BUILD TYPE: ${CMAKE_BUILD_TYPE}")
  set(TARGET coincidence_algorithm.exe)
  add_executable(${TARGET} ${SOURCES})  
  message("-- TARGET: ${TARGET}")
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

#add_subdirectory(interfacebuilder/backend/external/pybind11)
target_link_libraries(${TARGET} PUBLIC ${Python3_LIBRARY} OpenMP::OpenMP_CXX symspg)

if(CMAKE_BUILD_TYPE MATCHES Release)
  #pybind11_add_module(backend ${SOURCES})
  set_target_properties(coincidence_algorithm
      PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/heterostructure_builder"
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/heterostructure_builder"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/heterostructure_builder"
  )
endif()