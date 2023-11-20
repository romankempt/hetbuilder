cmake_minimum_required(VERSION 3.21.3)

# set the project name
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

message(STATUS "Project source directory: ${PROJECT_SOURCE_DIR}")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# specify sources
file(GLOB SOURCES    
    ${CMAKE_CURRENT_SOURCE_DIR}/backend/*.cpp
)

# target name
set(TARGET hetbuilder_backend)

if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE MATCHES RELEASE))
    set(CMAKE_BUILD_TYPE RELEASE)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=gnu++17")
endif()

if((CMAKE_BUILD_TYPE MATCHES DEBUG) OR (CMAKE_BUILD_TYPE MATCHES Debug))
  set(CMAKE_BUILD_TYPE DEBUG)
  set(CMAKE_CXX_FLAGS_RELEASE "-O0 -std=gnu++17")
endif()

# Set -fPIC on x86_64
if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC"  )
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC"  )
endif("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

# looking for python
find_package(Python REQUIRED COMPONENTS Interpreter Development)
include_directories(${Python_INCLUDE_DIRS})
message(STATUS "Found Python header files ${Python_INCLUDE_DIRS}")

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

#looking for spglib
set(spglib_FOUND FALSE)
if (NOT ${CMAKE_BUILD_TYPE} STREQUAL DEBUG)
    find_library(spglib symspg libsymspg libsymspg.a ${Python3_LIBRARY})
    set(spglib_FOUND TRUE)
endif()
if(NOT spglib_FOUND)
    message(STATUS "Could not find an installation of the space group library (spglib).")
    message(STATUS "Building spglib from submodule.")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/spglib)
else()
    message(STATUS "Found spglib ${spglib}")    
endif()

# looking for pybind11
set(pybind11_FOUND FALSE)
if (NOT ${CMAKE_BUILD_TYPE} STREQUAL DEBUG)
    find_package(pybind11 CONFIG REQUIRED)
endif()
if(NOT pybind11_FOUND)
    message(STATUS "Could not find an installation of pybind11.")
    message(STATUS "Linking to submodule instead.")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11)
endif()

Python_add_library(${TARGET} MODULE ${CMAKE_CURRENT_SOURCE_DIR}/backend/pybindings.cpp ${SOURCES} WITH_SOABI)
#pybind11_add_module(${TARGET} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/backend/pybindings.cpp ${SOURCES})

#adding XtalComp
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/XtalComp)
target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/XtalComp)

# looking for OpenMP
option(USE_OPENMP "Build with OpenMP support" ON)
if(USE_OPENMP)
    find_package(OpenMP 4.5)
    if(OPENMP_FOUND)
        message(STATUS "Found OpenMP ${OPENMP}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")    
        target_link_libraries(${TARGET} 
                                PRIVATE symspg
                                PRIVATE XtalComp
                                PRIVATE OpenMP::OpenMP_CXX  
                                PRIVATE pybind11::module pybind11::lto pybind11::windows_extras                                                              
                                )
    else()
        message(STATUS "Could not find OpenMP")
    endif()
else()
    message(STATUS "Building without OpenMP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    target_link_libraries(${TARGET} 
                            PRIVATE symspg
                            PRIVATE XtalComp
                            PRIVATE pybind11::module pybind11::lto pybind11::windows_extras
                            )
endif()

install(TARGETS ${TARGET} DESTINATION hetbuilder)

# Quiet a warning, since this project is only valid with SKBUILD
# set(ignoreMe "${SKBUILD}")