#Minimum version of CMake
cmake_minimum_required(VERSION 2.8)
#Project's name
project(Swarming)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

FIND_PACKAGE(OpenMP REQUIRED)
message(STATUS $ENV{OPENMP_FOUND})
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#SOURCES
add_subdirectory(src)

#TESTS
add_subdirectory(test)
