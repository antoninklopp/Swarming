#Minimum version of CMake
cmake_minimum_required(VERSION 2.8)
#Project's name
project(Swarming)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

#SOURCES
add_subdirectory(src)

#TESTS
# add_subdirectory(test)
