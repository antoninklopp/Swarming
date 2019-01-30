#Minimum version of CMake
cmake_minimum_required(VERSION 2.8)

include_directories(${PROJECT_SOURCE_DIR}/src)
message(${PROJECT_SOURCE_DIR})

#EXECUTABLES

add_executable(test_time test_time.cpp)

# LINKING LIBRARIES

target_link_libraries(test_time lib_boids)
