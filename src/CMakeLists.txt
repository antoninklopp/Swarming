cmake_minimum_required(VERSION 2.8)

include_directories(library)
add_library(lib_boids SHARED agent.cxx vector.cxx workspace.cxx)

add_executable(main main.cxx)
target_link_libraries(main lib_boids)
