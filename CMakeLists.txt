cmake_minimum_required(VERSION 3.25)
project(lab_4_opp)
add_executable(lab_4_opp parallel.cpp)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(lab_4_opp ${MPI_LIBRARIES})
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS} --std=c++11 -O0")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")



