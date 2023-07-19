# Ian Wark July 6 2017

cmake_minimum_required(VERSION 3.1)
project(DinoKnot_multimodel)

# enable c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \ -std=c++17 \ -g")

# create simfold library
add_subdirectory(simfold)

set(DinoKnot_multimodel_SOURCE
    W_final.cpp 
    pseudo_loop.cpp 
    V_final.cpp 
    VM_final.cpp
    h_common.cpp 
    hfold_validation.cpp
    DinoKnot.cpp
    Result.cpp
    cmdline.cc
)

# pass SIMFOLD_HOME to binary
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# get pthread library
find_package (Threads)

# set path for executable to be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executables
add_executable(DinoKnot ${DinoKnot_multimodel_SOURCE})

# link to simfold
target_link_libraries(DinoKnot LINK_PUBLIC simfold)

# link to pthread
target_link_libraries(DinoKnot LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT})
