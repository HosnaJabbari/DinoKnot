# Ian Wark July 6 2017

cmake_minimum_required(VERSION 3.1)
project(DinoKnot_multimodel)

# enable c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# create simfold library
add_subdirectory(simfold)

set(DinoKnot_multimodel_SOURCE
    W_final.h W_final.cpp 
    pseudo_loop.h pseudo_loop.cpp 
    V_final.h V_final.cpp 
    VM_final.h VM_final.cpp
    h_common.h h_common.cpp 
    hfold_validation.h hfold_validation.cpp
    dinoknot.h DinoKnot.cpp
    Result.h Result.cpp
)

# pass SIMFOLD_HOME to binary
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# get pthread library
find_package (Threads)

# set path for executable to be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executables
add_executable(DinoKnot_multimodel ${DinoKnot_multimodel_SOURCE})

# link to simfold
target_link_libraries(DinoKnot_multimodel LINK_PUBLIC simfold)

# link to pthread
target_link_libraries(DinoKnot_multimodel LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT})
