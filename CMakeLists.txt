# Ian Wark July 6 2017

cmake_minimum_required(VERSION 3.1)
project(DinoKnot_multimodel)

# enable c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -Wno-deprecated-declarations")

# create simfold library
add_subdirectory(simfold)

set(DinoKnot_multimodel_SOURCE
    src/W_final.cpp 
    src/pseudo_loop.cpp 
    src/V_final.cpp 
    src/VM_final.cpp
    src/h_common.cpp 
    src/hfold_validation.cpp
    src/DinoKnot.cpp
    src/Result.cpp
    src/cmdline.cpp
    src/param_path.cpp
)

# pass SIMFOLD_HOME to binary
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# get pthread library
find_package (Threads)

# set path for executable to be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# create executables
add_executable(DinoKnot ${DinoKnot_multimodel_SOURCE})

# Install DinoKnot into the bin directory (Needed for conda)
install(TARGETS DinoKnot RUNTIME DESTINATION bin)

# link to simfold
target_link_libraries(DinoKnot LINK_PUBLIC simfold)

# link to pthread
target_link_libraries(DinoKnot LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT})
