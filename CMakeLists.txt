# Ian Wark July 6 2017

cmake_minimum_required(VERSION 3.1)
project(HFold_interacting_multimodel)

# enable c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# only continue if we know we already have boost
# create simfold library
add_subdirectory(simfold)

set(HFold_interacting_multimodel_SOURCE
    W_final.h W_final.cpp 
    pseudo_loop.h pseudo_loop.cpp 
    V_final.h V_final.cpp 
    VM_final.h VM_final.cpp
    h_common.h h_common.cpp 
    hfold_validation.h hfold_validation.cpp
    hfold_interacting.h HFold_interacting.cpp
)

# get pthread library
find_package (Threads)

# set path for executable to be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executables
add_executable(HFold_interacting_multimodel ${HFold_interacting_multimodel_SOURCE})

# link to simfold
target_link_libraries(HFold_interacting_multimodel LINK_PUBLIC simfold)

# link to pthread
target_link_libraries(HFold_interacting_multimodel LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT})
