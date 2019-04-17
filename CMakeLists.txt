# Specify the minimum version for CMake
project(dfire_rna)
cmake_minimum_required(VERSION 2.8) 

set(CMAKE_BINARY_DIR ${CMAKE_HOME_DIRECTORY}/bin)
set(CMAKE_SOURCE_DIR ${CMAKE_HOME_DIRECTORY}/src)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(SRCPDB 
    ${CMAKE_SOURCE_DIR}/PDB.cc  
    ${CMAKE_SOURCE_DIR}/pdb_utils.cc 
    )

set(SRCDFIRE 
    ${CMAKE_SOURCE_DIR}/dfire_calculator.cc
    ${CMAKE_SOURCE_DIR}/dfire_PDB.cc
    )

include_directories(${CMAKE_HOME_DIRECTORY}/include)


# executables

add_executable(DFIRE_RNA  ${SRCPDB} ${SRCDFIRE} ${CMAKE_SOURCE_DIR}/main.cc )
