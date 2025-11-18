cmake_minimum_required(VERSION 3.15...4.0)

project(diaglib VERSION 1.0
        DESCRIPTION "A collection of iterative, matrix-free eigensolvers"
        LANGUAGES Fortran)

# Put all the .mod files in the build/include directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include")
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# By now find a whatever system BLAS
find_package(BLAS REQUIRED)

# Set explicitly the build type at runtime
if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release")
endif()

#message(STATUS "SONO QUI")

# Add appropriate subdirectories
add_subdirectory(src/fortran)
add_subdirectory(test/test_fortran)

