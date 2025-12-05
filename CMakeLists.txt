cmake_minimum_required(VERSION 3.15...4.0)

project(diaglib VERSION 1.0
        DESCRIPTION "A collection of iterative, matrix-free eigensolvers"
        LANGUAGES Fortran)

# Put all the .mod files in the build/include directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include")
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Find required BLAS
set(BLA_SIZEOF_INTEGER 4) # DiagLib only uses kind 4 integers
set(BLA_THREAD CACHE STRING "Choice between SEQ or OMP BLAS")
set(BLA_VENDOR CACHE STRING "Choice of BLAS vendor")

if(NOT BLA_VENDOR)
find_package(MKL)
        if (NOT MKL_FOUND)
                find_package(BLAS REQUIRED)
        endif()
else()
        find_package(BLAS REQUIRED)
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)

# Add appropriate subdirectories
add_subdirectory(src/fortran)
add_subdirectory(test/test_fortran)

