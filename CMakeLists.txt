cmake_minimum_required(VERSION 3.15...4.0)

project(diaglib VERSION 1.0
        DESCRIPTION "A collection of iterative, matrix-free eigensolvers"
        LANGUAGES Fortran)

find_package(BLAS REQUIRED)

add_subdirectory(src/fortran)
add_subdirectory(test/test_fortran)
