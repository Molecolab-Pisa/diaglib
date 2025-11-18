cmake_minimum_required(VERSION 3.15...4.0)

project(diaglib VERSION 1.0
        DESCRIPTION "A collection of iterative, matrix-free eigensolvers"
        LANGUAGES Fortran)

# Put all the .mod files in the build/include directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include")
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# By now find a whatever system BLAS
find_package(BLAS REQUIRED)

# Set Fortran flags depending on build type
if(NOT CMAKE_BUILD_TYPE OR
   CMAKE_BUILD_TYPE STREQUAL "Release")

        set(CMAKE_Fortran_FLAGS "-O3 -march=native")

elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")

        if(${CMAKE_Fortran_COMPILER} MATCHES "usr/bin/gfortran")

                set (CMAKE_Fortran_FLAGS "-fbacktrace -fcheck=all")

        elseif(${CMAKE_Fortran_COMPILER} MATCHES ".*ifort.*" OR
               ${CMAKE_Fortran_COMPILER} MATCHES ".*ifx.*")

                set (CMAKE_Fortran_FLAGS "-traceback -check all")

        endif()
endif()

# Some recap
message(STATUS "Fortran compiler is set to: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Build type is set to: ${CMAKE_BUILD}")
message(STATUS "Fortran flags are set to: ${CMAKE_Fortran_FLAGS}")

#message(STATUS "SONO QUIH")

# Add appropriate subdirectories
add_subdirectory(src/fortran)
add_subdirectory(test/test_fortran)
