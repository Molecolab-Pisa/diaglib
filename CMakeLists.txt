cmake_minimum_required(VERSION 3.15...4.0)

project(diaglib VERSION 1.0
        DESCRIPTION "A collection of iterative, matrix-free eigensolvers"
        LANGUAGES Fortran)

#Include CPack to install the package
include(CPack)
#Include this to get the standard variables for installation paths
include(GNUInstallDirs)
#Include this to automatically make the configuration file
include(CMakePackageConfigHelpers)


# If wanted force double precision integers
set(DGL_INTEGER_KIND 4 CACHE STRING "Kind of integers used in DiagLib")
if (DGL_INTEGER_KIND EQUAL 8)
        if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
                add_compile_options(-fdefault-integer-8)
        elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
                add_compile_options(-i8)
        endif()
endif()

# Find required BLAS
set(BLA_SIZEOF_INTEGER ${DGL_INTEGER_KIND})
set(BLA_THREAD CACHE STRING "Choice between SEQ or OMP BLAS")
set(BLA_VENDOR CACHE STRING "Choice of BLAS vendor")

if(NOT BLA_VENDOR)
        find_package(MKL)
        if (NOT MKL_FOUND)
                find_package(BLAS REQUIRED)
        endif()
elseif(BLA_VENDOR MATCHES "Intel")
        find_package(MKL REQUIRED)
else()
        find_package(BLAS REQUIRED)
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)


# Add appropriate subdirectories
add_subdirectory(src/fortran)
add_subdirectory(test/test_fortran)

# uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
