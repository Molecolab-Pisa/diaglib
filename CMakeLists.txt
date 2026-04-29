cmake_minimum_required(VERSION 3.15...4.0)

project(diaglib VERSION 1.0
        DESCRIPTION "A collection of iterative, matrix-free eigensolvers"
        LANGUAGES Fortran C)

#Include CPack to install the package
include(CPack)
#Include this to get the standard variables for installation paths
include(GNUInstallDirs)
#Include this to automatically make the configuration file
include(CMakePackageConfigHelpers)
#Include Ctest for testing
include(CTest)

# default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# set Fortran flags for build types
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-std=f2003 -pedantic")
    set(CMAKE_Fortran_FLAGS_DEBUG "-ggdb -Wall -Wextra -fcheck=all -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "-stand 2003")
    set(CMAKE_Fortran_FLAGS_DEBUG "-ggdb -warn all -check all -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost")
endif()
set(CMAKE_Fortran_PREPROCESS ON)

if(CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS "-pedantic")
    set(CMAKE_C_FLAGS_DEBUG "-ggdb -Wall -Wextra")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
    #set(CMAKE_C_FLAGS "-stand 2003")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -ggdb -traceback")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -xHost")
endif()

# If wanted force double precision integers
set(DGL_INTEGER_KIND 4 CACHE STRING "Kind of integers used in DiagLib")

# Set required size of integers for BLAS, MKL and internal integers
set(BLA_SIZEOF_INTEGER ${DGL_INTEGER_KIND})
if (DGL_INTEGER_KIND EQUAL 8)
        set(MKL_INTERFACE_FULL gf_ilp64)
        add_compile_definitions(DGL_INT_KIND_8)
        message(STATUS "Using 64-bit integers")
else()
        set(MKL_INTERFACE_FULL gf_lp64)
        add_compile_definitions(DGL_INT_KIND_4)
        message(STATUS "Using 32-bit integers")
endif()

# Set required parameters for BLAS
set(BLA_THREAD CACHE STRING "Choice between SEQ or OMP BLAS")
set(BLA_VENDOR CACHE STRING "Choice of BLAS vendor")


if(NOT BLA_VENDOR)
        find_package(MKL QUIET)
        if (NOT MKL_FOUND)
                find_package(BLAS REQUIRED)
                find_package(LAPACK REQUIRED)
        else()
                message(STATUS "MKL_VERSION: ${MKL_VERSION}")
                message(STATUS "MKL_ROOT: ${MKL_ROOT}")
                message(STATUS "MKL_ARCH: ${MKL_ARCH}")
                message(STATUS "MKL_LINK: ${MKL_LINK}")
                message(STATUS "MKL_INTERFACE_FULL: ${MKL_INTERFACE_FULL}")
                message(STATUS "MKL_THREADING: ${MKL_THREADING}")
                message(STATUS "MKL_MPI: ${MKL_MPI}")
        endif()
elseif(BLA_VENDOR MATCHES "Intel")
        find_package(MKL REQUIRED)
else()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)


# Add appropriate subdirectories
add_subdirectory(src/fortran)
add_subdirectory(src/fortran_include)
add_subdirectory(src/c_interface)
add_subdirectory(test/test_fortran)
add_subdirectory(test/test_c)
add_subdirectory(test/test_python)

# uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
