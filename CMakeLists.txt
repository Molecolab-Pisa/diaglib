cmake_minimum_required(VERSION 3.15...4.1)

project(diaglib VERSION 1.0
        DESCRIPTION "A collection of iterative, matrix-free eigensolvers"
        LANGUAGES Fortran C)

#Include CPack to install the package
include(CPack)
set(CPACK_SOURCE_GENERATOR "SOURCE_TGZ")
set(CPACK_SOURCE_IGNORE_FILES /*.git* /*ipynb* /*build* /*vscode*)

#Include this to get the standard variables for installation paths
include(GNUInstallDirs)

#Include this to automatically make the configuration file
include(CMakePackageConfigHelpers)

#Include Ctest for testing
include(CTest)

#Include project specific macros
include(${CMAKE_SOURCE_DIR}/cmake_macros.cmake)


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
set(BLA_SIZEOF_INTEGER ${DGL_INTEGER_KIND})
set(BLA_VENDOR CACHE STRING "Vendor of BLAS/LAPACK found")

# Set required size of integers for BLAS, LAPACK and internal integers
if (DGL_INTEGER_KIND EQUAL 8)
        add_compile_definitions(DGL_INT_KIND_8)
        message(STATUS "Using 64-bit integers")
        #set(CMAKE_Fortran_FLAGS "-fdefault-integer-8 ${CMAKE_Fortran_FLAGS}")
        add_compile_options(-fdefault-integer-8)
else()
        add_compile_definitions(DGL_INT_KIND_4)
        message(STATUS "Using 32-bit integers")
endif()

# Find BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Find OpenMP
find_package(OpenMP REQUIRED)


set(BUILD_PYTHON OFF CACHE STRING "Build python interface")

# Add appropriate subdirectories
add_subdirectory(src/fortran)
add_subdirectory(src/c_interface)
if (BUILD_TESTING)
        message(STATUS "Tests will be built")
        add_subdirectory(test/test_fortran)
        add_subdirectory(test/test_c)
endif()
if (BUILD_PYTHON)
        find_package(Python3 3.9 REQUIRED)
        create_venv()
        pipinstall(numpy)
        add_subdirectory(src/python_interface)
        if (BUILD_TESTING)
                add_subdirectory(test/test_python)
        endif()
endif()

# uninstall target
if(NOT TARGET uninstall)
        configure_file(
                "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY)

        add_custom_target(uninstall
                COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
                )
endif()
