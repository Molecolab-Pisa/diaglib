cmake_minimum_required(VERSION 3.15...4.1)

project(diaglib VERSION 2.0
        DESCRIPTION "A collection of iterative, matrix-free eigensolvers"
        LANGUAGES Fortran C)

# Is DiagLib the main project, or is it included in another one
# (add_subdirectory, FetchContent)? Global settings are only changed in the first case.
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
        set(DGL_TOP_LEVEL ON)
else()
        set(DGL_TOP_LEVEL OFF)
endif()

#Include this to get the standard variables for installation paths
include(GNUInstallDirs)

#Include this to automatically make the configuration file
include(CMakePackageConfigHelpers)

#Include Ctest for testing
if(DGL_TOP_LEVEL)
        include(CTest)
endif()

# default build type
if(DGL_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# The Fortran sources need the preprocessor (integer kinds)
set(CMAKE_Fortran_PREPROCESS ON)

# Compiler options for the DiagLib targets. They are added to the targets, so that the
# flags chosen by the user (CMAKE_<LANG>_FLAGS, CMAKE_<LANG>_FLAGS_<CONFIG>) are preserved.
# DGL_Fortran_F2008/DGL_Fortran_F2018 select the language standard.
option(DGL_NATIVE_ARCH "Optimize for the build machine (the binaries will not be portable)" OFF)

set(DGL_Fortran_OPTIONS "")
set(DGL_C_OPTIONS "")
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(DGL_Fortran_F2008 -std=f2008)
    set(DGL_Fortran_F2018 -std=f2018)
    list(APPEND DGL_Fortran_OPTIONS -pedantic "$<$<CONFIG:Debug>:-Wall;-Wextra;-fcheck=all;-fbacktrace>")
    if(DGL_NATIVE_ARCH)
        list(APPEND DGL_Fortran_OPTIONS -march=native)
    endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(DGL_Fortran_F2008 "SHELL:-stand f08")
    set(DGL_Fortran_F2018 "SHELL:-stand f18")
    list(APPEND DGL_Fortran_OPTIONS "$<$<CONFIG:Debug>:-warn;all;-check;all;-traceback>")
    if(DGL_NATIVE_ARCH)
        list(APPEND DGL_Fortran_OPTIONS -xHost)
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU")
    list(APPEND DGL_C_OPTIONS -pedantic "$<$<CONFIG:Debug>:-Wall;-Wextra>")
    if(DGL_NATIVE_ARCH)
        list(APPEND DGL_C_OPTIONS -march=native)
    endif()
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
    list(APPEND DGL_C_OPTIONS "$<$<CONFIG:Debug>:-traceback>")
    if(DGL_NATIVE_ARCH)
        list(APPEND DGL_C_OPTIONS -xHost)
    endif()
endif()

# Extra warnings for the DiagLib libraries, treated as errors (used by the CI pipeline).
# They are only added to the library targets: the test programs call BLAS/LAPACK routines
# through implicit interfaces on purpose, to check them as a user would.
option(DGL_STRICT_WARNINGS "Compile the DiagLib libraries with extra warnings, treated as errors" OFF)

set(DGL_Fortran_STRICT "")
if(DGL_STRICT_WARNINGS)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        list(APPEND DGL_Fortran_STRICT -Wall -Wextra -Wconversion -Wimplicit-interface -Werror)
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        list(APPEND DGL_Fortran_STRICT "SHELL:-warn all" "SHELL:-warn errors")
    endif()
endif()

# Integer default kind is 4. If 8 is desired one should run `cmake -DDGL_INTEGER_KIND=8 ..`
set(DGL_INTEGER_KIND 4 CACHE STRING "Kind (size in bytes) of the integers used in DiagLib: 4 or 8")
set_property(CACHE DGL_INTEGER_KIND PROPERTY STRINGS 4 8)
if(NOT DGL_INTEGER_KIND MATCHES "^(4|8)$")
        message(FATAL_ERROR "DGL_INTEGER_KIND must be 4 or 8, not '${DGL_INTEGER_KIND}'")
endif()
set(BLA_SIZEOF_INTEGER ${DGL_INTEGER_KIND})

# To pic a specific vendor for Blas this varible should be input.
# Common options are "Intel10_64lp","Intel10_64ilp","OpenBLAS"
set(BLA_VENDOR CACHE STRING "Vendor of BLAS/LAPACK found")

# FindBLAS only takes BLA_SIZEOF_INTEGER into account since CMake 3.22: with older versions,
# a library with 32-bit integers could be selected for a build with 64-bit integers,
# which would silently give wrong results. In that case, the vendor has to be given explicitly.
if(DGL_INTEGER_KIND EQUAL 8 AND CMAKE_VERSION VERSION_LESS 3.22 AND NOT BLA_VENDOR)
        message(FATAL_ERROR "With CMake older than 3.22, a BLAS/LAPACK library with 64-bit integers "
                            "cannot be selected automatically. Set BLA_VENDOR to such a library "
                            "(e.g., -DBLA_VENDOR=Intel10_64ilp), or use CMake 3.22 or newer.")
endif()

# Set required size of integers for BLAS, LAPACK and internal integers.
# The default integer kind is not changed: all the integers passed to BLAS/LAPACK are
# integer(ip), which is checked by the compiler through the interfaces in dgl_lapack.f90.
# The C sources do not need the definitions: they use the generated diaglib_config.h
if (DGL_INTEGER_KIND EQUAL 8)
        add_compile_definitions($<$<COMPILE_LANGUAGE:Fortran>:DGL_INT_KIND_8>)
        message(STATUS "Using 64-bit integers")
else()
        add_compile_definitions($<$<COMPILE_LANGUAGE:Fortran>:DGL_INT_KIND_4>)
        message(STATUS "Using 32-bit integers")
endif()

# Find BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Find OpenMP (optional, only used for wall-clock timings)
find_package(OpenMP COMPONENTS Fortran)

# Build python interface?
option(BUILD_PYTHON "Install the python interface" OFF)

# Build the stress tests? They take a few minutes, and are meant for the CI pipeline and for
# checking the solvers after a change, rather than for every build.
option(DGL_STRESS_TESTS "Build the stress tests (ctest -L stress)" OFF)

# Environment for the tests: the test problems are small, and threaded BLAS/LAPACK
# libraries can be extremely slow on them (e.g., OpenBLAS with many threads)
set(DGL_TEST_ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")

# Add appropriate subdirectories
add_subdirectory(src/fortran)
add_subdirectory(src/c_interface)
if (DGL_TOP_LEVEL AND BUILD_TESTING)
        message(STATUS "Tests will be built")
        add_subdirectory(test/test_fortran)
        add_subdirectory(test/test_c)
        if (DGL_STRESS_TESTS)
                message(STATUS "Stress tests will be built")
                add_subdirectory(test/test_stress)
        endif()
endif()
if (BUILD_PYTHON)
        message(STATUS "Python interface will be installed")
        find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
        add_subdirectory(src/python_interface)
        if (DGL_TOP_LEVEL AND BUILD_TESTING)
                add_subdirectory(test/test_python)
        endif()
endif()

# Export the DiagLib targets and install the files needed by find_package(diaglib)
install(EXPORT diaglibTargets
        FILE DiaglibTargets.cmake
        NAMESPACE diaglib::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Diaglib
)

configure_package_config_file(${PROJECT_SOURCE_DIR}/src/fortran/Config.cmake.in
        "${PROJECT_BINARY_DIR}/diaglibConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Diaglib
)

write_basic_package_version_file("${PROJECT_BINARY_DIR}/diaglibConfigVersion.cmake"
        COMPATIBILITY SameMajorVersion
)

install(FILES
        "${PROJECT_BINARY_DIR}/diaglibConfig.cmake"
        "${PROJECT_BINARY_DIR}/diaglibConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Diaglib
)

# Costum target to uninstall the project form the install destination
if(NOT TARGET uninstall)
        configure_file(
                "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY)

        add_custom_target(uninstall
                COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
                )
endif()

# Include CPack to create the source package (`make package_source`).
# The CPACK_* variables must be set before including CPack.
if(DGL_TOP_LEVEL)
        set(CPACK_SOURCE_GENERATOR "TGZ")
        set(CPACK_SOURCE_IGNORE_FILES
                "/\\\\.git/"
                "/\\\\.gitignore$"
                "/\\\\.vscode/"
                "/\\\\.venv/"
                "/build[^/]*/"
                "\\\\.ipynb"
                "/__pycache__/"
        )
        include(CPack)
endif()
