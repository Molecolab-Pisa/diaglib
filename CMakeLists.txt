cmake_minimum_required(VERSION 3.15...4.0)

project(diaglib VERSION 1.0
        DESCRIPTION "A collection of iterative, matrix-free eigensolvers"
        LANGUAGES Fortran)

# Put all the .mod files in the build/include directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include")
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Find required BLAS
STRING(TOLOWER "${BLA_VENDOR}" BLA_VENDOR_l)
if(BLA_VENDOR_l MATCHES "intel")
        find_package(MKL REQUIRED)
elseif(NOT BLA_VENDOR_l)
        find_package(BLAS)
        if (NOT BLAS_FOUND)
                find_package(MKL REQUIRED)
        endif()
else()
        find_package(BLAS REQUIRED)
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)

# Add appropriate subdirectories
add_subdirectory(src/fortran)
add_subdirectory(test/test_fortran)

