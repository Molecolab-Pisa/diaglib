add_library(diaglib SHARED
        dgl_global_utils.f90
        dgl_orthogonalizations.f90
        dgl_minor_utils.f90
        dgl_external_interfaces.f90
        dgl_drivers_interfaces.f90
        dgl_init.f90
        davidson_driver.f90
        lobpcg_driver.f90
        )

# Linking proper BLAS
if(BLA_VENDOR_l MATCHES "intel")
        target_link_libraries(diaglib PUBLIC ${MKL_LIBRARIES})
        target_include_directories(diaglib PUBLIC ${MKL_INCLUDE_DIRS})
else()
        target_link_libraries(diaglib PUBLIC ${BLAS_LIBRARIES})
        target_include_directories(diaglib PUBLIC ${BLAS_INCLUDE_DIRS})
endif()

# Set Fortran flags depending on build type
STRING(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_l)
if(BUILD_TYPE_l MATCHES "release" OR NOT CMAKE_BUILD_TYPE)

        target_compile_options(diaglib PRIVATE
                -march=native)

elseif(BUILD_TYPE_l MATCHES "debug")

        if(${CMAKE_Fortran_COMPILER} MATCHES "gfortran" OR
           ${CMAKE_Fortran_COMPILER} MATCHES "f95")

                target_compile_options(diaglib PRIVATE
                        -fbacktrace -fcheck=all -Wconversion -Wuninitialized)

        elseif(${CMAKE_Fortran_COMPILER} MATCHES "ifort" OR
               ${CMAKE_Fortran_COMPILER} MATCHES "ifx")

                target_compile_options(diaglib PRIVATE
                        -traceback -check all -warn uninitialized)
        endif()
endif()

#Set fortran kind for integers and floats
STRING(TOLOWER "${DGL_REAL_KIND}" REAL_KIND)
STRING(TOLOWER "${DGL_INTEGER_KIND}" INTEGER_KIND)

if(REAL_KIND MATCHES single)
        target_compile_definitions(diaglib PRIVATE REAL_KIND=SINGLE)
elseif(REAL_KIND MATCHES double OR NOT DGL_REAL_KIND)
        target_compile_definitions(diaglib PRIVATE REAL_KIND=DOUBLE)
else()
        message(FATAL_ERROR "Unrecognized kind requested for floating point numbers")
endif()

if(INTEGER_KIND MATCHES single)
        target_compile_definitions(diaglib PRIVATE INTEGER_KIND=SINGLE)
elseif(INTEGER_KIND MATCHES double OR NOT DGL_INTEGER_KIND)
        target_compile_definitions(diaglib PRIVATE INTEGER_KIND=DOUBLE)
else()
        message(FATAL_ERROR "Unrecognized kind requested for integer numbers")
endif()

# Linking proper OpenMP
target_link_libraries(diaglib PRIVATE OpenMP::OpenMP_Fortran)

get_target_property(_ diaglib COMPILE_DEFINITIONS)
message(STATUS "Diaglib compile definitions are: ${_}")

get_target_property(_ diaglib COMPILE_OPTIONS)
message(STATUS "Diaglib compile options are: ${_}")
