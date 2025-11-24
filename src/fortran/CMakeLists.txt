add_library(diaglib SHARED
        dgl_global_utils.f90
        dgl_orthogonalizations.f90
        dgl_minor_utils.f90
        dgl_external_interfaces.f90
        dgl_drivers_interfaces.f90
        davidson_driver.f90
        lobpcg_driver.f90
        #gen_davidson_driver.f90
        )

target_link_libraries(diaglib PUBLIC ${BLAS_LIBRARIES})
target_include_directories(diaglib PUBLIC ${BLAS_INCLUDE_DIRS})

# Set Fortran flags depending on build type
STRING(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_l)
if(CMAKE_BUILD_TYPE_l MATCHES "release")

        target_compile_options(diaglib PUBLIC -march=native)

elseif(CMAKE_BUILD_TYPE_l MATCHES "debug")

        if(${CMAKE_Fortran_COMPILER} MATCHES "gfortran" OR
           ${CMAKE_Fortran_COMPILER} MATCHES "f95")

                target_compile_options(diaglib PUBLIC -fbacktrace -fcheck=all)

        elseif(${CMAKE_Fortran_COMPILER} MATCHES "ifort" OR
               ${CMAKE_Fortran_COMPILER} MATCHES "ifx")

       target_compile_options(diaglib PUBLIC -traceback -check all)

        endif()
endif()

get_target_property(_ diaglib COMPILE_OPTIONS)
message(STATUS "Diaglib compile options are: ${_}")
