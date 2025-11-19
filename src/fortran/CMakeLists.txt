add_library(diaglib SHARED 
        diaglib_global_utils.f90 
        diaglib_minor_utils.f90 
        orthogonalizations.f90 
        lobpcg_driver.f90 
        davidson_driver.f90
        gen_davidson_driver.f90
        )

target_link_libraries(diaglib PUBLIC ${BLAS_LIBRARIES})
target_include_directories(diaglib PUBLIC ${BLAS_INCLUDE_DIRS})

# Set Fortran flags depending on build type
if(CMAKE_BUILD_TYPE STREQUAL "Release")

        target_compile_options(diaglib PUBLIC -march=native)

elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")

        if(${CMAKE_Fortran_COMPILER} MATCHES "usr/bin/gfortran")

                target_compile_options(diaglib PUBLIC -fbacktrace -fcheck=all)

        elseif(${CMAKE_Fortran_COMPILER} MATCHES ".*ifort.*" OR
               ${CMAKE_Fortran_COMPILER} MATCHES ".*ifx.*")

       target_compile_options(diaglib PUBLIC -traceback -check all)

        endif()
endif()

get_target_property(_ diaglib COMPILE_OPTIONS)
message(STATUS "Diaglib compile options are: ${_}")
