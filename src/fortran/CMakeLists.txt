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
if(MKL_FOUND)
        target_link_libraries(diaglib PUBLIC MKL::MKL)
else()
        target_link_libraries(diaglib PUBLIC BLAS::BLAS)
endif()


# Set Fortran flags depending on build type
STRING (TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_l)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
   
        target_compile_options(diaglib PRIVATE -std=f2003 -pedantic)

        if(BUILD_TYPE_l MATCHES "debug")
                target_compile_options(diaglib PRIVATE
                        -fbacktrace -fcheck=all -Wall -Wextra)
        endif()

elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")

        target_compile_options(diaglib PRIVATE -stand f03)

        if(BUILD_TYPE_l MATCHES "debug")
                target_compile_options(diaglib PRIVATE
                        -traceback -check all -warn all)
        endif()

endif()

get_target_property(_ diaglib COMPILE_OPTIONS)
message(STATUS "Diaglib compile options are: ${_}")
        

# Linking proper OpenMP
target_link_libraries(diaglib PRIVATE OpenMP::OpenMP_Fortran)