add_library(diaglib SHARED 
        diaglib_global_utils.f90 
        diaglib_minor_utils.f90 
        orthogonalizations.f90 
        lobpcg_driver.f90 
        davidson_driver.f90)

target_link_libraries(diaglib PUBLIC ${BLAS_LIBRARIES})
target_include_directories(diaglib PUBLIC ${BLAS_INCLUDE_DIRS})
