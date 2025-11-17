add_library(diaglib SHARED diaglib.f90 1_diaglib_global_utils.f90)

target_link_libraries(diaglib PUBLIC ${BLAS_LIBRARIES})
target_include_directories(diaglib PUBLIC ${BLAS_INCLUDE_DIRS})
