# make a reference for comparison
add_executable(dgl_reference
        ${PROJECT_SOURCE_DIR}/src/fortran/dgl_interface.f90
        ${PROJECT_SOURCE_DIR}/src/fortran/dgl_lapack.f90
        utility.f90
        reference.f90
)
# test diaglib directly from fortran
add_executable(dgl_test
        utility.f90
        test_fortran.f90
        direct_matvecs.f90
)

target_compile_options(dgl_reference PRIVATE ${DGL_Fortran_F2008} ${DGL_Fortran_OPTIONS})
target_compile_options(dgl_test PRIVATE ${DGL_Fortran_F2008} ${DGL_Fortran_OPTIONS})

# Keep the .mod files of the two executables separate, as they share utility.f90
set_target_properties(dgl_reference PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_reference)
set_target_properties(dgl_test PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_test)

# Link DiagLib
target_link_libraries(dgl_test PRIVATE diaglib)

# Linking BLAS and LAPACK
target_link_libraries(dgl_reference PRIVATE BLAS::BLAS LAPACK::LAPACK)

# Run test (the reference has to be computed first)
add_test(NAME dgl_reference COMMAND dgl_reference)
add_test(NAME dgl_test COMMAND dgl_test)
set_tests_properties(dgl_reference dgl_test PROPERTIES ENVIRONMENT "${DGL_TEST_ENVIRONMENT}")
set_tests_properties(dgl_test PROPERTIES DEPENDS dgl_reference)
