# Stress tests: every driver is run many times from random guesses on problems with different
# spectra, metrics and numbers of roots, and the results are compared with dense LAPACK
# references. They take a few minutes, so they are only built if DGL_STRESS_TESTS is on.
add_executable(dgl_stress_sym stress_sym.f90)
add_executable(dgl_stress_nosym stress_nosym.f90)

foreach(target dgl_stress_sym dgl_stress_nosym)
        target_compile_options(${target} PRIVATE ${DGL_Fortran_F2008} ${DGL_Fortran_OPTIONS})
        set_target_properties(${target} PROPERTIES Fortran_MODULE_DIRECTORY
                ${CMAKE_CURRENT_BINARY_DIR}/mod_${target})
        # the references are computed with the BLAS/LAPACK library DiagLib is linked to
        target_link_libraries(${target} PRIVATE diaglib BLAS::BLAS LAPACK::LAPACK)
        add_test(NAME ${target} COMMAND ${target})
        set_tests_properties(${target} PROPERTIES ENVIRONMENT "${DGL_TEST_ENVIRONMENT}" LABELS stress)
endforeach()
