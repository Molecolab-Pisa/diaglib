# Extreme, long-running stress tests: randomized, time-boxed fuzzers for the symmetric/
# generalized drivers and for the non-symmetric/SMO-GD drivers (problems up to n=4000 and metric
# condition numbers up to 1e11, a fresh random configuration every trial), plus a sustained
# heavy-concurrency test that hammers the drivers from many OpenMP threads at once.
#
# These are meant to be run by hand before a release or after a change to the drivers, not as
# part of routine testing: they are NOT registered with ctest and NOT built by DGL_STRESS_TESTS,
# only by the separate DGL_EXTREME_STRESS_TESTS option (off by default, never turned on by CI).
# Each program takes a seed and a time budget (in minutes) on the command line, e.g.:
#
#   ./dgl_stress_extreme_sym 12345 30
#   ./dgl_thread_stress 20 64        # minutes, number of threads
#
# and reports genuine bugs (wrong result claimed converged, NaN, unexpected error code) separately
# from plain non-convergence within max_iter, which is expected for some of the harder random
# configurations and not itself a bug.
add_executable(dgl_stress_extreme_sym stress_extreme_sym.f90)
add_executable(dgl_stress_extreme_nosym stress_extreme_nosym.f90)

foreach(target dgl_stress_extreme_sym dgl_stress_extreme_nosym)
        target_compile_options(${target} PRIVATE ${DGL_Fortran_F2008} ${DGL_Fortran_OPTIONS})
        set_target_properties(${target} PROPERTIES Fortran_MODULE_DIRECTORY
                ${CMAKE_CURRENT_BINARY_DIR}/mod_${target})
        # the references are computed with the BLAS/LAPACK library DiagLib is linked to
        target_link_libraries(${target} PRIVATE diaglib BLAS::BLAS LAPACK::LAPACK)
endforeach()

add_executable(dgl_thread_stress thread_stress.c)
target_compile_options(dgl_thread_stress PRIVATE ${DGL_C_OPTIONS})
target_link_libraries(dgl_thread_stress PRIVATE diaglib_c)

find_library(DGL_MATH_LIBRARY m)
if(DGL_MATH_LIBRARY)
        target_link_libraries(dgl_thread_stress PRIVATE ${DGL_MATH_LIBRARY})
endif()

if(OpenMP_Fortran_FOUND)
        find_package(OpenMP COMPONENTS C)
        if(OpenMP_C_FOUND)
                target_link_libraries(dgl_thread_stress PRIVATE OpenMP::OpenMP_C)
        endif()
endif()
