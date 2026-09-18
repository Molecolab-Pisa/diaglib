# test diaglib directly from C
add_executable(dgl_test_c
        test_c.c
)

target_compile_options(dgl_test_c PRIVATE ${DGL_C_OPTIONS})

# Link DiagLib (this also provides the include directories of the C headers)
target_link_libraries(dgl_test_c PRIVATE diaglib_c)

# the math library has to be linked explicitly on some systems
find_library(DGL_MATH_LIBRARY m)
if(DGL_MATH_LIBRARY)
        target_link_libraries(dgl_test_c PRIVATE ${DGL_MATH_LIBRARY})
endif()

# Test concurrent calls from different threads, if DiagLib has been built with OpenMP
if(OpenMP_Fortran_FOUND)
        find_package(OpenMP COMPONENTS C)
        if(OpenMP_C_FOUND)
                target_link_libraries(dgl_test_c PRIVATE OpenMP::OpenMP_C)
        endif()
endif()

# Run test
add_test(NAME dgl_test_c COMMAND dgl_test_c)
set_tests_properties(dgl_test_c PROPERTIES ENVIRONMENT "${DGL_TEST_ENVIRONMENT}")

# check that the C header can be used from C++, if a C++ compiler is available
include(CheckLanguage)
check_language(CXX)
if(CMAKE_CXX_COMPILER)
        enable_language(CXX)
        add_executable(dgl_test_cxx test_cxx.cpp)
        target_link_libraries(dgl_test_cxx PRIVATE diaglib_c)
        add_test(NAME dgl_test_cxx COMMAND dgl_test_cxx)
        set_tests_properties(dgl_test_cxx PROPERTIES ENVIRONMENT "${DGL_TEST_ENVIRONMENT}")
else()
        message(STATUS "No C++ compiler found: the C++ header test will not be built")
endif()
