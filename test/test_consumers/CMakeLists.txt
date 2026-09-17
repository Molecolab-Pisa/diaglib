# Tests of an installed DiagLib, from the point of view of a user: this is a separate CMake
# project, which is configured against an installation (find_package(diaglib)) and not part of
# the DiagLib build. It checks that the installed package provides what the four interfaces
# need, and exercises the drivers on problems with known eigenvalues, nested calls, concurrent
# calls from several threads and the error codes.
#
#   cmake -S test/test_consumers -B <build> -DCMAKE_PREFIX_PATH=<prefix>
#   cmake --build <build>
#   ctest --test-dir <build>
#
cmake_minimum_required(VERSION 3.15)
project(dgl_consumers LANGUAGES Fortran C)

include(CTest)
find_package(diaglib REQUIRED)
message(STATUS "DiagLib ${diaglib_VERSION}, integers of ${diaglib_INTEGER_KIND} bytes")

# Fortran: uses the drivers directly, and OpenMP for the concurrent calls
find_package(OpenMP REQUIRED COMPONENTS Fortran)
add_executable(consumer_fortran consumer_fortran.f90)
target_link_libraries(consumer_fortran PRIVATE diaglib::diaglib OpenMP::OpenMP_Fortran)
add_test(NAME consumer_fortran COMMAND consumer_fortran)

# C: uses the C interface, with POSIX threads for the concurrent calls
find_package(Threads REQUIRED)
add_executable(consumer_c consumer_c.c)
target_link_libraries(consumer_c PRIVATE diaglib::diaglib_c Threads::Threads)
find_library(DGL_MATH_LIBRARY m)
if(DGL_MATH_LIBRARY)
        target_link_libraries(consumer_c PRIVATE ${DGL_MATH_LIBRARY})
endif()
add_test(NAME consumer_c COMMAND consumer_c)

# C++: includes the same header and uses std::thread
include(CheckLanguage)
check_language(CXX)
if(CMAKE_CXX_COMPILER)
        enable_language(CXX)
        add_executable(consumer_cxx consumer_cxx.cpp)
        set_target_properties(consumer_cxx PROPERTIES CXX_STANDARD 17)
        target_link_libraries(consumer_cxx PRIVATE diaglib::diaglib_c Threads::Threads)
        add_test(NAME consumer_cxx COMMAND consumer_cxx)
else()
        message(STATUS "No C++ compiler found: the C++ consumer will not be built")
endif()

# python: needs the installed pyDiaglib module (in the prefix, or in PYTHONPATH) and numpy
find_package(Python3 3.9 COMPONENTS Interpreter)
if(Python3_FOUND)
        execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy"
                RESULT_VARIABLE DGL_NUMPY_MISSING OUTPUT_QUIET ERROR_QUIET)
        file(GLOB DGL_PYTHON_DIRS "${diaglib_DIR}/../../python*/site-packages")
        if(DGL_NUMPY_MISSING)
                message(WARNING "numpy is not available: the python consumer will not be run")
        elseif(NOT DGL_PYTHON_DIRS)
                message(STATUS "pyDiaglib is not installed in this prefix: the python consumer will not be run")
        else()
                add_test(NAME consumer_python COMMAND ${Python3_EXECUTABLE}
                        ${CMAKE_CURRENT_SOURCE_DIR}/consumer_python.py $<TARGET_FILE:diaglib::diaglib_c>)
                set_tests_properties(consumer_python PROPERTIES ENVIRONMENT "PYTHONPATH=${DGL_PYTHON_DIRS}")
        endif()
endif()

# the test problems are small: threaded BLAS/LAPACK libraries can be very slow on them
set_tests_properties(consumer_fortran consumer_c PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
