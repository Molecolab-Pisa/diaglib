# The python test needs numpy
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy"
        RESULT_VARIABLE DGL_NUMPY_MISSING OUTPUT_QUIET ERROR_QUIET)
if(DGL_NUMPY_MISSING)
        message(WARNING "numpy is not available for ${Python3_EXECUTABLE}: the python test will not be run")
        return()
endif()

add_test(NAME dgl_test_python COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_python.py
        $<TARGET_FILE:diaglib_c>
)
set_tests_properties(dgl_test_python PROPERTIES
        ENVIRONMENT "PYTHONPATH=${PROJECT_SOURCE_DIR}/src/python_interface/src;${DGL_TEST_ENVIRONMENT}")
