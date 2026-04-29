get_target_property(C_DIAGLIB_BINDIR diaglib_c BINARY_DIR)
add_test(NAME dgl_test_python COMMAND ${CMAKE_SOURCE_DIR}/.venv/bin/python
        ${CMAKE_CURRENT_SOURCE_DIR}/test_python.py
        ${C_DIAGLIB_BINDIR}/libdiaglib_c.so
        ${DGL_INTEGER_KIND}
)