# Install the pyDiaglib module. A relative path is relative to CMAKE_INSTALL_PREFIX.
# Alternatively, the module can be installed with `pip install src/python_interface`.
set(DGL_PYTHON_INSTALL_DIR
        "${CMAKE_INSTALL_LIBDIR}/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages"
        CACHE PATH "Installation directory of the pyDiaglib module")

install(FILES src/pyDiaglib.py
        DESTINATION ${DGL_PYTHON_INSTALL_DIR})
