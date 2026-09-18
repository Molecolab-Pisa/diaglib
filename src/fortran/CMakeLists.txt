add_library(diaglib SHARED
        dgl_interface.f90
        dgl_lapack.f90
        dgl_global_utils.f90
        dgl_orthogonalizations.f90
        dgl_minor_utils.f90
        davidson_driver.f90
        davidson_nosym_driver.f90
        lobpcg_driver.f90
        smogd_driver.f90
)

target_compile_options(diaglib PRIVATE ${DGL_Fortran_F2008} ${DGL_Fortran_OPTIONS} ${DGL_Fortran_STRICT})

# Set where the .mod files will be placed, and let the installed library
# find its dependencies (libraries in the same directory, BLAS/LAPACK)
set_target_properties(diaglib PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}
        INSTALL_RPATH "$ORIGIN"
        INSTALL_RPATH_USE_LINK_PATH TRUE)

# Include the directory with the .mod files
target_include_directories(diaglib PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Linking BLAS and LAPACK. They are private dependencies of the shared library,
# so that users of DiagLib do not need to find them.
target_link_libraries(diaglib PRIVATE BLAS::BLAS LAPACK::LAPACK)

# Linking OpenMP, if available
if(OpenMP_Fortran_FOUND)
        target_link_libraries(diaglib PRIVATE OpenMP::OpenMP_Fortran)
endif()

# Install libray
install(TARGETS diaglib
        EXPORT diaglibTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install the .mod file of the public module. The drivers are submodules of dgl_interface,
# so that it does not depend on the internal modules, which are not installed.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/dgl_interface.mod
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
