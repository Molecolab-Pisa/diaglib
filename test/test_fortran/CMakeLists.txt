# make a reference for comparison
add_executable(dgl_reference
        utility.f90
        reference.f90
)
# test diaglib directly from fortran
add_executable(dgl_test
        utility.f90
        test_fortran.f90
        direct_matvecs.f90
)


# Link DiagLib
target_link_libraries(dgl_test PRIVATE diaglib)

# Linking BLAS and LAPACK
target_link_libraries(dgl_reference PUBLIC BLAS::BLAS)
target_link_libraries(dgl_reference PUBLIC LAPACK::LAPACK)


# Ensure proper RPATH is set for installation
set_target_properties(dgl_test PROPERTIES
        #BUILD_WITH_INSTALL_RPATH ON    # Use install RPATH during build
        INSTALL_RPATH_USE_LINK_PATH OFF # Use link paths
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"  # Point to install lib dir
)
set_target_properties(dgl_reference PROPERTIES
        #BUILD_WITH_INSTALL_RPATH ON    # Use install RPATH during build
        INSTALL_RPATH_USE_LINK_PATH OFF # Use link paths
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"  # Point to install lib dir
)

# Install the executables
install(TARGETS dgl_test dgl_reference
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Run test
add_test(NAME dgl_reference COMMAND dgl_reference)
add_test(NAME dgl_test COMMAND dgl_test)