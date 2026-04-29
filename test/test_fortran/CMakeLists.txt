# make a reference for comparison
add_executable(dgl_reference
        reference.f90
        utility.f90
)
# test diaglib directly from fortran
add_executable(dgl_test
        test_fortran.f90
        utility.f90
        direct_matvecs.f90
)


# Link DiagLib
target_link_libraries(dgl_test PRIVATE diaglib)

# Linking BLAS to reference
if(MKL_FOUND)
        target_link_libraries(dgl_reference PUBLIC MKL::MKL)
else()
        target_link_libraries(dgl_reference PUBLIC BLAS::BLAS)
        target_link_libraries(dgl_reference PUBLIC LAPACK::LAPACK)
endif()


# Ensure proper RPATH is set for installation
set_target_properties(dgl_test PROPERTIES
        #BUILD_WITH_INSTALL_RPATH ON    # Use install RPATH during build
        INSTALL_RPATH_USE_LINK_PATH OFF # Use link paths
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"  # Point to install lib dir
)

# Install the executables
install(TARGETS dgl_test dgl_reference
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Run test
#add_test(NAME dgl_reference COMMAND dgl_reference)
add_test(NAME dgl_test COMMAND dgl_test)