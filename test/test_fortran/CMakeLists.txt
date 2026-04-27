# test diaglib directly from fortran
add_executable(test_diaglib
        test_fortran.f90
        utility.f90
        direct_matvecs.f90
)

# make a reference for comparison
add_executable(reference_diaglib
        reference.f90
        utility.f90
)


# Link DiagLib
target_link_libraries(test_diaglib PRIVATE diaglib)

# Linking BLAS to reference
if(MKL_FOUND)
        target_link_libraries(reference_diaglib PUBLIC MKL::MKL)
else()
        target_link_libraries(reference_diaglib PUBLIC BLAS::BLAS)
        target_link_libraries(reference_diaglib PUBLIC LAPACK::LAPACK)
endif()


# Ensure proper RPATH is set for installation
set_target_properties(test_diaglib PROPERTIES
        #BUILD_WITH_INSTALL_RPATH ON    # Use install RPATH during build
        INSTALL_RPATH_USE_LINK_PATH OFF # Use link paths
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"  # Point to install lib dir
)

# Install the executables
install(TARGETS test_diaglib reference_diaglib
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Run test
add_test(NAME reference_diaglib COMMAND reference_diaglib)
add_test(NAME test_diaglib COMMAND test_diaglib)