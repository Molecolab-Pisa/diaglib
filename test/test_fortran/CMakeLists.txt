# test diaglib directly from fortran
add_executable(test_fortran
        test_fortran.f90
        utility.f90
        direct_matvecs.f90
)

# make a reference for comparison
add_executable(reference
        reference.f90
        utility.f90
)


# Link DiagLib
target_link_libraries(test_fortran PRIVATE diaglib)

# Linking BLAS to reference
if(MKL_FOUND)
        target_link_libraries(reference PUBLIC MKL::MKL)
else()
        target_link_libraries(reference PUBLIC BLAS::BLAS)
        target_link_libraries(reference PUBLIC LAPACK::LAPACK)
endif()


# Ensure proper RPATH is set for installation
set_target_properties(test_fortran PROPERTIES
#        BUILD_WITH_INSTALL_RPATH ON    # Use install RPATH during build
        INSTALL_RPATH_USE_LINK_PATH OFF # Use link paths
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib64"  # Point to install lib dir
)

# Install fortran tests
install(TARGETS test_fortran
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install the reference
install(TARGETS test_fortran
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)