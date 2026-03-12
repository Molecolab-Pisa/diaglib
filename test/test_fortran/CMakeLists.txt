add_executable(test_fortran
        test_fortran.f90
        utility.f90
        direct_matvecs.f90
)
 
# Set the name of the executable
set_target_properties(test_fortran PROPERTIES
        OUTPUT_NAME test_fortran.exe)

# Link DiagLib
target_link_libraries(test_fortran PRIVATE diaglib)

# Ensure proper RPATH is set for installation
set_target_properties(test_fortran PROPERTIES
#        BUILD_WITH_INSTALL_RPATH ON    # Use install RPATH during build
        INSTALL_RPATH_USE_LINK_PATH OFF # Use link paths
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib64"  # Point to install lib dir
)

# Install main
install(TARGETS test_fortran
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)