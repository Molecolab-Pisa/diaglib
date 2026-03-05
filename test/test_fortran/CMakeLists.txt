add_executable(main.exe
        main.f90
        utility.f90
        direct_matvecs.f90
        direct_tests.f90
        sparsematrices_tests.f90
)
        
# Link DiagLib
target_link_libraries(main.exe PRIVATE diaglib)

# Ensure proper RPATH is set for installation
set_target_properties(main.exe PROPERTIES
    BUILD_WITH_INSTALL_RPATH ON    # Use install RPATH during build
    INSTALL_RPATH_USE_LINK_PATH OFF # Use link paths
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib64"  # Point to install lib dir
)


# Install matrices
install(DIRECTORY matrices/
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        FILES_MATCHING PATTERN "*.mtx"
)

# Install main
install(TARGETS main.exe
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)