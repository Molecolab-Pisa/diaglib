# test diaglib directly from C
add_executable(dgl_test_c
        test_c.c
)

# Include header file
target_include_directories(dgl_test_c PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../c_include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Link DiagLib
target_link_libraries(dgl_test_c PRIVATE diaglib_c)

# Ensure proper RPATH is set for installation
set_target_properties(dgl_test_c PROPERTIES
#        BUILD_WITH_INSTALL_RPATH ON    # Use install RPATH during build
        INSTALL_RPATH_USE_LINK_PATH OFF # Use link paths
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"  # Point to install lib dir
)

# Install C tests
install(TARGETS dgl_test_c
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES diaglib.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Run test
add_test(NAME dgl_test_c COMMAND dgl_test_c)