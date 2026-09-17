add_library(diaglib_c SHARED
        dgl_utils_c.f90
        davidson_driver_c.f90
        lobpcg_driver_c.f90
        davidson_nosym_driver_c.f90
        smogd_driver_c.f90
)

# Generate the header that tells C users the size of DiagLib integers
configure_file(diaglib_config.h.in
        "${CMAKE_CURRENT_BINARY_DIR}/diaglib_config.h" @ONLY)

# Include the directories with the C headers
target_include_directories(diaglib_c PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# The interface needs a more recent standard
target_compile_options(diaglib_c PRIVATE ${DGL_Fortran_F2018} ${DGL_Fortran_OPTIONS})

# Link actual diaglib
target_link_libraries(diaglib_c PRIVATE diaglib)

# OpenMP makes the stored callbacks threadprivate, so that the drivers can be called
# from different threads at the same time
if(OpenMP_Fortran_FOUND)
        target_link_libraries(diaglib_c PRIVATE OpenMP::OpenMP_Fortran)
endif()

# Let the installed library find libdiaglib in the same directory
set_target_properties(diaglib_c PROPERTIES
        INSTALL_RPATH "$ORIGIN"
        INSTALL_RPATH_USE_LINK_PATH TRUE)

# Install libray
install(TARGETS diaglib_c
        EXPORT diaglibTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install the C headers
install(FILES
        diaglib.h
        "${CMAKE_CURRENT_BINARY_DIR}/diaglib_config.h"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
