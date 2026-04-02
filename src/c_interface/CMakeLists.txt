add_library(diaglib_c SHARED
        dgl_utils_c.f90
        davidson_driver_c.f90
        dgl_interface_c.f90
)

target_include_directories(diaglib PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(diaglib_c PRIVATE diaglib)