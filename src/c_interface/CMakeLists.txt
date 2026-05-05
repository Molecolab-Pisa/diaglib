add_library(diaglib_c SHARED
        dgl_utils_c.f90
        davidson_driver_c.f90
        lobpcg_driver_c.f90
        davidson_nosym_driver_c.f90
        smogd_driver_c.f90
)

target_include_directories(diaglib_c PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)


if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-std=f2018 -pedantic")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "-stand 2018")
endif()

target_link_libraries(diaglib_c PRIVATE diaglib)

install(FILES diaglib.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)