add_library(diaglib_c SHARED
        dgl_utils_c.f90
        dgl_interface_c.f90
        davidson_driver_c.f90
        lobpcg_driver_c.f90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-std=f2018 -pedantic")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "-stand 2018")
endif()

target_link_libraries(diaglib_c PRIVATE diaglib)