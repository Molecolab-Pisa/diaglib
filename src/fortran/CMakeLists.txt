add_library(diaglib SHARED
        dgl_global_utils.f90
        dgl_orthogonalizations.f90
        dgl_minor_utils.f90
        dgl_external_interfaces.f90
        davidson_driver.f90
        davidson_nosym_driver.f90
        lobpcg_driver.f90
        smogd_driver.f90
        dgl_interface.f90
)
        
target_include_directories(diaglib PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Linking proper BLAS
if(MKL_FOUND)
        target_link_libraries(diaglib PUBLIC MKL::MKL)
else()
        target_link_libraries(diaglib PUBLIC BLAS::BLAS)
endif()


# Set Fortran flags depending on build type
STRING (TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_l)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
   
        target_compile_options(diaglib PRIVATE -std=f2003 -pedantic)

        if(BUILD_TYPE_l MATCHES "debug")
                target_compile_options(diaglib PRIVATE
                        -fbacktrace -fcheck=all -Wall -Wextra)
        endif()

elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")

        target_compile_options(diaglib PRIVATE -stand f03)

        if(BUILD_TYPE_l MATCHES "debug")
                target_compile_options(diaglib PRIVATE
                        -traceback -check all -warn all)
        endif()

endif()

# Print out flags that will be used for diaglib compilation
get_target_property(_ diaglib COMPILE_OPTIONS)
message(STATUS "Diaglib compile options are: ${_}")
        

# Linking proper OpenMP
target_link_libraries(diaglib PUBLIC OpenMP::OpenMP_Fortran)


# Install libray
install(TARGETS diaglib
        EXPORT diaglibTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Export DiagLib as a cmake package 
install(EXPORT diaglibTargets
    FILE DiaglibTargets.cmake
    NAMESPACE diaglib::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Diaglib
)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/diaglibConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Diaglib
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/diaglibConfig.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Diaglib
)