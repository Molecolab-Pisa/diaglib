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

set_target_properties(diaglib PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_INCLUDEDIR})

target_include_directories(diaglib PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Linking BLAS and LAPACK
target_link_libraries(diaglib PUBLIC BLAS::BLAS)
target_link_libraries(diaglib PUBLIC LAPACK::LAPACK)


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

# Export DiagLib target
install(EXPORT diaglibTargets
        FILE DiaglibTargets.cmake
        NAMESPACE diaglib::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Diaglib
)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/diaglibConfig.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Diaglib
)

install(DIRECTORY 
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/..)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/diaglibConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Diaglib
)