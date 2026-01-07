add_executable(main.exe
        test_fortran.f90
        ) 

# Put the test into the bin directory
set_property(TARGET main.exe
        PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
        
# Link DiagLib
target_link_libraries(main.exe PUBLIC diaglib)
