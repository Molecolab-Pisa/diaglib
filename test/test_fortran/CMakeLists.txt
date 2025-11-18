add_executable(main.exe test_fortran.f90
        ${CMAKE_SOURCE_DIR}/src/fortran/diaglib_global_utils.f90) 

# Link DiagLib and carry every dependency with it
target_link_libraries(main.exe PUBLIC diaglib)

