add_executable(main.exe
        test_fortran.f90
        ) 

# Link DiagLib and carry every dependency with it
target_link_libraries(main.exe PUBLIC diaglib)
#target_include_directories(main.exe PUBLIC ${CMAKE_SOURCE_DIR}/src/fortran)