add_executable(main.exe
        test_fortran.f90
        ) 

        
# Link DiagLib
target_link_libraries(main.exe PUBLIC diaglib)
