add_executable(main.exe
        test_solvers.f90
        matvecs.f90
        test_fortran.f90
)

        
# Link DiagLib
target_link_libraries(main.exe PRIVATE diaglib)