add_executable(main.exe test_fortran.f90 ${CMAKE_SOURCE_DIR}/src/fortran/1_diaglib_global_utils.f90 ${CMAKE_SOURCE_DIR}/src/fortran/diaglib.f90)

target_link_libraries(main.exe PUBLIC diaglib)
target_include_directories(main.exe PUBLIC ${CMAKE_SOURCE_DIR}/src/fortran)
