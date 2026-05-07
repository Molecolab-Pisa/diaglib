if(NOT Python3_EXECUTABLE)
        message(FATAL_ERROR "No python interpreter found for installation")
endif()

# Install pyDiaglib as a pip package in the previously created venv
add_custom_target(pyDiaglib ALL COMMAND ${Python3_EXECUTABLE} -m pip install ${CMAKE_CURRENT_SOURCE_DIR})

# Install the necessary stuff
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src"
        DESTINATION ${CMAKE_INSTALL_PREFIX})

install(FILES pyproject.toml
        DESTINATION ${CMAKE_INSTALL_PREFIX})