if(NOT Python3_EXECUTABLE)
        message(FATAL_ERROR "No python interpreter found for installation")
endif()

#pipinstall(build)
#add_custom_target(pyDiaglib ALL
#                COMMAND ${Python3_EXECUTABLE} -m build ${CMAKE_CURRENT_SOURCE_DIR})
pipinstall(${CMAKE_CURRENT_SOURCE_DIR})