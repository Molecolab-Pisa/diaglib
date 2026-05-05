# Source - https://stackoverflow.com/q/78281273
# Posted by JuliusCaesar, modified by community. See post 'Timeline' for change history
# Retrieved 2026-04-29, License - CC BY-SA 4.0

macro(create_venv)
  find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter Development)
  set(SYSTEM_PYTHON_EXE_PATH ${Python3_EXECUTABLE})
  message(STATUS "+++++++++++++++++  Found system python located at ${SYSTEM_PYTHON_EXE_PATH}++++++++++++++++")
  # Taken from https://discourse.cmake.org/t/possible-to-create-a-python-virtual-env-from-cmake-and-then-find-it-with-findpython3/1132
  execute_process(COMMAND "${Python3_EXECUTABLE}" -m venv "${PROJECT_BINARY_DIR}/.venv" --upgrade-deps COMMAND_ERROR_IS_FATAL ANY)

  # Here is the trick update the environment with VIRTUAL_ENV variable (mimic the activate script)
  set(VENV_PATH "${PROJECT_BINARY_DIR}/.venv")
  set(ENV{VIRTUAL_ENV} "${VENV_PATH}")

  set(Python3_FIND_VIRTUALENV FIRST)
  # unset Python3_EXECUTABLE because it is also an input variable (see documentation, Artifacts Specification section)
  unset(Python3_EXECUTABLE)
  # Launch a new search
  find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter Development)
  message(STATUS "+++++++++++++++++  Found venv python located at ${Python3_EXECUTABLE}++++++++++++++++")
  if(SYSTEM_PYTHON_EXE_PATH STREQUAL Python3_EXECUTABLE)
    message(FATAL_ERROR "Python3 executable is the same as the system one, this is not expected")
  endif()

endmacro()

macro(pipinstall package)

  # taken from https://www.scivision.dev/cmake-install-python-package/
  #set(REQUIREMENTS_ARG "-r ${TOOLBOX_SOURCE_DIR}/python3/requirements_generic.txt")
  message(STATUS "Trying to install ${package}")
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -m pip install ${package} COMMAND_ERROR_IS_FATAL ANY
  )

endmacro()