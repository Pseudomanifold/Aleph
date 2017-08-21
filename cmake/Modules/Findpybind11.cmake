# Module: Findpybind11.cmake
# Author: Bastian Rieck <bastian.rieck@iwr.uni-heidelberg.de>
#
# CMake find module for pybind11.

INCLUDE( FindPackageHandleStandardArgs )

SET_IF_EMPTY( PYBIND11_DIR "$ENV{PYBIND11_DIR}" )

FIND_PATH(
  PYBIND11_INCLUDE_DIR
    pybind11/pybind11.h
  HINTS
    ${PYBIND11_DIR}
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS( PYBIND11 DEFAULT_MSG
  PYBIND11_INCLUDE_DIR
)

IF( PYBIND11_FOUND )
  SET( PYBIND11_INCLUDE_DIRS ${PYBIND11_INCLUDE_DIR} )

  MARK_AS_ADVANCED(
    PYBIND11_INCLUDE_DIR
    PYBIND11_DIR
  )
ELSE()
  SET( PYBIND11_DIR "" CACHE STRING
    "An optional hint to a pybind11 directory"
  )
ENDIF()
