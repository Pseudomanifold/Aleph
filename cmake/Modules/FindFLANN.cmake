# Module: FindFLANN.cmake
# Author: Bastian Rieck <bastian.rieck@iwr.uni-heidelberg.de>
#
# CMake find module for FLANN (Fast Library for Approximate Nearest Neighbors).
# If not present in the local packaging system, the libraries may be obtained
# from:
#
#   http://www.cs.ubc.ca/~mariusm/index.php/FLANN/FLANN
#
# Some package names for FLANN:
#
# * libflann-dev
# * aur/lib/flann

INCLUDE( FindPackageHandleStandardArgs )

SET_IF_EMPTY( FLANN_DIR "$ENV{FLANN_DIR}" )

FIND_PATH(
  FLANN_INCLUDE_DIR
    flann/flann.hpp
  HINTS
    ${FLANN_DIR}
)

FIND_LIBRARY( FLANN_LIBRARY
  NAMES flann_cpp
  HINTS ${FLANN_DIR}
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS( FLANN DEFAULT_MSG
  FLANN_INCLUDE_DIR
  FLANN_LIBRARY
)

IF( FLANN_FOUND )
  SET( FLANN_LIBRARIES ${FLANN_LIBRARY} )
  SET( FLANN_INCLUDE_DIRS ${FLANN_INCLUDE_DIR} )

  FILE( STRINGS "${FLANN_INCLUDE_DIR}/flann/config.h" FLANN_CONFIG )
  FOREACH( LINE ${FLANN_CONFIG} )
    STRING( REGEX MATCH
      "FLANN_VERSION_ \"([0-9\\.]+)\""
      LINE_MATCHES
      ${LINE}
    )

    # We have a match, i.e. there is some sort of FLANN version
    # identification available.
    IF( NOT "${LINE_MATCHES}" STREQUAL "" )
      SET(
        FLANN_VERSION
          ${CMAKE_MATCH_1}
        CACHE STRING
          "Detected version of the FLANN library"
      )
    ENDIF()
  ENDFOREACH()

  MARK_AS_ADVANCED(
    FLANN_LIBRARY
    FLANN_INCLUDE_DIR
    FLANN_VERSION
    FLANN_DIR
  )
ELSE()
  SET( FLANN_DIR "" CACHE STRING
    "An optional hint to a FLANN directory"
  )
ENDIF()
