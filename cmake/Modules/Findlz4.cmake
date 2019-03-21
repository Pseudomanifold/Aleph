# Module: Findlz4.cmake
# Author: Bastian Rieck <bastian.rieck@bsse.ethz.ch>
#
# CMake find module for `lz4` (a package for compressing data). If not
# present in your local systems, the libraries may be installed as one
# of the following names:
#
# - lz4 (Homebrew, Ubuntu, ...)
# - aur/lz4-git

INCLUDE( FindPackageHandleStandardArgs )

SET_IF_EMPTY( FLANN_DIR "$ENV{LZ4_DIR}" )

FIND_PATH(
  LZ4_INCLUDE_DIR
    lz4.h
  HINTS
    ${LZ4_DIR}
)

FIND_LIBRARY( LZ4_LIBRARY
  NAMES lz4
  HINTS ${LZ4_DIR}
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS( LZ4 DEFAULT_MSG
  LZ4_INCLUDE_DIR
  LZ4_LIBRARY
)

IF( LZ4_FOUND )
  SET( LZ4_LIBRARIES ${LZ4_LIBRARY} )
  SET( LZ4_INCLUDE_DIRS ${LZ4_INCLUDE_DIR} )

  MARK_AS_ADVANCED(
    LZ4_LIBRARY
    LZ4_INCLUDE_DIR
    LZ4_DIR
  )
ELSE()
  SET( LZ4_DIR "" CACHE STRING
    "An optional hint to a LZ4 directory"
  )
ENDIF()
