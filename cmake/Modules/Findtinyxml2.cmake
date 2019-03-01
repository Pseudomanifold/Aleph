# Module: Findtinyxml2.cmake
# Author: Asem Alaa <asem.a.abdelaziz@ieee.org>
#
# CMake find module for tinyxml2 with Conan support.

INCLUDE( FindPackageHandleStandardArgs )

FIND_PATH( tinyxml2_INCLUDE_DIR
  NAMES
    tinyxml2.h
  PATHS
    ${CONAN_INCLUDE_DIRS_TINYXML2}
)

FIND_LIBRARY( tinyxml2_LIBRARY
  NAMES
    ${CONAN_LIBS_TINYXML2}
  PATHS
    ${CONAN_LIB_DIRS_TINYXML2}
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS( tinyxml2 DEFAULT_MSG
  tinyxml2_INCLUDE_DIR
)

IF( tinyxml2_FOUND )
  SET( tinyxml2_INCLUDE_DIRS ${tinyxml2_INCLUDE_DIR} )
  SET( tinyxml2_LIBRARIES ${tinyxml2_LIBRARY} )

  # Due to an immature recipe file of tinyxml2, the generated tinyxml2 files are not exported from
  # the build directory to the package directory. Issue reported here:
  # https://github.com/nicolastagliani/conan-tinyxml2/issues/3
  # So currenty, as a workaround, we will retrieve the tinyxml2Config.cmake file from the build directory!
  get_filename_component( tinyxml2_CONFIG_PATH ${CONAN_TINYXML2_ROOT} DIRECTORY )
  get_filename_component( tinyxml2_HASH ${CONAN_TINYXML2_ROOT} NAME )
  get_filename_component( tinyxml2_CONFIG_PATH ${tinyxml2_CONFIG_PATH} DIRECTORY )
  set( tinyxml2_CONFIG_PATH  ${tinyxml2_CONFIG_PATH}/build/${tinyxml2_HASH} )
  set( tinyxml2_CONFIG_FILENAME tinyxml2Config.cmake )

  find_file( tinyxml2_CONFIG_DIR
      ${tinyxml2_CONFIG_FILENAME}
    HINTS
      ${tinyxml2_CONFIG_PATH}
  )

  IF( tinyxml2_CONFIG_DIR-NOTFOUND )
    set( tinyxml2_CONFIG "" )
  ELSE()
    set( tinyxml2_CONFIG ${tinyxml2_CONFIG_DIR} )
  ENDIF()

  MARK_AS_ADVANCED(
    tinyxml2_INCLUDE_DIR
    tinyxml2_LIBRARY
    tinyxml2_DIR
    tinyxml2_CONFIG
  )
ELSE()
  SET( tinyxml2_DIR "" CACHE STRING
    "An optional hint to a tinyxml2 directory"
  )
ENDIF()
