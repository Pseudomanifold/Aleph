CMAKE_MINIMUM_REQUIRED( VERSION 2.8.7 )
PROJECT( Aleph CXX )

IF( CMAKE_MAJOR_VERSION GREATER 3 OR ( CMAKE_MAJOR_VERSION EQUAL 3 AND CMAKE_MINOR_VERSION GREATER 0 ) )
  SET( CMAKE_CXX_STANDARD 11 )
  SET( CMAKE_CXX_STANDARD_REQUIRED ON )
ELSE()
  # This is not the nicest way of activating C++11, but it is guaranteed to
  # work with older versions
  ADD_DEFINITIONS( "-std=c++11" )
ENDIF()

LIST( APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/Modules
)

SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )

MESSAGE( STATUS "CMAKE_ROOT:           ${CMAKE_ROOT}" )
MESSAGE( STATUS "CMAKE_SYSTEM_VERSION: ${CMAKE_SYSTEM_VERSION}" )

########################################################################
# Versioning
########################################################################

SET( ALEPH_MAJOR_VERSION 0 )
SET( ALEPH_MINOR_VERSION 1 )
SET( ALEPH_PATCH_VERSION 0 )

SET( ALEPH_VERSION ${ALEPH_MAJOR_VERSION}.${ALEPH_MINOR_VERSION}.${ALEPH_PATCH_VERSION} )

########################################################################
# Macros and additional flags
########################################################################

INCLUDE( CheckCXXCompilerFlag )

FILE( GLOB filenames "cmake/Macros/*.cmake" )
FOREACH( filename ${filenames} )
  INCLUDE( ${filename} )
ENDFOREACH()

ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wall" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wconversion" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wextra" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wnon-virtual-dtor" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wold-style-cast" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Woverloaded-virtual" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wself-init" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wunsafe-loop-optimization" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-pedantic" )

########################################################################
# Compile options
########################################################################

SET( BUILD_PYTHON_BINDINGS
  "ON"
  CACHE
  BOOL
  "Build with Python bindings"
)

SET( BUILD_EXAMPLES
  "ON"
  CACHE
  BOOL
  "Build with examples"
)

SET( BUILD_TOOLS
  "ON"
  CACHE
  BOOL
  "Build with tools"
)

########################################################################
# Additional packages
########################################################################

FIND_PACKAGE( Boost REQUIRED )
FIND_PACKAGE( OpenMP )

# This is somewhat brute-force: I am adding OpenMP flags regardless of
# the target. However, this makes it easier to use code for which this
# framework may or may not be used.
#
# I do not want to have too many of these declarations for different
# targets flying around.
#
IF( OpenMP_CXX_FOUND )
  # Workaround: OpenMP with Clang breaks the build for Travis CI. Hence,
  # disabling this feature for now, at least on older versions of Clang.
  IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.0.0" ) )
    ADD_FLAGS( CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS} )
  ENDIF()
ENDIF()

FIND_PACKAGE( FLANN )

IF( FLANN_FOUND )
  SET( ALEPH_WITH_FLANN TRUE )
ENDIF()

FIND_PACKAGE( HDF5 COMPONENTS CXX )

IF( HDF5_FOUND )
  SET( ALEPH_WITH_HDF5 TRUE )
ENDIF()

FIND_PACKAGE( PkgConfig )

IF( PKG_CONFIG_FOUND )
  PKG_SEARCH_MODULE( JSON RapidJSON )

  IF( JSON_FOUND )
    # RapidJSON is now available in the normal include path, under the
    # prefix of RapidJSON.
    INCLUDE_DIRECTORIES( SYSTEM ${JSON_INCLUDE_DIRS} )

    SET( ALEPH_WITH_RAPID_JSON TRUE )
  ENDIF()

  PKG_SEARCH_MODULE( EIGEN eigen3 )

  IF( EIGEN_FOUND )
    INCLUDE_DIRECTORIES( SYSTEM ${EIGEN_INCLUDE_DIRS} )

    SET( ALEPH_WITH_EIGEN TRUE )
  ENDIF()

  # Not used here, but in subordinate directories, so it is nicer to
  # look for the package in a *single* location.
  PKG_SEARCH_MODULE( Python3 python3 )
ENDIF()

FIND_PACKAGE( tinyxml2 )

IF( tinyxml2_FOUND )
  SET( ALEPH_WITH_TINYXML2 TRUE )
ENDIF()

########################################################################
# Configuration files
########################################################################

CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/include/aleph/config/Eigen.hh.in ${CMAKE_SOURCE_DIR}/include/aleph/config/Eigen.hh )
CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/include/aleph/config/FLANN.hh.in ${CMAKE_SOURCE_DIR}/include/aleph/config/FLANN.hh )
CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/include/aleph/config/HDF5.hh.in ${CMAKE_SOURCE_DIR}/include/aleph/config/HDF5.hh )
CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/include/aleph/config/RapidJSON.hh.in ${CMAKE_SOURCE_DIR}/include/aleph/config/RapidJSON.hh )
CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/include/aleph/config/TinyXML2.hh.in ${CMAKE_SOURCE_DIR}/include/aleph/config/TinyXML2.hh )

########################################################################
# Include directories
########################################################################

INCLUDE_DIRECTORIES(
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/include
)

ADD_SUBDIRECTORY( bindings )
ADD_SUBDIRECTORY( include )
ADD_SUBDIRECTORY( src )
ADD_SUBDIRECTORY( examples )

########################################################################
# Tests
########################################################################

# Enable testing globally. Individual tests may be found in the `tests`
# subdirectory of the repository.
ENABLE_TESTING()

ADD_SUBDIRECTORY( tests )

########################################################################
# Place compile commands in the source directory
########################################################################


IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()

########################################################################
# Add static analysis targets: `cppcheck`
########################################################################

ADD_CUSTOM_TARGET( cppcheck
  COMMAND cppcheck
            --enable=all
            --std=c++11
            --template "{file},{line},{severity},{message}"
            --verbose
            --project=${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
             2> /tmp/cppcheck.csv
)

########################################################################
# Package setup
########################################################################

EXPORT( PACKAGE Aleph )

SET( ALEPH_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include )

INCLUDE( CMakePackageConfigHelpers )

CONFIGURE_PACKAGE_CONFIG_FILE(
  ${CMAKE_SOURCE_DIR}/cmake/AlephConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/AlephConfig.cmake
  INSTALL_DESTINATION lib/cmake/
  PATH_VARS ALEPH_INCLUDE_DIR
)

WRITE_BASIC_PACKAGE_VERSION_FILE( ${CMAKE_CURRENT_BINARY_DIR}/AlephConfigVersion.cmake
  VERSION       ${ALEPH_VERSION}
  COMPATIBILITY AnyNewerVersion
)

INSTALL(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/AlephConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/AlephConfigVersion.cmake
  DESTINATION
    lib/cmake
  COMPONENT
    library
)
