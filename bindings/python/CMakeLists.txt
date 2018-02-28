IF( BUILD_PYTHON_BINDINGS )
  FIND_PACKAGE( pybind11 )
  IF( Python3_FOUND AND PYBIND11_FOUND )
    INCLUDE_DIRECTORIES(
      ${PYBIND11_INCLUDE_DIRS}
      ${Python3_INCLUDE_DIRS}
    )

    ADD_LIBRARY( aleph SHARED aleph.cc )

    SET_TARGET_PROPERTIES( aleph
      PROPERTIES
        PREFIX ""
    )

    IF( APPLE )
      SET_TARGET_PROPERTIES( aleph
        PROPERTIES
          LINK_FLAGS "-undefined dynamic_lookup"
          SUFFIX     ".so"
      )
    ENDIF()

    IF( HDF5_FOUND )
      TARGET_INCLUDE_DIRECTORIES( aleph SYSTEM PUBLIC ${HDF5_CXX_INCLUDE_DIRS} )
      TARGET_LINK_LIBRARIES(      aleph               ${HDF5_CXX_LIBRARIES} )
    ENDIF()

    IF( FLANN_FOUND AND FLANN_VERSION VERSION_GREATER "1.8.4" )
      TARGET_LINK_LIBRARIES( aleph lz4 )
    ENDIF()
  ELSE()
    MESSAGE( STATUS "Cannot build Python bindings because no suitable version of Python was found" )
  ENDIF()
ENDIF()
