IF( BUILD_PYTHON_BINDINGS )
  IF( Python3_FOUND )
    INCLUDE_DIRECTORIES(
      ${PYBIND11_INCLUDE_DIRS}
      ${Python3_INCLUDE_DIRS}
    )

    ADD_LIBRARY( aleph SHARED aleph.cc )

    SET_TARGET_PROPERTIES( aleph
      PROPERTIES
        PREFIX ""
    )
  ELSE()
    MESSAGE( STATUS "Cannot build Python bindings because no suitable version of Python was found" )
  ENDIF()
ENDIF()
