IF( PKG_CONFIG_FOUND )
  PKG_SEARCH_MODULE( Python3 python3 )

  IF( Python3_FOUND )
    INCLUDE_DIRECTORIES(
      ${PYBIND11_INCLUDE_DIRS}
      ${Python3_INCLUDE_DIRS}
    )

    ADD_LIBRARY( aleph aleph.cc )
  ENDIF()
ENDIF()
