# Workaround for CMake incompatability with a recent Clang compiler when
# it comes to the 'automoc' feature.
IF(   ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.0.0" )
   OR ( CMAKE_CXX_COMPILER_ID MATCHES "GNU" ) )

  SET( CMAKE_AUTOMOC ON )

  FIND_PACKAGE( Qt5 COMPONENTS Core Charts Widgets )

  IF( Qt5_FOUND )
    ADD_EXECUTABLE( Perdita
      DataSetItem.cc
      DataSetModel.cc
      ExpandSimplicialComplexDialog.cc
      LoadDataSetHelpers.cc
      MainWindow.cc
      Perdita.cc
      PersistenceDiagramHelpers.cc
      PersistenceDiagramView.cc
      PersistenceDiagramNormDialog.cc
      SimplicialComplexHelpers.cc
      UnsignedValidator.cc
    )
    TARGET_LINK_LIBRARIES( Perdita Qt5::Charts Qt5::Widgets )
  ENDIF()

ENDIF()
