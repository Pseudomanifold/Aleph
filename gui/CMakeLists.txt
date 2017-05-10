SET( CMAKE_AUTOMOC ON )

FIND_PACKAGE( Qt5 COMPONENTS Core Charts Widgets )

IF( Qt5_FOUND )
  ADD_EXECUTABLE( Perdita
    Perdita.cc
    MainWindow.cc
  )
  TARGET_LINK_LIBRARIES( Perdita Qt5::Charts Qt5::Widgets )
ENDIF()
