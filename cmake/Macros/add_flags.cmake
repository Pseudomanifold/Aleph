# Macro : add_flags.cmake
# Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
#         Bastian Rieck <bastian.rieck@iwr.uni-heidelberg.de>
#
# Provides macro ADD_FLAGS, which appends a string "${flags}" to a given
# variable "${variable}".
#
# Usage:
#   ADD_FLAGS( variable flags )

MACRO( ADD_FLAGS variable flags )
  STRING( STRIP "${flags}" flagsStripped )

  IF( NOT "${flagsStripped}" STREQUAL "" )
    SET( ${variable} "${${variable}} ${flags}" )
    STRING( STRIP "${${variable}}" ${variable} )
  ENDIF()
ENDMACRO()
