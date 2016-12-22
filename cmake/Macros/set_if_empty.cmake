# Module: set_if_empty
# Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
#         Bastian Rieck <bastian.rieck@iwr.uni-heidelberg.de>
#
# Provides macro SET_IF_EMPTY. Given a variable and a value, the variable is
# set to the value if it is still empty.

MACRO( SET_IF_EMPTY variable value )
  IF( "${${variable}}" STREQUAL "" )
    SET( ${variable} ${value} )
  ENDIF()
ENDMACRO()
