# Macro : enable_if_supported
# Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
#         Bastian Rieck <bastian.rieck@iwr.uni-heidelberg.de>
#
# Provides macro ENABLE_IF_SUPPORTED, which checks whether the CXX compiler
# understands a given compiler flag. If so, it is added to the given variable.
#
# Usage:
#   ENABLE_IF_SUPPORTED( variable flag )

MACRO( ENABLE_IF_SUPPORTED variable flag )
  STRING( REGEX REPLACE "^-" "" flagName "${flag}" )
  CHECK_CXX_COMPILER_FLAG( "${flag}"
    ALEPH_HAVE_FLAG_${flagName}
    )
  IF( ALEPH_HAVE_FLAG_${flagName} )
    SET( ${variable} "${${variable}} ${flag}")
  ENDIF()
ENDMACRO()
