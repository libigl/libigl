# - Try to find the SVD3X3 library
# Once done this will define
#
#  SVD3X3_FOUND - system has SVD3X3
#  SVD3X3_INCLUDE_DIR - the SVD3X3 include directory

FIND_PATH(SVD3X3_INCLUDE_DIR Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.hpp
   /usr/include
   /usr/local/include
   ${PROJECT_SOURCE_DIR}/../libigl/external/Singular_Value_Decomposition/
   ${PROJECT_SOURCE_DIR}/../../external/Singular_Value_Decomposition/
   NO_DEFAULT_PATH
)

SET(SVD3X3_FOUND "NO")
IF (SVD3X3_INCLUDE_DIR)
	SET(SVD3X3_FOUND "YES")
ENDIF (SVD3X3_INCLUDE_DIR)

if(SVD3X3_INCLUDE_DIR)
   message(STATUS "Found SVD3X3: ${SVD3X3_INCLUDE_DIR}")
else(SVD3X3_INCLUDE_DIR)
  if (NOT SVD3X3_FIND_QUIETLY)
   message(FATAL_ERROR "could NOT find SVD3X3")
 endif(NOT SVD3X3_FIND_QUIETLY)
endif(SVD3X3_INCLUDE_DIR)

MARK_AS_ADVANCED(SVD3X3_INCLUDE_DIR SVD3X3_LIBRARIES)
