# - Try to find the LIM library
# Once done this will define
#
#  LIM_FOUND - system has LIM
#  LIM_INCLUDE_DIR - the LIM include directory
#  LIM_SOURCES - the LIM source files

FIND_PATH(LIM_INCLUDE_DIR LIMSolverInterface.h
   /usr/include
   /usr/local/include
   ${PROJECT_SOURCE_DIR}/../libigl/external/lim/
   ${PROJECT_SOURCE_DIR}/../../external/lim/
   NO_DEFAULT_PATH
)

set(
  LIM_SOURCES
  ${LIM_INCLUDE_DIR}/NMSolver.cpp
  ${LIM_INCLUDE_DIR}/LIMSolver.cpp
  ${LIM_INCLUDE_DIR}/LIMSolver2D.cpp
  ${LIM_INCLUDE_DIR}/LIMSolver3D.cpp
  ${LIM_INCLUDE_DIR}/TriangleMesh.cpp
  ${LIM_INCLUDE_DIR}/TetrahedronMesh.cpp
  ${LIM_INCLUDE_DIR}/Dirichlet_LIMSolver2D.cpp
  ${LIM_INCLUDE_DIR}/Dirichlet_LIMSolver3D.cpp
  ${LIM_INCLUDE_DIR}/UniformLaplacian_LIMSolver2D.cpp
  ${LIM_INCLUDE_DIR}/UniformLaplacian_LIMSolver3D.cpp
  ${LIM_INCLUDE_DIR}/Laplacian_LIMSolver2D.cpp
  ${LIM_INCLUDE_DIR}/Laplacian_LIMSolver3D.cpp
  ${LIM_INCLUDE_DIR}/LGARAP_LIMSolver2D.cpp
  ${LIM_INCLUDE_DIR}/LGARAP_LIMSolver3D.cpp
  ${LIM_INCLUDE_DIR}/GreenStrain_LIMSolver2D.cpp
  ${LIM_INCLUDE_DIR}/GreenStrain_LIMSolver3D.cpp
  ${LIM_INCLUDE_DIR}/LSConformal_LIMSolver2D.cpp
  ${LIM_INCLUDE_DIR}/Poisson_LIMSolver2D.cpp
  )

SET(LIM_FOUND "NO")
IF (LIM_INCLUDE_DIR)
	SET(LIM_FOUND "YES")
ENDIF (LIM_INCLUDE_DIR)

if(LIM_INCLUDE_DIR)
   message(STATUS "Found LIM: ${LIM_INCLUDE_DIR}")
else(LIM_INCLUDE_DIR)
  if (NOT LIM_FIND_QUIETLY)
   message(FATAL_ERROR "could NOT find LIM")
 endif(NOT LIM_FIND_QUIETLY)
endif(LIM_INCLUDE_DIR)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DANSI_DECLARATORS")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DANSI_DECLARATORS")

MARK_AS_ADVANCED(LIM_INCLUDE_DIR LIM_LIBRARIES LIM_SOURCES)
