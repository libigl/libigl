cmake_minimum_required(VERSION 2.6)
project(403_BoundedBiharmonicWeights)

include("../CMakeLists.shared")

#find_package(EMBREE REQUIRED)
#include_directories(${EMBREE_INCLUDE_DIRS})

if(NOT MOSEK_FOUND)
  add_definitions(-DIGL_NO_MOSEK)
  if(LIBIGL_USE_STATIC_LIBRARY)
    set(LIBIGLMOSEK_LIBRARY "")
  endif(LIBIGL_USE_STATIC_LIBRARY)
  set(MOSEK_LIBRARIES "")
else(NOT MOSEK_FOUND)
  include_directories( ${MOSEK_INCLUDE_DIR} )
endif(NOT MOSEK_FOUND)

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME}_bin ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES} ${MOSEK_LIBRARIES})
#target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES} ${MOSEK_LIBRARIES} ${EMBREE_LIBRARIES})
