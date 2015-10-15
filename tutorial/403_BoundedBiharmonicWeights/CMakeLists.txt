cmake_minimum_required(VERSION 2.6)
project(403_BoundedBiharmonicWeights)

include("../CMakeLists.shared")

find_package(MOSEK QUIET)
if(MOSEK_FOUND)
  include_directories( ${MOSEK_INCLUDE_DIR} )
else(MOSEK_FOUND)
  add_definitions(-DIGL_NO_MOSEK)
  if(LIBIGL_USE_STATIC_LIBRARY)
    set(LIBIGLMOSEK_LIBRARY "")
  endif(LIBIGL_USE_STATIC_LIBRARY)
  set(MOSEK_LIBRARIES "")
endif(MOSEK_FOUND)

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

add_executable(${PROJECT_NAME}_bin ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES} ${MOSEK_LIBRARIES})
