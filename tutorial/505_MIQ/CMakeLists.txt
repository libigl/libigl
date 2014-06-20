cmake_minimum_required(VERSION 2.6)
project(XXX_MIQ)

include("../CMakeLists.shared")
include_directories("/Users/olkido/Documents/igl/MIQ/src/CoMISo/Solver")
include_directories("/Users/olkido/Documents/igl/MIQ/src/CoMISo/gmm/include")
include_directories("/Users/olkido/Documents/igl/MIQ/src")
set(SOURCES
${PROJECT_SOURCE_DIR}/main.cpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DINCLUDE_TEMPLATES")

add_executable(XXX_MIQ ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(XXX_MIQ ${SHARED_LIBRARIES} "/Users/olkido/Documents/igl/MIQ/src/CoMISo/build/libCoMISo.dylib")
