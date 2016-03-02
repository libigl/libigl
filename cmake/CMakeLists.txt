cmake_minimum_required(VERSION 2.6)
project(libigl)

### Available options ###
option(LIBIGL_USE_STATIC_LIBRARY "Use libIGL as static library" OFF)
option(LIBIGL_WITH_BBW         "Use BBW"            OFF)
option(LIBIGL_WITH_BOOLEAN     "Use Cork boolean"   OFF)
option(LIBIGL_WITH_CGAL        "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO      "Use CoMiso"         OFF)
option(LIBIGL_WITH_EMBREE      "Use Embree"         OFF)
option(LIBIGL_WITH_LIM         "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB      "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK       "Use MOSEK"          OFF)
option(LIBIGL_WITH_OPENGL      "Use OpenGL"         OFF)
option(LIBIGL_WITH_PNG         "Use PNG"            OFF)
option(LIBIGL_WITH_TETGEN      "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE    "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER      "Use Nanogui Viewer" OFF)
option(LIBIGL_WITH_XML         "Use XML"            OFF)

### Compilation configuration ###
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

#if(MSVC)
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj") ### Enable parallel compilation for Visual Studio
#endif()


### Compiling libraries based on chosen options ###
set(LIBIGL_INCLUDE_DIRS "")
set(LIBIGL_LIBRARIES "")
set(LIBIGL_EXTRA_LIBRARIES "")
set(LIBIGL_EXTRA_SOURCES "")
set(LIBIGL_DEFINITIONS "")

if(LIBIGL_USE_STATIC_LIBRARY)
  add_definitions(-DIGL_STATIC_LIBRARY)
endif()

#if(APPLE)
#  set(CMAKE_CXX_LINK_FLAGS "-framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
#endif()

#if(UNIX AND NOT APPLE)
#  set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -lGL -lGLU -lrt -lX11 -lXxf86vm -lXrandr -lpthread -lXi  -lXcursor -lXinerama ")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
#endif()

### macro definition ###
set(LIBIGL_ROOT "${PROJECT_SOURCE_DIR}/..")
set(LIBIGL_SOURCE_DIR "${LIBIGL_ROOT}/include")
set(LIBIGL_EXTERNAL "${LIBIGL_ROOT}/external")

macro(CompileIGL_Module module)
  file(GLOB SOURCES_${module}
  "${LIBIGL_SOURCE_DIR}/igl/${module}/*.cpp")
  add_library(igl${module} STATIC ${SOURCES_${module}})
  list(APPEND LIBIGL_LIBRARIES "igl${module}")
endmacro()


### Eigen ###
set(EIGEN_INCLUDE_DIR "${LIBIGL_EXTERNAL}/nanogui/ext/eigen")
include_directories(${EIGEN_INCLUDE_DIR})
list(APPEND LIBIGL_INCLUDE_DIRS "${EIGEN_INCLUDE_DIR}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations") # remove some annoying warning for Eigen (at least with gcc)

### OpenMP ###
find_package(OpenMP)
if(OPENMP_FOUND AND NOT WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

### IGL Common ###
list(APPEND LIBIGL_INCLUDE_DIRS "${LIBIGL_SOURCE_DIR}")
if(LIBIGL_USE_STATIC_LIBRARY)
  include_directories("${LIBIGL_SOURCE_DIR}")
  file(GLOB SOURCES_IGL "${LIBIGL_SOURCE_DIR}/igl/*.cpp")
  add_library(igl STATIC ${SOURCES_IGL})
  list(APPEND LIBIGL_LIBRARIES "igl")
endif()


### Compile the mosek part ###
if(LIBIGL_WITH_MOSEK)
  find_package(MOSEK REQUIRED)
  if(NOT MOSEK_FOUND)
    message(SEND_ERROR "Could NOT find MOSEK")
  endif()
  list(APPEND LIBIGL_INCLUDE_DIRS ${MOSEK_INCLUDE_DIR})
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${MOSEK_LIBRARIES})

  if(LIBIGL_USE_STATIC_LIBRARY)
    include_directories(${MOSEK_INCLUDE_DIR})
    CompileIGL_Module("mosek")
  endif()
else()
  add_definitions(-DIGL_NO_MOSEK)
  list(APPEND LIBIGL_DEFINITIONS "-DIGL_NO_MOSEK")
endif()

### Compile the BBW part ###
if(LIBIGL_WITH_BBW)
  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("bbw")
  endif()
endif()

### Compile the comiso part ###
if(LIBIGL_WITH_COMISO)
  find_package(LIBCOMISOH REQUIRED)
  if(NOT LIBCOMISOH_FOUND)
    message(SEND_ERROR "Could NOT find COMISO")
  endif()
  #include_directories( ${LIBCOMISO_INCLUDE_DIRS})
  list(GET LIBCOMISO_INCLUDE_DIRS 0 COMISO_ROOT)

  add_subdirectory("${COMISO_ROOT}/CoMISo/" "CoMISo")

  set(COMISO_INCLUDE_DIRS
    "${COMISO_ROOT}/CoMISo/ext/gmm-4.2/include"
    "${COMISO_ROOT}/")

  list(APPEND LIBIGL_INCLUDE_DIRS ${COMISO_INCLUDE_DIRS})

  if(APPLE)
    find_library(accelerate_library Accelerate)
    list(APPEND LIBIGL_EXTRA_LIBRARIES "CoMISo" ${accelerate_library})
  else()
    list(APPEND LIBIGL_EXTRA_LIBRARIES "CoMISo")
  endif()

  if(LIBIGL_USE_STATIC_LIBRARY)
    include_directories(${COMISO_INCLUDE_DIRS})
    CompileIGL_Module("comiso")
  endif()
endif()

### Compile the embree part ###
if(LIBIGL_WITH_EMBREE)
  set(EMBREE_DIR "${LIBIGL_EXTERNAL}/embree")

  option(ENABLE_ISPC_SUPPORT OFF)
  option(RTCORE_TASKING_SYSTEM "INTERNAL")
  option(ENABLE_TUTORIALS OFF)
  option(ENABLE_STATIC_LIB ON)

  if(MSVC)
    option(USE_STATIC_RUNTIME " " OFF)
  endif()

  add_subdirectory("${EMBREE_DIR}" "embree")
  list(APPEND LIBIGL_INCLUDE_DIRS "${EMBREE_DIR}/include")
  list(APPEND LIBIGL_EXTRA_LIBRARIES "embree")

  if(LIBIGL_USE_STATIC_LIBRARY)
    include_directories(${EMBREE_DIR}/include)
    CompileIGL_Module("embree")
  endif()
endif()

### Compile the lim part ###
if(LIBIGL_WITH_LIM)
  set(LIM_DIR "${LIBIGL_EXTERNAL}/lim")
  add_subdirectory("${LIM_DIR}" "lim")
  list(APPEND LIBIGL_INCLUDE_DIRS ${LIM_DIR})
  list(APPEND LIBIGL_LIBRARIES "lim") # it depends on ligigl, so placing it here solve linking problems

  if(LIBIGL_USE_STATIC_LIBRARY)
    include_directories(${LIM_DIR})
    CompileIGL_Module("lim")
  endif()
endif()

### Compile the matlab part ###
if(LIBIGL_WITH_MATLAB)
  find_package(MATLAB REQUIRED)
  if(NOT MATLAB_FOUND)
    message(SEND_ERROR "Could NOT find MATLAB")
  endif()
  list(APPEND LIBIGL_INCLUDE_DIRS ${MATLAB_INCLUDE_DIR})
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${MATLAB_LIBRARIES})

  if(LIBIGL_USE_STATIC_LIBRARY)
    include_directories(${MATLAB_INCLUDE_DIR})
    CompileIGL_Module("matlab")
  endif()
endif()


### Compile the opengl and png parts ###
if(LIBIGL_WITH_OPENGL)
  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("opengl")
    CompileIGL_Module("opengl2")
  endif()

  if(LIBIGL_WITH_PNG)
    set(STB_IMAGE_DIR "${LIBIGL_EXTERNAL}/stb_image")
    add_subdirectory("${STB_IMAGE_DIR}" "stb_image")
    list(APPEND LIBIGL_INCLUDE_DIRS ${STB_IMAGE_DIR})
    list(APPEND LIBIGL_EXTRA_LIBRARIES "stb_image")

    if(LIBIGL_USE_STATIC_LIBRARY)
      include_directories(${STB_IMAGE_DIR})
      CompileIGL_Module("png")
    endif()
  endif()
endif()

### Compile the tetgen part ###
if(LIBIGL_WITH_TETGEN)
  set(TETGEN_DIR "${LIBIGL_EXTERNAL}/tetgen")
  add_subdirectory("${TETGEN_DIR}" "tetgen")
  list(APPEND LIBIGL_INCLUDE_DIRS ${TETGEN_DIR})
  list(APPEND LIBIGL_EXTRA_LIBRARIES "tetgen")

  if(LIBIGL_USE_STATIC_LIBRARY)
    include_directories(${TETGEN_DIR})
    CompileIGL_Module("tetgen")
  endif()
endif()

### Compile the triangle part ###
if(LIBIGL_WITH_TRIANGLE)
  set(TRIANGLE_DIR "${LIBIGL_EXTERNAL}/triangle")
  add_subdirectory("${TRIANGLE_DIR}" "triangle")
  list(APPEND LIBIGL_INCLUDE_DIRS ${TRIANGLE_DIR})
  list(APPEND LIBIGL_EXTRA_LIBRARIES "triangle")

  if(LIBIGL_USE_STATIC_LIBRARY)
    include_directories(${TRIANGLE_DIR})
    CompileIGL_Module("triangle")
  endif()
endif()

### Compile the viewer ###
if(LIBIGL_WITH_VIEWER) # to finish cleaning
  set(NANOGUI_DIR "${LIBIGL_EXTERNAL}/nanogui")
  set(NANOGUI_INCLUDE_DIRS
    "${NANOGUI_DIR}/include"
    "${NANOGUI_DIR}/ext/nanovg/src"
    "${NANOGUI_DIR}/ext/glfw/include")

  option(GLFW_INSTALL OFF)
  option(NANOGUI_BUILD_PYTHON OFF)
  option(NANOGUI_BUILD_EXAMPLE OFF)
  option(NANOGUI_BUILD_SHARED  OFF)
  add_subdirectory("${NANOGUI_DIR}" "nanogui")
  list(APPEND LIBIGL_INCLUDE_DIRS ${NANOGUI_INCLUDE_DIRS})
  list(APPEND LIBIGL_EXTRA_LIBRARIES "nanogui" "glfw" ${GLFW_LIBRARIES} pthread)

  if(WIN32 OR UNIX) # GLEW should not be needed on Linux
    include_directories("${NANOGUI_DIR}/ext/glew/include")
    list(APPEND LIBIGL_INCLUDE_DIRS "${NANOGUI_DIR}/ext/glew/include")
    list(APPEND LIBIGL_EXTRA_SOURCES "${NANOGUI_DIR}/ext/glew/src/glew.c")
  endif()

  if(LIBIGL_USE_STATIC_LIBRARY)
    include_directories(${NANOGUI_INCLUDE_DIRS})
    CompileIGL_Module("viewer")
  endif()
endif()

### Compile the xml part ###
if(LIBIGL_WITH_XML)
  set(TINYXML2_DIR "${LIBIGL_EXTERNAL}/tinyxml2")

  add_library(tinyxml2 STATIC ${TINYXML2_DIR}/tinyxml2.cpp ${TINYXML2_DIR}/tinyxml2.h)
  set_target_properties(tinyxml2 PROPERTIES
          COMPILE_DEFINITONS "TINYXML2_EXPORT"
          VERSION "3.0.0"
          SOVERSION "3")
  list(APPEND LIBIGL_INCLUDE_DIRS ${TINYXML2_DIR})
  list(APPEND LIBIGL_EXTRA_LIBRARIES "tinyxml2")

  if(LIBIGL_USE_STATIC_LIBRARY)
    include_directories(${TINYXML2_DIR})
    CompileIGL_Module("xml")
  endif()
endif()


### Compile the boolean and cgal parts ###
if(LIBIGL_WITH_CGAL) # to be cleaned
  set(CGAL_DIR "${LIBIGL_EXTERNAL}/cgal/build")
  #option(CGAL_BRANCH_BUILD ON)
  #option(BUILD_SHARED_LIBS OFF)
  #option(WITH_CGAL_ImageIO OFF)
  #option(WITH_CGAL_QT3 OFF)
  #option(WITH_CGAL_QT5 OFF)
  set(CMAKELISTS_SHARED_INCLUDED FALSE)
  find_package(CGAL REQUIRED)
  message(${CGAL_USE_FILE})
  include(${CGAL_USE_FILE})
  #add_subdirectory("${LIBIGL_EXTERNAL}/cgal" "cgal")

#  find_package(CGAL REQUIRED)
#  if(NOT CGAL_FOUND)
#    message(SEND_ERROR "Could NOT find CGAL")
#  endif()

  if(LIBIGL_WITH_CORK)
#    find_package(CORK REQUIRED)
#    if(NOT MOSEK_FOUND)
#      message(SEND_ERROR "Could NOT find CORK")
#    endif()

#    if(LIBIGL_USE_STATIC_LIBRARY)
#      include_directories(${CORK_INCLUDE_DIR})
#      CompileIGL_Module("boolean")
#    endif()
  else()
    add_definitions(-DIGL_NO_CORK)
    list(APPEND LIBIGL_DEFINITIONS "-DIGL_NO_CORK")
  endif()

  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("cgal")
  endif()
endif()

# Pass the list of compiled libraries to the parent if there is one
if(NOT ${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
  list(REVERSE LIBIGL_LIBRARIES)

  set(LIBIGL_INCLUDE_DIRS ${LIBIGL_INCLUDE_DIRS} PARENT_SCOPE)
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} PARENT_SCOPE)
  set(LIBIGL_EXTRA_LIBRARIES ${LIBIGL_EXTRA_LIBRARIES} PARENT_SCOPE)
  set(LIBIGL_EXTRA_SOURCES ${LIBIGL_EXTRA_SOURCES} PARENT_SCOPE)
  set(LIBIGL_DEFINITIONS ${LIBIGL_DEFINITIONS} PARENT_SCOPE)
endif()
