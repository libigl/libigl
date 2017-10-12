cmake_minimum_required(VERSION 3.1)

### Available options ###
option(LIBIGL_USE_STATIC_LIBRARY    "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR      "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL             "Use CGAL"           ON)
option(LIBIGL_WITH_COMISO           "Use CoMiso"         ON)
option(LIBIGL_WITH_CORK             "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE           "Use Embree"         OFF)
option(LIBIGL_WITH_LIM              "Use LIM"            ON)
option(LIBIGL_WITH_MATLAB           "Use Matlab"         ON)
option(LIBIGL_WITH_MOSEK            "Use MOSEK"          ON)
option(LIBIGL_WITH_NANOGUI          "Use Nanogui menu"   OFF)
option(LIBIGL_WITH_OPENGL           "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW      "Use GLFW"           ON)
option(LIBIGL_WITH_PNG              "Use PNG"            ON)
option(LIBIGL_WITH_TETGEN           "Use Tetgen"         ON)
option(LIBIGL_WITH_TRIANGLE         "Use Triangle"       ON)
option(LIBIGL_WITH_VIEWER           "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_XML              "Use XML"            ON)
option(LIBIGL_WITH_PYTHON           "Use Python"         OFF)

if(LIBIGL_WITH_VIEWER AND (NOT LIBIGL_WITH_OPENGL_GLFW OR NOT LIBIGL_WITH_OPENGL) )
  message(FATAL_ERROR "LIBIGL_WITH_VIEWER=ON requires LIBIGL_WITH_OPENGL_GLFW=ON and LIBIGL_WITH_OPENGL=ON")
endif()

################################################################################

### Configuration
set(LIBIGL_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
set(LIBIGL_SOURCE_DIR "${LIBIGL_ROOT}/include")
set(LIBIGL_EXTERNAL "${LIBIGL_ROOT}/external")

### Multiple dependencies are buried in Nanogui
set(NANOGUI_DIR "${LIBIGL_EXTERNAL}/nanogui")

# Dependencies are linked as INTERFACE targets unless libigl is compiled as a static library
if(LIBIGL_USE_STATIC_LIBRARY)
  set(IGL_SCOPE PUBLIC)
else()
  set(IGL_SCOPE INTERFACE)
endif()

################################################################################
### IGL Common
################################################################################

add_library(igl_common INTERFACE)
target_include_directories(igl_common SYSTEM INTERFACE ${LIBIGL_SOURCE_DIR})
if(LIBIGL_USE_STATIC_LIBRARY)
  target_compile_definitions(igl_common INTERFACE -DIGL_STATIC_LIBRARY)
endif()

# Transitive C++11 flags
include(CXXFeatures)
target_compile_features(igl_common INTERFACE ${CXX11_FEATURES})

# Other compilation flags
if(MSVC)
  # Enable parallel compilation for Visual Studio
  target_compile_options(igl_common INTERFACE /MP /bigobj)
endif()

if(BUILD_SHARED_LIBS)
  # Generate position independent code
  set_target_properties(igl_common PROPERTIES INTERFACE_POSITION_INDEPENDENT_CODE ON)
endif()

# Eigen
if(TARGET Eigen3::Eigen)
  # If an imported target already exists, use it
  target_link_libraries(igl_common INTERFACE Eigen3::Eigen)
else()
  target_include_directories(igl_common SYSTEM INTERFACE ${NANOGUI_DIR}/ext/eigen)
endif()

################################################################################

function(compile_igl_module module_dir prefix)
  string(REPLACE "/" "_" module_name "${module_dir}")
  if(LIBIGL_USE_STATIC_LIBRARY)
    file(GLOB SOURCES_IGL_${module_name}
      "${LIBIGL_SOURCE_DIR}/igl/${prefix}/${module_dir}/*.cpp")
    add_library(igl_${module_name} STATIC ${SOURCES_IGL_${module_name}} ${ARGN})
    if(MSVC)
      target_compile_options(igl_${module_name} PRIVATE /w) # disable all warnings (not ideal but...)
    else()
      #target_compile_options(igl_${module_name} PRIVATE -w) # disable all warnings (not ideal but...)
    endif()
  else()
    add_library(igl_${module_name} INTERFACE)
  endif()

  target_link_libraries(igl_${module_name} ${IGL_SCOPE} igl_common)
  if(NOT module_name STREQUAL "core")
	  target_link_libraries(igl_${module_name} ${IGL_SCOPE} igl_core)
  endif()

  # Alias target because it looks nicer
  message(STATUS "Creating target: igl::${module_name}")
  add_library(igl::${module_name} ALIAS igl_${module_name})
endfunction()

################################################################################
### IGL Core
################################################################################

if(LIBIGL_USE_STATIC_LIBRARY)
  file(GLOB SOURCES_IGL
    "${LIBIGL_SOURCE_DIR}/igl/*.cpp"
    "${LIBIGL_SOURCE_DIR}/igl/copyleft/*.cpp")
endif()
compile_igl_module("core" "" ${SOURCES_IGL})

################################################################################
## Compile the AntTweakBar part ###
if(LIBIGL_WITH_ANTTWEAKBAR)
  set(ANTTWEAKBAR_DIR "${LIBIGL_EXTERNAL}/AntTweakBar")
  if(NOT TARGET AntTweakBar)
    add_subdirectory("${ANTTWEAKBAR_DIR}" AntTweakBar)
  endif()
  compile_igl_module("anttweakbar" "")
  target_link_libraries(igl_anttweakbar ${IGL_SCOPE} AntTweakBar)
endif()

################################################################################
### Compile the cgal parts ###
if(LIBIGL_WITH_CGAL)
  # CGAL Core is needed for
  # `Exact_predicates_exact_constructions_kernel_with_sqrt`
  find_package(CGAL COMPONENTS Core)
  if(CGAL_FOUND)
    compile_igl_module("cgal" "copyleft/")
    find_package(Boost 1.48 REQUIRED thread system)
    target_include_directories(igl_cgal ${IGL_SCOPE} ${CGAL_INCLUDE_DIRS})
    target_link_libraries(igl_cgal ${IGL_SCOPE} CGAL::CGAL CGAL::CGAL_Core ${Boost_LIBRARIES})
  else()
    set(LIBIGL_WITH_CGAL OFF CACHE BOOL "" FORCE)
  endif()
endif()

################################################################################
# Compile CoMISo
# NOTE: this cmakefile works only with the
# comiso available here: https://github.com/libigl/CoMISo
if(LIBIGL_WITH_COMISO)
  compile_igl_module("comiso" "copyleft/")
  if(NOT TARGET CoMISo)
    add_subdirectory("${LIBIGL_EXTERNAL}/CoMISo" CoMISo)
  endif()
  target_link_libraries(igl_comiso ${IGL_SCOPE} CoMISo)
endif()

################################################################################
### Compile the cork parts ###
if(LIBIGL_WITH_CORK)
  set(CORK_DIR "${LIBIGL_EXTERNAL}/cork")
  if(NOT TARGET cork)
    # call this "lib-cork" instead of "cork", otherwise cmake gets confused about
    # "cork" executable
    add_subdirectory("${CORK_DIR}" "lib-cork")
  endif()
  compile_igl_module("cork" "copyleft/")
  target_include_directories(igl_cork ${IGL_SCOPE} cork)
  target_include_directories(igl_cork ${IGL_SCOPE} "${CORK_DIR}/src")
endif()

################################################################################
### Compile the embree part ###
if(LIBIGL_WITH_EMBREE)
  set(EMBREE_DIR "${LIBIGL_EXTERNAL}/embree")

  set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL " " FORCE)
  set(EMBREE_TASKING_SYSTEM "INTERNAL" CACHE BOOL " " FORCE)
  set(EMBREE_TUTORIALS OFF CACHE BOOL " " FORCE)
  set(EMBREE_MAX_ISA NONE CACHE STRINGS " " FORCE)

  # set(ENABLE_INSTALLER OFF CACHE BOOL " " FORCE)
  if(MSVC)
    # set(EMBREE_STATIC_RUNTIME OFF CACHE BOOL " " FORCE)
    set(EMBREE_STATIC_LIB OFF CACHE BOOL " " FORCE)
  else()
    set(EMBREE_STATIC_LIB ON CACHE BOOL " " FORCE)
  endif()

  if(NOT TARGET embree)
    add_subdirectory("${EMBREE_DIR}" "embree")
  endif()

  if(MSVC)
    add_custom_target(Copy-Embree-DLL ALL # Adds a post-build event
        COMMAND ${CMAKE_COMMAND} -E copy_if_different # which executes "cmake - E
        $<TARGET_FILE:embree> # <--this is in-file
        "${CMAKE_BINARY_DIR}" # <--this is out-file path
        DEPENDS embree) # Execute after embree target has been built
  endif()

  compile_igl_module("embree" "")
  target_link_libraries(igl_embree ${IGL_SCOPE} embree)
  target_include_directories(igl_embree ${IGL_SCOPE} ${EMBREE_DIR}/include)
  if(NOT MSVC)
    target_compile_definitions(igl_embree ${IGL_SCOPE} -DENABLE_STATIC_LIB)
  endif()
endif()

################################################################################
### Compile the lim part ###
if(LIBIGL_WITH_LIM)
  set(LIM_DIR "${LIBIGL_EXTERNAL}/lim")
  if(NOT TARGET lim)
    add_subdirectory("${LIM_DIR}" "lim")
  endif()
  compile_igl_module("lim" "")
  target_link_libraries(igl_lim ${IGL_SCOPE} lim)
  target_include_directories(igl_lim ${IGL_SCOPE} ${LIM_DIR})
endif()

################################################################################
### Compile the matlab part ###
if(LIBIGL_WITH_MATLAB)
  find_package(MATLAB)
  if(MATLAB_FOUND)
    compile_igl_module("matlab" "")
    target_link_libraries(igl_matlab ${IGL_SCOPE} ${MATLAB_LIBRARIES})
    target_include_directories(igl_matlab ${IGL_SCOPE} ${MATLAB_INCLUDE_DIR})
  else()
    set(LIBIGL_WITH_MATLAB OFF CACHE BOOL "" FORCE)
  endif()
endif()

################################################################################
### Compile the mosek part ###
if(LIBIGL_WITH_MOSEK)
  find_package(MOSEK)
  if(MOSEK_FOUND)
    compile_igl_module("mosek" "")
    target_link_libraries(igl_mosek ${IGL_SCOPE} ${MOSEK_LIBRARIES})
    target_include_directories(igl_mosek ${IGL_SCOPE} ${MOSEK_INCLUDE_DIRS})
    target_compile_definitions(igl_mosek ${IGL_SCOPE} -DLIBIGL_WITH_MOSEK)
  else()
    set(LIBIGL_WITH_MOSEK OFF CACHE BOOL "" FORCE)
  endif()
endif()

################################################################################
### Compile the opengl parts ###

if(LIBIGL_WITH_OPENGL)
  # OpenGL modules
  find_package(OpenGL REQUIRED)
  compile_igl_module("opengl" "")
  compile_igl_module("opengl2" "")
  target_link_libraries(igl_opengl ${IGL_SCOPE} ${OPENGL_gl_LIBRARY})
  target_link_libraries(igl_opengl2 ${IGL_SCOPE} ${OPENGL_gl_LIBRARY})
  target_include_directories(igl_opengl SYSTEM ${IGL_SCOPE} ${OPENGL_INCLUDE_DIR})
  target_include_directories(igl_opengl2 SYSTEM ${IGL_SCOPE} ${OPENGL_INCLUDE_DIR})

  # GLEW for linux and windows
  if(NOT TARGET glew)
    add_library(glew STATIC ${NANOGUI_DIR}/ext/glew/src/glew.c)
    target_include_directories(glew SYSTEM PUBLIC ${NANOGUI_DIR}/ext/glew/include)
    target_compile_definitions(glew PUBLIC -DGLEW_BUILD -DGLEW_NO_GLU)
  endif()
  target_link_libraries(igl_opengl ${IGL_SCOPE} glew)
  target_link_libraries(igl_opengl2 ${IGL_SCOPE} glew)

  # Nanogui
  if(LIBIGL_WITH_NANOGUI)
    if(LIBIGL_WITH_PYTHON)
      set(NANOGUI_BUILD_PYTHON ON CACHE BOOL " " FORCE)
    else()
      set(NANOGUI_BUILD_PYTHON OFF CACHE BOOL " " FORCE)
    endif()
    set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
    set(NANOGUI_BUILD_SHARED  OFF CACHE BOOL " " FORCE)
    add_subdirectory(${NANOGUI_DIR} nanogui)
    target_include_directories(nanogui PUBLIC
      "${NANOGUI_DIR}/include"
      "${NANOGUI_DIR}/ext/nanovg/src")
  endif()

  # GLFW module
  if(LIBIGL_WITH_OPENGL_GLFW)
    compile_igl_module("opengl/glfw" "")
    if(NOT TARGET glfw)
      set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
      set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
      set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
      set(GLFW_BUILD_INSTALL OFF CACHE BOOL " " FORCE)
      add_subdirectory(${NANOGUI_DIR}/ext/glfw glfw)
    endif()
    target_include_directories(glfw ${IGL_SCOPE} ${NANOGUI_DIR}/ext/glfw/include)
    target_link_libraries(igl_opengl_glfw ${IGL_SCOPE} igl_opengl glfw)
  endif()

  # Viewer module
  if(LIBIGL_WITH_VIEWER)
    compile_igl_module("viewer" "")
    target_link_libraries(igl_viewer ${IGL_SCOPE} igl_core glfw glew ${OPENGL_gl_LIBRARY})
    target_include_directories(igl_viewer SYSTEM ${IGL_SCOPE} ${OPENGL_INCLUDE_DIR})
    if(TARGET nanogui)
      target_link_libraries(igl_viewer ${IGL_SCOPE} nanogui)
      target_compile_definitions(igl_viewer ${IGL_SCOPE} -DIGL_VIEWER_WITH_NANOGUI)
    endif()
  endif()

endif()

################################################################################
### Compile the png parts ###
if(LIBIGL_WITH_PNG)
  if(TARGET igl_opengl)
    set(STB_IMAGE_DIR "${LIBIGL_EXTERNAL}/stb_image")
    if(NOT TARGET stb_image)
      add_subdirectory("${STB_IMAGE_DIR}" "stb_image")
    endif()
    compile_igl_module("png" "")
    target_link_libraries(igl_png ${IGL_SCOPE} igl_stb_image igl_opengl)
  else()
    set(LIBIGL_WITH_PNG OFF CACHE BOOL "" FORCE)
  endif()
endif()

################################################################################
### Compile the tetgen part ###
if(LIBIGL_WITH_TETGEN)
  set(TETGEN_DIR "${LIBIGL_EXTERNAL}/tetgen")
  if(NOT TARGET tetgen)
    add_subdirectory("${TETGEN_DIR}" "tetgen")
  endif()
  compile_igl_module("tetgen" "copyleft/")
  target_link_libraries(igl_tetgen ${IGL_SCOPE} tetgen)
  target_include_directories(igl_tetgen ${IGL_SCOPE} ${TETGEN_DIR})
endif()

################################################################################
### Compile the triangle part ###
if(LIBIGL_WITH_TRIANGLE)
  set(TRIANGLE_DIR "${LIBIGL_EXTERNAL}/triangle")
  if(NOT TARGET triangle)
    add_subdirectory("${TRIANGLE_DIR}" "triangle")
  endif()
  compile_igl_module("triangle" "")
  target_link_libraries(igl_triangle ${IGL_SCOPE} triangle)
  target_include_directories(igl_triangle ${IGL_SCOPE} ${TRIANGLE_DIR})
endif()

################################################################################
### Compile the xml part ###
if(LIBIGL_WITH_XML)
  set(TINYXML2_DIR "${LIBIGL_EXTERNAL}/tinyxml2")
  if(NOT TARGET tinyxml2)
    add_library(tinyxml2 STATIC ${TINYXML2_DIR}/tinyxml2.cpp ${TINYXML2_DIR}/tinyxml2.h)
    target_include_directories(tinyxml2 PUBLIC ${TINYXML2_DIR})
    set_target_properties(tinyxml2 PROPERTIES
            COMPILE_DEFINITIONS "TINYXML2_EXPORT"
            VERSION "3.0.0"
            SOVERSION "3")
  endif()
  compile_igl_module("xml" "")
  target_link_libraries(igl_xml ${IGL_SCOPE} tinyxml2)
endif()
