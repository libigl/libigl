cmake_minimum_required(VERSION 3.1)

### Available options ###
option(LIBIGL_USE_STATIC_LIBRARY    "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR      "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL             "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO           "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK             "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE           "Use Embree"         OFF)
option(LIBIGL_WITH_LIM              "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB           "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK            "Use MOSEK"          OFF)
option(LIBIGL_WITH_NANOGUI          "Use Nanogui menu"   OFF)
option(LIBIGL_WITH_OPENGL           "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW      "Use GLFW"           OFF)
option(LIBIGL_WITH_PNG              "Use PNG"            OFF)
option(LIBIGL_WITH_TETGEN           "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE         "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER           "Use OpenGL viewer"  OFF)
option(LIBIGL_WITH_XML              "Use XML"            OFF)
option(LIBIGL_WITH_PYTHON           "Use Python"         OFF)

if(LIBIGL_WITH_VIEWER AND (NOT LIBIGL_WITH_OPENGL_GLFW OR NOT LIBIGL_WITH_OPENGL) )
  message(FATAL_ERROR "LIBIGL_WITH_VIEWER=ON requires LIBIGL_WITH_OPENGL_GLFW=ON and LIBIGL_WITH_OPENGL=ON")
endif()

### Compilation configuration ###
function(set_compilation_flags target_name)
  if(MSVC)
    ### Enable parallel compilation for Visual Studio
    target_compile_options(${target_name} PRIVATE /MP /bigobj)
    #target_compile_options(${target_name} PRIVATE /w) # disable all warnings (not ideal but...)
  else()
    #target_compile_options(${target_name} PRIVATE -w) # disable all warnings (not ideal but...)
  endif()
endfunction()

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

if(BUILD_SHARED_LIBS)
  # Generate position independent code
  set_target_properties(igl_${module_name} PROPERTIES INTERFACE_POSITION_INDEPENDENT_CODE ON)
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
    target_link_libraries(igl_${module_name} PUBLIC igl_common)
    set_compilation_flags(igl_${module_name})

  else()
    add_library(igl_${module_name} INTERFACE)
    target_link_libraries(igl_${module_name} INTERFACE igl_common)
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
### Compile the AntTweakBar part ###
# if(LIBIGL_WITH_ANTTWEAKBAR)
#   set(ANTTWEAKBAR_DIR "${LIBIGL_EXTERNAL}/AntTweakBar")
#   set(ANTTWEAKBAR_INCLUDE_DIR "${ANTTWEAKBAR_DIR}/include")
#   set(ANTTWEAKBAR_C_SRC_FILES
#     "${ANTTWEAKBAR_DIR}/src/TwEventGLFW.c"
#     "${ANTTWEAKBAR_DIR}/src/TwEventGLUT.c"
#     "${ANTTWEAKBAR_DIR}/src/TwEventSDL.c"
#     "${ANTTWEAKBAR_DIR}/src/TwEventSDL12.c"
#     "${ANTTWEAKBAR_DIR}/src/TwEventSDL13.c"
#     )
#   set(ANTTWEAKBAR_CPP_SRC_FILES
#     "${ANTTWEAKBAR_DIR}/src/LoadOGL.cpp"
#     "${ANTTWEAKBAR_DIR}/src/LoadOGLCore.cpp"
#     "${ANTTWEAKBAR_DIR}/src/TwBar.cpp"
#     "${ANTTWEAKBAR_DIR}/src/TwColors.cpp"
#     "${ANTTWEAKBAR_DIR}/src/TwEventSFML.cpp"
#     "${ANTTWEAKBAR_DIR}/src/TwFonts.cpp"
#     "${ANTTWEAKBAR_DIR}/src/TwMgr.cpp"
#     "${ANTTWEAKBAR_DIR}/src/TwOpenGL.cpp"
#     "${ANTTWEAKBAR_DIR}/src/TwOpenGLCore.cpp"
#     "${ANTTWEAKBAR_DIR}/src/TwPrecomp.cpp"
#     )
#     # These are probably needed for windows/Linux, should append if
#     # windows/Linux
#     #"${ANTTWEAKBAR_DIR}/src/TwEventWin.c"
#     #"${ANTTWEAKBAR_DIR}/src/TwEventX11.c"
#     #"${ANTTWEAKBAR_DIR}/src/TwDirect3D10.cpp"
#     #"${ANTTWEAKBAR_DIR}/src/TwDirect3D11.cpp"
#     #"${ANTTWEAKBAR_DIR}/src/TwDirect3D9.cpp"
#   list(
#     APPEND
#     ANTTWEAKBAR_SRC_FILES
#     "${ANTTWEAKBAR_C_SRC_FILES}"
#     "${ANTTWEAKBAR_CPP_SRC_FILES}")
#   add_library(AntTweakBar STATIC "${ANTTWEAKBAR_SRC_FILES}")
#   target_include_directories(AntTweakBar PUBLIC "${ANTTWEAKBAR_INCLUDE_DIR}")
#   if(APPLE)
#     set_target_properties(
#       AntTweakBar
#       PROPERTIES
#       COMPILE_FLAGS
#       "-fPIC -fno-strict-aliasing -x objective-c++")
#     target_compile_definitions(
#       AntTweakBar PUBLIC _MACOSX __PLACEMENT_NEW_INLINE)
#   endif()
#   list(APPEND LIBIGL_INCLUDE_DIRS "${ANTTWEAKBAR_INCLUDE_DIR}")
#   set(LIBIGL_ANTTWEAKBAR_EXTRA_LIBRARIES "AntTweakBar")
#   list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_ANTTWEAKBAR_EXTRA_LIBRARIES})
#   if(LIBIGL_USE_STATIC_LIBRARY)
#     compile_igl_module("anttweakbar" "")
#     target_include_directories(igl_anttweakbar PRIVATE ${ANTTWEAKBAR_INCLUDE_DIR})
#   endif()
# endif()

################################################################################
### Compile the cgal parts ###
if(LIBIGL_WITH_CGAL)
  # CGAL Core is needed for
  # `Exact_predicates_exact_constructions_kernel_with_sqrt`
  compile_igl_module("cgal" "copyleft/")
  find_package(CGAL REQUIRED COMPONENTS Core)
  find_package(Boost 1.48 REQUIRED thread system)
  target_link_libraries(igl_cgal ${IGL_SCOPE} CGAL::CGAL ${Boost_LIBRARIES})
endif()

################################################################################
# Compile CoMISo
# NOTE: this cmakefile works only with the
# comiso available here: https://github.com/libigl/CoMISo
if(LIBIGL_WITH_COMISO)
  compile_igl_module("comiso" "copyleft/")
  add_subdirectory(${LIBIGL_EXTERNAL}/CoMISo CoMISo)
  target_link_libraries(igl_comiso ${IGL_SCOPE} CoMISo)
endif()

################################################################################
### Compile the cork parts ###
# if(LIBIGL_WITH_CORK)
#   set(CORK_DIR "${LIBIGL_EXTERNAL}/cork")
#   set(CORK_INCLUDE_DIR "${CORK_DIR}/src")
#   # call this "lib-cork" instead of "cork", otherwise cmake gets confused about
#   # "cork" executable
#   add_subdirectory("${CORK_DIR}" "lib-cork")
#   list(APPEND LIBIGL_INCLUDE_DIRS "${CORK_INCLUDE_DIR}")
#   list(APPEND LIBIGL_CORK_EXTRA_LIBRARIES "cork")
#   list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_CORK_EXTRA_LIBRARIES})
#   if(LIBIGL_USE_STATIC_LIBRARY)
#     compile_igl_module("cork" "copyleft/")
#     target_include_directories(igl_cork PRIVATE ${CORK_INCLUDE_DIR})
#   endif()
# endif()

################################################################################
### Compile the embree part ###
# if(LIBIGL_WITH_EMBREE)
#   set(EMBREE_DIR "${LIBIGL_EXTERNAL}/embree")

#   set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL " " FORCE)
#   set(EMBREE_TASKING_SYSTEM "INTERNAL" CACHE BOOL " " FORCE)
#   set(EMBREE_TUTORIALS OFF CACHE BOOL " " FORCE)
#   set(EMBREE_MAX_ISA NONE CACHE STRINGS " " FORCE)

#   # set(ENABLE_INSTALLER OFF CACHE BOOL " " FORCE)
#   if(MSVC)
#     # set(EMBREE_STATIC_RUNTIME OFF CACHE BOOL " " FORCE)
#     set(EMBREE_STATIC_LIB OFF CACHE BOOL " " FORCE)
#   else()
#     set(EMBREE_STATIC_LIB ON CACHE BOOL " " FORCE)
#   endif()

#   add_subdirectory("${EMBREE_DIR}" "embree")
#   list(APPEND LIBIGL_INCLUDE_DIRS "${EMBREE_DIR}/include")
#   list(APPEND LIBIGL_EMBREE_EXTRA_LIBRARIES "embree")
#   list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_EMBREE_EXTRA_LIBRARIES})

#   if(NOT MSVC)
#     list(APPEND LIBIGL_DEFINITIONS "-DENABLE_STATIC_LIB")
#   endif()

#   if(MSVC)
#     add_custom_target(Copy-Embree-DLL ALL        # Adds a post-build event to MyTest
#         COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
#             "${CMAKE_BINARY_DIR}/libigl/embree/$<CONFIGURATION>/embree.dll"      # <--this is in-file
#           "${CMAKE_BINARY_DIR}/embree.dll")                 # <--this is out-file path  endif()
#   endif()

#   if(LIBIGL_USE_STATIC_LIBRARY)
#     compile_igl_module("embree" "")
#     target_include_directories(igl_embree PRIVATE ${EMBREE_DIR}/include)
#   if(NOT MSVC)
#     target_compile_definitions(igl_embree PRIVATE -DENABLE_STATIC_LIB)
#   endif()
#   endif()
# endif()

################################################################################
### Compile the lim part ###
# if(LIBIGL_WITH_LIM)
#   set(LIM_DIR "${LIBIGL_EXTERNAL}/lim")
#   add_subdirectory("${LIM_DIR}" "lim")
#   list(APPEND LIBIGL_INCLUDE_DIRS ${LIM_DIR})
#   ## it depends on ligigl, so placing it here solve linking problems
#   #list(APPEND LIBIGL_LIBRARIES "lim")
#   # ^--- Alec: I don't understand this comment. Does lim need to come before
#   # libigl libraries? Why can't lim be placed where it belongs in
#   # LIBIGL_EXTRA_LIBRARIES?
#   set(LIBIGL_LIM_EXTRA_LIBRARIES "lim")
#   list(APPEND LIBIGL_EXTRA_LIBRARIES "${LIBIGL_LIM_EXTRA_LIBRARIES}")

#   if(LIBIGL_USE_STATIC_LIBRARY)
#     compile_igl_module("lim" "")
#     target_include_directories(igl_lim PRIVATE ${LIM_DIR})
#   endif()
# endif()

################################################################################
### Compile the matlab part ###
# if(LIBIGL_WITH_MATLAB)
#   find_package(MATLAB REQUIRED)
#   list(APPEND LIBIGL_INCLUDE_DIRS ${MATLAB_INCLUDE_DIR})
#   list(APPEND LIBIGL_MATLAB_EXTRA_LIBRARIES ${MATLAB_LIBRARIES})
#   list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_MATLAB_EXTRA_LIBRARIES})

#   if(LIBIGL_USE_STATIC_LIBRARY)
#     compile_igl_module("matlab" "")
#     target_include_directories(igl_matlab PRIVATE ${MATLAB_INCLUDE_DIR})
#   endif()
# endif()

################################################################################
### Compile the mosek part ###
# if(LIBIGL_WITH_MOSEK)
#   find_package(MOSEK REQUIRED)
#   list(APPEND LIBIGL_INCLUDE_DIRS ${MOSEK_INCLUDE_DIR})
#   list(APPEND LIBIGL_MOSEK_EXTRA_LIBRARIES ${MOSEK_LIBRARIES})
#   list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_MOSEK_EXTRA_LIBRARIES})
#   if(LIBIGL_USE_STATIC_LIBRARY)
#     compile_igl_module("mosek" "")
#     target_include_directories(igl_mosek PRIVATE ${MOSEK_INCLUDE_DIR})
#   endif()
# else()
#   list(APPEND LIBIGL_DEFINITIONS "-DIGL_NO_MOSEK")
# endif()

################################################################################
### Compile the opengl parts ###

if(LIBIGL_WITH_OPENGL)
  compile_igl_module("opengl" "")
  compile_igl_module("opengl2" "")

  find_package(OpenGL REQUIRED)
  target_link_libraries(igl_opengl ${IGL_SCOPE} ${OPENGL_gl_LIBRARY})
  target_link_libraries(igl_opengl2 ${IGL_SCOPE} ${OPENGL_gl_LIBRARY})
  target_include_directories(igl_opengl SYSTEM ${IGL_SCOPE} ${OPENGL_INCLUDE_DIR})
  target_include_directories(igl_opengl2 SYSTEM ${IGL_SCOPE} ${OPENGL_INCLUDE_DIR})

  ### GLEW for linux and windows
  if(NOT TARGET glew)
    add_library(glew STATIC ${NANOGUI_DIR}/ext/glew/src/glew.c)
    target_include_directories(glew SYSTEM PUBLIC ${NANOGUI_DIR}/ext/glew/include)
    target_compile_definitions(glew PUBLIC -DGLEW_BUILD -DGLEW_NO_GLU)
  endif()
  target_link_libraries(igl_opengl ${IGL_SCOPE} glew)
  target_link_libraries(igl_opengl2 ${IGL_SCOPE} glew)

  if(LIBIGL_WITH_OPENGL_GLFW)
    # GLFW
    compile_igl_module("opengl/glfw" "")
    if(NOT TARGET glfw)
      set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
      set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
      set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
      set(GLFW_BUILD_INSTALL OFF CACHE BOOL " " FORCE)
      add_subdirectory(${NANOGUI_DIR}/ext/glfw glfw)
      target_include_directories(glfw ${IGL_SCOPE} ${NANOGUI_DIR}/ext/glfw/include)
    endif()
    target_link_libraries(igl_opengl_glfw ${IGL_SCOPE} igl_opengl glfw)

    ### Compile the viewer
    if(LIBIGL_WITH_VIEWER)
      compile_igl_module("viewer" "")
      target_link_libraries(igl_viewer ${IGL_SCOPE} igl_core glfw glew ${OPENGL_gl_LIBRARY})
      target_include_directories(igl_opengl2 SYSTEM ${IGL_SCOPE} ${OPENGL_INCLUDE_DIR})

      if(LIBIGL_WITH_NANOGUI)
        target_compile_definitions(igl_viewer PUBLIC -DIGL_VIEWER_WITH_NANOGUI)
        if(LIBIGL_WITH_PYTHON)
          set(NANOGUI_BUILD_PYTHON ON CACHE BOOL " " FORCE)
        else()
          set(NANOGUI_BUILD_PYTHON OFF CACHE BOOL " " FORCE)
        endif()
        set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
        set(NANOGUI_BUILD_SHARED  OFF CACHE BOOL " " FORCE)
        add_subdirectory(${NANOGUI_DIR} nanogui)
        target_link_libraries(igl_viewer ${IGL_SCOPE} nanogui)
      endif()
    endif()

  endif()

endif()

################################################################################
### Compile the png parts ###
# if(LIBIGL_WITH_PNG)
#   if(LIBIGL_WITH_NANOGUI)
#     set(STBI_LOAD OFF CACHE BOOL " " FORCE)
#   endif()
#   set(STB_IMAGE_DIR "${LIBIGL_EXTERNAL}/stb_image")
#   add_subdirectory("${STB_IMAGE_DIR}" "stb_image")
#   list(APPEND LIBIGL_INCLUDE_DIRS ${STB_IMAGE_DIR})
#   list(APPEND LIBIGL_PNG_EXTRA_LIBRARIES "stb_image")
#   list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_PNG_EXTRA_LIBRARIES})
#   if(LIBIGL_USE_STATIC_LIBRARY)
#     compile_igl_module("png" "")
#     target_include_directories(igl_png PRIVATE ${STB_IMAGE_DIR})
#     if(NOT APPLE)
#       target_include_directories(igl_png PRIVATE "${NANOGUI_DIR}/ext/glew/include")
#     endif()
#   endif()
# endif()

################################################################################
### Compile the tetgen part ###
# if(LIBIGL_WITH_TETGEN)
#   set(TETGEN_DIR "${LIBIGL_EXTERNAL}/tetgen")
#   add_subdirectory("${TETGEN_DIR}" "tetgen")
#   list(APPEND LIBIGL_INCLUDE_DIRS ${TETGEN_DIR})
#   list(APPEND LIBIGL_TETGEN_EXTRA_LIBRARIES "tetgen")
#   list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_TETGEN_EXTRA_LIBRARIES})

#   if(LIBIGL_USE_STATIC_LIBRARY)
#     compile_igl_module("tetgen" "copyleft/")
#     target_include_directories(igl_tetgen PRIVATE ${TETGEN_DIR})
#   endif()
# endif()

################################################################################
### Compile the triangle part ###
# if(LIBIGL_WITH_TRIANGLE)
#   set(TRIANGLE_DIR "${LIBIGL_EXTERNAL}/triangle")
#   add_subdirectory("${TRIANGLE_DIR}" "triangle")
#   list(APPEND LIBIGL_INCLUDE_DIRS ${TRIANGLE_DIR})
#   list(APPEND LIBIGL_TRIANGLE_EXTRA_LIBRARIES "triangle")
#   list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_TRIANGLE_EXTRA_LIBRARIES})

#   if(LIBIGL_USE_STATIC_LIBRARY)
#     compile_igl_module("triangle" "")
#     target_include_directories(igl_triangle PRIVATE ${TRIANGLE_DIR})
#   endif()
# endif()

################################################################################
### Compile the xml part ###
# if(LIBIGL_WITH_XML)
#   set(TINYXML2_DIR "${LIBIGL_EXTERNAL}/tinyxml2")
#   add_library(tinyxml2 STATIC ${TINYXML2_DIR}/tinyxml2.cpp ${TINYXML2_DIR}/tinyxml2.h)
#   set_target_properties(tinyxml2 PROPERTIES
#           COMPILE_DEFINITIONS "TINYXML2_EXPORT"
#           VERSION "3.0.0"
#           SOVERSION "3")
#   list(APPEND LIBIGL_INCLUDE_DIRS ${TINYXML2_DIR})
#   list(APPEND LIBIGL_XML_EXTRA_LIBRARIES "tinyxml2")
#   list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_XML_EXTRA_LIBRARIES})
#   if(LIBIGL_USE_STATIC_LIBRARY)
#     compile_igl_module("xml" "")
#     target_include_directories(igl_xml PRIVATE ${TINYXML2_DIR})
#   endif()
# endif()

# Function to print list nicely
# function(print_list title list)
#   message("-- ${title}:")
#   foreach(elt ${list})
#     message("\t ${elt}")
#   endforeach()
# endfunction()
