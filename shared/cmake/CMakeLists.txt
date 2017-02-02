cmake_minimum_required(VERSION 2.8.12)
project(libigl)

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
if(MSVC)
  ### Enable parallel compilation for Visual Studio
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /w") # disable all warnings (not ideal but...)
else()
  #### Libigl requires a modern C++ compiler that supports c++11
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w") # disable all warnings (not ideal but...)
endif()
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

### OpenMP ### (OpenMP is disable for now)
#find_package(OpenMP)
#if(OPENMP_FOUND AND NOT WIN32)
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#  list(APPEND LIBIGL_DEFINITIONS "${OpenMP_CXX_FLAGS}")
#endif()


### Compiling libraries based on chosen options ###
set(LIBIGL_INCLUDE_DIRS "")
set(LIBIGL_LIBRARIES "")
set(LIBIGL_EXTRA_LIBRARIES "")
set(LIBIGL_DEFINITIONS "")

if(LIBIGL_USE_STATIC_LIBRARY)
  list(APPEND LIBIGL_DEFINITIONS "-DIGL_STATIC_LIBRARY")
endif()

### macro definition ###
set(LIBIGL_ROOT "${PROJECT_SOURCE_DIR}/../..")
set(LIBIGL_SOURCE_DIR "${LIBIGL_ROOT}/include")
set(LIBIGL_EXTERNAL "${LIBIGL_ROOT}/external")

### Multiple dependencies are buried in Nanogui
set(NANOGUI_DIR "${LIBIGL_EXTERNAL}/nanogui")

### Eigen ###
set(EIGEN_INCLUDE_DIR "${NANOGUI_DIR}/ext/eigen")
list(APPEND LIBIGL_INCLUDE_DIRS "${EIGEN_INCLUDE_DIR}")

macro(CompileIGL_Module module_dir prefix)
  string(REPLACE "/" "_" module_name "${module_dir}")
  file(GLOB SOURCES_IGL_${module_name}
  "${LIBIGL_SOURCE_DIR}/igl/${prefix}/${module_dir}/*.cpp"
  )
  set(module_lib_name igl_${module_name})
  add_library(${module_lib_name} STATIC ${SOURCES_IGL_${module_name}})
  target_include_directories(${module_lib_name} PRIVATE ${EIGEN_INCLUDE_DIR})
  target_include_directories(${module_lib_name} PRIVATE ${LIBIGL_SOURCE_DIR})
  target_compile_definitions(${module_lib_name} PRIVATE -DIGL_STATIC_LIBRARY)
  list(APPEND LIBIGL_LIBRARIES "${module_lib_name}")
endmacro()

### IGL Common ###
list(APPEND LIBIGL_INCLUDE_DIRS "${LIBIGL_SOURCE_DIR}")
if(LIBIGL_USE_STATIC_LIBRARY)
  file(GLOB SOURCES_IGL
  "${LIBIGL_SOURCE_DIR}/igl/*.cpp"
  "${LIBIGL_SOURCE_DIR}/igl/copyleft/*.cpp")
  add_library(igl STATIC ${SOURCES_IGL})
  target_include_directories(igl PRIVATE ${EIGEN_INCLUDE_DIR})
  target_include_directories(igl PRIVATE ${LIBIGL_SOURCE_DIR})
  target_compile_definitions(igl PRIVATE -DIGL_STATIC_LIBRARY)
  list(APPEND LIBIGL_LIBRARIES "igl")
endif()


### Compile the AntTweakBar part ###
if(LIBIGL_WITH_ANTTWEAKBAR)
  set(ANTTWEAKBAR_DIR "${LIBIGL_EXTERNAL}/anttweakbar")
  set(ANTTWEAKBAR_INCLUDE_DIR "${ANTTWEAKBAR_DIR}/include")
  set(ANTTWEAKBAR_C_SRC_FILES
    "${ANTTWEAKBAR_DIR}/src/TwEventGLFW.c"
    "${ANTTWEAKBAR_DIR}/src/TwEventGLUT.c"
    "${ANTTWEAKBAR_DIR}/src/TwEventSDL.c"
    "${ANTTWEAKBAR_DIR}/src/TwEventSDL12.c"
    "${ANTTWEAKBAR_DIR}/src/TwEventSDL13.c"
    )
  set(ANTTWEAKBAR_CPP_SRC_FILES
    "${ANTTWEAKBAR_DIR}/src/LoadOGL.cpp"
    "${ANTTWEAKBAR_DIR}/src/LoadOGLCore.cpp"
    "${ANTTWEAKBAR_DIR}/src/TwBar.cpp"
    "${ANTTWEAKBAR_DIR}/src/TwColors.cpp"
    "${ANTTWEAKBAR_DIR}/src/TwEventSFML.cpp"
    "${ANTTWEAKBAR_DIR}/src/TwFonts.cpp"
    "${ANTTWEAKBAR_DIR}/src/TwMgr.cpp"
    "${ANTTWEAKBAR_DIR}/src/TwOpenGL.cpp"
    "${ANTTWEAKBAR_DIR}/src/TwOpenGLCore.cpp"
    "${ANTTWEAKBAR_DIR}/src/TwPrecomp.cpp"
    )
    # These are probably needed for windows/Linux, should append if
    # windows/Linux
    #"${ANTTWEAKBAR_DIR}/src/TwEventWin.c"
    #"${ANTTWEAKBAR_DIR}/src/TwEventX11.c"
    #"${ANTTWEAKBAR_DIR}/src/TwDirect3D10.cpp"
    #"${ANTTWEAKBAR_DIR}/src/TwDirect3D11.cpp"
    #"${ANTTWEAKBAR_DIR}/src/TwDirect3D9.cpp"
  list(
    APPEND
    ANTTWEAKBAR_SRC_FILES
    "${ANTTWEAKBAR_C_SRC_FILES}"
    "${ANTTWEAKBAR_CPP_SRC_FILES}")
  add_library(AntTweakBar STATIC "${ANTTWEAKBAR_SRC_FILES}")
  target_include_directories(AntTweakBar PUBLIC "${ANTTWEAKBAR_INCLUDE_DIR}")
  if(APPLE)
    set_target_properties(
      AntTweakBar
      PROPERTIES
      COMPILE_FLAGS
      "-fPIC -fno-strict-aliasing -x objective-c++")
    target_compile_definitions(
      AntTweakBar PUBLIC _MACOSX __PLACEMENT_NEW_INLINE)
  endif()
  list(APPEND LIBIGL_INCLUDE_DIRS "${ANTTWEAKBAR_INCLUDE_DIR}")
  set(LIBIGL_ANTTWEAKBAR_EXTRA_LIBRARIES "AntTweakBar")
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_ANTTWEAKBAR_EXTRA_LIBRARIES})
  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("anttweakbar" "")
    target_include_directories(igl_anttweakbar PRIVATE ${ANTTWEAKBAR_INCLUDE_DIR})
  endif()
endif()

### Compile the cgal parts ###
if(LIBIGL_WITH_CGAL) # to be cleaned
  find_package(CGAL REQUIRED)
  # set(Boost_USE_MULTITHREADED      ON)
  # set(Boost_USE_STATIC_LIBS      ON)
  #
  # find_package(BOOST REQUIRED)
  include(${CGAL_USE_FILE})
  list(APPEND LIBIGL_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS})
  list(APPEND LIBIGL_INCLUDE_DIRS ${CGAL_INCLUDE_DIRS})
  list(APPEND LIBIGL_CGAL_EXTRA_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES})
  list(APPEND LIBIGL_CGAL_EXTRA_LIBRARIES ${CGAL_LIBRARIES})
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_CGAL_EXTRA_LIBRARIES})
  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("cgal" "copyleft/")
    target_include_directories(igl_cgal PRIVATE
      ${CGAL_3RD_PARTY_INCLUDE_DIRS}
      ${CGAL_INCLUDE_DIRS})
  endif()
endif()

#Compile CoMISo
# NOTE: this cmakefile works only with the
# comiso available here: https://github.com/libigl/CoMISo
if(LIBIGL_WITH_COMISO)
  set(COMISO_DIR "${LIBIGL_EXTERNAL}/CoMISo")
  set(COMISO_INCLUDE_DIRS
    "${COMISO_DIR}/ext/gmm-4.2/include"
    "${COMISO_DIR}/../")
  list(APPEND LIBIGL_INCLUDE_DIRS ${COMISO_INCLUDE_DIRS})
  #add_definitions(-DINCLUDE_TEMPLATES) (what need this?)
  list(APPEND LIBIGL_DEFINITIONS "-DINCLUDE_TEMPLATES")
  if(APPLE)
    find_library(accelerate_library Accelerate)
    list(APPEND LIBIGL_COMISO_EXTRA_LIBRARIES "CoMISo" ${accelerate_library})
  elseif(UNIX)
    find_package(BLAS REQUIRED)
    list(APPEND LIBIGL_COMISO_EXTRA_LIBRARIES "CoMISo" ${BLAS_LIBRARIES})
  endif(APPLE)
  if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_SCL_SECURE_NO_DEPRECATE")
    #link_directories("${COMISO_ROOT}/CoMISo/ext/OpenBLAS-v0.2.14-Win64-int64/lib/")
    list(APPEND LIBIGL_COMISO_EXTRA_LIBRARIES "CoMISo" "${COMISO_DIR}/ext/OpenBLAS-v0.2.14-Win64-int64/lib/libopenblas.dll.a.lib")
  endif()
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_COMISO_EXTRA_LIBRARIES})
  add_subdirectory("${COMISO_DIR}" "CoMISo")
  if(MSVC)
    # Copy the dll
    add_custom_target(Copy-CoMISo-DLL ALL       # Adds a post-build event to MyTest
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${COMISO_DIR}/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libopenblas.dll"
        "${CMAKE_CURRENT_BINARY_DIR}/../libopenblas.dll"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${COMISO_DIR}/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libgcc_s_seh-1.dll"
        "${CMAKE_CURRENT_BINARY_DIR}/../libgcc_s_seh-1.dll"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${COMISO_DIR}/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libgfortran-3.dll"
        "${CMAKE_CURRENT_BINARY_DIR}/../libgfortran-3.dll"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${COMISO_DIR}/ext/OpenBLAS-v0.2.14-Win64-int64/bin/libquadmath-0.dll"
        "${CMAKE_CURRENT_BINARY_DIR}/../libquadmath-0.dll")
  endif()
  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("comiso" "copyleft/")
    target_include_directories(igl_comiso PRIVATE ${COMISO_INCLUDE_DIRS})
    target_compile_definitions(igl_comiso PRIVATE -DINCLUDE_TEMPLATES)
  endif()
endif()

### Compile the cork parts ###
if(LIBIGL_WITH_CORK)
  set(CORK_DIR "${LIBIGL_EXTERNAL}/cork")
  set(CORK_INCLUDE_DIR "${CORK_DIR}/src")
  # call this "lib-cork" instead of "cork", otherwise cmake gets confused about
  # "cork" executable
  add_subdirectory("${CORK_DIR}" "lib-cork")
  list(APPEND LIBIGL_INCLUDE_DIRS "${CORK_INCLUDE_DIR}")
  list(APPEND LIBIGL_CORK_EXTRA_LIBRARIES "cork")
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_CORK_EXTRA_LIBRARIES})
  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("cork" "copyleft/")
    target_include_directories(igl_cork PRIVATE ${CORK_INCLUDE_DIR})
  endif()
endif()

### Compile the embree part ###
if(LIBIGL_WITH_EMBREE)
  set(EMBREE_DIR "${LIBIGL_EXTERNAL}/embree")

  set(ENABLE_ISPC_SUPPORT OFF CACHE BOOL " " FORCE)
  set(RTCORE_TASKING_SYSTEM "INTERNAL" CACHE BOOL " " FORCE)
  set(ENABLE_TUTORIALS OFF CACHE BOOL " " FORCE)
  set(ENABLE_STATIC_LIB ON CACHE BOOL " " FORCE)
  set(ENABLE_INSTALLER OFF CACHE BOOL " " FORCE)
  if(MSVC)
    set(USE_STATIC_RUNTIME OFF CACHE BOOL " " FORCE)
  endif()

  add_subdirectory("${EMBREE_DIR}" "embree")
  list(APPEND LIBIGL_INCLUDE_DIRS "${EMBREE_DIR}/include")
  list(APPEND LIBIGL_EMBREE_EXTRA_LIBRARIES "embree")
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_EMBREE_EXTRA_LIBRARIES})
  list(APPEND LIBIGL_DEFINITIONS "-DENABLE_STATIC_LIB")

  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("embree" "")
    target_include_directories(igl_embree PRIVATE ${EMBREE_DIR}/include)
    target_compile_definitions(igl_embree PRIVATE -DENABLE_STATIC_LIB)
  endif()
endif()

### Compile the lim part ###
if(LIBIGL_WITH_LIM)
  set(LIM_DIR "${LIBIGL_EXTERNAL}/lim")
  add_subdirectory("${LIM_DIR}" "lim")
  list(APPEND LIBIGL_INCLUDE_DIRS ${LIM_DIR})
  ## it depends on ligigl, so placing it here solve linking problems
  #list(APPEND LIBIGL_LIBRARIES "lim")
  # ^--- Alec: I don't understand this comment. Does lim need to come before
  # libigl libraries? Why can't lim be placed where it belongs in
  # LIBIGL_EXTRA_LIBRARIES?
  set(LIBIGL_LIM_EXTRA_LIBRARIES "lim")
  list(APPEND LIBIGL_EXTRA_LIBRARIES "${LIBIGL_LIM_EXTRA_LIBRARIES}")

  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("lim" "")
    target_include_directories(igl_lim PRIVATE ${LIM_DIR})
  endif()
endif()

### Compile the matlab part ###
if(LIBIGL_WITH_MATLAB)
  find_package(MATLAB REQUIRED)
  list(APPEND LIBIGL_INCLUDE_DIRS ${MATLAB_INCLUDE_DIR})
  list(APPEND LIBIGL_MATLAB_EXTRA_LIBRARIES ${MATLAB_LIBRARIES})
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_MATLAB_EXTRA_LIBRARIES})

  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("matlab" "")
    target_include_directories(igl_matlab PRIVATE ${MATLAB_INCLUDE_DIR})
  endif()
endif()

### Compile the mosek part ###
if(LIBIGL_WITH_MOSEK)
  find_package(MOSEK REQUIRED)
  list(APPEND LIBIGL_INCLUDE_DIRS ${MOSEK_INCLUDE_DIR})
  list(APPEND LIBIGL_MOSEK_EXTRA_LIBRARIES ${MOSEK_LIBRARIES})
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_MOSEK_EXTRA_LIBRARIES})
  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("mosek" "")
    target_include_directories(igl_mosek PRIVATE ${MOSEK_INCLUDE_DIR})
  endif()
else()
  list(APPEND LIBIGL_DEFINITIONS "-DIGL_NO_MOSEK")
endif()

### Compile the opengl parts ###
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_INSTALL OFF CACHE BOOL " " FORCE)
if(LIBIGL_WITH_OPENGL)
  find_package(OpenGL REQUIRED)

  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("opengl" "")
    CompileIGL_Module("opengl2" "")
    if(NOT APPLE)
      target_include_directories(igl_opengl PRIVATE "${NANOGUI_DIR}/ext/glew/include")
      target_include_directories(igl_opengl2 PRIVATE "${NANOGUI_DIR}/ext/glew/include")
    endif()
  endif()
  set(LIBIGL_OPENGL_EXTRA_LIBRARIES ${OPENGL_gl_LIBRARY})

  ### GLEW for linux and windows
  if((UNIX AND NOT APPLE) OR WIN32) ### Compile glew if needed
    set(GLEW_INSTALL OFF CACHE BOOL " " FORCE)
    add_subdirectory("${NANOGUI_DIR}/ext/glew" "glew")
  endif()
  if(NOT APPLE)
    list(APPEND LIBIGL_INCLUDE_DIRS "${NANOGUI_DIR}/ext/glew/include")
    list(APPEND LIBIGL_OPENGL_EXTRA_LIBRARIES "glew")
  endif()

  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_OPENGL_EXTRA_LIBRARIES})

  if(LIBIGL_WITH_OPENGL_GLFW)
    set(GLFW_INCLUDE_DIRS "${NANOGUI_DIR}/ext/glfw/include")
    if(NOT APPLE)
      list(APPEND GLFW_INCLUDE_DIRS "${NANOGUI_DIR}/ext/glew/include")
    endif()

    # Note: if add_subdirectory("${NANOGUI_DIR}" "nanogui") runs below it will
    # add GLFW as a side-effect; in this case, CMake will complain about
    # duplicates if we add them here.
    if (NOT (LIBIGL_WITH_VIEWER AND LIBIGL_WITH_NANOGUI))
        add_subdirectory("${NANOGUI_DIR}/ext/glfw" "glfw")
    endif()
    set(LIBIGL_OPENGL_GLFW_EXTRA_LIBRARIES "glfw" ${GLFW_LIBRARIES})
    if(LIBIGL_USE_STATIC_LIBRARY)
      CompileIGL_Module("opengl/glfw" "")
      target_include_directories(igl_opengl_glfw PRIVATE ${GLFW_INCLUDE_DIRS})
      target_include_directories(igl_opengl_glfw PRIVATE ${OPENGL_INCLUDE_DIR})
      if(NOT APPLE)
        target_include_directories(igl_opengl PRIVATE "${NANOGUI_DIR}/ext/glew/include")
      endif()
    endif()
    list(APPEND LIBIGL_INCLUDE_DIRS ${GLFW_INCLUDE_DIRS})
    list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_OPENGL_GLFW_EXTRA_LIBRARIES})

    ### Compile the viewer ###
    if(LIBIGL_WITH_VIEWER)
      if(LIBIGL_WITH_NANOGUI)
        list(APPEND LIBIGL_DEFINITIONS "-DIGL_VIEWER_WITH_NANOGUI")
        if (LIBIGL_WITH_PYTHON)
          set(NANOGUI_BUILD_PYTHON ON CACHE BOOL " " FORCE)
        else()
          set(NANOGUI_BUILD_PYTHON OFF CACHE BOOL " " FORCE)
        endif()
        set(NANOGUI_BUILD_EXAMPLE OFF CACHE BOOL " " FORCE)
        set(NANOGUI_BUILD_SHARED  OFF CACHE BOOL " " FORCE)
        add_subdirectory("${NANOGUI_DIR}" "nanogui")
        set(VIEWER_INCLUDE_DIRS
          "${NANOGUI_DIR}/include"
          "${NANOGUI_DIR}/ext/nanovg/src")
        set(LIBIGL_VIEWER_EXTRA_LIBRARIES "nanogui" ${NANOGUI_EXTRA_LIBS})
      else()
        set(VIEWER_INCLUDE_DIRS "${NANOGUI_DIR}/ext/glfw/include")
      endif()
      list(APPEND VIEWER_INCLUDE_DIRS ${GLFW_INCLUDE_DIRS})
      list(APPEND LIBIGL_INCLUDE_DIRS ${VIEWER_INCLUDE_DIRS})
      list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_VIEWER_EXTRA_LIBRARIES})
      if(LIBIGL_USE_STATIC_LIBRARY)
        CompileIGL_Module("viewer" "")
        if(LIBIGL_WITH_NANOGUI)
          target_compile_definitions(igl_viewer PRIVATE -DIGL_VIEWER_WITH_NANOGUI)
        endif()
        target_include_directories(igl_viewer PRIVATE ${VIEWER_INCLUDE_DIRS})
        if(NOT APPLE)
          target_include_directories(igl_viewer PRIVATE "${NANOGUI_DIR}/ext/glew/include")
        endif()
      endif()
    endif()
  endif()

endif()

### Compile the png parts ###
if(LIBIGL_WITH_PNG)
  if(LIBIGL_WITH_NANOGUI)
    set(STBI_LOAD OFF CACHE BOOL " " FORCE)
  endif()
  set(STB_IMAGE_DIR "${LIBIGL_EXTERNAL}/stb_image")
  add_subdirectory("${STB_IMAGE_DIR}" "stb_image")
  list(APPEND LIBIGL_INCLUDE_DIRS ${STB_IMAGE_DIR})
  list(APPEND LIBIGL_PNG_EXTRA_LIBRARIES "stb_image")
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_PNG_EXTRA_LIBRARIES})
  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("png" "")
    target_include_directories(igl_png PRIVATE ${STB_IMAGE_DIR})
    if(NOT APPLE)
      target_include_directories(igl_png PRIVATE "${NANOGUI_DIR}/ext/glew/include")
    endif()
  endif()
endif()

### Compile the tetgen part ###
if(LIBIGL_WITH_TETGEN)
  set(TETGEN_DIR "${LIBIGL_EXTERNAL}/tetgen")
  add_subdirectory("${TETGEN_DIR}" "tetgen")
  list(APPEND LIBIGL_INCLUDE_DIRS ${TETGEN_DIR})
  list(APPEND LIBIGL_TETGEN_EXTRA_LIBRARIES "tetgen")
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_TETGEN_EXTRA_LIBRARIES})

  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("tetgen" "copyleft/")
    target_include_directories(igl_tetgen PRIVATE ${TETGEN_DIR})
  endif()
endif()

### Compile the triangle part ###
if(LIBIGL_WITH_TRIANGLE)
  set(TRIANGLE_DIR "${LIBIGL_EXTERNAL}/triangle")
  add_subdirectory("${TRIANGLE_DIR}" "triangle")
  list(APPEND LIBIGL_INCLUDE_DIRS ${TRIANGLE_DIR})
  list(APPEND LIBIGL_TRIANGLE_EXTRA_LIBRARIES "triangle")
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_TRIANGLE_EXTRA_LIBRARIES})

  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("triangle" "")
    target_include_directories(igl_triangle PRIVATE ${TRIANGLE_DIR})
  endif()
endif()

### Compile the xml part ###
if(LIBIGL_WITH_XML)
  set(TINYXML2_DIR "${LIBIGL_EXTERNAL}/tinyxml2")
  add_library(tinyxml2 STATIC ${TINYXML2_DIR}/tinyxml2.cpp ${TINYXML2_DIR}/tinyxml2.h)
  set_target_properties(tinyxml2 PROPERTIES
          COMPILE_DEFINITIONS "TINYXML2_EXPORT"
          VERSION "3.0.0"
          SOVERSION "3")
  list(APPEND LIBIGL_INCLUDE_DIRS ${TINYXML2_DIR})
  list(APPEND LIBIGL_XML_EXTRA_LIBRARIES "tinyxml2")
  list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_XML_EXTRA_LIBRARIES})
  if(LIBIGL_USE_STATIC_LIBRARY)
    CompileIGL_Module("xml" "")
    target_include_directories(igl_xml PRIVATE ${TINYXML2_DIR})
  endif()
endif()

# Function to print list nicely
function(print_list title list)
  message("-- ${title}:")
  foreach(elt ${list})
    message("\t ${elt}")
  endforeach()
endfunction()

# Pass the list of compiled libraries to the parent if there is one
if(NOT ${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
  list(REVERSE LIBIGL_LIBRARIES)

  set(LIBIGL_INCLUDE_DIRS ${LIBIGL_INCLUDE_DIRS} PARENT_SCOPE)
  set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} PARENT_SCOPE)
  set(LIBIGL_ANTTWEAKBAR_EXTRA_LIBRARIES ${LIBIGL_ANTTWEAKBAR_EXTRA_LIBRARIES} PARENT_SCOPE)
  set(LIBIGL_CGAL_EXTRA_LIBRARIES        ${LIBIGL_CGAL_EXTRA_LIBRARIES}        PARENT_SCOPE)
  set(LIBIGL_COMISO_EXTRA_LIBRARIES      ${LIBIGL_COMISO_EXTRA_LIBRARIES}      PARENT_SCOPE)
  set(LIBIGL_CORK_EXTRA_LIBRARIES        ${LIBIGL_CORK_EXTRA_LIBRARIES}        PARENT_SCOPE)
  set(LIBIGL_EMBREE_EXTRA_LIBRARIES      ${LIBIGL_EMBREE_EXTRA_LIBRARIES}      PARENT_SCOPE)
  set(LIBIGL_LIM_EXTRA_LIBRARIES         ${LIBIGL_LIM_EXTRA_LIBRARIES}         PARENT_SCOPE)
  set(LIBIGL_MATLAB_EXTRA_LIBRARIES      ${LIBIGL_MATLAB_EXTRA_LIBRARIES}      PARENT_SCOPE)
  set(LIBIGL_MOSEK_EXTRA_LIBRARIES       ${LIBIGL_MOSEK_EXTRA_LIBRARIES}       PARENT_SCOPE)
  set(LIBIGL_NANOGUI_EXTRA_LIBRARIES     ${LIBIGL_NANOGUI_EXTRA_LIBRARIES}     PARENT_SCOPE)
  set(LIBIGL_OPENGL_EXTRA_LIBRARIES      ${LIBIGL_OPENGL_EXTRA_LIBRARIES}     PARENT_SCOPE)
  set(LIBIGL_OPENGL_GLFW_EXTRA_LIBRARIES ${LIBIGL_OPENGL_GLFW_EXTRA_LIBRARIES} PARENT_SCOPE)
  set(LIBIGL_PNG_EXTRA_LIBRARIES         ${LIBIGL_PNG_EXTRA_LIBRARIES}         PARENT_SCOPE)
  set(LIBIGL_TETGEN_EXTRA_LIBRARIES      ${LIBIGL_TETGEN_EXTRA_LIBRARIES}      PARENT_SCOPE)
  set(LIBIGL_TRIANGLE_EXTRA_LIBRARIES    ${LIBIGL_TRIANGLE_EXTRA_LIBRARIES}    PARENT_SCOPE)
  set(LIBIGL_VIEWER_EXTRA_LIBRARIES      ${LIBIGL_VIEWER_EXTRA_LIBRARIES}      PARENT_SCOPE)
  set(LIBIGL_XML_EXTRA_LIBRARIES         ${LIBIGL_XML_EXTRA_LIBRARIES}         PARENT_SCOPE)
  set(LIBIGL_EXTRA_LIBRARIES ${LIBIGL_EXTRA_LIBRARIES} PARENT_SCOPE)
  set(LIBIGL_DEFINITIONS ${LIBIGL_DEFINITIONS} PARENT_SCOPE)

  ### ligIGL information ###
  print_list("libigl includes" "${LIBIGL_INCLUDE_DIRS}")
  print_list("libigl libraries" "${LIBIGL_LIBRARIES}")
  print_list("libigl extra libraries" "${LIBIGL_EXTRA_LIBRARIES}")
  print_list("libigl definitions" "${LIBIGL_DEFINITIONS}")
endif()
