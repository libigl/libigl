cmake_minimum_required(VERSION 3.1)

# https://github.com/libigl/libigl/issues/751
# http://lists.llvm.org/pipermail/llvm-commits/Week-of-Mon-20160425/351643.html
if(APPLE)
  if(NOT CMAKE_LIBTOOL)
    find_program(CMAKE_LIBTOOL NAMES libtool)
  endif()
  if(CMAKE_LIBTOOL)
    set(CMAKE_LIBTOOL ${CMAKE_LIBTOOL} CACHE PATH "libtool executable")
    message(STATUS "Found libtool - ${CMAKE_LIBTOOL}")
    get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
    foreach(lang ${languages})
      # Added -c 
      set(CMAKE_${lang}_CREATE_STATIC_LIBRARY
        "${CMAKE_LIBTOOL} -c -static -o <TARGET> <LINK_FLAGS> <OBJECTS> ")
    endforeach()
  endif()
endif()

### Find packages to populate default options ###
#
# COMPONENTS should match subsequent calls
find_package(Matlab COMPONENTS MEX_COMPILER MX_LIBRARY ENG_LIBRARY) # --> Matlab_FOUND
find_package(MOSEK) # --> MOSEK_FOUND
find_package(OpenGL) # --> OPENGL_FOUND

### Available options ###
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" ON)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           ON)
option(LIBIGL_WITH_COMISO            "Use CoMiso"         ON)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            ON)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         "${Matlab_FOUND}")
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          "${MOSEK_FOUND}")
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         "${OPENGL_FOUND}")
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           "${OPENGL_FOUND}")
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          OFF)
option(LIBIGL_WITH_PNG               "Use PNG"            ON)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"         ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"       ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  "${OPENGL_FOUND}")
option(LIBIGL_WITH_XML               "Use XML"            ON)
option(LIBIGL_WITH_PYTHON            "Use Python"         OFF)

if(LIBIGL_WITH_VIEWER AND (NOT LIBIGL_WITH_OPENGL_GLFW OR NOT LIBIGL_WITH_OPENGL) )
  message(FATAL_ERROR "LIBIGL_WITH_VIEWER=ON requires LIBIGL_WITH_OPENGL_GLFW=ON and LIBIGL_WITH_OPENGL=ON")
endif()

################################################################################

### Configuration
set(LIBIGL_ROOT "${CMAKE_CURRENT_LIST_DIR}/../..")
set(LIBIGL_SOURCE_DIR "${LIBIGL_ROOT}/include")
set(LIBIGL_EXTERNAL "${LIBIGL_ROOT}/external")

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
  target_include_directories(igl_common SYSTEM INTERFACE ${LIBIGL_EXTERNAL}/eigen)
endif()

# C++11 Thread library
find_package(Threads REQUIRED)
target_link_libraries(igl_common INTERFACE ${CMAKE_THREAD_LIBS_INIT})

################################################################################

include(DownloadProject)

# Shortcut function
function(igl_download_project name)
  download_project(
    PROJ         ${name}
    SOURCE_DIR   ${LIBIGL_EXTERNAL}/${name}
    DOWNLOAD_DIR ${LIBIGL_EXTERNAL}/.cache/${name}
    ${ARGN}
  )
endfunction()

################################################################################

## CGAL dependencies on Windows: GMP & MPFR
function(igl_download_cgal_deps)
  if(WIN32)
    igl_download_project(gmp
        URL     https://cgal.geometryfactory.com/CGAL/precompiled_libs/auxiliary/x64/GMP/5.0.1/gmp-all-CGAL-3.9.zip
        URL_MD5 508c1292319c832609329116a8234c9f
    )
    igl_download_project(mpfr
        URL https://cgal.geometryfactory.com/CGAL/precompiled_libs/auxiliary/x64/MPFR/3.0.0/mpfr-all-CGAL-3.9.zip
        URL_MD5 48840454eef0ff18730050c05028734b
    )
    set(ENV{GMP_DIR} "${LIBIGL_EXTERNAL}/gmp")
    set(ENV{MPFR_DIR} "${LIBIGL_EXTERNAL}/mpfr")
  endif()
endfunction()

################################################################################

function(compile_igl_module module_dir)
  string(REPLACE "/" "_" module_name "${module_dir}")
  if(module_name STREQUAL "core")
    set(module_libname "igl")
  else()
    set(module_libname "igl_${module_name}")
  endif()
  if(LIBIGL_USE_STATIC_LIBRARY)
    file(GLOB SOURCES_IGL_${module_name}
      "${LIBIGL_SOURCE_DIR}/igl/${module_dir}/*.cpp"
      "${LIBIGL_SOURCE_DIR}/igl/copyleft/${module_dir}/*.cpp")
    add_library(${module_libname} STATIC ${SOURCES_IGL_${module_name}} ${ARGN})
    if(MSVC)
      target_compile_options(${module_libname} PRIVATE /w) # disable all warnings (not ideal but...)
    else()
      #target_compile_options(${module_libname} PRIVATE -w) # disable all warnings (not ideal but...)
    endif()
  else()
    add_library(${module_libname} INTERFACE)
  endif()

  target_link_libraries(${module_libname} ${IGL_SCOPE} igl_common)
  if(NOT module_name STREQUAL "core")
    target_link_libraries(${module_libname} ${IGL_SCOPE} igl)
  endif()

  # Alias target because it looks nicer
  message(STATUS "Creating target: igl::${module_name} (${module_libname})")
  add_library(igl::${module_name} ALIAS ${module_libname})
endfunction()

################################################################################
### IGL Core
################################################################################

if(LIBIGL_USE_STATIC_LIBRARY)
  file(GLOB SOURCES_IGL
    "${LIBIGL_SOURCE_DIR}/igl/*.cpp"
    "${LIBIGL_SOURCE_DIR}/igl/copyleft/*.cpp")
endif()
compile_igl_module("core" ${SOURCES_IGL})

################################################################################
## Compile the AntTweakBar part ###
if(LIBIGL_WITH_ANTTWEAKBAR)
  set(ANTTWEAKBAR_DIR "${LIBIGL_EXTERNAL}/AntTweakBar")
  if(NOT TARGET AntTweakBar)
    add_subdirectory("${ANTTWEAKBAR_DIR}" AntTweakBar)
  endif()
  compile_igl_module("anttweakbar")
  target_link_libraries(igl_anttweakbar ${IGL_SCOPE} AntTweakBar)
  target_include_directories(igl_anttweakbar ${IGL_SCOPE} "${ANTTWEAKBAR_DIR}/include")
endif()

################################################################################
### Compile the CGAL part ###
if(LIBIGL_WITH_CGAL)
  # Try to find the CGAL library
  # CGAL Core is needed for
  # `Exact_predicates_exact_constructions_kernel_with_sqrt`
  if(NOT TARGET CGAL::CGAL)
    set(CGAL_DIR "${LIBIGL_EXTERNAL}/cgal")
    igl_download_cgal_deps()
    if(EXISTS ${LIBIGL_EXTERNAL}/boost)
      set(BOOST_ROOT "${LIBIGL_EXTERNAL}/boost")
    endif()
    set(CGAL_Boost_USE_STATIC_LIBS ON CACHE BOOL "" FORCE)
    find_package(CGAL CONFIG COMPONENTS Core PATHS ${CGAL_DIR} NO_DEFAULT_PATH)
  endif()

  # If CGAL has been found, then build the libigl module
  if(TARGET CGAL::CGAL AND TARGET CGAL::CGAL_Core)
    compile_igl_module("cgal")
    target_link_libraries(igl_cgal ${IGL_SCOPE} CGAL::CGAL CGAL::CGAL_Core)
  else()
    set(LIBIGL_WITH_CGAL OFF CACHE BOOL "" FORCE)
  endif()
endif()

# Helper function for `igl_copy_cgal_dll()`
function(igl_copy_imported_dll src_target dst_target)
  get_target_property(other_libs ${src_target} INTERFACE_LINK_LIBRARIES)
  set(locations)
  list(APPEND locations ${main_lib} ${other_libs})
  foreach(location ${locations})
    string(REGEX MATCH "^(.*)\\.[^.]*$" dummy ${location})
    set(location "${CMAKE_MATCH_1}.dll")
    if(EXISTS "${location}" AND location MATCHES "^.*\\.dll$")
      add_custom_command(TARGET ${dst_target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different "${location}" $<TARGET_FILE_DIR:${dst_target}>)
    endif()
  endforeach()
endfunction()

# Convenient functions to copy CGAL dlls into a target (executable) destination folder (for Windows)
function(igl_copy_cgal_dll target)
  if(WIN32 AND LIBIGL_WITH_CGAL)
    igl_copy_imported_dll(CGAL::CGAL ${target})
    igl_copy_imported_dll(CGAL::CGAL_Core ${target})
  endif()
endfunction()

################################################################################
### Compile the CoMISo part ###
# NOTE: this cmakefile works only with the
# comiso available here: https://github.com/libigl/CoMISo
if(LIBIGL_WITH_COMISO)
  compile_igl_module("comiso")
  if(NOT TARGET CoMISo)
    add_subdirectory("${LIBIGL_EXTERNAL}/CoMISo" CoMISo)
  endif()
  target_link_libraries(igl_comiso ${IGL_SCOPE} CoMISo)
endif()

################################################################################
### Compile the cork part ###
if(LIBIGL_WITH_CORK)
  set(CORK_DIR "${LIBIGL_EXTERNAL}/cork")
  if(NOT TARGET cork)
    # call this "lib-cork" instead of "cork", otherwise cmake gets confused about
    # "cork" executable
    add_subdirectory("${CORK_DIR}" "lib-cork")
  endif()
  compile_igl_module("cork")
  target_include_directories(igl_cork ${IGL_SCOPE} cork)
  target_include_directories(igl_cork ${IGL_SCOPE} "${CORK_DIR}/src")
  target_link_libraries(igl_cork ${IGL_SCOPE} cork)
endif()

################################################################################
### Compile the embree part ###
if(LIBIGL_WITH_EMBREE)
  set(EMBREE_DIR "${LIBIGL_EXTERNAL}/embree")

  set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL " " FORCE)
  set(EMBREE_TASKING_SYSTEM "INTERNAL" CACHE BOOL " " FORCE)
  set(EMBREE_TUTORIALS OFF CACHE BOOL " " FORCE)
  set(EMBREE_MAX_ISA NONE CACHE STRINGS " " FORCE)
  set(BUILD_TESTING OFF CACHE BOOL " " FORCE)

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

  compile_igl_module("embree")
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
  compile_igl_module("lim")
  target_link_libraries(igl_lim ${IGL_SCOPE} lim)
  target_include_directories(igl_lim ${IGL_SCOPE} ${LIM_DIR})
endif()

################################################################################
### Compile the matlab part ###
if(LIBIGL_WITH_MATLAB)
  find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER MX_LIBRARY ENG_LIBRARY)
  compile_igl_module("matlab")
  target_link_libraries(igl_matlab ${IGL_SCOPE} ${Matlab_LIBRARIES})
  target_include_directories(igl_matlab ${IGL_SCOPE} ${Matlab_INCLUDE_DIRS})
endif()

################################################################################
### Compile the mosek part ###
if(LIBIGL_WITH_MOSEK)
  find_package(MOSEK REQUIRED)
  compile_igl_module("mosek")
  target_link_libraries(igl_mosek ${IGL_SCOPE} ${MOSEK_LIBRARIES})
  target_include_directories(igl_mosek ${IGL_SCOPE} ${MOSEK_INCLUDE_DIRS})
  target_compile_definitions(igl_mosek ${IGL_SCOPE} -DLIBIGL_WITH_MOSEK)
endif()

################################################################################
### Compile the opengl part ###
if(LIBIGL_WITH_OPENGL)
  # OpenGL module
  find_package(OpenGL REQUIRED)
  compile_igl_module("opengl")
  target_link_libraries(igl_opengl ${IGL_SCOPE} ${OPENGL_gl_LIBRARY})
  target_include_directories(igl_opengl SYSTEM ${IGL_SCOPE} ${OPENGL_INCLUDE_DIR})

  # glad module
  if(NOT TARGET glad)
    add_subdirectory(${LIBIGL_EXTERNAL}/glad glad)
  endif()
  target_link_libraries(igl_opengl ${IGL_SCOPE} glad)
endif()

################################################################################
### Compile the GLFW part ###
if(LIBIGL_WITH_OPENGL_GLFW)
  if(TARGET igl::opengl)
    # GLFW module
    compile_igl_module("opengl/glfw")
    if(NOT TARGET glfw)
      set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
      set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
      set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
      set(GLFW_INSTALL OFF CACHE BOOL " " FORCE)
      add_subdirectory(${LIBIGL_EXTERNAL}/glfw glfw)
    endif()
    target_link_libraries(igl_opengl_glfw ${IGL_SCOPE} igl_opengl glfw)
  endif()
endif()

################################################################################
### Compile the ImGui part ###
if(LIBIGL_WITH_OPENGL_GLFW_IMGUI)
  if(TARGET igl::opengl_glfw)
    # ImGui module
    compile_igl_module("opengl/glfw/imgui")
    if(NOT TARGET imgui)
      add_subdirectory(${LIBIGL_EXTERNAL}/imgui imgui)
    endif()
    target_link_libraries(igl_opengl_glfw_imgui ${IGL_SCOPE} igl_opengl_glfw imgui)
  endif()
endif()

################################################################################
### Compile the png part ###
if(LIBIGL_WITH_PNG)
  # png/ module is anomalous because it also depends on opengl it really should
  # be moved into the opengl/ directory and namespace ...
  if(TARGET igl_opengl)
    set(STB_IMAGE_DIR "${LIBIGL_EXTERNAL}/stb_image")
    if(NOT TARGET stb_image)
      add_subdirectory("${STB_IMAGE_DIR}" "stb_image")
    endif()
    compile_igl_module("png" "")
    target_link_libraries(igl_png ${IGL_SCOPE} igl_stb_image igl_opengl)
  endif()
endif()

################################################################################
### Compile the tetgen part ###
if(LIBIGL_WITH_TETGEN)
  set(TETGEN_DIR "${LIBIGL_EXTERNAL}/tetgen")
  if(NOT TARGET tetgen)
    add_subdirectory("${TETGEN_DIR}" "tetgen")
  endif()
  compile_igl_module("tetgen")
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
  compile_igl_module("triangle")
  target_link_libraries(igl_triangle ${IGL_SCOPE} triangle)
  target_include_directories(igl_triangle ${IGL_SCOPE} ${TRIANGLE_DIR})
endif()

################################################################################
### Compile the xml part ###
if(LIBIGL_WITH_XML)
  set(TINYXML2_DIR "${LIBIGL_EXTERNAL}/tinyxml2")
  if(NOT TARGET tinyxml2)
    add_library(tinyxml2 STATIC ${TINYXML2_DIR}/tinyxml2.cpp ${TINYXML2_DIR}/tinyxml2.h)
    set_target_properties(tinyxml2 PROPERTIES
            COMPILE_DEFINITIONS "TINYXML2_EXPORT"
            VERSION "3.0.0"
            SOVERSION "3")
  endif()
  compile_igl_module("xml")
  target_link_libraries(igl_xml ${IGL_SCOPE} tinyxml2)
  target_include_directories(igl_xml ${IGL_SCOPE} ${TINYXML2_DIR})
endif()
