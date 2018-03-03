cmake_minimum_required(VERSION 3.1)

### Find packages to populate default options ###
#
# COMPONENTS should match subsequent calls
find_package(CGAL COMPONENTS Core) # --> CGAL_FOUND
find_package(Boost 1.48 COMPONENTS thread system) # --> BOOST_FOUND
if(CGAL_FOUND AND BOOST_FOUND)
  set(CGAL_AND_BOOST_FOUND TRUE)
endif()
find_package(Matlab COMPONENTS MEX_COMPILER MX_LIBRARY ENG_LIBRARY) # --> Matlab_FOUND
find_package(MOSEK) # --> MOSEK_FOUND
find_package(OpenGL) # --> OPENGL_FOUND

### Available options ###
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"           "${CGAL_AND_BOOST_FOUND}")
option(LIBIGL_WITH_COMISO            "Use CoMiso"         ON)
option(LIBIGL_WITH_CORK              "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"         OFF)
option(LIBIGL_WITH_LIM               "Use LIM"            ON)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         "${Matlab_FOUND}")
option(LIBIGL_WITH_MOSEK             "Use MOSEK"          "${MOSEK_FOUND}")
option(LIBIGL_WITH_NANOGUI           "Use Nanogui menu"   OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         "${OPENGL_FOUND}")
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           "${OPENGL_FOUND}")
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use IMGUI"          OFF)
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
  if(LIBIGL_WITH_CGAL)
    target_compile_options(igl_common INTERFACE "/MD$<$<CONFIG:Debug>:d>")
  endif()
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

function(compile_igl_module module_dir)
  string(REPLACE "/" "_" module_name "${module_dir}")
  if(LIBIGL_USE_STATIC_LIBRARY)
    file(GLOB SOURCES_IGL_${module_name}
      "${LIBIGL_SOURCE_DIR}/igl/${module_dir}/*.cpp"
      "${LIBIGL_SOURCE_DIR}/igl/copyleft/${module_dir}/*.cpp")
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
endif()

################################################################################
### Compile the cgal parts ###
if(LIBIGL_WITH_CGAL)
  # CGAL Core is needed for
  # `Exact_predicates_exact_constructions_kernel_with_sqrt`
  if(EXISTS ${LIBIGL_EXTERNAL}/boost)
    set(BOOST_ROOT "${LIBIGL_EXTERNAL}/boost")
  endif()
  find_package(CGAL COMPONENTS Core)
  if(CGAL_FOUND)
    compile_igl_module("cgal")
    if(WIN32)
      set(Boost_USE_STATIC_LIBS ON) # Favor static Boost libs on Windows
    endif()
    target_include_directories(igl_cgal ${IGL_SCOPE} "${GMP_INCLUDE_DIR}" "${MPFR_INCLUDE_DIR}")
    find_package(Boost 1.48 REQUIRED thread system)
    target_include_directories(igl_cgal ${IGL_SCOPE} ${CGAL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
    target_link_libraries(igl_cgal ${IGL_SCOPE} CGAL::CGAL CGAL::CGAL_Core ${Boost_LIBRARIES})
  else()
    set(LIBIGL_WITH_CGAL OFF CACHE BOOL "" FORCE)
  endif()
endif()

# Helper function for `igl_copy_cgal_dll()`
function(igl_copy_imported_dll src_target dst_target)
  get_target_property(configurations ${src_target} IMPORTED_CONFIGURATIONS)
  foreach(config ${configurations})
    get_target_property(main_lib ${src_target} IMPORTED_LOCATION_${config})
    get_target_property(other_libs ${src_target} IMPORTED_LINK_INTERFACE_LIBRARIES_${config})
    set(locations)
    list(APPEND locations ${main_lib} ${other_libs})
    foreach(location ${locations})
      string(REGEX MATCH "^(.*)\\.[^.]*$" dummy ${location})
      set(location "${CMAKE_MATCH_1}.dll")
      if(EXISTS "${location}" AND location MATCHES "^.*\\.dll$")
        add_custom_command(TARGET ${dst_target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different "${location}" $<TARGET_FILE_DIR:${dst_target}>)
      endif()
    endforeach()
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
# Compile CoMISo
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
### Compile the cork parts ###
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
### Compile the opengl parts ###

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
    target_include_directories(tinyxml2 PUBLIC ${TINYXML2_DIR})
    set_target_properties(tinyxml2 PROPERTIES
            COMPILE_DEFINITIONS "TINYXML2_EXPORT"
            VERSION "3.0.0"
            SOVERSION "3")
  endif()
  compile_igl_module("xml")
  target_link_libraries(igl_xml ${IGL_SCOPE} tinyxml2)
endif()
