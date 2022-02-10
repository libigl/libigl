# Sanity check for backward compatibility
get_filename_component(PATH_X "${PROJECT_SOURCE_DIR}" REALPATH)
get_filename_component(PATH_Y "${CMAKE_CURRENT_LIST_DIR}/.." REALPATH)
if(NOT PATH_X STREQUAL PATH_Y)
    message(FATAL_ERROR "You included libigl.cmake directly from your own project. This behavior "
                        "is not supported anymore. Please add libigl to your project via "
                        "add_subdirectory(<path_to_libigl>). See the libigl example project for "
                        "more information: https://github.com/libigl/libigl-example-project/")
endif()

# Dependencies are linked as INTERFACE targets unless libigl is compiled as a static library
if(LIBIGL_USE_STATIC_LIBRARY)
    set(IGL_SCOPE PUBLIC)
else()
    set(IGL_SCOPE INTERFACE)
endif()

# Global options
include(igl_windows)

# Libigl permissive modules
igl_include(core)
igl_include_optional(embree)
igl_include_optional(opengl)
igl_include_optional(glfw)
igl_include_optional(imgui)
igl_include_optional(predicates)
igl_include_optional(png)
igl_include_optional(xml)

# Libigl copyleft modules
igl_include_optional(copyleft core)
igl_include_optional(copyleft cgal)
igl_include_optional(copyleft comiso)
igl_include_optional(copyleft tetgen)

# Libigl restricted modules
igl_include_optional(restricted matlab)
igl_include_optional(restricted mosek)
igl_include_optional(restricted triangle)
