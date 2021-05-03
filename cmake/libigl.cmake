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

# Libigl modules
include(modules/core)
# include(modules/embree)
# include(modules/opengl)
# include(modules/glfw)
# include(modules/imgui)
# include(modules/predicates)
# include(modules/png)
# include(modules/xml)

# Libigl nonfree modules
# include(modules/nonfree/matlab)
# include(modules/nonfree/mosek)
# include(modules/nonfree/triangle)

# Libigl copyleft modules
# include(modules/copyleft/core)
# include(modules/copyleft/cgal)
# include(modules/copyleft/comiso)
# include(modules/copyleft/cork)
# include(modules/copyleft/tetgen)
