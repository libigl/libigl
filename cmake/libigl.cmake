# Sanity check for backward compatibility
get_filename_component(PATH_X "${libigl_SOURCE_DIR}" REALPATH)
get_filename_component(PATH_Y "${CMAKE_CURRENT_LIST_DIR}/.." REALPATH)
if(NOT PATH_X STREQUAL PATH_Y)
    message(FATAL_ERROR "You included libigl.cmake directly from your own project. This behavior "
                        "is not supported anymore. Please add libigl to your project via "
                        "add_subdirectory(<path_to_libigl>). See the libigl example project for "
                        "more information: https://github.com/libigl/libigl-example-project/")
endif()

# Global options
include(igl_windows)

# Libigl permissive modules
igl_include(core)
