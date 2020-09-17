get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

set(IMGUIZMO_PATH ${LIBIGL_SOURCE_DIR}/igl/opengl/glfw/imgui/imguizmo)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries( ${PROJECT_NAME} igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui tutorials)
