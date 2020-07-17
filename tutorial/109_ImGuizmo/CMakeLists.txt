get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

set(IMGUIZMO_PATH ${LIBIGL_SOURCE_DIR}/igl/opengl/glfw/imgui/imguizmo)
add_executable(${PROJECT_NAME}_bin main.cpp)
add_library( imguizmo ${IMGUIZMO_PATH}/ImGuizmo.cpp ${IMGUIZMO_PATH}/ImGuizmo.h )
target_link_libraries( imguizmo igl::opengl_glfw_imgui )
target_link_libraries( ${PROJECT_NAME}_bin igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui imguizmo tutorials)
