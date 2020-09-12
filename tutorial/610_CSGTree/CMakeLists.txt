get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} igl::core igl::cgal igl::opengl igl::opengl_glfw tutorials)
igl_copy_cgal_dll(${PROJECT_NAME})
