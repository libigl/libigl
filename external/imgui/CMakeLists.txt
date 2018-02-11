cmake_minimum_required(VERSION 3.1)

add_library(imgui
		imgui/imgui.cpp
		imgui/imgui_demo.cpp
		imgui/imgui_draw.cpp
		imgui_impl_glfw_gl3.cpp
)

target_include_directories(imgui PRIVATE imgui)
target_include_directories(imgui SYSTEM INTERFACE .)

target_link_libraries(imgui PRIVATE glad glfw)
