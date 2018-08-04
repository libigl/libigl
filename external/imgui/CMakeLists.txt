cmake_minimum_required(VERSION 3.1)

add_library(imgui
		imgui/imgui.cpp
		imgui/imgui_demo.cpp
		imgui/imgui_draw.cpp
		imgui/examples/imgui_impl_glfw.cpp
		imgui_impl_opengl3.cpp
)

target_include_directories(imgui PRIVATE imgui)
target_include_directories(imgui PUBLIC . imgui/examples)

target_link_libraries(imgui PRIVATE glad glfw)
