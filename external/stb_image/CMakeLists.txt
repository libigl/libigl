cmake_minimum_required(VERSION 3.1)
project(igl_stb_image)

if(NOT CMAKE_VERSION VERSION_LESS "3.3")
	# https://cmake.org/cmake/help/latest/policy/CMP0063.html
	cmake_policy(SET CMP0063 NEW)
else()
	message(WARNING "CMake version is < 3.3, visibility properties for static libraries will be ignored.")
endif()

add_library(${PROJECT_NAME} igl_stb_image.cpp stb_image.h stb_image_write.h)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Export interface symbols and hide internal implementation
include(GenerateExportHeader)
generate_export_header(${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_VISIBILITY_PRESET hidden)
