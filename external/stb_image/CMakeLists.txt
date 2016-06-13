cmake_minimum_required(VERSION 2.6)
project(stb_image)

option(STBI_LOAD              "Compile stbi"                         ON)
option(STBI_WRITE             "Compile stbi_writeUse CoMiso"         ON)

if (STBI_LOAD)
	add_definitions(-DSTBI_LOAD)
endif()

if (STBI_WRITE)
	add_definitions(-DSTBI_WRITE)
endif()

add_library(stb_image STATIC stb_image.cpp stb_image.h stb_image_write.h)
