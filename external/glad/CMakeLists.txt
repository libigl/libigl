cmake_minimum_required(VERSION 3.1)
project(glad)

add_library(glad src/glad.c)
target_include_directories(glad PUBLIC include)

if(NOT WIN32)
  target_link_libraries(glad PUBLIC ${CMAKE_DL_LIBS})
endif()
