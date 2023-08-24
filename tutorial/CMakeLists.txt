# Choose which chapters to compile
option(LIBIGL_TUTORIALS_CHAPTER1 "Compile libigl tutorials chapter 1" ON)
option(LIBIGL_TUTORIALS_CHAPTER2 "Compile libigl tutorials chapter 2" ON)
option(LIBIGL_TUTORIALS_CHAPTER3 "Compile libigl tutorials chapter 3" ON)
option(LIBIGL_TUTORIALS_CHAPTER4 "Compile libigl tutorials chapter 4" ON)
option(LIBIGL_TUTORIALS_CHAPTER5 "Compile libigl tutorials chapter 5" ON)
option(LIBIGL_TUTORIALS_CHAPTER6 "Compile libigl tutorials chapter 6" ON)
option(LIBIGL_TUTORIALS_CHAPTER7 "Compile libigl tutorials chapter 7" ON)
option(LIBIGL_TUTORIALS_CHAPTER8 "Compile libigl tutorials chapter 8" ON)
option(LIBIGL_TUTORIALS_CHAPTER9 "Compile libigl tutorials chapter 9" ON)

# Download tutorial data
include(libigl_tutorial_data)

if(NOT TARGET igl::glfw)
    message(WARNING "Most tutorial executables depend on OpenGL and glfw. Use `cmake -DLIBIGL_GLFW=ON ..`")
endif()

# Chapter 1
if(LIBIGL_TUTORIALS_CHAPTER1)
    igl_add_tutorial(101_FileIO)
endif()
