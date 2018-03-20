cmake_minimum_required(VERSION 3.1)

add_library(AntTweakBar STATIC
  # C source files
  src/TwEventGLFW.c
  src/TwEventGLUT.c
  src/TwEventSDL.c
  src/TwEventSDL12.c
  src/TwEventSDL13.c

  # C++ source files
  src/LoadOGL.cpp
  src/LoadOGLCore.cpp
  src/TwBar.cpp
  src/TwColors.cpp
  src/TwEventSFML.cpp
  src/TwFonts.cpp
  src/TwMgr.cpp
  src/TwOpenGL.cpp
  src/TwOpenGLCore.cpp
  src/TwPrecomp.cpp
)

# These are probably needed for windows/Linux, should append if windows/Linux
# src/TwEventWin.c
# src/TwEventX11.c
# src/TwDirect3D10.cpp
# src/TwDirect3D11.cpp
# src/TwDirect3D9.cpp

target_include_directories(AntTweakBar PUBLIC include)

if(APPLE)
  set_target_properties(
    AntTweakBar
    PROPERTIES
    COMPILE_FLAGS
    "-fPIC -fno-strict-aliasing -x objective-c++")
  target_compile_definitions(
    AntTweakBar PUBLIC _MACOSX __PLACEMENT_NEW_INLINE)
endif()
