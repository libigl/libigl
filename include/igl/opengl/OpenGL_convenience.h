// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_OPENGL_OPENGL_CONVENIENCE_H
#define IGL_OPENGL_OPENGL_CONVENIENCE_H

// Always use this:
//     #include "OpenGL_convenience.h"
// Convenience includer for opengl.

// For now this includes glu, glew and glext (perhaps these should be
// separated)
/*#if __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  include <OpenGL/glext.h>
#elif defined(_WIN32)
#    define NOMINMAX
#    include <Windows.h>
#    undef NOMINMAX
#    include <GL/glew.h>
#    include <GL/gl.h>
#else
#  define GL_GLEXT_PROTOTYPES
#  include <GL/gl.h>
#  include <GL/glext.h>
#  include <GL/glu.h>
#endif*/

// From Nanogui
#if defined(__APPLE__)
    #define GLFW_INCLUDE_GLCOREARB
#elif defined(WIN32)
    #define GLEW_STATIC
    #include <GL/glew.h>
#else
    #define GLEW_STATIC
    #include <GL/glew.h>
    #define GL_GLEXT_PROTOTYPES
#endif

//#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#endif
