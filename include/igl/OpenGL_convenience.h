#ifndef OPENGL_CONVENIENCE_H
#define OPENGL_CONVENIENCE_H
#ifndef IGL_NO_OPENGL

// Always use this:
//     #include "OpenGL_convenience.h"
// Convenience includer for opengl.

// For now this includes glu, glew and glext (perhaps these should be
// separated)
#if __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
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
#endif

#endif
#endif

