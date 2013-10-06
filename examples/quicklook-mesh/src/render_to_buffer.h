#define GL_GLEXT_PROTOTYPES
#include <GL/osmesa.h>

// Inputs:
//   filename   path to mesh
//   background  background color
//   width  width of image buffer
//   height  height of image buffer
// Outputs:
//   buffer   width*height*4 buffer of bytes (should already be allocated) RGBA
// Returns true only only upon success.
bool render_to_buffer(
  const char * filename,
  const float *background_color,
  const int width,
  const int height,
  GLubyte * buffer);
