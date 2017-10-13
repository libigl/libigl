#ifndef IGL_STB_IMAGE_H
#define IGL_STB_IMAGE_H

#include <igl_stb_image_export.h>

namespace igl {
  IGL_STB_IMAGE_EXPORT unsigned char * stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp);
  IGL_STB_IMAGE_EXPORT void stbi_image_free(void *retval_from_stbi_load);
  IGL_STB_IMAGE_EXPORT int stbi_write_png(char const *filename, int w, int h, int comp, const void *data, int stride_in_bytes);
} // namespace igl

#endif
