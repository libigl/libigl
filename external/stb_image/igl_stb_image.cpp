#include "igl_stb_image.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


unsigned char * igl::stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp)
{
  return ::stbi_load(filename, x, y, comp, req_comp);
}

void igl::stbi_image_free(void *retval_from_stbi_load)
{
  ::stbi_image_free(retval_from_stbi_load);
}

int igl::stbi_write_png(
  char const *filename, int w, int h, int comp, const void *data, int stride_in_bytes)
{
  return ::stbi_write_png(filename, w, h, comp, data, stride_in_bytes);
}
