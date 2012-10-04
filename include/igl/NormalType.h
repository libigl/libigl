#ifndef IGL_NORMALTYPE_H
#define IGL_NORMALTYPE_H

namespace igl
{
  enum NormalType
  {
    PER_VERTEX_NORMALS,
    PER_FACE_NORMALS,
    PER_CORNER_NORMALS,
  };
#  define NUM_NORMAL_TYPE 3
}

#endif

