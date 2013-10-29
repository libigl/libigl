#ifndef IGL_HIT_H
#define IGL_HIT_H

namespace igl
{
  // Reimplementation of the embree::Hit struct from embree1.0
  struct Hit
  {
    int id; // primitive id
    float u,v; // barycentric coordinates
    float t; // distance = direction*t to intersection
  };
}
#endif 
