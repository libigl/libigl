#ifndef IGL_PARAMS_H
#define IGL_PARAMS_H
#ifdef IGL

// Defines pMesh
#include "grafic.h"

class IGLParams
{
  public:
    // Show lines on cap (not just surface)
    bool lines_on_cap;
    // Show hot-dog when clipping
    bool hot_dog_view;
    // Number of hot-dog slices
    int num_hot_dog_slices;
  public:
    IGLParams():
      lines_on_cap(false),
      hot_dog_view(false),
      num_hot_dog_slices(10)
    {};
    // width of one hot dog slice
    double width(pMesh mesh) const
    {
      // Ridiculous way to get largest length
      return
        (mesh->xmax-mesh->xmin > 
          mesh->ymax-mesh->ymin ?
            (mesh->xmax-mesh->xmin > 
              mesh->zmax-mesh->zmin ?
                mesh->xmax-mesh->xmin :
                  mesh->zmax-mesh->zmin) :
            (mesh->ymax-mesh->ymin > 
              mesh->zmax-mesh->zmin ?
                mesh->ymax-mesh->ymin :
                  mesh->zmax-mesh->zmin)) /
        (double)num_hot_dog_slices;
    };
};
#endif
#endif
