#ifndef IGL_PARAMS_H
#define IGL_PARAMS_H
#ifdef IGL

// Defines pMesh
#include "grafic.h"
enum ColorMap {
  COLOR_MAP_DEFAULT = 0,
  COLOR_MAP_JET = 1,
  COLOR_MAP_EASTER = 2,
  COLOR_MAP_WINDING_THEN_EASTER = 3,
  NUM_COLOR_MAP = 4,
};

class IGLParams
{
  public:
    // Show lines on cap (not just surface)
    bool lines_on_cap;
    // Show hot-dog when clipping
    bool hot_dog_view;
    // Number of hot-dog slices
    int num_hot_dog_slices;
    // ratio of in to out hot dog slice widths
    double hot_dog_ratio;
    ColorMap color_map;
    float easter_red[3];
    float easter_s;
    float easter_v;
    float tet_color[4];
    float open_color[4];
    float nme_color[4];
    bool render_on_C_UPDATE;
    bool render_on_next;
    float dot_size;
    // Place holder alpha value used in cutTriangle
    float alpha_holder;
    // Fade parameters
    double fade_flip;
    double fade_max_s;
    double fade_min_s;
    double fade_max_v;
    double fade_min_v;
  public:
    IGLParams():
      lines_on_cap(false),
      hot_dog_view(false),
      num_hot_dog_slices(10),
      hot_dog_ratio(0.5),
      color_map(COLOR_MAP_DEFAULT),
      render_on_C_UPDATE(false),
      render_on_next(false),
      alpha_holder(1.0),
      fade_flip(false),
      fade_max_s(0.5),
      fade_min_s(0.25),
      fade_max_v(0.5),
      fade_min_v(0.01)
    {
      easter_red[0] = 0.8;
      easter_red[1] = 0.1;
      easter_red[2] = 0.8;
      easter_red[3] = 1.0;
      easter_s = 0.1;
      easter_v = 0.8;
      tet_color[0] = 1.0;
      tet_color[1] = 0.0;
      tet_color[2] = 0.0;
      tet_color[3] = 1.0;
      nme_color[0] = 0.22;
      nme_color[1] = 1.0;
      nme_color[2] = 0.2;
      nme_color[3] = 1.0;
      open_color[0] = 0.2;
      open_color[1] = 0.22;
      open_color[2] = 1.0;
      open_color[3] = 1.0;
      dot_size = 10;
    };
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
    // get color from colormap, mlist.c
    //
    // Inputs:
    //   x  value between 0 and 1
    void rgb(double x, double * rgb);
    ;
};
#endif
#endif
