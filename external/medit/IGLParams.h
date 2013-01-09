#ifndef IGL_PARAMS_H
#define IGL_PARAMS_H
#ifdef IGL

struct IGLParams
{
  // Show lines on cap (not just surface)
  bool lines_on_cap;
  // Show hot-dog when clipping
  bool hot_dog_view;
  // Number of hot-dog slices
  int num_hot_dog_slices;
  IGLParams():
    lines_on_cap(false),
    hot_dog_view(false),
    num_hot_dog_slices(10)
  {};
};
#endif
#endif
