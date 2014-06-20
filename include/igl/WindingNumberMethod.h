#ifndef IGL_WINDINGNUMBERMETHOD_H
#define IGL_WINDINGNUMBERMETHOD_H
namespace igl
{
  enum WindingNumberMethod
  {
    EXACT_WINDING_NUMBER_METHOD = 0, // Exact hierarchical evaluation
    APPROX_SIMPLE_WINDING_NUMBER_METHOD = 1,
    APPROX_CACHE_WINDING_NUMBER_METHOD = 2, // Approximate hierarchical evaluation
    NUM_WINDING_NUMBER_METHODS = 3
  };
}
#endif
