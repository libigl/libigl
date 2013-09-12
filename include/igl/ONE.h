#ifndef IGL_ONE_H
#define IGL_ONE_H
// Often one needs a reference to a dummy variable containing one as its
// value, for example when using AntTweakBar's
// TwSetParam( "3D View", "opened", TW_PARAM_INT32, 1, &INT_ONE);
namespace igl
{
  const char CHAR_ONE = 1;
  const int INT_ONE = 1;
  const unsigned int UNSIGNED_INT_ONE = 1;
  const double DOUBLE_ONE = 1;
  const float FLOAT_ONE = 1;
}
#endif

