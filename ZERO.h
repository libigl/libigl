#ifndef ZERO_H
#define ZERO_H
// Often one needs a reference to a dummy variable containing zero as its
// value, for example when using AntTweakBar's
// TwSetParam( "3D View", "opened", TW_PARAM_INT32, 1, &INT_ZERO);
namespace igl
{
  const char CHAR_ZERO = 0;
  const int INT_ZERO = 0;
  const unsigned int UNSIGNED_INT_ZERO = 0;
  const double DOUBLE_ZERO = 0;
  const float FLOAT_ZERO = 0;
}
#endif
