#ifndef IGL_PI_H
#define IGL_PI_H
namespace igl
{
  // Use standard mathematical constants' M_PI if available
#ifdef M_PI
  const double PI = M_PI;
#else
  const double PI = 3.1415926535897932384626433832795;
#endif
}
#endif
