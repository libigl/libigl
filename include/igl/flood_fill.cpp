// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "flood_fill.h"
#include <limits>

IGL_INLINE void igl::flood_fill(
  const Eigen::RowVector3i & res,
  Eigen::VectorXd & S)
{
  using namespace Eigen;
  using namespace std;
  const auto flood = [&res,&S] (
     const int xi, 
     const int yi, 
     const int zi,
     const int signed_xi, 
     const int signed_yi, 
     const int signed_zi,
     const double s)
    {
      // flood fill this value back on this row
      for(int bxi = xi;signed_xi<--bxi;)
      {
        S(bxi+res(0)*(yi + res(1)*zi)) = s;
      }
      // flood fill this value back on any previous rows
      for(int byi = yi;signed_yi<--byi;)
      {
        for(int xi = 0;xi<res(0);xi++)
        {
          S(xi+res(0)*(byi + res(1)*zi)) = s;
        }
      }
      // flood fill this value back on any previous "sheets"
      for(int bzi = zi;signed_zi<--bzi;)
      {
        for(int yi = 0;yi<res(1);yi++)
        {
          for(int xi = 0;xi<res(0);xi++)
          {
            S(xi+res(0)*(yi + res(1)*bzi)) = s;
          }
        }
      }
    };
  int signed_zi = -1;
  double s = numeric_limits<double>::quiet_NaN();
  for(int zi = 0;zi<res(2);zi++)
  {
    int signed_yi = -1;
    if(zi != 0)
    {
      s = S(0+res(0)*(0 + res(1)*(zi-1)));
    }
    for(int yi = 0;yi<res(1);yi++)
    {
      // index of last signed item on this row
      int signed_xi = -1;
      if(yi != 0)
      {
        s = S(0+res(0)*(yi-1 + res(1)*zi));
      }
      for(int xi = 0;xi<res(0);xi++)
      {
        int i = xi+res(0)*(yi + res(1)*zi);
        if(S(i)!=S(i))
        {
          if(s == s)
          {
            S(i) = s;
          }
          continue;
        }
        s = S(i);
        flood(xi,yi,zi,signed_xi,signed_yi,signed_zi,s);
        // remember this position
        signed_xi = xi;
        signed_yi = yi;
        signed_zi = zi;
      }
    }
  }
}
