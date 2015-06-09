// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "median.h"
#include "matrix_to_list.h"

#include <vector>
#include <algorithm>

IGL_INLINE bool igl::median(const Eigen::VectorXd & V, double & m)
{
  using namespace std;
  if(V.size() == 0)
  {
    return false;
  }
  vector<double> vV;
  matrix_to_list(V,vV);
  // http://stackoverflow.com/a/1719155/148668
  size_t n = vV.size()/2;
  nth_element(vV.begin(),vV.begin()+n,vV.end());
  if(vV.size()%2==0)
  {
    nth_element(vV.begin(),vV.begin()+n-1,vV.end());
    m = 0.5*(vV[n]+vV[n-1]);
  }else
  {
    m = vV[n];
  }
  return true;
}
