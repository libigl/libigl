// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Gustavo Segovia <gustavo.segovia@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "writeCSV.h"

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>

#include <vector>

template <typename Scalar>
IGL_INLINE bool igl::writeCSV(
  const std::string str, 
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& M)
{
  using namespace std;
  using namespace igl;
  std::ofstream s(str.c_str());
  
  if(!s.is_open())
  {
    fprintf(stderr,"IOError: writeCSV() could not open %s\n",str.c_str());
    return false;
  }
  
  s << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  for(int i=0;i<(int)M.rows();++i) {
    for(int j=0;j<(int)M.cols();++j) {
      if (j != 0) s << ", ";
      s << M(i,j);
    }
    s << std::endl;
  }
  s.close();
  return true;
}
