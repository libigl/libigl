// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "writeWRL.h"
#include <iostream>
#include <fstream>
template <typename DerivedV, typename DerivedF>
IGL_INLINE bool igl::writeWRL(
  const std::string & str,
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F)
{
  using namespace std;
  using namespace Eigen;
  assert(V.cols() == 3 && "V should have 3 columns");
  assert(F.cols() == 3 && "F should have 3 columns");
  ofstream s(str);
  if(!s.is_open())
  {
    cerr<<"IOError: writeWRL() could not open "<<str<<endl;
    return false;
  }
  // Append column of -1 to F
  Matrix<typename DerivedF::Scalar,Dynamic,4> FF(F.rows(),4);
  FF.leftCols(3) = F;
  FF.col(3).setConstant(-1);

  s<<R"(
DEF default Transform {
translation 0 0 0
children [
Shape {
geometry DEF default-FACES IndexedFaceSet {
ccw TRUE
)"<<
    V.format(
      IOFormat(
        FullPrecision,
        DontAlignCols,
        " ",",\n","","",
        "coord DEF default-COORD Coordinate { point [ \n","]\n}\n"))<<
    FF.format(
      IOFormat(
        FullPrecision,
        DontAlignCols,
        ",","\n","","",
        "coordIndex [ \n"," ]\n"))<<
    "}\n}\n]\n}\n";
  return true;
}

#ifdef IGL_STATIC_LIBRARY
template bool igl::writeWRL<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
#endif
