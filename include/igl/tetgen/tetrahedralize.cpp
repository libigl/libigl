// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "tetrahedralize.h"
#include "mesh_to_tetgenio.h"
#include "tetgenio_to_tetmesh.h"

// IGL includes 
#include <igl/matrix_to_list.h>
#include <igl/list_to_matrix.h>
#include <igl/boundary_facets.h>

// STL includes
#include <cassert>
#include <iostream>

IGL_INLINE int igl::tetgen::tetrahedralize(
  const std::vector<std::vector<REAL > > & V, 
  const std::vector<std::vector<int> > & F, 
  const std::string switches,
  std::vector<std::vector<REAL > > & TV, 
  std::vector<std::vector<int > > & TT, 
  std::vector<std::vector<int> > & TF)
{
  using namespace std;
  tetgenio in,out;
  bool success;
  success = mesh_to_tetgenio(V,F,in);
  if(!success)
  {
    return -1;
  }
  try
  {
    char * cswitches = new char[switches.size() + 1];
    std::strcpy(cswitches,switches.c_str());
    ::tetrahedralize(cswitches,&in, &out);
    delete[] cswitches;
  }catch(int e)
  {
    cerr<<"^"<<__FUNCTION__<<": TETGEN CRASHED... KABOOOM!!!"<<endl;
    return 1;
  }
  if(out.numberoftetrahedra == 0)
  {
    cerr<<"^"<<__FUNCTION__<<": Tetgen failed to create tets"<<endl;
    return 2;
  }
  success = tetgenio_to_tetmesh(out,TV,TT,TF);
  if(!success)
  {
    return -1;
  }
  //boundary_facets(TT,TF);
  return 0;
}

template <
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedTV, 
  typename DerivedTT, 
  typename DerivedTF>
IGL_INLINE int igl::tetgen::tetrahedralize(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  const std::string switches,
  Eigen::PlainObjectBase<DerivedTV>& TV,
  Eigen::PlainObjectBase<DerivedTT>& TT,
  Eigen::PlainObjectBase<DerivedTF>& TF)
{
  using namespace igl;
  using namespace std;
  vector<vector<REAL> > vV,vTV;
  vector<vector<int> > vF,vTT,vTF;
  matrix_to_list(V,vV);
  matrix_to_list(F,vF);
  int e = tetrahedralize(vV,vF,switches,vTV,vTT,vTF);
  if(e == 0)
  {
    bool TV_rect = list_to_matrix(vTV,TV);
    if(!TV_rect)
    {
      return false;
    }
    bool TT_rect = list_to_matrix(vTT,TT);
    if(!TT_rect)
    {
      return false;
    }
    bool TF_rect = list_to_matrix(vTF,TF);
    if(!TF_rect)
    {
      return false;
    }
  }
  return e;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template int igl::tetgen::tetrahedralize<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
