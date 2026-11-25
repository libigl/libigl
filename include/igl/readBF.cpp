// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "readBF.h"
#include "list_to_matrix.h"
#include <vector>
#include <cstdio>
#include <fstream>
#include <cassert>
#include <functional>
template < 
  typename DerivedWI,
  typename DerivedP,
  typename DerivedO>
IGL_INLINE bool igl::readBF(
  const std::string & filename,
  Eigen::PlainObjectBase<DerivedWI> & WI,
  Eigen::PlainObjectBase<DerivedP> & P,
  Eigen::PlainObjectBase<DerivedO> & O)
{
  std::ifstream is(filename);
  if(!is.is_open())
  {
    return false;
  }
  std::string line;
  std::vector<typename DerivedWI::Scalar> vWI;
  std::vector<typename DerivedP::Scalar> vP;
  std::vector<std::vector<typename DerivedO::Scalar> > vO;
  while(getline(is, line))
  {
    int wi,p;
    double cx,cy,cz;
    if(sscanf(line.c_str(), "%d %d %lg %lg %lg",&wi,&p,&cx,&cy,&cz) != 5)
    {
      return false;
    }
    vWI.push_back(wi);
    vP.push_back(p);
    vO.push_back({cx,cy,cz});
  }
  list_to_matrix(vWI,WI);
  list_to_matrix(vP,P);
  list_to_matrix(vO,O);
  return true;
}

template < 
  typename DerivedWI,
  typename DerivedbfP,
  typename DerivedO,
  typename DerivedC,
  typename DerivedBE,
  typename DerivedP>
IGL_INLINE bool igl::readBF(
  const std::string & filename,
  Eigen::PlainObjectBase<DerivedWI> & WI,
  Eigen::PlainObjectBase<DerivedbfP> & bfP,
  Eigen::PlainObjectBase<DerivedO> & offsets,
  Eigen::PlainObjectBase<DerivedC> & C,
  Eigen::PlainObjectBase<DerivedBE> & BE,
  Eigen::PlainObjectBase<DerivedP> & P)
{
  if(!readBF(filename,WI,bfP,offsets))
  {
    return false;
  }

  C.resize(WI.rows(),3);
  std::vector<bool> computed(C.rows(),false);
  // better not be cycles in bfP
  std::function<Eigen::RowVector3d(const int)> locate_tip;
  locate_tip = 
    [&offsets,&computed,&bfP,&locate_tip,&C](const int w)->Eigen::RowVector3d
  {
    if(w<0) return Eigen::RowVector3d(0,0,0);
    if(computed[w]) return C.row(w);
    computed[w] = true;
    return C.row(w) = locate_tip(bfP(w)) + offsets.row(w);
  };
  int num_roots = (bfP.array() == -1).count();
  BE.resize(WI.rows()-num_roots,2);
  P.resize(BE.rows());
  // The parent indices bfP and weight indices WI come straight from the file
  // and are used below to index into C/offsets (via locate_tip), WI, and the
  // bone arrays P and BE. Validate their range first: a root has parent -1,
  // any other parent must reference an existing joint, and a weight index is
  // either -1 (no weight) or an index into P/BE. Roots cannot have weights:
  // a weighted bone's parent is used as an index into WI.
  for(int c = 0;c<WI.rows();c++)
  {
    if(bfP(c) < -1 || bfP(c) >= WI.rows())
    {
      fprintf(stderr,"IOError: readBF() parent index %d out of range\n",(int)bfP(c));
      return false;
    }
    if(WI(c) >= BE.rows())
    {
      fprintf(stderr,"IOError: readBF() weight index %d out of range\n",(int)WI(c));
      return false;
    }
    if(WI(c) >= 0 && bfP(c) == -1)
    {
      fprintf(stderr,"IOError: readBF() root %d has weight index %d\n",c,(int)WI(c));
      return false;
    }
  }
  for(int c = 0;c<C.rows();c++)
  {
    locate_tip(c);
    assert(c>=0);
    // weight associated with this bone
    const int wi = WI(c);
    if(wi >= 0)
    {
      // index into C
      const int p = bfP(c);
      assert(p >= 0 && "No weights for roots allowed");
      // index into BE
      const int pwi = WI(p);
      P(wi) = pwi;
      BE(wi,0) = p;
      BE(wi,1) = c;
    }
  }
  return true;
}

#ifdef IGL_STATIC_LIBRARY
template bool igl::readBF<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
