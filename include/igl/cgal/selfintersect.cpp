// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "selfintersect.h"
#include "SelfIntersectMesh.h"
#include <igl/C_STR.h>
#include <list>

IGL_INLINE void igl::selfintersect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const SelfintersectParam & params,
  Eigen::MatrixXd & VV,
  Eigen::MatrixXi & FF,
  Eigen::MatrixXi & IF,
  Eigen::VectorXi & J)
{
  using namespace std;
  if(params.detect_only)
  {
    //// This is probably a terrible idea, but CGAL is throwing floating point
    //// exceptions.

//#ifdef __APPLE__
//#define IGL_THROW_FPE 11
//    const auto & throw_fpe = [](int e)
//    {
//      throw "IGL_THROW_FPE";
//    };
//    signal(SIGFPE,throw_fpe);
//#endif

    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    SelfIntersectMesh<Kernel> SIM = SelfIntersectMesh<Kernel>(V,F,params,VV,FF,IF,J);

//#ifdef __APPLE__
//    signal(SIGFPE,SIG_DFL);
//#endif

  }else
  {
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    SelfIntersectMesh<Kernel> SIM = SelfIntersectMesh<Kernel>(V,F,params,VV,FF,IF,J);
  }
}
