// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "pseudonormal_test.h"
#include "barycentric_coordinates.h"
#include "doublearea.h"
#include "project_to_line_segment.h"
#include <cassert>

IGL_INLINE void igl::pseudonormal_test(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & FN,
  const Eigen::MatrixXd & VN,
  const Eigen::MatrixXd & EN,
  const Eigen::VectorXi & EMAP,
  const Eigen::RowVector3d & q,
  const int f,
  const Eigen::RowVector3d & c,
  double & s,
  Eigen::RowVector3d & n)
{
  using namespace Eigen;
  const auto & qc = q-c;
  RowVector3d b;
  // Using barycentric coorindates to determine whether close to a vertex/edge
  // seems prone to error when dealing with nearly degenerate triangles: Even
  // the barycenter (1/3,1/3,1/3) can be made arbitrarily close to an
  // edge/vertex
  //
  const RowVector3d A = V.row(F(f,0));
  const RowVector3d B = V.row(F(f,1));
  const RowVector3d C = V.row(F(f,2));

  const double area = [&A,&B,&C]()
  {
    Matrix<double,1,1> area;
    doublearea(A,B,C,area);
    return area(0);
  }();
  // These were chosen arbitrarily. In a floating point scenario, I'm not sure
  // the best way to determine if c is on a vertex/edge or in the middle of the
  // face: specifically, I'm worrying about degenerate triangles where
  // barycentric coordinates are error-prone.
  const double MIN_DOUBLE_AREA = 1e-4;
  const double epsilon = 1e-12;
  if(area>MIN_DOUBLE_AREA)
  {
    barycentric_coordinates( c,A,B,C,b);
    // Determine which normal to use
    const int type = (b.array()<=epsilon).cast<int>().sum();
    switch(type)
    {
      case 2:
        // Find vertex
        for(int x = 0;x<3;x++)
        {
          if(b(x)>epsilon)
          {
            n = VN.row(F(f,x));
            break;
          }
        }
        break;
      case 1:
        // Find edge
        for(int x = 0;x<3;x++)
        {
          if(b(x)<=epsilon)
          {
            n = EN.row(EMAP(F.rows()*x+f));
            break;
          }
        }
        break;
      default:
        assert(false && "all barycentric coords zero.");
      case 0:
        n = FN.row(f);
        break;
    }
  }else
  {
    // Check each vertex
    bool found = false;
    for(int v = 0;v<3 && !found;v++)
    {
      if( (c-V.row(F(f,v))).norm() < epsilon)
      {
        found = true;
        n = VN.row(F(f,v));
      }
    }
    // Check each edge
    for(int e = 0;e<3 && !found;e++)
    {
      const RowVector3d s = V.row(F(f,(e+1)%3));
      const RowVector3d d = V.row(F(f,(e+2)%3));
      Matrix<double,1,1> sqr_d_j_x(1,1);
      Matrix<double,1,1> t(1,1);
      project_to_line_segment(c,s,d,t,sqr_d_j_x);
      if(sqrt(sqr_d_j_x(0)) < epsilon)
      {
        n = EN.row(EMAP(F.rows()*e+f));
        found = true;
      }
    }
    // Finally just use face
    if(!found)
    {
      n = FN.row(f);
    }
  }
  s = (qc.dot(n) >= 0 ? 1. : -1.);
}

IGL_INLINE void igl::pseudonormal_test(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const Eigen::MatrixXd & EN,
  const Eigen::MatrixXd & VN,
  const Eigen::RowVector2d & q,
  const int e,
  const Eigen::RowVector2d & c,
  double & s,
  Eigen::RowVector2d & n)
{
  using namespace Eigen;
  const auto & qc = q-c;
  const double len = (V.row(E(e,1))-V.row(E(e,0))).norm();
  // barycentric coordinates
  RowVector2d b((c-V.row(E(e,1))).norm()/len,(c-V.row(E(e,0))).norm()/len);
  // Determine which normal to use
  const double epsilon = 1e-12;
  const int type = (b.array()<=epsilon).cast<int>().sum();
  switch(type)
  {
    case 1:
      // Find vertex
      for(int x = 0;x<2;x++)
      {
        if(b(x)>epsilon)
        {
          n = VN.row(E(e,x));
          break;
        }
      }
      break;
    default:
      assert(false && "all barycentric coords zero.");
    case 0:
      n = EN.row(e);
      break;
  }
  s = (qc.dot(n) >= 0 ? 1. : -1.);
}
