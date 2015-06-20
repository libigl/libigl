// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "winding_number.h"
#include "WindingNumberAABB.h"

#include <igl/PI.h>
#include <cmath>

// IF THIS IS EVER TEMPLATED BE SURE THAT V IS COLMAJOR
IGL_INLINE void igl::winding_number(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & O,
  Eigen::VectorXd & W)
{
  using namespace Eigen;
  // make room for output
  W.resize(O.rows(),1);
  switch(F.cols())
  {
    case 2:
      return winding_number_2(
        V.data(),
        V.rows(),
        F.data(),
        F.rows(),
        O.data(),
        O.rows(),
        W.data());
    case 3:
    {
      WindingNumberAABB<Vector3d> hier(V,F);
      hier.grow();
      // loop over origins
      const int no = O.rows();
#   pragma omp parallel for if (no>IGL_WINDING_NUMBER_OMP_MIN_VALUE)
      for(int o = 0;o<no;o++)
      {
        Vector3d p = O.row(o);
        W(o) = hier.winding_number(p);
      }
      break;
    }
    default: assert(false && "Bad simplex size"); break;
  }
}

template <typename Scalar, typename DerivedF>
IGL_INLINE void igl::winding_number_3(
  const Scalar * V,
  const int n,
  const DerivedF * F,
  const int m,
  const Scalar * O,
  const int no,
  Scalar * S)
{
  // Initialize output
  // loop over origins
#pragma omp parallel for if (no>IGL_WINDING_NUMBER_OMP_MIN_VALUE)
  for(int o = 0;o<no;o++)
  {
    S[o] = 0;
  }
  // Only use parallel for if there are many facets and more than one origin.
  // Assumes that if there is exactly one origin then this is being called
  // within an outer for loop which may be parallel
#pragma omp parallel for if (m>IGL_WINDING_NUMBER_OMP_MIN_VALUE && no>1)
  // loop over faces
  for(int f = 0;f<m;f++)
  {
    // Gather corners 
    Scalar C[3][3];
    // loop around triangle
    for(int t=0;t<3;t++)
    {
      // loop over dimensions
      for(int d = 0;d<3;d++)
      {
        // Indices are offset by 1
        int Ff = F[m*t + f];
        C[t][d] = V[d*n + Ff];
      }
    }
    // loop over origins
    for(int o = 0;o<no;o++)
    {
      // Gather vectors to corners
      Scalar v[3][3];
      Scalar vl[3];
      // loop around triangle
      for(int t=0;t<3;t++)
      {
        vl[t] = 0;
        // loop over dimensions
        for(int d = 0;d<3;d++)
        {
          v[t][d] = C[t][d] - O[d*no + o];
          // compute edge length contribution
          vl[t] += v[t][d]*v[t][d];
        }
        // finish edge length computation
        // Matlab crashes on NaN
        if(vl[t]!=0)
        {
          vl[t] = sqrt(vl[t]);
        }
      }
      //printf("\n");
      // Compute determinant
      Scalar detf = 
        v[0][0]*v[1][1]*v[2][2]+
        v[1][0]*v[2][1]*v[0][2]+
        v[2][0]*v[0][1]*v[1][2]-
        v[2][0]*v[1][1]*v[0][2]-
        v[1][0]*v[0][1]*v[2][2]-
        v[0][0]*v[2][1]*v[1][2];
      // Compute pairwise dotproducts
      Scalar dp[3];
      dp[0] = v[1][0]*v[2][0];
      dp[0] += v[1][1]*v[2][1];
      dp[0] += v[1][2]*v[2][2];
      dp[1] = v[2][0]*v[0][0];
      dp[1] += v[2][1]*v[0][1];
      dp[1] += v[2][2]*v[0][2];
      dp[2] = v[0][0]*v[1][0];
      dp[2] += v[0][1]*v[1][1];
      dp[2] += v[0][2]*v[1][2];
      // Compute winding number
      // Only divide by TWO_PI instead of 4*pi because there was a 2 out front
      Scalar val = atan2(detf,
        vl[0]*vl[1]*vl[2] + 
        dp[0]*vl[0] +
        dp[1]*vl[1] +
        dp[2]*vl[2]) / (2.*igl::PI);
#pragma omp atomic
      S[o] += val;
    }
  }
}

template <typename DerivedF>
IGL_INLINE void igl::winding_number_2(
  const double * V,
  const int n,
  const DerivedF * F,
  const int m,
  const double * O,
  const int no,
  double * S)
{
  // Initialize output
  // loop over origins
#pragma omp parallel for if (no>IGL_WINDING_NUMBER_OMP_MIN_VALUE)
  for(int o = 0;o<no;o++)
  {
    S[o] = 0;
  }
#pragma omp parallel for if (m>IGL_WINDING_NUMBER_OMP_MIN_VALUE && no>1)
  // loop over faces
  for(int f = 0;f<m;f++)
  {
    // Index of source and destination
    int s = F[m*0 + f];
    int d = F[m*1 + f];
    // Positions of source and destination
    double vs[2];
    double vd[2];
    vs[0] = V[0*n + s];
    vs[1] = V[1*n + s];
    vd[0] = V[0*n + d];
    vd[1] = V[1*n + d];
    
    // loop over origins
    for(int o = 0;o<no;o++)
    {
      // Gather vectors to source and destination
      double o2vs[2];
      double o2vd[2];
      // and lengths
      double o2vsl = 0;
      double o2vdl = 0;
      for(int i = 0;i<2;i++)
      {
        o2vs[i] = O[i*no + o] - vs[i];
        o2vd[i] = O[i*no + o] - vd[i];
        o2vsl += o2vs[i]*o2vs[i];
        o2vdl += o2vd[i]*o2vd[i];
      }
      o2vsl = sqrt(o2vsl);
      o2vdl = sqrt(o2vdl);
      // Normalize
      for(int i = 0;i<2;i++)
      {
        // Matlab crashes on NaN
        if(o2vsl!=0)
        {
          o2vs[i] /= o2vsl;
        }
        if(o2vdl!=0)
        {
          o2vd[i] /= o2vdl;
        }
      }
      double val =
        atan2(o2vd[0]*o2vs[1]-o2vd[1]*o2vs[0],o2vd[0]*o2vs[0]+o2vd[1]*o2vs[1])/
        (2.*igl::PI);
#pragma omp atomic
      S[o] += val;
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::winding_number_2<double>(double const*, int, double const*, int, double const*, int, double*);
template void igl::winding_number_3<double>(double const*, int, double const*, int, double const*, int, double*);
template void igl::winding_number_3<double, int>(double const*, int, int const*, int, double const*, int, double*);
#endif
