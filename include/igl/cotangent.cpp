// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "cotangent.h"
#include "doublearea.h"
#include "edge_lengths.h"
#include "face_areas.h"
#include "volume.h"
#include "dihedral_angles.h"

#include "verbose.h"


template <typename DerivedV, typename DerivedF, typename DerivedC>
IGL_INLINE void igl::cotangent(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedC>& C)
{
  using namespace igl;
  using namespace std;
  using namespace Eigen;
  // simplex size (3: triangles, 4: tetrahedra)
  int simplex_size = F.cols();
  // Number of elements
  int m = F.rows();

  // Law of cosines + law of sines
  switch(simplex_size)
  {
    case 3:
    {
      // Triangles
      //Matrix<typename DerivedC::Scalar,Dynamic,3> l;
      //edge_lengths(V,F,l);
      // edge lengths numbered same as opposite vertices
      Matrix<typename DerivedC::Scalar,Dynamic,3> l;
      igl::edge_lengths(V,F,l);
      // double area
      Matrix<typename DerivedC::Scalar,Dynamic,1> dblA;
      doublearea(l,dblA);
      // cotangents and diagonal entries for element matrices
      // correctly divided by 4 (alec 2010)
      C.resize(m,3);
      for(int i = 0;i<m;i++)
      {
        C(i,0) = (l(i,1)*l(i,1) + l(i,2)*l(i,2) - l(i,0)*l(i,0))/dblA(i)/4.0;
        C(i,1) = (l(i,2)*l(i,2) + l(i,0)*l(i,0) - l(i,1)*l(i,1))/dblA(i)/4.0;
        C(i,2) = (l(i,0)*l(i,0) + l(i,1)*l(i,1) - l(i,2)*l(i,2))/dblA(i)/4.0;
      }
      break;
    }
    case 4:
    {

      // edge lengths numbered same as opposite vertices
      Matrix<typename DerivedC::Scalar,Dynamic,6> l;
      edge_lengths(V,F,l);
      Matrix<typename DerivedC::Scalar,Dynamic,4> s;
      face_areas(l,s);
      Matrix<typename DerivedC::Scalar,Dynamic,6> cos_theta,theta;
      dihedral_angles_intrinsic(l,s,theta,cos_theta);

      // volume
      Matrix<typename DerivedC::Scalar,Dynamic,1> vol;
      volume(l,vol);


      // Law of sines
      // http://mathworld.wolfram.com/Tetrahedron.html
      Matrix<typename DerivedC::Scalar,Dynamic,6> sin_theta(m,6);
      sin_theta.col(0) = vol.array() / ((2./(3.*l.col(0).array())).array() * s.col(1).array() * s.col(2).array());
      sin_theta.col(1) = vol.array() / ((2./(3.*l.col(1).array())).array() * s.col(2).array() * s.col(0).array());
      sin_theta.col(2) = vol.array() / ((2./(3.*l.col(2).array())).array() * s.col(0).array() * s.col(1).array());
      sin_theta.col(3) = vol.array() / ((2./(3.*l.col(3).array())).array() * s.col(3).array() * s.col(0).array());
      sin_theta.col(4) = vol.array() / ((2./(3.*l.col(4).array())).array() * s.col(3).array() * s.col(1).array());
      sin_theta.col(5) = vol.array() / ((2./(3.*l.col(5).array())).array() * s.col(3).array() * s.col(2).array());


      // http://arxiv.org/pdf/1208.0354.pdf Page 18
      C = (1./6.) * l.array() * cos_theta.array() / sin_theta.array();

// LEGACY
//      // Tetrahedra
//      typedef Matrix<typename MatV::Scalar,3,1> Vec3;
//      typedef Matrix<typename MatV::Scalar,3,3> Mat3;
//      typedef Matrix<typename MatV::Scalar,3,4> Mat3x4;
//      typedef Matrix<typename MatV::Scalar,4,4> Mat4x4;
//
//      // preassemble right hand side
//      // COLUMN-MAJOR ORDER FOR LAPACK
//      Mat3x4 rhs;
//      rhs <<
//        1,0,0,-1,
//        0,1,0,-1,
//        0,0,1,-1;
//
//      bool diag_all_pos = true;
//      C.resize(m,6);
//
//      // loop over tetrahedra
//      for(int j = 0;j<F.rows();j++)
//      {
//        // points a,b,c,d make up the tetrahedra
//        size_t a = F(j,0);
//        size_t b = F(j,1);
//        size_t c = F(j,2);
//        size_t d = F(j,3);
//        //const std::vector<double> & pa = vertices[a];
//        //const std::vector<double> & pb = vertices[b];
//        //const std::vector<double> & pc = vertices[c];
//        //const std::vector<double> & pd = vertices[d];
//        Vec3 pa = V.row(a);
//        Vec3 pb = V.row(b);
//        Vec3 pc = V.row(c);
//        Vec3 pd = V.row(d);
//
//        // Following definition that appears in the appendix of: ``Interactive
//        // Topology-aware Surface Reconstruction,'' by Sharf, A. et al
//        // http://www.cs.bgu.ac.il/~asharf/Projects/InSuRe/Insure_siggraph_final.pdf
//
//        // compute transpose of jacobian Jj
//        Mat3 JTj;
//        JTj.row(0) = pa-pd;
//        JTj.row(1) = pb-pd;
//        JTj.row(2) = pc-pd;
//
//        // compute abs(determinant of JTj)/6 (volume of tet)
//        // determinant of transpose of A equals determinant of A
//        double volume = fabs(JTj.determinant())/6.0;
//        //printf("volume[%d] = %g\n",j+1,volume);
//
//        // solve Jj' * Ej = [-I -1], for Ej
//        // in other words solve JTj * Ej = [-I -1], for Ej
//        Mat3x4 Ej = JTj.inverse() * rhs;
//        // compute Ej'*Ej
//        Mat4x4 EjTEj = Ej.transpose() * Ej;
//
//        // Kj =  det(JTj)/6 * Ej'Ej 
//        Mat4x4 Kj = EjTEj*volume;
//        diag_all_pos &= ((Kj(0,0)>0) & (Kj(1,1)>0)) & ((Kj(2,2)>0) & (Kj(3,3)>0));
//        C(j,0) = Kj(1,2);
//        C(j,1) = Kj(2,0);
//        C(j,2) = Kj(0,1);
//        C(j,3) = Kj(3,0);
//        C(j,4) = Kj(3,1);
//        C(j,5) = Kj(3,2);
//      }
//      if(diag_all_pos)
//      {
//#ifdef VERBOSE 
//        verbose("cotangent.h: Flipping sign of cotangent, so that cots are positive\n");
//#endif
//        C *= -1.0;
//      }
      break;
    }
    default:
    {
      fprintf(stderr,
          "cotangent.h: Error: Simplex size (%d) not supported\n", simplex_size);
      assert(false);
    }
  }
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
template void igl::cotangent<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
