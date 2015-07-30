// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "triangulate.h"
#ifdef ANSI_DECLARATORS
#  define IGL_PREVIOUSLY_DEFINED_ANSI_DECLARATORS ANSI_DECLARATORS
#  undef ANSI_DECLARATORS
#endif
#ifdef REAL
#  define IGL_PREVIOUSLY_DEFINED_REAL REAL
#  undef REAL
#endif
#ifdef VOID
#  define IGL_PREVIOUSLY_DEFINED_VOID VOID
#  undef VOID
#endif
#define ANSI_DECLARATORS
#define REAL double
#define VOID int

extern "C"
{
#include <triangle.h>
}

#undef ANSI_DECLARATORS
#ifdef IGL_PREVIOUSLY_DEFINED_ANSI_DECLARATORS
#  define ANSI_DECLARATORS IGL_PREVIOUSLY_DEFINED_ANSI_DECLARATORS
#endif

#undef REAL
#ifdef IGL_PREVIOUSLY_DEFINED_REAL
#  define REAL IGL_PREVIOUSLY_DEFINED_REAL
#endif

#undef VOID
#ifdef IGL_PREVIOUSLY_DEFINED_VOID
#  define VOID IGL_PREVIOUSLY_DEFINED_VOID
#endif

IGL_INLINE void igl::triangle::triangulate(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& E,
  const Eigen::MatrixXd& H,
  const std::string flags,
  Eigen::MatrixXd& V2,
  Eigen::MatrixXi& F2)
{
  using namespace std;
  using namespace Eigen;

  // Prepare the flags
  string full_flags = flags + "pzBV";

  // Prepare the input struct
  triangulateio in;

  assert(V.cols() == 2);

  in.numberofpoints = V.rows();
  in.pointlist = (double*)calloc(V.rows()*2,sizeof(double));
  for (unsigned i=0;i<V.rows();++i)
    for (unsigned j=0;j<2;++j)
      in.pointlist[i*2+j] = V(i,j);

  in.numberofpointattributes = 0;
  in.pointmarkerlist = (int*)calloc(V.rows(),sizeof(int));
   for (unsigned i=0;i<V.rows();++i)
     in.pointmarkerlist[i] = 1;

  in.trianglelist = NULL;
  in.numberoftriangles = 0;
  in.numberofcorners = 0;
  in.numberoftriangleattributes = 0;
  in.triangleattributelist = NULL;

  in.numberofsegments = E.rows();
  in.segmentlist = (int*)calloc(E.rows()*2,sizeof(int));
  for (unsigned i=0;i<E.rows();++i)
    for (unsigned j=0;j<2;++j)
      in.segmentlist[i*2+j] = E(i,j);
  in.segmentmarkerlist = (int*)calloc(E.rows(),sizeof(int));
  for (unsigned i=0;i<E.rows();++i)
    in.segmentmarkerlist[i] = 1;

  in.numberofholes = H.rows();
  in.holelist = (double*)calloc(H.rows()*2,sizeof(double));
  for (unsigned i=0;i<H.rows();++i)
    for (unsigned j=0;j<2;++j)
      in.holelist[i*2+j] = H(i,j);
  in.numberofregions = 0;

  // Prepare the output struct
  triangulateio out;

  out.pointlist = NULL;
  out.trianglelist = NULL;
  out.segmentlist = NULL;

  // Call triangle
  ::triangulate(const_cast<char*>(full_flags.c_str()), &in, &out, 0);

  // Return the mesh
  V2.resize(out.numberofpoints,2);
  for (unsigned i=0;i<V2.rows();++i)
    for (unsigned j=0;j<2;++j)
      V2(i,j) = out.pointlist[i*2+j];

  F2.resize(out.numberoftriangles,3);
  for (unsigned i=0;i<F2.rows();++i)
    for (unsigned j=0;j<3;++j)
      F2(i,j) = out.trianglelist[i*3+j];

  // Cleanup in
  free(in.pointlist);
  free(in.pointmarkerlist);
  free(in.segmentlist);
  free(in.segmentmarkerlist);
  free(in.holelist);

  // Cleanup out
  free(out.pointlist);
  free(out.trianglelist);
  free(out.segmentlist);

}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif
