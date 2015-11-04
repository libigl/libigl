// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "remesh_intersections.h"
#include "SelfIntersectMesh.h"
#include "assign_scalar.h"
#include "projected_delaunay.h"
#include <iostream>
#include <cassert>

template <
  typename DerivedV,
  typename DerivedF,
  typename Kernel,
  typename DerivedVV,
  typename DerivedFF,
  typename DerivedJ,
  typename DerivedIM>
IGL_INLINE void igl::cgal::remesh_intersections(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const std::vector<CGAL::Triangle_3<Kernel> > & T,
  const std::map<
    typename DerivedF::Index,
    std::pair<typename DerivedF::Index,
      std::vector<CGAL::Object> > > & offending,
  const std::map<
    std::pair<typename DerivedF::Index,typename DerivedF::Index>,
    std::vector<typename DerivedF::Index> > & edge2faces,
  Eigen::PlainObjectBase<DerivedVV> & VV,
  Eigen::PlainObjectBase<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedJ> & J,
  Eigen::PlainObjectBase<DerivedIM> & IM)
{
  using namespace std;
  using namespace Eigen;
  typedef typename DerivedF::Index          Index;
  typedef CGAL::Point_3<Kernel>    Point_3;
  //typedef CGAL::Segment_3<Kernel>  Segment_3; 
  //typedef CGAL::Triangle_3<Kernel> Triangle_3; 
  typedef CGAL::Plane_3<Kernel>    Plane_3;
  //typedef CGAL::Point_2<Kernel>    Point_2;
  //typedef CGAL::Segment_2<Kernel>  Segment_2; 
  //typedef CGAL::Triangle_2<Kernel> Triangle_2; 
  typedef CGAL::Triangulation_vertex_base_2<Kernel>  TVB_2;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel> CTFB_2;
  typedef CGAL::Triangulation_data_structure_2<TVB_2,CTFB_2> TDS_2;
  typedef CGAL::Exact_intersections_tag Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,TDS_2,Itag> 
    CDT_2;
  typedef CGAL::Constrained_triangulation_plus_2<CDT_2> CDT_plus_2;
  typedef std::pair<Index,Index> EMK;
  typedef std::vector<Index> EMV;
  //typedef std::map<EMK,EMV> EdgeMap;
  typedef std::pair<Index,Index> EMK;
  //typedef std::vector<CGAL::Object> ObjectList;
  typedef std::vector<Index> IndexList;

  int NF_count = 0;
  // list of new faces, we'll fix F later
  vector<
    typename Eigen::Matrix<typename DerivedFF::Scalar,Dynamic,Dynamic>
    > NF(offending.size());
  // list of new vertices
  typedef vector<Point_3> Point_3List;
  Point_3List NV;
  Index NV_count = 0;
  vector<CDT_plus_2> cdt(offending.size());
  vector<Plane_3> P(offending.size());
  // Use map for *all* faces
  map<typename CDT_plus_2::Vertex_handle,Index> v2i;
#ifdef IGL_SELFINTERSECTMESH_DEBUG
  double t_proj_del = 0;
#endif
  // Unfortunately it looks like CGAL has trouble allocating memory when
  // multiple openmp threads are running. Crashes durring CDT...
  //// Loop over offending triangles
  //const size_t noff = offending.size();
//# pragma omp parallel for if (noff>1000)
  for(const auto & okv : offending)
  {
    // index in F
    const Index f = okv.first;
    const Index o = okv.second.first;
    {
#ifdef IGL_SELFINTERSECTMESH_DEBUG
      const double t_before = get_seconds();
#endif
      CDT_plus_2 cdt_o;
      projected_delaunay(T[f],okv.second.second,cdt_o);
      cdt[o] = cdt_o;
#ifdef IGL_SELFINTERSECTMESH_DEBUG
      t_proj_del += (get_seconds()-t_before);
#endif
    }
    // Q: Is this also delaunay in 3D?
    // A: No, because the projection is affine and delaunay is not affine
    // invariant.
    // Q: Then, can't we first get the 2D delaunay triangulation, then lift it
    // to 3D and flip any offending edges?
    // Plane of projection (also used by projected_delaunay)
    P[o] = Plane_3(T[f].vertex(0),T[f].vertex(1),T[f].vertex(2));
    // Build index map
    {
      Index i=0;
      for(
        typename CDT_plus_2::Finite_vertices_iterator vit = cdt[o].finite_vertices_begin();
        vit != cdt[o].finite_vertices_end();
        ++vit)
      {
        if(i<3)
        {
          //cout<<T[f].vertex(i)<<
          //  (T[f].vertex(i) == P[o].to_3d(vit->point())?" == ":" != ")<<
          //  P[o].to_3d(vit->point())<<endl;
#ifndef NDEBUG
          // I want to be sure that the original corners really show up as the
          // original corners of the CDT. I.e. I don't trust CGAL to maintain
          // the order
          assert(T[f].vertex(i) == P[o].to_3d(vit->point()));
#endif
          // For first three, use original index in F
//#   pragma omp critical
          v2i[vit] = F(f,i);
        }else
        {
          const Point_3 vit_point_3 = P[o].to_3d(vit->point());
          // First look up each edge's neighbors to see if exact point has
          // already been added (This makes everything a bit quadratic)
          bool found = false;
          for(int e = 0; e<3 && !found;e++)
          {
            // Index of F's eth edge in V
            Index i = F(f,(e+1)%3);
            Index j = F(f,(e+2)%3);
            // Be sure that i<j
            if(i>j)
            {
              swap(i,j);
            }
            assert(edge2faces.count(EMK(i,j))==1);
            const EMV & facesij = edge2faces.find(EMK(i,j))->second;
            // loop over neighbors
            for(
              typename IndexList::const_iterator nit = facesij.begin();
              nit != facesij.end() && !found;
              nit++)
            {
              // don't consider self
              if(*nit == f)
              {
                continue;
              }
              // index of neighbor in offending (to find its cdt)
              assert(offending.count(*nit) == 1);
              Index no = offending.find(*nit)->second.first;
              // Loop over vertices of that neighbor's cdt (might not have been
              // processed yet, but then it's OK because it'll just be empty)
              for(
                  typename CDT_plus_2::Finite_vertices_iterator uit = cdt[no].finite_vertices_begin();
                  uit != cdt[no].finite_vertices_end() && !found;
                  ++uit)
              {
                if(vit_point_3 == P[no].to_3d(uit->point()))
                {
                  assert(v2i.count(uit) == 1);
//#   pragma omp critical
                  v2i[vit] = v2i[uit];
                  found = true;
                }
              }
            }
          }
          if(!found)
          {
//#   pragma omp critical
            {
              v2i[vit] = V.rows()+NV_count;
              NV.push_back(vit_point_3);
              NV_count++;
            }
          }
        }
        i++;
      }
    }
    {
      Index i = 0;
      // Resize to fit new number of triangles
      NF[o].resize(cdt[o].number_of_faces(),3);
//#   pragma omp atomic
      NF_count+=NF[o].rows();
      // Append new faces to NF
      for(
        typename CDT_plus_2::Finite_faces_iterator fit = cdt[o].finite_faces_begin();
        fit != cdt[o].finite_faces_end();
        ++fit)
      {
        NF[o](i,0) = v2i[fit->vertex(0)];
        NF[o](i,1) = v2i[fit->vertex(1)];
        NF[o](i,2) = v2i[fit->vertex(2)];
        i++;
      }
    }
  }
#ifdef IGL_SELFINTERSECTMESH_DEBUG
  cout<<"CDT: "<<tictoc()<<"  "<<t_proj_del<<endl;
#endif

  assert(NV_count == (Index)NV.size());
  // Build output
#ifndef NDEBUG
  //{
  //  Index off_count = 0;
  //  for(Index f = 0;f<F.rows();f++)
  //  {
  //    off_count+= (offensive[f]?1:0);
  //  }
  //  assert(off_count==(Index)offending.size());
  //}
#endif
  // Append faces
  FF.resize(F.rows()-offending.size()+NF_count,3);
  J.resize(FF.rows());
  // First append non-offending original faces
  // There's an Eigen way to do this in one line but I forget
  Index off = 0;
  for(Index f = 0;f<F.rows();f++)
  {
    if(!offending.count(f))
    {
      FF.row(off) = F.row(f);
      J(off) = f;
      off++;
    }
  }
  assert(off == (Index)(F.rows()-offending.size()));
  // Now append replacement faces for offending faces
  for(const auto & okv : offending)
  {
    // index in F
    const Index f = okv.first;
    const Index o = okv.second.first;
    FF.block(off,0,NF[o].rows(),3) = NF[o];
    J.block(off,0,NF[o].rows(),1).setConstant(f);
    off += NF[o].rows();
  }
  // Append vertices
  VV.resize(V.rows()+NV_count,3);
  VV.block(0,0,V.rows(),3) = V.template cast<typename DerivedVV::Scalar>();
  {
    Index i = 0;
    for(
      typename Point_3List::const_iterator nvit = NV.begin();
      nvit != NV.end();
      nvit++)
    {
      for(Index d = 0;d<3;d++)
      {
        const Point_3 & p = *nvit;
        // Don't convert via double if output type is same as Kernel
        assign_scalar(p[d], VV(V.rows()+i,d));
      }
      i++;
    }
  }
  IM.resize(VV.rows(),1);
  map<Point_3,Index> vv2i;
  // Safe to check for duplicates using double for original vertices: if
  // incoming reps are different then the points are unique.
  for(Index v = 0;v<V.rows();v++)
  {
    typename Kernel::FT p0,p1,p2;
    assign_scalar(V(v,0),p0);
    assign_scalar(V(v,1),p1);
    assign_scalar(V(v,2),p2);
    const Point_3 p(p0,p1,p2);
    if(vv2i.count(p)==0)
    {
      vv2i[p] = v;
    }
    assert(vv2i.count(p) == 1);
    IM(v) = vv2i[p];
  }
  // Must check for duplicates of new vertices using exact.
  {
    Index v = V.rows();
    for(
      typename Point_3List::const_iterator nvit = NV.begin();
      nvit != NV.end();
      nvit++)
    {
      const Point_3 & p = *nvit;
      if(vv2i.count(p)==0)
      {
        vv2i[p] = v;
      }
      assert(vv2i.count(p) == 1);
      IM(v) = vv2i[p];
      v++;
    }
  }
#ifdef IGL_SELFINTERSECTMESH_DEBUG
  cout<<"Output + dupes: "<<tictoc()<<endl;
#endif
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epeck, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, CGAL::Epick, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index const, std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index>, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> const, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epeck, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epeck>, std::allocator<CGAL::Triangle_3<CGAL::Epeck> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
template void igl::cgal::remesh_intersections<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, CGAL::Epick, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, std::vector<CGAL::Triangle_3<CGAL::Epick>, std::allocator<CGAL::Triangle_3<CGAL::Epick> > > const&, std::map<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > >, std::less<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::allocator<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index const, std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::vector<CGAL::Object, std::allocator<CGAL::Object> > > > > > const&, std::map<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index>, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::less<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> >, std::allocator<std::pair<std::pair<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> const, std::vector<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index, std::allocator<Eigen::Matrix<int, -1, 3, 0, -1, 3>::Index> > > > > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
#endif
