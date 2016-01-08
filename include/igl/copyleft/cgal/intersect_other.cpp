// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "intersect_other.h"
#include "CGAL_includes.hpp"
#include "mesh_to_cgal_triangle_list.h"
#include "remesh_intersections.h"

#ifndef IGL_FIRST_HIT_EXCEPTION
#define IGL_FIRST_HIT_EXCEPTION 10
#endif

// Un-exposed helper functions
namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      template <typename DerivedF>
      static IGL_INLINE void push_result(
        const Eigen::PlainObjectBase<DerivedF> & F,
        const int f,
        const int f_other,
        const CGAL::Object & result,
        std::map<
          typename DerivedF::Index,
          std::vector<std::pair<typename DerivedF::Index, CGAL::Object> > > &
          offending,
        std::map<
          std::pair<typename DerivedF::Index,typename DerivedF::Index>,
          std::vector<typename DerivedF::Index> > & edge2faces)
      {
        typedef typename DerivedF::Index Index;
        typedef std::pair<Index,Index> EMK;
        if(offending.count(f) == 0)
        {
          // first time marking, initialize with new id and empty list
          offending[f] = {};
          for(Index e = 0; e<3;e++)
          {
            // append face to edge's list
            Index i = F(f,(e+1)%3) < F(f,(e+2)%3) ? F(f,(e+1)%3) : F(f,(e+2)%3);
            Index j = F(f,(e+1)%3) < F(f,(e+2)%3) ? F(f,(e+2)%3) : F(f,(e+1)%3);
            edge2faces[EMK(i,j)].push_back(f);
          }
        }
        offending[f].push_back({f_other,result});
      }
      template <
        typename Kernel,
        typename DerivedVA,
        typename DerivedFA,
        typename DerivedVB,
        typename DerivedFB,
        typename DerivedIF,
        typename DerivedVVA,
        typename DerivedFFA,
        typename DerivedJA,
        typename DerivedIMA,
        typename DerivedVVB,
        typename DerivedFFB,
        typename DerivedJB,
        typename DerivedIMB>
      static IGL_INLINE bool intersect_other_helper(
        const Eigen::PlainObjectBase<DerivedVA> & VA,
        const Eigen::PlainObjectBase<DerivedFA> & FA,
        const Eigen::PlainObjectBase<DerivedVB> & VB,
        const Eigen::PlainObjectBase<DerivedFB> & FB,
        const RemeshSelfIntersectionsParam & params,
        Eigen::PlainObjectBase<DerivedIF> & IF,
        Eigen::PlainObjectBase<DerivedVVA> & VVA,
        Eigen::PlainObjectBase<DerivedFFA> & FFA,
        Eigen::PlainObjectBase<DerivedJA>  & JA,
        Eigen::PlainObjectBase<DerivedIMA> & IMA,
        Eigen::PlainObjectBase<DerivedVVB> & VVB,
        Eigen::PlainObjectBase<DerivedFFB> & FFB,
        Eigen::PlainObjectBase<DerivedJB>  & JB,
        Eigen::PlainObjectBase<DerivedIMB> & IMB)
      {

        using namespace std;
        using namespace Eigen;

        typedef typename DerivedFA::Index Index;
        // 3D Primitives
        typedef CGAL::Point_3<Kernel>    Point_3;
        typedef CGAL::Segment_3<Kernel>  Segment_3; 
        typedef CGAL::Triangle_3<Kernel> Triangle_3; 
        typedef CGAL::Plane_3<Kernel>    Plane_3;
        typedef CGAL::Tetrahedron_3<Kernel> Tetrahedron_3; 
        // 2D Primitives
        typedef CGAL::Point_2<Kernel>    Point_2;
        typedef CGAL::Segment_2<Kernel>  Segment_2; 
        typedef CGAL::Triangle_2<Kernel> Triangle_2; 
        // 2D Constrained Delaunay Triangulation types
        typedef CGAL::Triangulation_vertex_base_2<Kernel>  TVB_2;
        typedef CGAL::Constrained_triangulation_face_base_2<Kernel> CTAB_2;
        typedef CGAL::Triangulation_data_structure_2<TVB_2,CTAB_2> TDS_2;
        typedef CGAL::Exact_intersections_tag Itag;
        typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,TDS_2,Itag> 
          CDT_2;
        typedef CGAL::Constrained_triangulation_plus_2<CDT_2> CDT_plus_2;
        // Axis-align boxes for all-pairs self-intersection detection
        typedef std::vector<Triangle_3> Triangles;
        typedef typename Triangles::iterator TrianglesIterator;
        typedef typename Triangles::const_iterator TrianglesConstIterator;
        typedef 
          CGAL::Box_intersection_d::Box_with_handle_d<double,3,TrianglesIterator> 
          Box;
        typedef 
          std::map<Index,std::vector<std::pair<Index,CGAL::Object> > > 
          OffendingMap;
        typedef std::map<std::pair<Index,Index>,std::vector<Index> >  EdgeMap;
        typedef std::pair<Index,Index> EMK;

        Triangles TA,TB;
        // Compute and process self intersections
        mesh_to_cgal_triangle_list(VA,FA,TA);
        mesh_to_cgal_triangle_list(VB,FB,TB);
        // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Box_intersection_d/Chapter_main.html#Section_63.5 
        // Create the corresponding vector of bounding boxes
        std::vector<Box> A_boxes,B_boxes;
        const auto box_up = [](Triangles & T, std::vector<Box> & boxes) -> void
        {
          boxes.reserve(T.size());
          for ( 
            TrianglesIterator tit = T.begin(); 
            tit != T.end(); 
            ++tit)
          {
            boxes.push_back(Box(tit->bbox(), tit));
          }
        };
        box_up(TA,A_boxes);
        box_up(TB,B_boxes);
        OffendingMap offendingA,offendingB;
        EdgeMap edge2facesA,edge2facesB;

        std::list<int> lIF;
        const auto cb = [&](const Box &a, const Box &b) -> void
        {
          using namespace std;
          // index in F and T
          int fa = a.handle()-TA.begin();
          int fb = b.handle()-TB.begin();
          const Triangle_3 & A = *a.handle();
          const Triangle_3 & B = *b.handle();
          if(CGAL::do_intersect(A,B))
          {
            // There was an intersection
            lIF.push_back(fa);
            lIF.push_back(fb);
            if(params.first_only)
            {
              throw IGL_FIRST_HIT_EXCEPTION;
            }
            if(!params.detect_only)
            {
              CGAL::Object result = CGAL::intersection(A,B);

              push_result(FA,fa,fb,result,offendingA,edge2facesA);
              push_result(FB,fb,fa,result,offendingB,edge2facesB);
            }
          }
        };
        try{
          CGAL::box_intersection_d(
            A_boxes.begin(), A_boxes.end(),
            B_boxes.begin(), B_boxes.end(),
            cb);
        }catch(int e)
        {
          // Rethrow if not FIRST_HIT_EXCEPTION
          if(e != IGL_FIRST_HIT_EXCEPTION)
          {
            throw e;
          }
          // Otherwise just fall through
        }

        // Convert lIF to Eigen matrix
        assert(lIF.size()%2 == 0);
        IF.resize(lIF.size()/2,2);
        {
          int i=0;
          for(
            list<int>::const_iterator ifit = lIF.begin();
            ifit!=lIF.end();
            )
          {
            IF(i,0) = (*ifit);
            ifit++; 
            IF(i,1) = (*ifit);
            ifit++;
            i++;
          }
        }
        if(!params.detect_only)
        {
          remesh_intersections(VA,FA,TA,offendingA,edge2facesA,VVA,FFA,JA,IMA);
          remesh_intersections(VB,FB,TB,offendingB,edge2facesB,VVB,FFB,JB,IMB);
        }

        return IF.rows() > 0;
      }
    }
  }
}

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedIF,
  typename DerivedVVA,
  typename DerivedFFA,
  typename DerivedJA,
  typename DerivedIMA,
  typename DerivedVVB,
  typename DerivedFFB,
  typename DerivedJB,
  typename DerivedIMB>
IGL_INLINE bool igl::copyleft::cgal::intersect_other(
  const Eigen::PlainObjectBase<DerivedVA> & VA,
  const Eigen::PlainObjectBase<DerivedFA> & FA,
  const Eigen::PlainObjectBase<DerivedVB> & VB,
  const Eigen::PlainObjectBase<DerivedFB> & FB,
  const RemeshSelfIntersectionsParam & params,
  Eigen::PlainObjectBase<DerivedIF> & IF,
  Eigen::PlainObjectBase<DerivedVVA> & VVA,
  Eigen::PlainObjectBase<DerivedFFA> & FFA,
  Eigen::PlainObjectBase<DerivedJA>  & JA,
  Eigen::PlainObjectBase<DerivedIMA> & IMA,
  Eigen::PlainObjectBase<DerivedVVB> & VVB,
  Eigen::PlainObjectBase<DerivedFFB> & FFB,
  Eigen::PlainObjectBase<DerivedJB>  & JB,
  Eigen::PlainObjectBase<DerivedIMB> & IMB)
{
  if(params.detect_only)
  {
    return intersect_other_helper<CGAL::Epick>
      (VA,FA,VB,FB,params,IF,VVA,FFA,JA,IMA,VVB,FFB,JB,IMB);
  }else
  {
    return intersect_other_helper<CGAL::Epeck>
      (VA,FA,VB,FB,params,IF,VVA,FFA,JA,IMA,VVB,FFB,JB,IMB);
  }
}

IGL_INLINE bool igl::copyleft::cgal::intersect_other(
  const Eigen::MatrixXd & VA,
  const Eigen::MatrixXi & FA,
  const Eigen::MatrixXd & VB,
  const Eigen::MatrixXi & FB,
  const bool first_only,
  Eigen::MatrixXi & IF)
{
  Eigen::MatrixXd VVA,VVB;
  Eigen::MatrixXi FFA,FFB;
  Eigen::VectorXi JA,JB,IMA,IMB;
  return intersect_other(
    VA,FA,VB,FB,{true,first_only},IF,VVA,FFA,JA,IMA,VVB,FFB,JB,IMB);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif
