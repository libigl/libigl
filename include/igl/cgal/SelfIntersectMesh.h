#ifndef IGL_SELFINTERSECTMESH_H
#define IGL_SELFINTERSECTMESH_H

#include "CGAL_includes.hpp"
#include "selfintersect.h"

#include <Eigen/Dense>
#include <list>
#include <map>
#include <vector>

#ifndef IGL_FIRST_HIT_EXCEPTION
#define IGL_FIRST_HIT_EXCEPTION 10
#endif

// The easiest way to keep track of everything is to use a class

namespace igl
{
  // Kernel is a CGAL kernel like:
  //     CGAL::Exact_predicates_inexact_constructions_kernel
  // or 
  //     CGAL::Exact_predicates_exact_constructions_kernel

  template <typename Kernel>
  class SelfIntersectMesh
  {
    public:
      // 3D Primitives
      typedef CGAL::Point_3<Kernel>    Point_3;
      typedef CGAL::Segment_3<Kernel>  Segment_3; 
      typedef CGAL::Triangle_3<Kernel> Triangle_3; 
      typedef CGAL::Plane_3<Kernel>    Plane_3;
      typedef CGAL::Tetrahedron_3<Kernel> Tetrahedron_3; 
      typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3; 
      typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3; 
      // 2D Primitives
      typedef CGAL::Point_2<Kernel>    Point_2;
      typedef CGAL::Segment_2<Kernel>  Segment_2; 
      typedef CGAL::Triangle_2<Kernel> Triangle_2; 
      // 2D Constrained Delaunay Triangulation types
      typedef CGAL::Triangulation_vertex_base_2<Kernel>  TVB_2;
      typedef CGAL::Constrained_triangulation_face_base_2<Kernel> CTFB_2;
      typedef CGAL::Triangulation_data_structure_2<TVB_2,CTFB_2> TDS_2;
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

      // Input mesh
      const Eigen::MatrixXd & V;
      const Eigen::MatrixXi & F;
      // Number of self-intersecting triangle pairs
      int count;
      std::vector<std::list<CGAL::Object> > F_objects;
      Triangles T;
      std::list<int> lIF;
      std::vector<bool> offensive;
      std::vector<int> offending_index;
      std::vector<int> offending;
      // Make a short name for the edge map's key
      typedef std::pair<int,int> EMK;
      // Make a short name for the type stored at each edge, the edge map's
      // value
      typedef std::list<int> EMV;
      // Make a short name for the edge map
      typedef std::map<EMK,EMV> EdgeMap;
      EdgeMap edge2faces;
    public:
      SelfintersectParam params;
    public:
      // Constructs (VV,FF) a new mesh with self-intersections of (V,F)
      // subdivided
      inline SelfIntersectMesh(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const SelfintersectParam & params,
        Eigen::MatrixXd & VV,
        Eigen::MatrixXi & FF,
        Eigen::MatrixXi & IF);
    private:
      // Helper function to mark a face as offensive
      //
      // Inputs:
      //   f  index of face in F
      inline void mark_offensive(const int f);
      // Helper function to count intersections between faces
      //
      // Input:
      //   fa  index of face A in F
      //   fb  index of face B in F
      inline void count_intersection(const int fa,const int fb);
      // Helper function for box_intersect. Intersect two triangles A and B,
      // append the intersection object (point,segment,triangle) to a running
      // list for A and B
      //
      // Inputs:
      //   A  triangle in 3D
      //   B  triangle in 3D
      //   fa  index of A in F (and F_objects)
      //   fb  index of A in F (and F_objects)
      // Returns true only if A intersects B
      //
      inline bool intersect(
          const Triangle_3 & A, 
          const Triangle_3 & B, 
          const int fa,
          const int fb);
      // Helper function for box_intersect. In the case where A and B have
      // already been identified to share a vertex, then we only want to add
      // possible segment intersections. Assumes truly duplicate triangles are
      // not given as input
      //
      // Inputs:
      //   A  triangle in 3D
      //   B  triangle in 3D
      //   fa  index of A in F (and F_objects)
      //   fb  index of B in F (and F_objects)
      //   va  index of shared vertex in A (and F_objects)
      //   vb  index of shared vertex in B (and F_objects)
      //// Returns object of intersection (should be Segment or point)
      //   Returns true if intersection (besides shared point)
      //
      inline bool single_shared_vertex(
          const Triangle_3 & A,
          const Triangle_3 & B,
          const int fa,
          const int fb,
          const int va,
          const int vb);
      // Helper handling one direction
      inline bool single_shared_vertex(
          const Triangle_3 & A,
          const Triangle_3 & B,
          const int fa,
          const int fb,
          const int va);
      // Helper function for box_intersect. In the case where A and B have
      // already been identified to share two vertices, then we only want to add
      // a possible coplanar (Triangle) intersection. Assumes truly degenerate
      // facets are not givine as input.
      inline bool double_shared_vertex(
          const Triangle_3 & A,
          const Triangle_3 & B,
          const int fa,
          const int fb);

    public:
      // Callback function called during box self intersections test. Means
      // boxes a and b intersect. This method then checks if the triangles in
      // each box intersect and if so, then processes the intersections
      //
      // Inputs:
      //   a  box containing a triangle
      //   b  box containing a triangle
      inline void box_intersect(const Box& a, const Box& b);
    private:
      // Compute 2D delaunay triangulation of a given 3d triangle and a list of
      // intersection objects (points,segments,triangles). CGAL uses an affine
      // projection rather than an isometric projection, so we're not
      // guaranteed that the 2D delaunay triangulation here will be a delaunay
      // triangulation in 3D.
      //
      // Inputs:
      //   A  triangle in 3D
      //   A_objects_3  updated list of intersection objects for A
      // Outputs:
      //   cdt  Contrained delaunay triangulation in projected 2D plane
      inline void projected_delaunay(
          const Triangle_3 & A,
          const std::list<CGAL::Object> & A_objects_3,
          CDT_plus_2 & cdt);
      // Getters:
    public:
      //const std::list<int>& get_lIF() const{ return lIF;}
      static inline void box_intersect(
        SelfIntersectMesh * SIM, 
        const SelfIntersectMesh::Box &a, 
        const SelfIntersectMesh::Box &b);
  };
}

// Implementation

#include "mesh_to_cgal_triangle_list.h"

#include <igl/REDRUM.h>
#include <igl/C_STR.h>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <algorithm>
#include <exception>
#include <cassert>
#include <iostream>

// References:
// http://minregret.googlecode.com/svn/trunk/skyline/src/extern/CGAL-3.3.1/examples/Polyhedron/polyhedron_self_intersection.cpp
// http://www.cgal.org/Manual/3.9/examples/Boolean_set_operations_2/do_intersect.cpp

// Q: Should we be using CGAL::Polyhedron_3?
// A: No! Input is just a list of unoriented triangles. Polyhedron_3 requires
// a 2-manifold.
// A: But! It seems we could use CGAL::Triangulation_3. Though it won't be easy
// to take advantage of functions like insert_in_facet because we want to
// constrain segments. Hmmm. Actualy Triangulation_3 doesn't look right...

//static void box_intersect(SelfIntersectMesh * SIM,const Box & A, const Box & B)
//{
//  return SIM->box_intersect(A,B);
//}


// CGAL's box_self_intersection_d uses C-style function callbacks without
// userdata. This is a leapfrog method for calling a member function. It should
// be bound as if the prototype was:
//   static void box_intersect(const Box &a, const Box &b)
// using boost:
//  boost::function<void(const Box &a,const Box &b)> cb
//    = boost::bind(&::box_intersect, this, _1,_2);
//   
template <typename Kernel>
inline void igl::SelfIntersectMesh<Kernel>::box_intersect(
  igl::SelfIntersectMesh<Kernel> * SIM, 
  const typename igl::SelfIntersectMesh<Kernel>::Box &a, 
  const typename igl::SelfIntersectMesh<Kernel>::Box &b)
{
  SIM->box_intersect(a,b);
}

template <typename Kernel>
inline igl::SelfIntersectMesh<Kernel>::SelfIntersectMesh(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const SelfintersectParam & params,
  Eigen::MatrixXd & VV,
  Eigen::MatrixXi & FF,
  Eigen::MatrixXi & IF):
  V(V),
  F(F),
  count(0),
  F_objects(F.rows()),
  T(),
  lIF(),
  offensive(F.rows(),false),
  offending_index(F.rows(),-1),
  offending(),
  edge2faces(),
  params(params)
{
  using namespace std;
  using namespace Eigen;
  // Compute and process self intersections
  mesh_to_cgal_triangle_list(V,F,T);
  // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Box_intersection_d/Chapter_main.html#Section_63.5 
  // Create the corresponding vector of bounding boxes
  std::vector<Box> boxes;
  boxes.reserve(T.size());
  for ( 
    TrianglesIterator tit = T.begin(); 
    tit != T.end(); 
    ++tit)
  {
    boxes.push_back(Box(tit->bbox(), tit));
  }
  // Leapfrog callback
  boost::function<void(const Box &a,const Box &b)> cb
    = boost::bind(&box_intersect, this, _1,_2);
  // Run the self intersection algorithm with all defaults
  try{
    CGAL::box_self_intersection_d(boxes.begin(), boxes.end(),cb);
  }catch(int e)
  {
    // Rethrow if not IGL_FIRST_HIT_EXCEPTION
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
      typename list<int>::const_iterator ifit = lIF.begin();
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

  if(params.detect_only)
  {
    return;
  }

  int NF_count = 0;
  // list of new faces, we'll fix F later
  vector<MatrixXi> NF(offending.size());
  // list of new vertices
  list<Point_3> NV;
  int NV_count = 0;
  vector<CDT_plus_2> cdt(offending.size());
  vector<Plane_3> P(offending.size());
  // Use map for *all* faces
  map<typename CDT_plus_2::Vertex_handle,int> v2i;
  // Loop over offending triangles
  for(int o = 0;o<(int)offending.size();o++)
  {
    // index in F
    const int f = offending[o];
    projected_delaunay(T[f],F_objects[f],cdt[o]);
    // Q: Is this also delaunay in 3D?
    // A: No, because the projection is affine and delaunay is not affine
    // invariant.
    // Q: Then, can't we first get the 2D delaunay triangulation, then lift it
    // to 3D and flip any offending edges?
    // Plane of projection (also used by projected_delaunay)
    P[o] = Plane_3(T[f].vertex(0),T[f].vertex(1),T[f].vertex(2));
    // Build index map
    {
      int i=0;
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
          // original corners of the CDT. I.e. I don't trust CGAL to maintain the
          // order
          assert(T[f].vertex(i) == P[o].to_3d(vit->point()));
#endif
          // For first three, use original index in F
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
            int i = F(f,(e+1)%3);
            int j = F(f,(e+2)%3);
            // Be sure that i<j
            if(i>j)
            {
              swap(i,j);
            }
            assert(edge2faces.count(EMK(i,j))==1);
            // loop over neighbors
            for(
              list<int>::const_iterator nit = edge2faces[EMK(i,j)].begin();
              nit != edge2faces[EMK(i,j)].end() && !found;
              nit++)
            {
              // don't consider self
              if(*nit == f)
              {
                continue;
              }
              // index of neighbor in offending (to find its cdt)
              int no = offending_index[*nit];
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
                  v2i[vit] = v2i[uit];
                  found = true;
                }
              }
            }
          }
          if(!found)
          {
            v2i[vit] = V.rows()+NV_count;
            NV.push_back(vit_point_3);
            NV_count++;
          }
        }
        i++;
      }
    }
    {
      int i = 0;
      // Resize to fit new number of triangles
      NF[o].resize(cdt[o].number_of_faces(),3);
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
  assert(NV_count == (int)NV.size());
  // Build output
#ifndef NDEBUG
  {
    int off_count = 0;
    for(int f = 0;f<F.rows();f++)
    {
      off_count+= (offensive[f]?1:0);
    }
    assert(off_count==(int)offending.size());
  }
#endif
  // Append faces
  FF.resize(F.rows()-offending.size()+NF_count,3);
  // First append non-offending original faces
  // There's an Eigen way to do this in one line but I forget
  int off = 0;
  for(int f = 0;f<F.rows();f++)
  {
    if(!offensive[f])
    {
      FF.row(off++) = F.row(f);
    }
  }
  assert(off == (int)(F.rows()-offending.size()));
  // Now append replacement faces for offending faces
  for(int o = 0;o<(int)offending.size();o++)
  {
    FF.block(off,0,NF[o].rows(),3) = NF[o];
    off += NF[o].rows();
  }
  // Append vertices
  VV.resize(V.rows()+NV_count,3);
  VV.block(0,0,V.rows(),3) = V;
  {
    int i = 0;
    for(
      typename list<Point_3>::const_iterator nvit = NV.begin();
      nvit != NV.end();
      nvit++)
    {
      for(int d = 0;d<3;d++)
      {
        const Point_3 & p = *nvit;
        VV(V.rows()+i,d) = CGAL::to_double(p[d]);
        // This distinction does not seem necessary:
//#ifdef INEXACT_CONSTRUCTION
//        VV(V.rows()+i,d) = CGAL::to_double(p[d]);
//#else
//        VV(V.rows()+i,d) = CGAL::to_double(p[d].exact());
//#endif
      }
      i++;
    }
  }

  // Q: Does this give the same result as TETGEN?
  // A: For the cow and beast, yes.

  // Q: Is tetgen faster than this CGAL implementation?
  // A: Well, yes. But Tetgen is only solving the detection (predicates)
  // problem. This is also appending the intersection objects (construction).
  // But maybe tetgen is still faster for the detection part. For example, this
  // CGAL implementation on the beast takes 98 seconds but tetgen detection
  // takes 14 seconds

}


template <typename Kernel>
inline void igl::SelfIntersectMesh<Kernel>::mark_offensive(const int f)
{
  using namespace std;
  lIF.push_back(f);
  if(!offensive[f])
  {
    offensive[f]=true;
    offending_index[f]=offending.size();
    offending.push_back(f);
    // Add to edge map
    for(int e = 0; e<3;e++)
    {
      // Index of F's eth edge in V
      int i = F(f,(e+1)%3);
      int j = F(f,(e+2)%3);
      // Be sure that i<j
      if(i>j)
      {
        swap(i,j);
      }
      // Create new list if there is no entry
      if(edge2faces.count(EMK(i,j))==0)
      {
        edge2faces[EMK(i,j)] = list<int>();
      }
      // append to list
      edge2faces[EMK(i,j)].push_back(f);
    }
  }
}

template <typename Kernel>
inline void igl::SelfIntersectMesh<Kernel>::count_intersection(
  const int fa,
  const int fb)
{
  mark_offensive(fa);
  mark_offensive(fb);
  this->count++;
  // We found the first intersection
  if(params.first_only && this->count >= 1)
  {
    throw IGL_FIRST_HIT_EXCEPTION;
  }
}

template <typename Kernel>
inline bool igl::SelfIntersectMesh<Kernel>::intersect(
  const Triangle_3 & A, 
  const Triangle_3 & B, 
  const int fa,
  const int fb)
{
  // Determine whether there is an intersection
  if(!CGAL::do_intersect(A,B))
  {
    return false;
  }
  if(!params.detect_only)
  {
    // Construct intersection
    CGAL::Object result = CGAL::intersection(A,B);
    F_objects[fa].push_back(result);
    F_objects[fb].push_back(result);
  }
  count_intersection(fa,fb);
  return true;
}

template <typename Kernel>
inline bool igl::SelfIntersectMesh<Kernel>::single_shared_vertex(
  const Triangle_3 & A,
  const Triangle_3 & B,
  const int fa,
  const int fb,
  const int va,
  const int vb)
{
  ////using namespace std;
  //CGAL::Object result = CGAL::intersection(A,B);
  //if(CGAL::object_cast<Segment_3 >(&result))
  //{
  //  // Append to each triangle's running list
  //  F_objects[fa].push_back(result);
  //  F_objects[fb].push_back(result);
  //  count_intersection(fa,fb);
  //}else
  //{
  //  // Then intersection must be at point
  //  // And point must be at shared vertex
  //  assert(CGAL::object_cast<Point_3>(&result));
  //}
  if(single_shared_vertex(A,B,fa,fb,va))
  {
    return true;
  }
  return single_shared_vertex(B,A,fb,fa,vb);
}

template <typename Kernel>
inline bool igl::SelfIntersectMesh<Kernel>::single_shared_vertex(
  const Triangle_3 & A,
  const Triangle_3 & B,
  const int fa,
  const int fb,
  const int va)
{
  // This was not a good idea. It will not handle coplanar triangles well.
  using namespace std;
  Segment_3 sa(
    A.vertex((va+1)%3),
    A.vertex((va+2)%3));

  if(CGAL::do_intersect(sa,B))
  {
    CGAL::Object result = CGAL::intersection(sa,B);
    if(const Point_3 * p = CGAL::object_cast<Point_3 >(&result))
    {
      if(!params.detect_only)
      {
        // Single intersection --> segment from shared point to intersection
        CGAL::Object seg = CGAL::make_object(Segment_3(
          A.vertex(va),
          *p));
        F_objects[fa].push_back(seg);
        F_objects[fb].push_back(seg);
      }
      count_intersection(fa,fb);
      return true;
    }else if(CGAL::object_cast<Segment_3 >(&result))
    {
      //cerr<<REDRUM("Coplanar at: "<<fa<<" & "<<fb<<" (single shared).")<<endl;
      // Must be coplanar
      if(params.detect_only)
      {
        count_intersection(fa,fb);
      }else
      {
        // WRONG:
        //// Segment intersection --> triangle from shared point to intersection
        //CGAL::Object tri = CGAL::make_object(Triangle_3(
        //  A.vertex(va),
        //  s->vertex(0),
        //  s->vertex(1)));
        //F_objects[fa].push_back(tri);
        //F_objects[fb].push_back(tri);
        //count_intersection(fa,fb);
        // Need to do full test. Intersection could be a general poly.
        bool test = intersect(A,B,fa,fb);
        ((void)test);
        assert(test);
      }
      return true;
    }else
    {
      cerr<<REDRUM("Segment ∩ triangle neither point nor segment?")<<endl;
      assert(false);
    }
  }

  return false;
}


template <typename Kernel>
inline bool igl::SelfIntersectMesh<Kernel>::double_shared_vertex(
  const Triangle_3 & A,
  const Triangle_3 & B,
  const int fa,
  const int fb)
{
  using namespace std;
  // Cheaper way to do this than calling do_intersect?
  if(
    // Can only have an intersection if co-planar
    (A.supporting_plane() == B.supporting_plane() ||
    A.supporting_plane() == B.supporting_plane().opposite()) &&
    CGAL::do_intersect(A,B))
  {
    // Construct intersection
    try
    {
      CGAL::Object result = CGAL::intersection(A,B);
      if(result)
      {
        if(CGAL::object_cast<Segment_3 >(&result))
        {
          // not coplanar
          return false;
        } else if(CGAL::object_cast<Point_3 >(&result))
        {
          // this "shouldn't" happen but does for inexact
          return false;
        } else
        {
          if(!params.detect_only)
          {
            // Triangle object
            F_objects[fa].push_back(result);
            F_objects[fb].push_back(result);
          }
          count_intersection(fa,fb);
          //cerr<<REDRUM("Coplanar at: "<<fa<<" & "<<fb<<" (double shared).")<<endl;
          return true;
        }
      }else
      {
        // CGAL::intersection is disagreeing with do_intersect
        return false;
      }
    }catch(...)
    {
      // This catches some cgal assertion:
      //     CGAL error: assertion violation!
      //     Expression : is_finite(d)
      //     File       : /opt/local/include/CGAL/GMP/Gmpq_type.h
      //     Line       : 132
      //     Explanation: 
      // But only if NDEBUG is not defined, otherwise there's an uncaught
      // "Floating point exception: 8" SIGFPE
      return false;
    }
  }
  // Shouldn't get here either
  assert(false);
  return false;
}

template <typename Kernel>
inline void igl::SelfIntersectMesh<Kernel>::box_intersect(
  const Box& a, 
  const Box& b)
{
  using namespace std;
  // index in F and T
  int fa = a.handle()-T.begin();
  int fb = b.handle()-T.begin();
  const Triangle_3 & A = *a.handle();
  const Triangle_3 & B = *b.handle();
  // I'm not going to deal with degenerate triangles, though at some point we
  // should
  assert(!a.handle()->is_degenerate());
  assert(!b.handle()->is_degenerate());
  // Number of combinatorially shared vertices
  int comb_shared_vertices = 0;
  // Number of geometrically shared vertices (*not* including combinatorially
  // shared)
  int geo_shared_vertices = 0;
  // Keep track of shared vertex indices (we only handles single shared
  // vertices as a special case, so just need last/first/only ones)
  int va=-1,vb=-1;
  int ea,eb;
  for(ea=0;ea<3;ea++)
  {
    for(eb=0;eb<3;eb++)
    {
      if(F(fa,ea) == F(fb,eb))
      {
        comb_shared_vertices++;
        va = ea;
        vb = eb;
      }else if(A.vertex(ea) == B.vertex(eb))
      {
        geo_shared_vertices++;
        va = ea;
        vb = eb;
      }
    }
  }
  const int total_shared_vertices = comb_shared_vertices + geo_shared_vertices;
  if(comb_shared_vertices== 3)
  {
    // Combinatorially duplicate face, these should be removed by preprocessing
    cerr<<REDRUM("Facets "<<fa<<" and "<<fb<<" are combinatorial duplicates")<<endl;
    goto done;
  }
  if(total_shared_vertices== 3)
  {
    // Geometrically duplicate face, these should be removed by preprocessing
    cerr<<REDRUM("Facets "<<fa<<" and "<<fb<<" are geometrical duplicates")<<endl;
    goto done;
  }
  //// SPECIAL CASES ARE BROKEN FOR COPLANAR TRIANGLES
  //if(total_shared_vertices > 0)
  //{
  //  bool coplanar = 
  //    CGAL::coplanar(A.vertex(0),A.vertex(1),A.vertex(2),B.vertex(0)) &&
  //    CGAL::coplanar(A.vertex(0),A.vertex(1),A.vertex(2),B.vertex(1)) &&
  //    CGAL::coplanar(A.vertex(0),A.vertex(1),A.vertex(2),B.vertex(2));
  //  if(coplanar)
  //  {
  //    cerr<<MAGENTAGIN("Facets "<<fa<<" and "<<fb<<
  //      " are coplanar and share vertices")<<endl;
  //    goto full;
  //  }
  //}

  if(total_shared_vertices == 2)
  {
    // Q: What about coplanar?
    //
    // o    o
    // |\  /|
    // | \/ |
    // | /\ |
    // |/  \|
    // o----o
    double_shared_vertex(A,B,fa,fb);

    goto done;
  }
  assert(total_shared_vertices<=1);
  if(total_shared_vertices==1)
  {
    assert(va>=0 && va<3);
    assert(vb>=0 && vb<3);
//#ifndef NDEBUG
//    CGAL::Object result =
//#endif
    single_shared_vertex(A,B,fa,fb,va,vb);
//#ifndef NDEBUG
//    if(!CGAL::object_cast<Segment_3 >(&result))
//    {
//      const Point_3 * p = CGAL::object_cast<Point_3 >(&result);
//      assert(p);
//      for(int ea=0;ea<3;ea++)
//      {
//        for(int eb=0;eb<3;eb++)
//        {
//          if(F(fa,ea) == F(fb,eb))
//          {
//            assert(*p==A.vertex(ea));
//            assert(*p==B.vertex(eb));
//          }
//        }
//      }
//    }
//#endif
  }else
  {
//full:
    // No geometrically shared vertices, do general intersect
    intersect(*a.handle(),*b.handle(),fa,fb);
  }
done:
  return;
}

// Compute 2D delaunay triangulation of a given 3d triangle and a list of
// intersection objects (points,segments,triangles). CGAL uses an affine
// projection rather than an isometric projection, so we're not guaranteed
// that the 2D delaunay triangulation here will be a delaunay triangulation
// in 3D.
//
// Inputs:
//   A  triangle in 3D
//   A_objects_3  updated list of intersection objects for A
// Outputs:
//   cdt  Contrained delaunay triangulation in projected 2D plane
template <typename Kernel>
inline void igl::SelfIntersectMesh<Kernel>::projected_delaunay(
  const Triangle_3 & A,
  const std::list<CGAL::Object> & A_objects_3,
  CDT_plus_2 & cdt)
{
  using namespace std;
  // http://www.cgal.org/Manual/3.2/doc_html/cgal_manual/Triangulation_2/Chapter_main.html#Section_2D_Triangulations_Constrained_Plus
  // Plane of triangle A
  Plane_3 P(A.vertex(0),A.vertex(1),A.vertex(2));
  // Insert triangle into vertices
  typename CDT_plus_2::Vertex_handle corners[3];
  for(int i = 0;i<3;i++)
  {
    corners[i] = cdt.insert(P.to_2d(A.vertex(i)));
  }
  // Insert triangle edges as constraints
  for(int i = 0;i<3;i++)
  {
    cdt.insert_constraint( corners[(i+1)%3], corners[(i+2)%3]);
  }
  // Insert constraints for intersection objects
  for(
    typename list<CGAL::Object>::const_iterator lit = A_objects_3.begin();
    lit != A_objects_3.end();
    lit++)
  {
    CGAL::Object obj = *lit;
    if(const Point_3 *ipoint = CGAL::object_cast<Point_3 >(&obj))
    {
      // Add point
      cdt.insert(P.to_2d(*ipoint));
    } else if(const Segment_3 *iseg = CGAL::object_cast<Segment_3 >(&obj))
    {
      // Add segment constraint
      cdt.insert_constraint(P.to_2d(iseg->vertex(0)),P.to_2d(iseg->vertex(1)));
    } else if(const Triangle_3 *itri = CGAL::object_cast<Triangle_3 >(&obj))
    {
      // Add 3 segment constraints
      cdt.insert_constraint(P.to_2d(itri->vertex(0)),P.to_2d(itri->vertex(1)));
      cdt.insert_constraint(P.to_2d(itri->vertex(1)),P.to_2d(itri->vertex(2)));
      cdt.insert_constraint(P.to_2d(itri->vertex(2)),P.to_2d(itri->vertex(0)));
    } else if(const std::vector<Point_3 > *polyp = 
        CGAL::object_cast< std::vector<Point_3 > >(&obj))
    {
      //cerr<<REDRUM("Poly...")<<endl;
      const std::vector<Point_3 > & poly = *polyp;
      const int m = poly.size();
      assert(m>=2);
      for(int p = 0;p<m;p++)
      {
        const int np = (p+1)%m;
        cdt.insert_constraint(P.to_2d(poly[p]),P.to_2d(poly[np]));
      }
    }else
    {
      cerr<<REDRUM("What is this object?!")<<endl;
      assert(false);
    }
  }
}

#endif

