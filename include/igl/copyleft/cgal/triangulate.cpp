#include "triangulate.h"
#include "assign_scalar.h"
#include "../../list_to_matrix.h"
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <queue>

template <
  typename Kernel,
  typename DerivedV,
  typename DerivedE,
  typename DerivedH,
  typename DerivedV2,
  typename DerivedF2>
IGL_INLINE void igl::copyleft::cgal::triangulate(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedE> & E,
  const Eigen::MatrixBase<DerivedH> & H,
  const bool retain_convex_hull,
  Eigen::PlainObjectBase<DerivedV2> & TV,
  Eigen::PlainObjectBase<DerivedF2> & TF)
{
  struct FaceInfo2
  {
    FaceInfo2(){}
    bool visited = false;
    bool in_domain = true;
  };

  typedef Eigen::Index Index;
  //typedef CGAL::Epeck Kernel;
  typedef CGAL::Point_2<Kernel>    Point_2;
  typedef CGAL::Segment_2<Kernel>  Segment_2;
  typedef CGAL::Triangle_2<Kernel> Triangle_2;
  typedef CGAL::Triangulation_vertex_base_2<Kernel>  TVB_2;
  typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,Kernel>    TFB_2;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel,TFB_2> CTFB_2;
  typedef CGAL::Triangulation_data_structure_2<TVB_2,CTFB_2> TDS_2;
  typedef CGAL::Exact_intersections_tag Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,TDS_2,Itag> CDT_2;
  typedef CGAL::Constrained_triangulation_plus_2<CDT_2> CDT_plus_2;
  typedef typename CDT_plus_2::Face_handle Face_handle;

  CDT_plus_2 cdt;
  // Insert all points
  std::vector<Point_2> points;points.reserve(V.rows());
  for(Eigen::Index i = 0;i<V.rows();i++)
  {
    Point_2 p(V(i,0),V(i,1));
    cdt.insert(p);
    points.emplace_back(p);
  }
  // Insert all edges 
  for(Eigen::Index e = 0;e<E.rows();e++)
  {
    cdt.insert_constraint(points[E(e,0)],points[E(e,1)]);
  }

  // https://doc.cgal.org/latest/Triangulation_2/index.html#title30
  //
  // "remove" connected component of face start by marking faces as
  // in_domain = false. start should not yet be visited.
  const auto remove_cc = [&](Face_handle start)
  {
    std::queue<Face_handle> queue;
    queue.push(start);
    while(!queue.empty())
    {
      Face_handle fh = queue.front();
      queue.pop();
      if(!fh->info().visited)
      {
        fh->info().visited = true;
        fh->info().in_domain = false;
        for(int i = 0; i < 3; i++)
        {
          typename CDT_plus_2::Edge e(fh,i);
          Face_handle n = fh->neighbor(i);
          if(!n->info().visited && !cdt.is_constrained(e))
          {
            queue.push(n);
          }
        }
      }
    }
  };

  // If _not_ meshsing convex hull, remove anything connected to infinite face
  if(!retain_convex_hull)
  {
    remove_cc(cdt.infinite_face());
  }
  // remove holes
  for(Eigen::Index h = 0;h<H.rows();h++)
  {
    remove_cc(cdt.locate({H(h,0),H(h,1)}));
  }

  std::vector<CGAL::Point_2<Kernel> > vertices;
  std::vector<std::vector<Index> > faces;
  // Read off vertices of the cdt, remembering index
  std::map<typename CDT_plus_2::Vertex_handle,Index> v2i;
  size_t count=0;
  for (
    auto itr = cdt.finite_vertices_begin();
    itr != cdt.finite_vertices_end();
    itr++)
  {
    vertices.push_back(itr->point());
    v2i[itr] = count;
    count++;
  }
  // Read off faces and store index triples
  for (
    auto itr = cdt.finite_faces_begin();
    itr != cdt.finite_faces_end();
    itr++)
  {
    if(itr->info().in_domain)
    {
      faces.push_back(
        { v2i[itr->vertex(0)], v2i[itr->vertex(1)], v2i[itr->vertex(2)] });
    }
  }
  TV.resize(vertices.size(),2);
  for(int v = 0;v<vertices.size();v++)
  {
    for(int d = 0;d<2;d++)
    {
      assign_scalar(vertices[v][d], TV(v,d));
    }
  }
  list_to_matrix(faces,TF);

}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::copyleft::cgal::triangulate<CGAL::Epeck, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, const bool , Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::copyleft::cgal::triangulate<CGAL::Epick, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&,       const bool , Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::copyleft::cgal::triangulate<CGAL::Epeck, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&,       const bool , Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
