#include "triangle_triangle_intersect.h"
#include "triangle_triangle_intersect_shared_edge.h"
#include "triangle_triangle_intersect_shared_vertex.h"
#include "tri_tri_intersect.h"
#include <Eigen/Geometry>
#include <stdio.h>

//#define IGL_TRIANGLE_TRIANGLE_INTERSECT_DEBUG
#ifdef IGL_TRIANGLE_TRIANGLE_INTERSECT_DEBUG
// CGAL::Epeck
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#warning "🐌🐌🐌🐌🐌🐌🐌🐌 Slow debug mode for igl::triangle_triangle_intersect"
#endif

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedE,
  typename DerivedEMAP,
  typename DerivedEF,
  typename Derivedp>
IGL_INLINE bool igl::triangle_triangle_intersect(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<DerivedE> & E,
  const Eigen::MatrixBase<DerivedEMAP> & EMAP,
  const Eigen::MatrixBase<DerivedEF> & EF,
  const int f,
  const int c,
  const Eigen::MatrixBase<Derivedp> & p,
  const int g)
{
  static_assert(
    std::is_same<typename DerivedV::Scalar,typename Derivedp::Scalar>::value, 
    "V and p should have same Scalar type");
  assert(V.cols() == 3);
  assert(p.cols() == 3);
#ifdef IGL_TRIANGLE_TRIANGLE_INTERSECT_DEBUG
  using Kernel = CGAL::Epeck;
  typedef CGAL::Point_3<Kernel>    Point_3;
  typedef CGAL::Segment_3<Kernel>  Segment_3;
  typedef CGAL::Triangle_3<Kernel> Triangle_3;
    bool cgal_found_intersection = false;
    Point_3 Vg[3];
    Point_3 Vf[3];
    for(int i = 0;i<3;i++)
    {
      Vg[i] = Point_3(V(F(g,i),0),V(F(g,i),1),V(F(g,i),2));
      if(i == c)
      {
        Vf[i] = Point_3(p(0),p(1),p(2));
      }else
      {
        Vf[i] = Point_3(V(F(f,i),0),V(F(f,i),1),V(F(f,i),2));
      }
    }
    Triangle_3 Tg(Vg[0],Vg[1],Vg[2]);
    Triangle_3 Tf(Vf[0],Vf[1],Vf[2]);
#endif

  // I'm leaving this debug printing stuff in for a bit until I trust this
  // better. 
  constexpr bool stinker = false;
  //const bool stinker = (f==1492 && g==1554);
  if(stinker) { printf("👀\n"); }
  bool found_intersection = false;
  // So edge opposite F(f,c) is the outer edge.
  const bool share_edge_opposite_c = [&]()
  {
    const int o = EMAP(f + c*F.rows());
    return (EF(o,0) == f && EF(o,1) == g) || (EF(o,1) == f && EF(o,0) == g);
  }();
  const int o = EMAP(f + c*F.rows());
  // Do they share the edge opposite c?
  if(share_edge_opposite_c)
  {
    if(stinker) { printf("⚠️ shares an edge\n"); }
    found_intersection = triangle_triangle_intersect_shared_edge(V,F,f,c,p,g,1e-8);
  }else
  {
    if(stinker) { printf("does not share an edge\n"); }
    // Do they share a vertex?
    int sf,sg;
    bool found_shared_vertex = false;
    for(sf = 0;sf<3;sf++)
    {
      if(sf == c){ continue;}
      for(sg = 0;sg<3;sg++)
      {
        if(F(f,sf) == F(g,sg))
        {
          found_shared_vertex = true;
          break;
        }
      }
      if(found_shared_vertex) { break;} 
    }
    if(found_shared_vertex)
    {
      if(stinker) { printf("⚠️ shares a vertex\n"); }
      found_intersection = 
        triangle_triangle_intersect_shared_vertex(V,F,f,sf,c,p,g,sg,1e-14);
    }else
    {
      bool coplanar;
      Eigen::RowVector3d i1,i2;
      found_intersection = igl::tri_tri_intersection_test_3d(
                V.row(F(g,0)).template cast<double>(), 
                V.row(F(g,1)).template cast<double>(), 
                V.row(F(g,2)).template cast<double>(),
                            p.template cast<double>(),
          V.row(F(f,(c+1)%3)).template cast<double>(),
          V.row(F(f,(c+2)%3)).template cast<double>(),
          coplanar,
          i1,i2);
      if(stinker) { printf("tri_tri_intersection_test_3d says %s\n",found_intersection?"☠️":"✅"); }
#ifdef IGL_TRIANGLE_TRIANGLE_INTERSECT_DEBUG
      if(CGAL::do_intersect(Tg,Tf))
      {
        cgal_found_intersection = true;
        printf("  ✅ sure it's anything\n");
      }
      assert(found_intersection == cgal_found_intersection);
#endif
    }
  }
  if(stinker) { printf("%s\n",found_intersection?"☠️":"✅"); }
  return found_intersection;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template bool igl::triangle_triangle_intersect<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 1, -1, false> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, int, int, Eigen::MatrixBase<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 1, -1, false> > const&, int);
template bool igl::triangle_triangle_intersect<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, int, int, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, int);
#endif
