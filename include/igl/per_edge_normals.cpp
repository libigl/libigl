#include "all_edges.h"
#include "doublearea.h"
#include "per_edge_normals.h"
#include "per_face_normals.h"
#include "unique_simplices.h"
#include <vector>

template <
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedN,
  typename DerivedE,
  typename DerivedEMAP>
IGL_INLINE void igl::per_edge_normals(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedN> & N,
  Eigen::PlainObjectBase<DerivedE> & E,
  Eigen::PlainObjectBase<DerivedEMAP> & EMAP)
{
  using namespace Eigen;
  using namespace std;
  assert(F.cols() == 3 && "Faces must be triangles");
  // number of faces
  const int m = F.rows();
  // All occurances of directed edges
  MatrixXi allE;
  all_edges(F,allE);
  // Find unique undirected edges and mapping
  VectorXi _;
  unique_simplices(allE,E,_,EMAP);
  // now sort(allE,2) == E(EMAP,:), that is, if EMAP(i) = j, then E.row(j) is
  // the undirected edge corresponding to the directed edge allE.row(i).
  MatrixXd FN;
  per_face_normals(V,F,FN);

  VectorXd dblA;
  doublearea(V,F,dblA);
  N.setConstant(E.rows(),3,0);
  for(int f = 0;f<m;f++)
  {
    for(int c = 0;c<3;c++)
    {
      N.row(EMAP(f+c*m)) += dblA(f) * FN.row(f);
    }
  }
  N.rowwise().normalize();
  
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instanciation
template void igl::per_edge_normals<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
