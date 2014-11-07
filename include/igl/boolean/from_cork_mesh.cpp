#ifndef IGL_NO_CORK
#include "from_cork_mesh.h"

template <
  typename DerivedV,
  typename DerivedF>
IGL_INLINE void igl::from_cork_mesh(
  const CorkTriMesh & mesh,
  Eigen::PlainObjectBase<DerivedV > & V,
  Eigen::PlainObjectBase<DerivedF > & F)
{
  using namespace std;
  F.resize(mesh.n_triangles,3);
  V.resize(mesh.n_vertices,3);
  for(size_t v = 0;v<mesh.n_vertices;v++)
  {
    for(size_t c = 0;c<3;c++)
    {
      V(v,c) = mesh.vertices[v*3+c];
    }
  }
  for(size_t f = 0;f<mesh.n_triangles;f++)
  {
    for(size_t c = 0;c<3;c++)
    {
      F(f,c) = mesh.triangles[f*3+c];
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::from_cork_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(CorkTriMesh const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::from_cork_mesh<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(CorkTriMesh const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
#endif

#endif
