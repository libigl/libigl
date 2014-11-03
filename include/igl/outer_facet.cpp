#include "outer_facet.h"
#include "sort.h"
#include "vertex_triangle_adjacency.h"

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedN,
  typename DerivedI,
  typename f_type>
IGL_INLINE void igl::outer_facet(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const Eigen::PlainObjectBase<DerivedN> & N,
  const Eigen::PlainObjectBase<DerivedI> & I,
  f_type & max_f,
  bool & flip)
{
  using namespace std;
  typedef typename DerivedV::Scalar Scalar;
  typedef typename DerivedV::Index Index;
  typedef 
    typename Eigen::Matrix<Index, DerivedV::RowsAtCompileTime,1> VectorXI;
  typedef 
    typename Eigen::Matrix<Scalar, DerivedV::RowsAtCompileTime,1> VectorXS;
  // "Direct repair of self-intersecting meshes" [Attene 14]
  const Index mi = I.size();
  Scalar max_x = -1e26;
  Index max_v = -1;
  Scalar max_nx = -1e26;
  for(Index i = 0;i<mi;i++)
  {
    const Index f = I(i);
    const Scalar nx = N(f,0);
    if(fabs(nx)>0)
    {
      for(Index c = 0;c<3;c++)
      {
        const Index v = F(f,c);
        if(v == max_v)
        {
          if(fabs(nx) > max_nx)
          {
            // Just update max face and normal
            max_f = f;
            max_nx = fabs(nx);
            flip = nx<0;
          }
        }else
        {
          const Scalar x = V(v);
          if(x>max_x)
          {
            // update max vertex, face and normal
            max_v = v;
            max_x = x;
            max_f = f;
            max_nx = fabs(nx);
            flip = nx<0;
          }
        }
      }
    }
  }
  assert(max_v >=0 && "Very degenerate case, no suitable face found.");
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::outer_facet<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> > const&, int&, bool&);
#endif
