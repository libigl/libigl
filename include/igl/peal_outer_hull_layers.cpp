#include "peal_outer_hull_layers.h"
#include "outer_hull.h"
#include <vector>

template <
  typename DerivedV,
  typename DerivedF,
  typename Derivedodd,
  typename Derivedflip>
IGL_INLINE void igl::peal_outer_hull_layers(
  const Eigen::PlainObjectBase<DerivedV > & V,
  const Eigen::PlainObjectBase<DerivedF > & F,
  Eigen::PlainObjectBase<Derivedodd > & odd,
  Eigen::PlainObjectBase<Derivedflip > & flip)
{
  using namespace Eigen;
  using namespace std;
  typedef typename DerivedF::Index Index;
  typedef Matrix<typename DerivedF::Scalar,Dynamic,DerivedF::ColsAtCompileTime> MatrixXF;
  typedef Matrix<Index,Dynamic,1> MatrixXI;
  typedef Matrix<typename Derivedflip::Scalar,Dynamic,Derivedflip::ColsAtCompileTime> MatrixXflip;
  const Index m = F.rows();

  // keep track of iteration parity and whether flipped in hull
  MatrixXF Fr = F;
  odd.resize(m,1);
  flip.resize(m,1);
  // Keep track of index map
  MatrixXI IM = MatrixXI::LinSpaced(m,0,m-1);
  // This is O(n * layers)
  bool odd_iter = true;
  MatrixXI P(m,1);
  Index iter = 0;
  while(Fr.size() > 0)
  {
    assert(Fr.rows() == IM.rows());
    // Compute outer hull of current Fr
    MatrixXF Fo;
    MatrixXI Jo;
    MatrixXflip flipr;
    outer_hull(V,Fr,Fo,Jo,flipr);
    assert(Fo.rows() == Jo.rows());
    // all faces in Fo of Fr
    vector<bool> in_outer(Fr.rows(),false);
    for(Index g = 0;g<Jo.rows();g++)
    {
      odd(IM(Jo(g))) = odd_iter;
      P(IM(Jo(g))) = iter;
      in_outer[Jo(g)] = true;
      flip(IM(Jo(g))) = flipr(Jo(g));
    }
    // Fr = Fr - Fo
    // update IM
    MatrixXF prev_Fr = Fr;
    MatrixXI prev_IM = IM;
    Fr.resize(prev_Fr.rows() - Fo.rows(),F.cols());
    IM.resize(Fr.rows());
    {
      Index g = 0;
      for(Index f = 0;f<prev_Fr.rows();f++)
      {
        if(!in_outer[f])
        {
          Fr.row(g) = prev_Fr.row(f);
          IM(g) = prev_IM(f);
          g++;
        }
      }
    }
    odd_iter = !odd_iter;
    iter++;
  }
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::peal_outer_hull_layers<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<bool, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&);
#endif
