#include "minkowski_sum.h"
#include "mesh_boolean.h"

#include "../../slice_mask.h"
#include "../../unique.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

template <
  typename DerivedVA,
  typename DerivedFA,
  typename Deriveds,
  typename Derivedd,
  typename DerivedW,
  typename DerivedG,
  typename DerivedJ>
IGL_INLINE void igl::copyleft::boolean::minkowski_sum(
  const Eigen::PlainObjectBase<DerivedVA> & VA,
  const Eigen::PlainObjectBase<DerivedFA> & FA,
  const Eigen::PlainObjectBase<Deriveds> & s,
  const Eigen::PlainObjectBase<Derivedd> & d,
  const bool resolve_overlaps, 
  Eigen::PlainObjectBase<DerivedW> & W,
  Eigen::PlainObjectBase<DerivedG> & G,
  Eigen::PlainObjectBase<DerivedJ> & J)
{
  using namespace Eigen;
  using namespace std;
  // silly base case
  if(FA.size() == 0)
  {
    W.resize(0,3);
    G.resize(0,3);
    return;
  }
  const int dim = VA.cols();
  assert(dim == 3 && "dim must be 3D");
  assert(s.size() == 3 && "s must be 3D point");
  assert(d.size() == 3 && "d must be 3D point");
  // segment vector
  const CGAL::Vector_3<CGAL::Epeck> v(d(0)-s(0),d(1)-s(1),d(2)-s(2));
  // number of vertices
  const int n = VA.rows();
  // duplicate vertices at s and d, we'll remove unreferernced later
  DerivedW WW(2*n,dim);
  for(int i = 0;i<n;i++)
  {
    for(int j = 0;j<dim;j++)
    {
      WW  (i,j) = VA(i,j) + s(j);
      WW(i+n,j) = VA(i,j) + d(j);
    }
  }
  // number of faces
  const int m = FA.rows();
  // Mask whether positive dot product, or negative: because of exactly zero,
  // these are not necessarily complementary
  Matrix<bool,Dynamic,1> P(m,1),N(m,1);
  // loop over faces
  int mp = 0,mn = 0;
  for(int f = 0;f<m;f++)
  {
    const CGAL::Plane_3<CGAL::Epeck> plane(
      CGAL::Point_3<CGAL::Epeck>(VA(FA(f,0),0),VA(FA(f,0),1),VA(FA(f,0),2)),
      CGAL::Point_3<CGAL::Epeck>(VA(FA(f,1),0),VA(FA(f,1),1),VA(FA(f,1),2)),
      CGAL::Point_3<CGAL::Epeck>(VA(FA(f,2),0),VA(FA(f,2),1),VA(FA(f,2),2)));
    const auto normal = plane.orthogonal_vector();
    const auto dt = normal * v;
    if(dt > 0)
    {
      P(f) = true;
      N(f) = false;
      mp++;
    }else if(dt < 0)
    {
      P(f) = false;
      N(f) = true;
      mn++;
    }else
    {
      P(f) = false;
      N(f) = false;
    }
  }

  typedef Matrix<typename DerivedG::Scalar,Dynamic,Dynamic> MatrixXI;
  typedef Matrix<typename DerivedG::Scalar,Dynamic,1> VectorXI;
  MatrixXI GT(mp+mn,3);
  GT<< slice_mask(FA,N,1), slice_mask((FA.array()+n).eval(),P,1);
  // J indexes FA for parts at s and m+FA for parts at d
  J = DerivedJ::LinSpaced(m,0,m-1);
  DerivedJ JT(mp+mn);
  JT << slice_mask(J,P,1), slice_mask(J,N,1);
  JT.block(mp,0,mn,1).array()+=m;

  // Original non-co-planar faces with positively oriented reversed
  MatrixXI BA(mp+mn,3);
  BA << slice_mask(FA,P,1).rowwise().reverse(), slice_mask(FA,N,1);
  // Quads along **all** sides
  MatrixXI GQ((mp+mn)*3,4);
  GQ<< 
    BA.col(1), BA.col(0), BA.col(0).array()+n, BA.col(1).array()+n,
    BA.col(2), BA.col(1), BA.col(1).array()+n, BA.col(2).array()+n,
    BA.col(0), BA.col(2), BA.col(2).array()+n, BA.col(0).array()+n;

  MatrixXI uGQ;
  VectorXI S,sI,sJ;
  //const auto & total_signed_distance = 
  [](
      const MatrixXI & F,
      VectorXI & S,
      MatrixXI & uF,
      VectorXI & I,
      VectorXI & J)
  {
    const int m = F.rows();
    const int d = F.cols();
    MatrixXI sF = F;
    const auto MN = sF.rowwise().minCoeff().eval();
    // rotate until smallest index is first
    for(int p = 0;p<d;p++)
    {
      for(int f = 0;f<m;f++)
      {
        if(sF(f,0) != MN(f))
        {
          for(int r = 0;r<d-1;r++)
          {
            std::swap(sF(f,r),sF(f,r+1));
          }
        }
      }
    }
    // swap orienation
    for(int f = 0;f<m;f++)
    {
      if(sF(f,d-1) < sF(f,1))
      {
        sF.block(f,1,1,d-1) = sF.block(f,1,1,d-1).reverse().eval();
      }
    }
    Matrix<bool,Dynamic,1> M = Matrix<bool,Dynamic,1>::Zero(m,1);
    {
      VectorXI P = VectorXI::LinSpaced(d,0,d-1);
      for(int p = 0;p<d;p++)
      {
        for(int f = 0;f<m;f++)
        {
          bool all = true;
          for(int r = 0;r<d;r++)
          {
            all = all && (sF(f,P(r)) == F(f,r));
          }
          M(f) = M(f) || all;
        }
        for(int r = 0;r<d-1;r++)
        {
          std::swap(P(r),P(r+1));
        }
      }
    }
    unique_rows(sF,uF,I,J);
    S = VectorXI::Zero(uF.rows(),1);
    assert(m == J.rows());
    for(int f = 0;f<m;f++)
    {
      S(J(f)) += M(f) ? 1 : -1;
    }
  }(MatrixXI(GQ),S,uGQ,sI,sJ);
  assert(S.rows() == uGQ.rows());
  const int nq = (S.array().abs()==2).count();
  GQ.resize(nq,4);
  {
    int k = 0;
    for(int q = 0;q<uGQ.rows();q++)
    {
      switch(S(q))
      {
        case -2:
          GQ.row(k++) = uGQ.row(q).reverse().eval();
          break;
        case 2:
          GQ.row(k++) = uGQ.row(q);
          break;
        default:
        // do not add
          break;
      }
    }
    assert(nq == k);
  }

  MatrixXI GG(GT.rows()+2*GQ.rows(),3);
  GG<< 
    GT,
    GQ.col(0), GQ.col(1), GQ.col(2), 
    GQ.col(0), GQ.col(2), GQ.col(3);
  J.resize(JT.rows()+2*GQ.rows(),1);
  J<<JT,DerivedJ::Constant(2*GQ.rows(),1,2*m+1);
  if(resolve_overlaps)
  {
    DerivedJ SJ;
    mesh_boolean(
      WW,GG,
      Matrix<typename DerivedVA::Scalar,Dynamic,Dynamic>(),MatrixXI(),
      MESH_BOOLEAN_TYPE_UNION,
      W,G,SJ);
    J = slice(DerivedJ(J),SJ,1);
  }
}

template <
  typename DerivedVA,
  typename DerivedFA,
  typename Deriveds,
  typename Derivedd,
  typename DerivedW,
  typename DerivedG,
  typename DerivedJ>
IGL_INLINE void igl::copyleft::boolean::minkowski_sum(
  const Eigen::PlainObjectBase<DerivedVA> & VA,
  const Eigen::PlainObjectBase<DerivedFA> & FA,
  const Eigen::PlainObjectBase<Deriveds> & s,
  const Eigen::PlainObjectBase<Derivedd> & d,
  Eigen::PlainObjectBase<DerivedW> & W,
  Eigen::PlainObjectBase<DerivedG> & G,
  Eigen::PlainObjectBase<DerivedJ> & J)
{
  return minkowski_sum(VA,FA,s,d,true,W,G,J);
}
