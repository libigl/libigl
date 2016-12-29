#include "ramer_douglas_peucker.h"

#include "find.h"
#include "project_to_line.h"
#include "EPS.h"
#include "slice_mask.h"
template <typename DerivedP, typename DerivedS, typename DerivedJ>
IGL_INLINE void igl::ramer_douglas_peucker(
  const Eigen::MatrixBase<DerivedP> & P,
  const typename DerivedP::Scalar tol,
  Eigen::PlainObjectBase<DerivedS> & S,
  Eigen::PlainObjectBase<DerivedJ> & J)
{
  typedef typename DerivedP::Scalar Scalar;
  // number of vertices
  const int n = P.rows();
  // Trivial base case
  if(n <= 1)
  {
    J = DerivedJ::Zero(n);
    S = P;
    return;
  }
  // number of dimensions
  const int m = P.cols();
  Eigen::Array<bool,Eigen::Dynamic,1> I =
    Eigen::Array<bool,Eigen::Dynamic,1>::Constant(n,1,true);
  const auto stol = tol*tol;
  std::function<void(const int,const int)> simplify;
  simplify = [&I,&P,&stol,&simplify](const int ixs, const int ixe)->void
  {
    assert(ixe>ixs);
    Scalar sdmax = 0;
    typename Eigen::Matrix<Scalar,Eigen::Dynamic,1>::Index ixc = -1;
    if((ixe-ixs)>1)
    {
      Scalar sdes = (P.row(ixe)-P.row(ixs)).squaredNorm();
      Eigen::Matrix<Scalar,Eigen::Dynamic,1> sD;
      const auto & Pblock = P.block(ixs+1,0,((ixe+1)-ixs)-2,P.cols());
      if(sdes<=EPS<Scalar>())
      {
        sD = (Pblock.rowwise()-P.row(ixs)).rowwise().squaredNorm();
      }else
      {
        Eigen::Matrix<Scalar,Eigen::Dynamic,1> T;
        project_to_line(Pblock,P.row(ixs).eval(),P.row(ixe).eval(),T,sD);
      }
      sdmax = sD.maxCoeff(&ixc);
      // Index full P
      ixc = ixc+(ixs+1);
    }
    if(sdmax <= stol)
    {
      if(ixs != ixe-1)
      {
        I.block(ixs+1,0,((ixe+1)-ixs)-2,1).setConstant(false);
      }
    }else
    {
      simplify(ixs,ixc);
      simplify(ixc,ixe);
    }
  };
  simplify(0,n-1);
  slice_mask(P,I,1,S);
  find(I,J);
}
