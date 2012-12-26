#include "remove_duplicate_vertices.h"
#include "round.h"
#include "unique.h"
#include "colon.h"
#include "slice.h"
#include <functional>

template <
  typename DerivedV, 
  typename DerivedSV, 
  typename DerivedSVI, 
  typename DerivedSVJ>
IGL_INLINE void igl::remove_duplicate_vertices(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const double epsilon,
  Eigen::PlainObjectBase<DerivedSV>& SV,
  Eigen::PlainObjectBase<DerivedSVI>& SVI,
  Eigen::PlainObjectBase<DerivedSVJ>& SVJ)
{
  using namespace igl;
  if(epsilon > 0)
  {
    Eigen::PlainObjectBase<DerivedV> rV,rSV;
    round((V/(10.0*epsilon)).eval(),rV);
    unique_rows(rV,rSV,SVI,SVJ);
    slice(V,SVI,colon<int>(0,V.cols()-1),SV);
  }else
  {
    unique_rows(V,SV,SVI,SVJ);
  }
}

template <
  typename DerivedV, 
  typename DerivedF,
  typename DerivedSV, 
  typename DerivedSVI, 
  typename DerivedSVJ,
  typename DerivedSF>
IGL_INLINE void igl::remove_duplicate_vertices(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  const double epsilon,
  Eigen::PlainObjectBase<DerivedSV>& SV,
  Eigen::PlainObjectBase<DerivedSVI>& SVI,
  Eigen::PlainObjectBase<DerivedSVJ>& SVJ,
  Eigen::PlainObjectBase<DerivedSF>& SF)
{
  using namespace Eigen;
  using namespace std;
  remove_duplicate_vertices(V,epsilon,SV,SVI,SVJ);
  // remap faces
  SF = F.unaryExpr(bind1st(mem_fun( 
    static_cast<VectorXi::Scalar&(VectorXi::*)(VectorXi::Index)>
      (&VectorXi::operator())), &SVJ)).eval();
}

#ifndef IGL_HEADER_ONLY
// Explicit instanciation
template void igl::remove_duplicate_vertices<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, double, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
