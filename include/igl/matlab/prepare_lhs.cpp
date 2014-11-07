#include "prepare_lhs.h"
#include <algorithm>
template <typename DerivedV>
IGL_INLINE void igl::prepare_lhs_double(
  const Eigen::PlainObjectBase<DerivedV> & V,
  mxArray *plhs[])
{
  using namespace std;
  plhs[0] = mxCreateDoubleMatrix(V.rows(),V.cols(), mxREAL);
  double * Vp = mxGetPr(plhs[0]);
  copy(&V.data()[0],&V.data()[0]+V.size(),Vp);
}

template <typename DerivedV>
IGL_INLINE void igl::prepare_lhs_logical(
  const Eigen::PlainObjectBase<DerivedV> & V,
  mxArray *plhs[])
{
  using namespace std;
  plhs[0] = mxCreateLogicalMatrix(V.rows(),V.cols());
  mxLogical * Vp = static_cast<mxLogical*>(mxGetData(plhs[0]));
  copy(&V.data()[0],&V.data()[0]+V.size(),Vp);
}

template <typename DerivedV>
IGL_INLINE void igl::prepare_lhs_index(
  const Eigen::PlainObjectBase<DerivedV> & V,
  mxArray *plhs[])
{
  // Treat indices as reals
  
  const auto Vd = (V.template cast<double>().array()+1).eval();
  return prepare_lhs_double(Vd,plhs);
}

#ifdef IGL_STATIC_LIBRARY
template void igl::prepare_lhs_index<Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> > const&, mxArray_tag**);
template void igl::prepare_lhs_index<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, mxArray_tag**);
template void igl::prepare_lhs_double<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, mxArray_tag**);
template void igl::prepare_lhs_index<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, mxArray_tag**);
#endif
