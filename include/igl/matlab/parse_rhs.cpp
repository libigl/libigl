#include "parse_rhs.h"
#include <algorithm>

template <typename DerivedV>
IGL_INLINE void igl::matlab::parse_rhs_double(
    const mxArray *prhs[], 
    Eigen::PlainObjectBase<DerivedV> & V)
{
  using namespace std;
  using namespace Eigen;
  // set number of mesh vertices
  const int n = mxGetM(prhs[0]);
  // set vertex position pointers
  double * Vp = mxGetPr(prhs[0]);
  const int dim = mxGetN(prhs[0]);

  typedef typename DerivedV::Scalar Scalar;
  Matrix<Scalar, DerivedV::ColsAtCompileTime, DerivedV::RowsAtCompileTime, RowMajor> VT;
  Scalar * V_data;
  if(DerivedV::IsRowMajor)
  {
    VT.resize(dim,n);
    V_data = VT.data();
  }else
  {
    V.resize(n,dim);
    V_data = V.data();
  }
  copy(Vp,Vp+n*dim,V_data);
  if(DerivedV::IsRowMajor)
  {
    V = VT.transpose();
  }
}

template <typename DerivedV>
IGL_INLINE void igl::matlab::parse_rhs_index(
    const mxArray *prhs[], 
    Eigen::PlainObjectBase<DerivedV> & V)
{
  parse_rhs_double(prhs,V);
  V.array() -= 1;
}

#ifdef IGL_STATIC_LIBRARY
template void igl::matlab::parse_rhs_index<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(mxArray_tag const**, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::matlab::parse_rhs_index<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(mxArray_tag const**, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::matlab::parse_rhs_double<Eigen::Matrix<double, -1, -1, 0, -1, -1> >(mxArray_tag const**, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
