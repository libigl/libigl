#include "setxor.h"
#include "setdiff.h"
#include "setunion.h"
#include "slice.h"
#include <iostream>
#include <igl/matlab_format.h>

template <
  typename DerivedA,
  typename DerivedB,
  typename DerivedC,
  typename DerivedIA,
  typename DerivedIB>
IGL_INLINE void igl::setxor(
  const Eigen::DenseBase<DerivedA> & A,
  const Eigen::DenseBase<DerivedB> & B,
  Eigen::PlainObjectBase<DerivedC> & C,
  Eigen::PlainObjectBase<DerivedIA> & IA,
  Eigen::PlainObjectBase<DerivedIB> & IB)
{
  DerivedC AB,BA;
  DerivedIA IAB,IBA;
  setdiff(A,B,AB,IAB);
  std::cout<<igl::matlab_format(AB.transpose(),"AB")<<std::endl;
  std::cout<<igl::matlab_format(IAB.transpose().array()+1,"IAB")<<std::endl;
  setdiff(B,A,BA,IBA);
  std::cout<<igl::matlab_format(BA.transpose(),"BA")<<std::endl;
  std::cout<<igl::matlab_format(IBA.transpose().array()+1,"IBA")<<std::endl;
  setunion(AB,BA,C,IA,IB);

  std::cout<<igl::matlab_format(C.transpose(),"C")<<std::endl;
  std::cout<<igl::matlab_format(IA.transpose().array()+1,"IA")<<std::endl;
  std::cout<<igl::matlab_format(IB.transpose().array()+1,"IB")<<std::endl;
  slice(IAB,DerivedIA(IA),IA);
  slice(IBA,DerivedIB(IB),IB);

  std::cout<<igl::matlab_format(IA.transpose().array()+1,"IA")<<std::endl;
  std::cout<<igl::matlab_format(IB.transpose().array()+1,"IB")<<std::endl;
}
