#ifndef BONE_VISIBLE_H
#define BONE_VISIBLE_H
#include <Eigen/Core>
//
// BONE_VISIBLE  test whether vertices of mesh are "visible" to a given bone,
// where "visible" is defined as in [Baran & Popovic 07].
//
// [flag] = bone_visible(V,F,s,d);
//
// Input:
//    s  row vector of position of start end point of bone
//    d  row vector of position of dest end point of bone
//    V  #V by 3 list of vertex positions
//    F  #F by 3 list of triangle indices
// Output:
//    flag  #V by 1 list of bools (true) visible, (false) obstructed
//
template <
  typename DerivedV, 
  typename DerivedF, 
  typename DerivedSD,
  typename Derivedflag>
void bone_visible(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const Eigen::PlainObjectBase<DerivedSD> & s,
  const Eigen::PlainObjectBase<DerivedSD> & d,
  Eigen::PlainObjectBase<Derivedflag>  & flag);
#endif
