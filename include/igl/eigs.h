#ifndef IGL_EIGS_H
#define IGL_EIGS_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace igl
{
  // Act like MATLAB's eigs function. Compute the first/last k eigen pairs of
  // the generalized eigen value problem:
  //
  //     A u = s B u
  //
  // Solutions are approximate and sorted. 
  //
  // Ideally one should use ARPACK and the Eigen unsupported ARPACK module.
  // This implementation does simple, naive power iterations.
  //
  // Inputs:
  //   A  #A by #A symmetric matrix
  //   B  #A by #A symmetric positive-definite matrix
  //   k  number of eigen pairs to compute
  // Outputs:
  //   sU  #A by k list of sorted eigen vectors (descending)
  //   sS  k list of sorted eigen values (descending)
  //
  // Known issues:
  //   - only one pair per eigen value is found (no multiples)
  //   - only the 'sm' small magnitude eigen values are well supported
  //   
  enum EigsType
  {
    EIGS_TYPE_SM = 0,
    EIGS_TYPE_LM = 1,
    NUM_EIGS_TYPES = 2
  };
  template <
    typename Atype,
    typename Btype,
    typename DerivedU,
    typename DerivedS>
  IGL_INLINE bool eigs(
    const Eigen::SparseMatrix<Atype> & A,
    const Eigen::SparseMatrix<Btype> & B,
    const size_t k,
    const EigsType type,
    Eigen::PlainObjectBase<DerivedU> & sU,
    Eigen::PlainObjectBase<DerivedS> & sS);
}

#ifndef IGL_STATIC_LIBRARY
#include "eigs.cpp"
#endif
#endif
