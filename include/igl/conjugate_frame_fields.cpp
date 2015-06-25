// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <igl/conjugate_frame_fields.h>
#include <igl/speye.h>
#include <igl/slice.h>
#include <igl/polyroots.h>
#include <Eigen/Sparse>

#include <iostream>

namespace igl {
  template <typename DerivedV, typename DerivedF, typename DerivedO>
  class ConjugateFFSolver
  {
  public:
    IGL_INLINE ConjugateFFSolver(const ConjugateFFSolverData<DerivedV, DerivedF> &_data,
                                 int _maxIter = 50,
                                 const typename DerivedV::Scalar &_lambdaOrtho = .1,
                                 const typename DerivedV::Scalar &_lambdaInit = 100,
                                 const typename DerivedV::Scalar &_lambdaMultFactor = 1.01,
                                 bool _doHardConstraints = true);
    IGL_INLINE bool solve(const Eigen::VectorXi &isConstrained,
                          const Eigen::PlainObjectBase<DerivedO> &initialSolution,
                          Eigen::PlainObjectBase<DerivedO> &output,
                          typename DerivedV::Scalar *lambdaOut = NULL);

  private:

    const ConjugateFFSolverData<DerivedV, DerivedF> &data;

    //polyVF data
    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> Acoeff, Bcoeff;
    Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 2> pvU, pvV;
    typename DerivedV::Scalar lambda;

    //parameters
    typename DerivedV::Scalar lambdaOrtho;
    typename DerivedV::Scalar lambdaInit,lambdaMultFactor;
    int maxIter;
    bool doHardConstraints;

    IGL_INLINE void localStep();
    IGL_INLINE void getPolyCoeffsForLocalSolve(const Eigen::Matrix<typename DerivedV::Scalar, 4, 1> &s,
                                               const Eigen::Matrix<typename DerivedV::Scalar, 4, 1> &z,
                                               Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> &polyCoeff);

    IGL_INLINE void globalStep(const Eigen::Matrix<int, Eigen::Dynamic, 1>  &isConstrained,
                               const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1>  &Ak,
                               const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1>  &Bk);
    IGL_INLINE void minQuadWithKnownMini(const Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &Q,
                         const Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &f,
                         const Eigen::VectorXi isConstrained,
                         const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> &xknown,
                                         Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> &x);
    IGL_INLINE void setFieldFromCoefficients();
    IGL_INLINE void setCoefficientsFromField();

  };
}

//Implementation
/***************************** Solver ***********************************/
template <typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE igl::ConjugateFFSolver<DerivedV, DerivedF, DerivedO>::
ConjugateFFSolver(const ConjugateFFSolverData<DerivedV, DerivedF> &_data,
                  int _maxIter,
                  const typename DerivedV::Scalar &_lambdaOrtho,
                  const typename DerivedV::Scalar &_lambdaInit,
                  const typename DerivedV::Scalar &_lambdaMultFactor,
                  bool _doHardConstraints):
data(_data),
lambdaOrtho(_lambdaOrtho),
lambdaInit(_lambdaInit),
maxIter(_maxIter),
lambdaMultFactor(_lambdaMultFactor),
doHardConstraints(_doHardConstraints)
{
  Acoeff.resize(data.numF,1);
  Bcoeff.resize(data.numF,1);
  pvU.setZero(data.numF, 2);
  pvV.setZero(data.numF, 2);
};



template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::ConjugateFFSolver<DerivedV, DerivedF, DerivedO>::
getPolyCoeffsForLocalSolve(const Eigen::Matrix<typename DerivedV::Scalar, 4, 1> &s,
                           const Eigen::Matrix<typename DerivedV::Scalar, 4, 1> &z,
                           Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> &polyCoeff)
{
  typename DerivedV::Scalar s0 = s(0);
  typename DerivedV::Scalar s1 = s(1);
  typename DerivedV::Scalar s2 = s(2);
  typename DerivedV::Scalar s3 = s(3);
  typename DerivedV::Scalar z0 = z(0);
  typename DerivedV::Scalar z1 = z(1);
  typename DerivedV::Scalar z2 = z(2);
  typename DerivedV::Scalar z3 = z(3);

  polyCoeff.resize(7,1);
  polyCoeff(0) =  s0*s0* s1*s1* s2*s2* s3* z3*z3 +  s0*s0* s1*s1* s2* s3*s3* z2*z2 +  s0*s0* s1* s2*s2* s3*s3* z1*z1 +  s0* s1*s1* s2*s2* s3*s3* z0*z0 ;
  polyCoeff(1) = 2* s0*s0* s1*s1* s2* s3* z2*z2 + 2* s0*s0* s1*s1* s2* s3* z3*z3 + 2* s0*s0* s1* s2*s2* s3* z1*z1 + 2* s0*s0* s1* s2*s2* s3* z3*z3 + 2* s0*s0* s1* s2* s3*s3* z1*z1 + 2* s0*s0* s1* s2* s3*s3* z2*z2 + 2* s0* s1*s1* s2*s2* s3* z0*z0 + 2* s0* s1*s1* s2*s2* s3* z3*z3 + 2* s0* s1*s1* s2* s3*s3* z0*z0 + 2* s0* s1*s1* s2* s3*s3* z2*z2 + 2* s0* s1* s2*s2* s3*s3* z0*z0 + 2* s0* s1* s2*s2* s3*s3* z1*z1 ;
  polyCoeff(2) =  s0*s0* s1*s1* s2* z2*z2 +  s0*s0* s1*s1* s3* z3*z3 +  s0*s0* s1* s2*s2* z1*z1 + 4* s0*s0* s1* s2* s3* z1*z1 + 4* s0*s0* s1* s2* s3* z2*z2 + 4* s0*s0* s1* s2* s3* z3*z3 +  s0*s0* s1* s3*s3* z1*z1 +  s0*s0* s2*s2* s3* z3*z3 +  s0*s0* s2* s3*s3* z2*z2 +  s0* s1*s1* s2*s2* z0*z0 + 4* s0* s1*s1* s2* s3* z0*z0 + 4* s0* s1*s1* s2* s3* z2*z2 + 4* s0* s1*s1* s2* s3* z3*z3 +  s0* s1*s1* s3*s3* z0*z0 + 4* s0* s1* s2*s2* s3* z0*z0 + 4* s0* s1* s2*s2* s3* z1*z1 + 4* s0* s1* s2*s2* s3* z3*z3 + 4* s0* s1* s2* s3*s3* z0*z0 + 4* s0* s1* s2* s3*s3* z1*z1 + 4* s0* s1* s2* s3*s3* z2*z2 +  s0* s2*s2* s3*s3* z0*z0 +  s1*s1* s2*s2* s3* z3*z3 +  s1*s1* s2* s3*s3* z2*z2 +  s1* s2*s2* s3*s3* z1*z1;
  polyCoeff(3) = 2* s0*s0* s1* s2* z1*z1 + 2* s0*s0* s1* s2* z2*z2 + 2* s0*s0* s1* s3* z1*z1 + 2* s0*s0* s1* s3* z3*z3 + 2* s0*s0* s2* s3* z2*z2 + 2* s0*s0* s2* s3* z3*z3 + 2* s0* s1*s1* s2* z0*z0 + 2* s0* s1*s1* s2* z2*z2 + 2* s0* s1*s1* s3* z0*z0 + 2* s0* s1*s1* s3* z3*z3 + 2* s0* s1* s2*s2* z0*z0 + 2* s0* s1* s2*s2* z1*z1 + 8* s0* s1* s2* s3* z0*z0 + 8* s0* s1* s2* s3* z1*z1 + 8* s0* s1* s2* s3* z2*z2 + 8* s0* s1* s2* s3* z3*z3 + 2* s0* s1* s3*s3* z0*z0 + 2* s0* s1* s3*s3* z1*z1 + 2* s0* s2*s2* s3* z0*z0 + 2* s0* s2*s2* s3* z3*z3 + 2* s0* s2* s3*s3* z0*z0 + 2* s0* s2* s3*s3* z2*z2 + 2* s1*s1* s2* s3* z2*z2 + 2* s1*s1* s2* s3* z3*z3 + 2* s1* s2*s2* s3* z1*z1 + 2* s1* s2*s2* s3* z3*z3 + 2* s1* s2* s3*s3* z1*z1 + 2* s1* s2* s3*s3* z2*z2 ;
  polyCoeff(4) =  s0*s0* s1* z1*z1 +  s0*s0* s2* z2*z2 +  s0*s0* s3* z3*z3 +  s0* s1*s1* z0*z0 + 4* s0* s1* s2* z0*z0 + 4* s0* s1* s2* z1*z1 + 4* s0* s1* s2* z2*z2 + 4* s0* s1* s3* z0*z0 + 4* s0* s1* s3* z1*z1 + 4* s0* s1* s3* z3*z3 +  s0* s2*s2* z0*z0 + 4* s0* s2* s3* z0*z0 + 4* s0* s2* s3* z2*z2 + 4* s0* s2* s3* z3*z3 +  s0* s3*s3* z0*z0 +  s1*s1* s2* z2*z2 +  s1*s1* s3* z3*z3 +  s1* s2*s2* z1*z1 + 4* s1* s2* s3* z1*z1 + 4* s1* s2* s3* z2*z2 + 4* s1* s2* s3* z3*z3 +  s1* s3*s3* z1*z1 +  s2*s2* s3* z3*z3 +  s2* s3*s3* z2*z2;
  polyCoeff(5) = 2* s0* s1* z0*z0 + 2* s0* s1* z1*z1 + 2* s0* s2* z0*z0 + 2* s0* s2* z2*z2 + 2* s0* s3* z0*z0 + 2* s0* s3* z3*z3 + 2* s1* s2* z1*z1 + 2* s1* s2* z2*z2 + 2* s1* s3* z1*z1 + 2* s1* s3* z3*z3 + 2* s2* s3* z2*z2 + 2* s2* s3* z3*z3 ;
  polyCoeff(6) =  s0* z0*z0 +  s1* z1*z1 +  s2* z2*z2 +  s3* z3*z3;

}


template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::ConjugateFFSolver<DerivedV, DerivedF, DerivedO>::
localStep()
{
  for (int j =0; j<data.numF; ++j)
  {
    Eigen::Matrix<typename DerivedV::Scalar, 4, 1> xproj; xproj << pvU.row(j).transpose(),pvV.row(j).transpose();
    Eigen::Matrix<typename DerivedV::Scalar, 4, 1> z = data.UH[j].transpose()*xproj;
    Eigen::Matrix<typename DerivedV::Scalar, 4, 1> x;

    Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> polyCoeff;
    getPolyCoeffsForLocalSolve(data.s[j], z, polyCoeff);
    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> roots;
    igl::polyRoots<typename DerivedV::Scalar, typename DerivedV::Scalar> (polyCoeff, roots );

    //  find closest real root to xproj
    typename DerivedV::Scalar minDist = 1e10;
    for (int i =0; i< 6; ++i)
    {
      if (fabs(imag(roots[i]))>1e-10)
        continue;
      Eigen::Matrix<typename DerivedV::Scalar, 4, 4> D = ((Eigen::Matrix<typename DerivedV::Scalar, 4, 1>::Ones()+real(roots(i))*data.s[j]).array().inverse()).matrix().asDiagonal();
      Eigen::Matrix<typename DerivedV::Scalar, 4, 1> candidate = data.UH[j]*D*z;
      typename DerivedV::Scalar dist = (candidate-xproj).norm();
      if (dist<minDist)
      {
        minDist = dist;
        x = candidate;
      }

    }

    pvU.row(j) << x(0),x(1);
    pvV.row(j) << x(2),x(3);
  }
}


template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::ConjugateFFSolver<DerivedV, DerivedF, DerivedO>::
setCoefficientsFromField()
{
  for (int i = 0; i <data.numF; ++i)
  {
    std::complex<typename DerivedV::Scalar> u(pvU(i,0),pvU(i,1));
    std::complex<typename DerivedV::Scalar> v(pvV(i,0),pvV(i,1));
    Acoeff(i) = u*u+v*v;
    Bcoeff(i) = u*u*v*v;
  }
}


template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::ConjugateFFSolver<DerivedV, DerivedF, DerivedO>::
globalStep(const Eigen::Matrix<int, Eigen::Dynamic, 1>  &isConstrained,
           const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1>  &Ak,
           const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1>  &Bk)
{
  setCoefficientsFromField();

  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > I;
  igl::speye(data.numF, data.numF, I);
  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > QA = data.DDA+lambda*data.planarityWeight+lambdaOrtho*I;
  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > fA = (-2*lambda*data.planarityWeight*Acoeff).sparseView();

  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > QB = data.DDB+lambda*data.planarityWeight;
  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > fB = (-2*lambda*data.planarityWeight*Bcoeff).sparseView();

  if(doHardConstraints)
  {
    minQuadWithKnownMini(QA, fA, isConstrained, Ak, Acoeff);
    minQuadWithKnownMini(QB, fB, isConstrained, Bk, Bcoeff);
  }
  else
  {
    Eigen::Matrix<int, Eigen::Dynamic, 1>isknown_; isknown_.setZero(data.numF,1);
    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> xknown_; xknown_.setZero(0,1);
    minQuadWithKnownMini(QA, fA, isknown_, xknown_, Acoeff);
    minQuadWithKnownMini(QB, fB, isknown_, xknown_, Bcoeff);
  }
  setFieldFromCoefficients();

}


template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::ConjugateFFSolver<DerivedV, DerivedF, DerivedO>::
setFieldFromCoefficients()
{
  for (int i = 0; i <data.numF; ++i)
  {
    //    poly coefficients: 1, 0, -Acoeff, 0, Bcoeff
    //    matlab code from roots (given there are no trailing zeros in the polynomial coefficients)
    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> polyCoeff(5,1);
    polyCoeff<<1., 0., -Acoeff(i), 0., Bcoeff(i);

    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> roots;
    polyRoots<std::complex<typename DerivedV::Scalar>>(polyCoeff,roots);

    std::complex<typename DerivedV::Scalar> u = roots[0];
    int maxi = -1;
    float maxd = -1;
    for (int k =1; k<4; ++k)
    {
      float dist = abs(roots[k]+u);
      if (dist>maxd)
      {
        maxd = dist;
        maxi = k;
      }
    }
    std::complex<typename DerivedV::Scalar> v = roots[maxi];
    pvU(i,0) = real(u); pvU(i,1) = imag(u);
    pvV(i,0) = real(v); pvV(i,1) = imag(v);
  }

}

template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::ConjugateFFSolver<DerivedV, DerivedF, DerivedO>::
minQuadWithKnownMini(const Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &Q,
                     const Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &f,
                     const Eigen::VectorXi isConstrained,
                     const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> &xknown,
                     Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> &x)
{
  int N = Q.rows();

  int nc = xknown.rows();
  Eigen::VectorXi known; known.setZero(nc,1);
  Eigen::VectorXi unknown; unknown.setZero(N-nc,1);

  int indk = 0, indu = 0;
  for (int i = 0; i<N; ++i)
    if (isConstrained[i])
    {
      known[indk] = i;
      indk++;
    }
    else
    {
      unknown[indu] = i;
      indu++;
    }

  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar>> Quu, Quk;

  igl::slice(Q,unknown, unknown, Quu);
  igl::slice(Q,unknown, known, Quk);


  std::vector<typename Eigen::Triplet<std::complex<typename DerivedV::Scalar> > > tripletList;

  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > fu(N-nc,1);

  igl::slice(f,unknown, Eigen::VectorXi::Zero(1,1), fu);

  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > rhs = (Quk*xknown).sparseView()+.5*fu;

  Eigen::SparseLU< Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar>>> solver;
  solver.compute(-Quu);
  if(solver.info()!=Eigen::Success)
  {
    std::cerr<<"Decomposition failed!"<<std::endl;
    return;
  }
  Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar>>  b  = solver.solve(rhs);
  if(solver.info()!=Eigen::Success)
  {
    std::cerr<<"Solving failed!"<<std::endl;
    return;
  }

  indk = 0, indu = 0;
  x.setZero(N,1);
  for (int i = 0; i<N; ++i)
    if (isConstrained[i])
      x[i] = xknown[indk++];
    else
      x[i] = b.coeff(indu++,0);

}


template<typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE bool igl::ConjugateFFSolver<DerivedV, DerivedF, DerivedO>::
solve(const Eigen::VectorXi &isConstrained,
      const Eigen::PlainObjectBase<DerivedO> &initialSolution,
      Eigen::PlainObjectBase<DerivedO> &output,
      typename DerivedV::Scalar *lambdaOut)
{
  int numConstrained = isConstrained.sum();
  // coefficient values
  Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> Ak, Bk;

  pvU.resize(data.numF,2);
  pvV.resize(data.numF,2);
  for (int fi = 0; fi <data.numF; ++fi)
  {
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b1 = data.B1.row(fi);
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b2 = data.B2.row(fi);
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &u3 = initialSolution.block(fi,0,1,3);
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &v3 = initialSolution.block(fi,3,1,3);
    pvU.row(fi)<< u3.dot(b1), u3.dot(b2);
    pvV.row(fi)<< v3.dot(b1), v3.dot(b2);
  }
  setCoefficientsFromField();
  Ak.resize(numConstrained,1);
  Bk.resize(numConstrained,1);
  int ind = 0;
  for (int i = 0; i <data.numF; ++i)
  {
    if(isConstrained[i])
    {
      Ak(ind) = Acoeff[i];
      Bk(ind) = Bcoeff[i];
      ind ++;
    }
  }



  typename DerivedV::Scalar smoothnessValue;
  Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> conjValues;
  typename DerivedV::Scalar meanConj;
  typename DerivedV::Scalar maxConj;

  data.evaluateConjugacy(pvU, pvV, conjValues);
  meanConj = conjValues.cwiseAbs().mean();
  maxConj = conjValues.cwiseAbs().maxCoeff();
  printf("Initial max non-conjugacy: %.5g\n",maxConj);

  smoothnessValue = (Acoeff.adjoint()*data.DDA*Acoeff + Bcoeff.adjoint()*data.DDB*Bcoeff).real()[0];
  printf("\n\nInitial smoothness: %.5g\n",smoothnessValue);

  lambda = lambdaInit;

  bool doit = false;
  for (int iter = 0; iter<maxIter; ++iter)
  {
    printf("\n\n--- Iteration %d ---\n",iter);

    typename DerivedV::Scalar oldMeanConj = meanConj;

    localStep();
    globalStep(isConstrained, Ak, Bk);


    smoothnessValue = (Acoeff.adjoint()*data.DDA*Acoeff + Bcoeff.adjoint()*data.DDB*Bcoeff).real()[0];

    printf("Smoothness: %.5g\n",smoothnessValue);

    data.evaluateConjugacy(pvU, pvV, conjValues);
    meanConj = conjValues.cwiseAbs().mean();
    maxConj = conjValues.cwiseAbs().maxCoeff();
    printf("Mean/Max non-conjugacy: %.5g, %.5g\n",meanConj,maxConj);
    typename DerivedV::Scalar diffMeanConj = fabs(oldMeanConj-meanConj);

    if (diffMeanConj<1e-4)
      doit = true;

    if (doit)
      lambda = lambda*lambdaMultFactor;
    printf(" %d %.5g %.5g\n",iter, smoothnessValue,maxConj);

  }

  output.setZero(data.numF,6);
  for (int fi=0; fi<data.numF; ++fi)
  {
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b1 = data.B1.row(fi);
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b2 = data.B2.row(fi);
    output.block(fi,0, 1, 3) = pvU(fi,0)*b1 + pvU(fi,1)*b2;
    output.block(fi,3, 1, 3) = pvV(fi,0)*b1 + pvV(fi,1)*b2;
  }

  if (lambdaOut)
    *lambdaOut = lambda;

  return true;
}



template <typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::conjugate_frame_fields(const Eigen::PlainObjectBase<DerivedV> &V,
                                            const Eigen::PlainObjectBase<DerivedF> &F,
                                            const Eigen::VectorXi &isConstrained,
                                            const Eigen::PlainObjectBase<DerivedO> &initialSolution,
                                            Eigen::PlainObjectBase<DerivedO> &output,
                                            int maxIter,
                                            const typename DerivedV::Scalar &lambdaOrtho,
                                            const typename DerivedV::Scalar &lambdaInit,
                                            const typename DerivedV::Scalar &lambdaMultFactor,
                                            bool doHardConstraints)
{
  igl::ConjugateFFSolverData<DerivedV, DerivedF> csdata(V, F);
  igl::ConjugateFFSolver<DerivedV, DerivedF, DerivedO> cs(csdata, maxIter, lambdaOrtho, lambdaInit, lambdaMultFactor, doHardConstraints);
  cs.solve(isConstrained, initialSolution, output);
}

template <typename DerivedV, typename DerivedF, typename DerivedO>
IGL_INLINE void igl::conjugate_frame_fields(const igl::ConjugateFFSolverData<DerivedV, DerivedF> &csdata,
                                            const Eigen::VectorXi &isConstrained,
                                            const Eigen::PlainObjectBase<DerivedO> &initialSolution,
                                            Eigen::PlainObjectBase<DerivedO> &output,
                                            int maxIter,
                                            const typename DerivedV::Scalar &lambdaOrtho,
                                            const typename DerivedV::Scalar &lambdaInit,
                                            const typename DerivedV::Scalar &lambdaMultFactor,
                                            bool doHardConstraints,
                                            typename DerivedV::Scalar *lambdaOut)
{
  igl::ConjugateFFSolver<DerivedV, DerivedF, DerivedO> cs(csdata, maxIter, lambdaOrtho, lambdaInit, lambdaMultFactor, doHardConstraints);
  cs.solve(isConstrained, initialSolution, output, lambdaOut);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::conjugate_frame_fields<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(igl::ConjugateFFSolverData<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, int, Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar const&, Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar const&, Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar const&, bool, Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar*);
#endif
