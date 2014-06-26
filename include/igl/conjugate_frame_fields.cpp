// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include <igl/conjugate_frame_fields.h>
#include <igl/edge_topology.h>
#include <igl/local_basis.h>
#include <igl/nchoosek.h>
#include <igl/sparse.h>
#include <igl/speye.h>
#include <igl/slice.h>
#include <igl/polyroots.h>
#include <igl/colon.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/principal_curvature.h>
#include <Eigen/Sparse>

#include <iostream>

namespace igl {

  template <typename DerivedV, typename DerivedF>
  class ConjugateFFSolverData
  {
    public:
      const Eigen::PlainObjectBase<DerivedV> &V; int numV;
      const Eigen::PlainObjectBase<DerivedF> &F; int numF;

      Eigen::MatrixXi EV; int numE;
      Eigen::MatrixXi F2E;
      Eigen::MatrixXi E2F;
      Eigen::VectorXd K;

      Eigen::VectorXi isBorderEdge;
      int numInteriorEdges;
      Eigen::Matrix<int,Eigen::Dynamic,2> E2F_int;
      Eigen::VectorXi indInteriorToFull;
      Eigen::VectorXi indFullToInterior;

      Eigen::PlainObjectBase<DerivedV> B1, B2, FN;


      Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic,1> kmin, kmax;
      Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic,2> dmin, dmax;
      Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic,3> dmin3, dmax3;

      Eigen::VectorXd nonPlanarityMeasure;
      Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > planarityWeight;

      //conjugacy matrix
      std::vector<Eigen::Matrix<typename DerivedV::Scalar, 4,4> > H;

      //conjugacy matrix eigenvectors and (scaled) eigenvalues
      std::vector<Eigen::Matrix<typename DerivedV::Scalar, 4,4> > UH;
      std::vector<Eigen::Matrix<typename DerivedV::Scalar, 4,1> > s;

      //laplacians
      Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar>> DDA, DDB;

  private:
    IGL_INLINE void computeCurvatureAndPrincipals();
    IGL_INLINE void precomputeConjugacyStuff();
    IGL_INLINE void computeLaplacians();
    IGL_INLINE void computek();
    IGL_INLINE void computeCoefficientLaplacian(int n, Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &D);
    
    IGL_INLINE void precomputeInteriorEdges();

public:
      IGL_INLINE ConjugateFFSolverData(const Eigen::PlainObjectBase<DerivedV> &_V,
                                   const Eigen::PlainObjectBase<DerivedF> &_F);
  };

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



    IGL_INLINE void evaluateConjugacy(Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> &conjValues);



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
/***************************** Data ***********************************/

template <typename DerivedV, typename DerivedF>
IGL_INLINE igl::ConjugateFFSolverData<DerivedV, DerivedF>::
ConjugateFFSolverData(const Eigen::PlainObjectBase<DerivedV> &_V,
                  const Eigen::PlainObjectBase<DerivedF> &_F):
V(_V),
numV(_V.rows()),
F(_F),
numF(_F.rows())
{
  igl::edge_topology(V,F,EV,F2E,E2F);
  numE = EV.rows();

  precomputeInteriorEdges();

  igl::local_basis(V,F,B1,B2,FN);

  computek();

  computeLaplacians();

  computeCurvatureAndPrincipals();
  precomputeConjugacyStuff();

};


template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::ConjugateFFSolverData<DerivedV, DerivedF>::computeCurvatureAndPrincipals()
{
  Eigen::MatrixXd VCBary;
  Eigen::MatrixXi FCBary;

  VCBary.setZero(numV+numF,3);
  FCBary.setZero(3*numF,3);
  igl::false_barycentric_subdivision(V, F, VCBary, FCBary);

  Eigen::MatrixXd dmax3_,dmin3_;
  igl::principal_curvature(VCBary, FCBary, dmax3_, dmin3_, kmax, kmin, 5,true);

  dmax3 = dmax3_.bottomRows(numF);
  dmin3 = dmin3_.bottomRows(numF);

  kmax = kmax.bottomRows(numF);
  kmin = kmin.bottomRows(numF);

  //  kmax = dmax3.rowwise().norm();
  //  kmin = dmin3.rowwise().norm();

  dmin3.rowwise().normalize();
  dmax3.rowwise().normalize();
  dmax.setZero(numF,2);
  dmin.setZero(numF,2);
  for (int i= 0; i <numF; ++i)
  {
    if(kmin[i] != kmin[i] || kmax[i] != kmax[i] || (dmin3.row(i).array() != dmin3.row(i).array()).any() || (dmax3.row(i).array() != dmax3.row(i).array()).any())
    {
      kmin[i] = 0;
      kmax[i] = 0;
      dmin3.row(i) = B1.row(i);
      dmax3.row(i) = B2.row(i);
    }
    else
    {
      dmax3.row(i) = (dmax3.row(i) - (dmax3.row(i).dot(FN.row(i)))*FN.row(i)).normalized();
      dmin3.row(i) = dmin3.row(i) - (dmin3.row(i).dot(FN.row(i)))*FN.row(i);
      dmin3.row(i) = (dmin3.row(i) - (dmin3.row(i).dot(dmax3.row(i)))*dmax3.row(i)).normalized();
      if ((dmin3.row(i).cross(dmax3.row(i))).dot(FN.row(i))<0)
        dmin3.row(i) = -dmin3.row(i);
    }
    dmax.row(i) << dmax3.row(i).dot(B1.row(i)), dmax3.row(i).dot(B2.row(i));
    dmax.row(i).normalize();
    dmin.row(i) << dmin3.row(i).dot(B1.row(i)), dmin3.row(i).dot(B2.row(i));
    dmin.row(i).normalize();

  }

  nonPlanarityMeasure = kmax.cwiseAbs().array()*kmin.cwiseAbs().array();
  typename DerivedV::Scalar minP = nonPlanarityMeasure.minCoeff();
  typename DerivedV::Scalar maxP = nonPlanarityMeasure.maxCoeff();
  nonPlanarityMeasure = (nonPlanarityMeasure.array()-minP)/(maxP-minP);
  Eigen::VectorXi I = igl::colon<typename DerivedF::Scalar>(0, numF-1);
  igl::sparse(I, I, nonPlanarityMeasure, numF, numF, planarityWeight);

}

template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::ConjugateFFSolverData<DerivedV, DerivedF>::precomputeConjugacyStuff()
{
  H.resize(numF);
  UH.resize(numF);
  s.resize(numF);

  for (int i = 0; i<numF; ++i)
  {
    //compute conjugacy matrix
    typename DerivedV::Scalar e1x = dmin(i,0), e1y = dmin(i,1), e2x = dmax(i,0), e2y = dmax(i,1), k1 = kmin[i], k2 = kmax[i];

    H[i]<<
    0,          0, k1*e1x*e1x, k1*e1x*e1y,
    0,          0, k1*e1x*e1y, k1*e1y*e1y,
    k2*e2x*e2x, k2*e2x*e2y,          0,          0,
    k2*e2x*e2y, k2*e2y*e2y,          0,          0;
    Eigen::Matrix<typename DerivedV::Scalar, 4, 4> Ht = H[i].transpose();
    H[i] = .5*(H[i]+Ht);

    Eigen::EigenSolver<Eigen::Matrix<typename DerivedV::Scalar, 4, 4> > es(H[i]);
    s[i] = es.eigenvalues().real();//ok to do this because H symmetric
    //scale
    s[i] = s[i]/(s[i].cwiseAbs().minCoeff());
    UH[i] = es.eigenvectors().real();


  }
}


template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::ConjugateFFSolverData<DerivedV, DerivedF>::computeLaplacians()
{
  computeCoefficientLaplacian(2, DDA);

  computeCoefficientLaplacian(4, DDB);
}

template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::ConjugateFFSolverData<DerivedV, DerivedF>::
precomputeInteriorEdges()
{
  // Flag border edges
  numInteriorEdges = 0;
  isBorderEdge.setZero(numE,1);
  indFullToInterior = -1.*Eigen::VectorXi::Ones(numE,1);

  for(unsigned i=0; i<numE; ++i)
  {
    if ((E2F(i,0) == -1) || ((E2F(i,1) == -1)))
      isBorderEdge[i] = 1;
    else
    {
      indFullToInterior[i] = numInteriorEdges;
      numInteriorEdges++;
    }
  }

  E2F_int.resize(numInteriorEdges, 2);
  indInteriorToFull.setZero(numInteriorEdges,1);
  int ii = 0;
  for (int k=0; k<numE; ++k)
  {
    if (isBorderEdge[k])
      continue;
    E2F_int.row(ii) = E2F.row(k);
    indInteriorToFull[ii] = k;
    ii++;
  }

}



template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::ConjugateFFSolverData<DerivedV, DerivedF>::
computeCoefficientLaplacian(int n, Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &D)
{
  std::vector<Eigen::Triplet<std::complex<typename DerivedV::Scalar> >> tripletList;

  // For every non-border edge
  for (unsigned eid=0; eid<numE; ++eid)
  {
    if (!isBorderEdge[eid])
    {
      int fid0 = E2F(eid,0);
      int fid1 = E2F(eid,1);

      tripletList.push_back(Eigen::Triplet<std::complex<typename DerivedV::Scalar> >(fid0,
                                                                                     fid0,
                                                                                     std::complex<typename DerivedV::Scalar>(1.)));
      tripletList.push_back(Eigen::Triplet<std::complex<typename DerivedV::Scalar> >(fid1,
                                                                                     fid1,
                                                                                     std::complex<typename DerivedV::Scalar>(1.)));
      tripletList.push_back(Eigen::Triplet<std::complex<typename DerivedV::Scalar> >(fid0,
                                                                                     fid1,
                                                                                     -1.*std::polar(1.,-1.*n*K[eid])));
      tripletList.push_back(Eigen::Triplet<std::complex<typename DerivedV::Scalar> >(fid1,
                                                                                     fid0,
                                                                                     -1.*std::polar(1.,1.*n*K[eid])));

    }
  }
  D.resize(numF,numF);
  D.setFromTriplets(tripletList.begin(), tripletList.end());


}

template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::ConjugateFFSolverData<DerivedV, DerivedF>::
computek()
{
  K.setZero(numE);
  // For every non-border edge
  for (unsigned eid=0; eid<numE; ++eid)
  {
    if (!isBorderEdge[eid])
    {
      int fid0 = E2F(eid,0);
      int fid1 = E2F(eid,1);

      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> N0 = FN.row(fid0);
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> N1 = FN.row(fid1);

      // find common edge on triangle 0 and 1
      int fid0_vc = -1;
      int fid1_vc = -1;
      for (unsigned i=0;i<3;++i)
      {
        if (F2E(fid0,i) == eid)
          fid0_vc = i;
        if (F2E(fid1,i) == eid)
          fid1_vc = i;
      }
      assert(fid0_vc != -1);
      assert(fid1_vc != -1);

      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> common_edge = V.row(F(fid0,(fid0_vc+1)%3)) - V.row(F(fid0,fid0_vc));
      common_edge.normalize();

      // Map the two triangles in a new space where the common edge is the x axis and the N0 the z axis
      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> P;
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> o = V.row(F(fid0,fid0_vc));
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> tmp = -N0.cross(common_edge);
      P << common_edge, tmp, N0;
      //      P.transposeInPlace();


      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> V0;
      V0.row(0) = V.row(F(fid0,0)) -o;
      V0.row(1) = V.row(F(fid0,1)) -o;
      V0.row(2) = V.row(F(fid0,2)) -o;

      V0 = (P*V0.transpose()).transpose();

      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> V1;
      V1.row(0) = V.row(F(fid1,0)) -o;
      V1.row(1) = V.row(F(fid1,1)) -o;
      V1.row(2) = V.row(F(fid1,2)) -o;
      V1 = (P*V1.transpose()).transpose();

      // compute rotation R such that R * N1 = N0
      // i.e. map both triangles to the same plane
      double alpha = -atan2(V1((fid1_vc+2)%3,2),V1((fid1_vc+2)%3,1));

      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> R;
      R << 1,          0,            0,
      0, cos(alpha), -sin(alpha) ,
      0, sin(alpha),  cos(alpha);
      V1 = (R*V1.transpose()).transpose();

      // measure the angle between the reference frames
      // k_ij is the angle between the triangle on the left and the one on the right
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> ref0 = V0.row(1) - V0.row(0);
      Eigen::Matrix<typename DerivedV::Scalar, 1, 3> ref1 = V1.row(1) - V1.row(0);

      ref0.normalize();
      ref1.normalize();

      double ktemp = atan2(ref1(1),ref1(0)) - atan2(ref0(1),ref0(0));

      // just to be sure, rotate ref0 using angle ktemp...
      Eigen::Matrix<typename DerivedV::Scalar, 2, 2> R2;
      R2 << cos(ktemp), -sin(ktemp), sin(ktemp), cos(ktemp);

      Eigen::Matrix<typename DerivedV::Scalar, 1, 2> tmp1 = R2*(ref0.head(2)).transpose();

      K[eid] = ktemp;
    }
  }

}


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
evaluateConjugacy(Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> &conjValues)
{
  conjValues.resize(data.numF,1);
  for (int j =0; j<data.numF; ++j)
  {
    Eigen::Matrix<typename DerivedV::Scalar, 4, 1> x; x<<pvU.row(j).transpose(), pvV.row(j).transpose();
    conjValues[j] = x.transpose()*data.H[j]*x;
  }
}


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

  evaluateConjugacy(conjValues);
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

    evaluateConjugacy(conjValues);
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
#endif
