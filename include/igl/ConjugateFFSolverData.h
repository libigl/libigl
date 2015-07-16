// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Olga Diamanti, 2015 Alec Jacobson
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_CONJUGATE_FF_SOLVER_DATA_H
#define IGL_CONJUGATE_FF_SOLVER_DATA_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

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
    IGL_INLINE void evaluateConjugacy(const Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 2> &pvU,
                                      const Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 2> &pvV,
                                      Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> &conjValues) const ;
  };
}

#include <igl/colon.h>
#include <igl/edge_topology.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/local_basis.h>
#include <igl/principal_curvature.h>
#include <igl/sparse.h>

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
  indFullToInterior = -1*Eigen::VectorXi::Ones(numE,1);

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

template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::ConjugateFFSolverData<DerivedV, DerivedF>::
evaluateConjugacy(const Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 2> &pvU,
                   const Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 2> &pvV,
                  Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 1> &conjValues) const 
{
  conjValues.resize(numF,1);
  for (int j =0; j<numF; ++j)
  {
    Eigen::Matrix<typename DerivedV::Scalar, 4, 1> x; x<<pvU.row(j).transpose(), pvV.row(j).transpose();
    conjValues[j] = x.transpose()*H[j]*x;
  }
}

#endif
