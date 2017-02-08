// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "LinSpaced.h"
#include "n_polyvector_general.h"
#include "edge_topology.h"
#include "local_basis.h"
#include "nchoosek.h"
#include "slice.h"
#include "polyroots.h"
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <iostream>
#include <iostream>

namespace igl {
  template <typename DerivedV, typename DerivedF>
  class GeneralPolyVectorFieldFinder
  {
  private:
    const Eigen::PlainObjectBase<DerivedV> &V;
    const Eigen::PlainObjectBase<DerivedF> &F; int numF;
    const int n;

    Eigen::MatrixXi EV; int numE;
    Eigen::MatrixXi F2E;
    Eigen::MatrixXi E2F;
    Eigen::VectorXd K;

    Eigen::VectorXi isBorderEdge;
    int numInteriorEdges;
    Eigen::Matrix<int,Eigen::Dynamic,2> E2F_int;
    Eigen::VectorXi indInteriorToFull;
    Eigen::VectorXi indFullToInterior;

    DerivedV B1, B2, FN;

    IGL_INLINE void computek();
    IGL_INLINE void setFieldFromGeneralCoefficients(const  std::vector<Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic,1>> &coeffs,
                                                    std::vector<Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 2> > &pv);
    IGL_INLINE void computeCoefficientLaplacian(int n, Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &D);
    IGL_INLINE void getGeneralCoeffConstraints(const Eigen::VectorXi &isConstrained,
                                    const Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic> &cfW,
                                    int k,
                                    const Eigen::VectorXi &rootsIndex,
                                    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic,1> &Ck);
    IGL_INLINE void precomputeInteriorEdges();


    IGL_INLINE void minQuadWithKnownMini(const Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &Q,
                                         const Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &f,
                                         const Eigen::VectorXi isConstrained,
                                         const Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> &xknown,
                                         Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic, 1> &x);

  public:
    IGL_INLINE GeneralPolyVectorFieldFinder(const Eigen::PlainObjectBase<DerivedV> &_V,
                                     const Eigen::PlainObjectBase<DerivedF> &_F,
                                     const int &_n);
    IGL_INLINE bool solve(const Eigen::VectorXi &isConstrained,
               const Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic> &cfW,
               const Eigen::VectorXi &rootsIndex,
               Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic> &output);

  };
}

template<typename DerivedV, typename DerivedF>
IGL_INLINE igl::GeneralPolyVectorFieldFinder<DerivedV, DerivedF>::
          GeneralPolyVectorFieldFinder(const Eigen::PlainObjectBase<DerivedV> &_V,
                                const Eigen::PlainObjectBase<DerivedF> &_F,
                                const int &_n):
V(_V),
F(_F),
numF(_F.rows()),
n(_n)
{

  igl::edge_topology(V,F,EV,F2E,E2F);
  numE = EV.rows();


  precomputeInteriorEdges();

  igl::local_basis(V,F,B1,B2,FN);

  computek();

};


template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::GeneralPolyVectorFieldFinder<DerivedV, DerivedF>::
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
IGL_INLINE void igl::GeneralPolyVectorFieldFinder<DerivedV, DerivedF>::
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



template<typename DerivedV, typename DerivedF>
IGL_INLINE bool igl::GeneralPolyVectorFieldFinder<DerivedV, DerivedF>::
                     solve(const Eigen::VectorXi &isConstrained,
                           const Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic> &cfW,
                           const Eigen::VectorXi &rootsIndex,
                           Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic> &output)
{
  // polynomial is of the form:
  // z^(2n) +
  // -c[0]z^(2n-1) +
  // c[1]z^(2n-2) +
  // -c[2]z^(2n-3) +
  // ... +
  // (-1)^n c[n-1]
  std::vector<Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic,1>> coeffs(n,Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic,1>::Zero(numF, 1));
  for (int i =0; i<n; ++i)
  {
    int degree = i+1;
    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic,1> Ck;
    getGeneralCoeffConstraints(isConstrained,
                               cfW,
                               i,
                               rootsIndex,
                               Ck);

    Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > DD;
    computeCoefficientLaplacian(degree, DD);
    Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > f; f.resize(numF,1);

    if (isConstrained.sum() == numF)
      coeffs[i] = Ck;
    else
      minQuadWithKnownMini(DD, f, isConstrained, Ck, coeffs[i]);
  }
  std::vector<Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 2> > pv;
  setFieldFromGeneralCoefficients(coeffs, pv);
  output.setZero(numF,3*n);
  for (int fi=0; fi<numF; ++fi)
  {
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b1 = B1.row(fi);
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b2 = B2.row(fi);
    for (int i=0; i<n; ++i)
      output.block(fi,3*i, 1, 3) = pv[i](fi,0)*b1 + pv[i](fi,1)*b2;
  }
  return true;
}

template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::GeneralPolyVectorFieldFinder<DerivedV, DerivedF>::setFieldFromGeneralCoefficients(const  std::vector<Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic,1>> &coeffs,
                                                            std::vector<Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 2>> &pv)
{
  pv.assign(n, Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 2>::Zero(numF, 2));
  for (int i = 0; i <numF; ++i)
  {

    //    poly coefficients: 1, 0, -Acoeff, 0, Bcoeff
    //    matlab code from roots (given there are no trailing zeros in the polynomial coefficients)
    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic,1> polyCoeff;
    polyCoeff.setZero(n+1,1);
    polyCoeff[0] = 1.;
    int sign = 1;
    for (int k =0; k<n; ++k)
    {
      sign = -sign;
      int degree = k+1;
      polyCoeff[degree] = (1.*sign)*coeffs[k](i);
    }

    Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic,1> roots;
    igl::polyRoots<std::complex<typename DerivedV::Scalar>, typename DerivedV::Scalar >(polyCoeff,roots);
    for (int k=0; k<n; ++k)
    {
      pv[k](i,0) = real(roots[k]);
      pv[k](i,1) = imag(roots[k]);
    }
  }

}


template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::GeneralPolyVectorFieldFinder<DerivedV, DerivedF>::computeCoefficientLaplacian(int n, Eigen::SparseMatrix<std::complex<typename DerivedV::Scalar> > &D)
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

//this gives the coefficients without the (-1)^k that multiplies them
template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::GeneralPolyVectorFieldFinder<DerivedV, DerivedF>::getGeneralCoeffConstraints(const Eigen::VectorXi &isConstrained,
                                                       const Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic> &cfW,
                                                       int k,
                                                       const Eigen::VectorXi &rootsIndex,
                                                       Eigen::Matrix<std::complex<typename DerivedV::Scalar>, Eigen::Dynamic,1> &Ck)
{
  int numConstrained = isConstrained.sum();
  Ck.resize(numConstrained,1);
  // int n = rootsIndex.cols();

  Eigen::MatrixXi allCombs;
  {
    Eigen::VectorXi V = igl::LinSpaced<Eigen::VectorXi >(n,0,n-1);
    igl::nchoosek(V,k+1,allCombs);
  }

  int ind = 0;
  for (int fi = 0; fi <numF; ++fi)
  {
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b1 = B1.row(fi);
    const Eigen::Matrix<typename DerivedV::Scalar, 1, 3> &b2 = B2.row(fi);
    if(isConstrained[fi])
    {
      std::complex<typename DerivedV::Scalar> ck(0);

      for (int j = 0; j < allCombs.rows(); ++j)
      {
        std::complex<typename DerivedV::Scalar> tk(1.);
        //collect products
        for (int i = 0; i < allCombs.cols(); ++i)
        {
          int index = allCombs(j,i);

          int ri = rootsIndex[index];
          Eigen::Matrix<typename DerivedV::Scalar, 1, 3> w;
          if (ri>0)
            w = cfW.block(fi,3*(ri-1),1,3);
          else
            w = -cfW.block(fi,3*(-ri-1),1,3);
          typename DerivedV::Scalar w0 = w.dot(b1);
          typename DerivedV::Scalar w1 = w.dot(b2);
          std::complex<typename DerivedV::Scalar> u(w0,w1);
          tk*= u;
        }
        //collect sum
        ck += tk;
      }
      Ck(ind) = ck;
      ind ++;
    }
  }


}

template<typename DerivedV, typename DerivedF>
IGL_INLINE void igl::GeneralPolyVectorFieldFinder<DerivedV, DerivedF>::computek()
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

//      assert(V0(0,2) < 1e-10);
//      assert(V0(1,2) < 1e-10);
//      assert(V0(2,2) < 1e-10);

      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> V1;
      V1.row(0) = V.row(F(fid1,0)) -o;
      V1.row(1) = V.row(F(fid1,1)) -o;
      V1.row(2) = V.row(F(fid1,2)) -o;
      V1 = (P*V1.transpose()).transpose();

//      assert(V1(fid1_vc,2) < 10e-10);
//      assert(V1((fid1_vc+1)%3,2) < 10e-10);

      // compute rotation R such that R * N1 = N0
      // i.e. map both triangles to the same plane
      double alpha = -atan2(V1((fid1_vc+2)%3,2),V1((fid1_vc+2)%3,1));

      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> R;
      R << 1,          0,            0,
      0, cos(alpha), -sin(alpha) ,
      0, sin(alpha),  cos(alpha);
      V1 = (R*V1.transpose()).transpose();

//      assert(V1(0,2) < 1e-10);
//      assert(V1(1,2) < 1e-10);
//      assert(V1(2,2) < 1e-10);

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

//      assert(tmp1(0) - ref1(0) < 1e-10);
//      assert(tmp1(1) - ref1(1) < 1e-10);

      K[eid] = ktemp;
    }
  }

}


IGL_INLINE void igl::n_polyvector_general(const Eigen::MatrixXd &V,
                             const Eigen::MatrixXi &F,
                             const Eigen::VectorXi& b,
                             const Eigen::MatrixXd& bc,
                             const Eigen::VectorXi &I,
                             Eigen::MatrixXd &output)
{
  Eigen::VectorXi isConstrained = Eigen::VectorXi::Constant(F.rows(),0);
  Eigen::MatrixXd cfW = Eigen::MatrixXd::Constant(F.rows(),bc.cols(),0);

  for(unsigned i=0; i<b.size();++i)
  {
    isConstrained(b(i)) = 1;
    cfW.row(b(i)) << bc.row(i);
  }
  int n = I.rows();
  igl::GeneralPolyVectorFieldFinder<Eigen::MatrixXd, Eigen::MatrixXi> pvff(V,F,n);
  pvff.solve(isConstrained, cfW, I, output);
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
#endif
