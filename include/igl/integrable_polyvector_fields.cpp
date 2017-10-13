// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include <igl/integrable_polyvector_fields.h>
#include <igl/field_local_global_conversions.h>
#include <igl/parallel_transport_angles.h>
#include <igl/local_basis.h>
#include <igl/edge_topology.h>
#include <igl/sparse.h>
#include <igl/sort.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/sort_vectors_ccw.h>
#include <iostream>

IGL_INLINE igl::integrable_polyvector_fields_parameters::integrable_polyvector_fields_parameters():
numIter(5),
wBarrier(0.1),
sBarrier(0.9),
wCurl(10),
wQuotCurl(10),
wSmooth(1.),
wCloseUnconstrained(1e-3),
wCloseConstrained(100),
redFactor_wsmooth(.8),
gamma(0.1),
tikh_gamma(1e-8)
{}




namespace igl {
  template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
  class IntegrableFieldSolver
  {
  private:

    IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC> &data;
    //Symbolic calculations
    IGL_INLINE void rj_barrier_face(const Eigen::RowVectorXd &vec2D_a,
                                    const double &s,
                                    Eigen::VectorXd &residuals,
                                    bool do_jac = false,
                                    // Alec: why use a reference if it can
                                    // point some undefined junk? This is asking
                                    // for trouble...
                                    Eigen::MatrixXd &J = *(Eigen::MatrixXd*)NULL);
    IGL_INLINE void rj_polycurl_edge(const Eigen::RowVectorXd &vec2D_a,
                                     const Eigen::RowVector2d &ea,
                                     const Eigen::RowVectorXd &vec2D_b,
                                     const Eigen::RowVector2d &eb,
                                     Eigen::VectorXd &residuals,
                                     bool do_jac = false,
                                     Eigen::MatrixXd &Jac = *(Eigen::MatrixXd*)NULL);
    IGL_INLINE void rj_quotcurl_edge_polyversion(const Eigen::RowVectorXd &vec2D_a,
                                                 const Eigen::RowVector2d &ea,
                                                 const Eigen::RowVectorXd &vec2D_b,
                                                 const Eigen::RowVector2d &eb,
                                                 Eigen::VectorXd &residuals,
                                                 bool do_jac = false,
                                                 Eigen::MatrixXd &Jac = *(Eigen::MatrixXd*)NULL);
    IGL_INLINE void rj_smoothness_edge(const Eigen::RowVectorXd &vec2D_a,
                                       const Eigen::RowVectorXd &vec2D_b,
                                       const double &k,
                                       const int nA,
                                       const int nB,
                                       Eigen::VectorXd &residuals,
                                       bool do_jac = false,
                                       Eigen::MatrixXd &Jac = *(Eigen::MatrixXd*)NULL);

  public:
    IGL_INLINE IntegrableFieldSolver(IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC> &cffsoldata);

    IGL_INLINE bool solve(integrable_polyvector_fields_parameters &params,
                          Eigen::PlainObjectBase<DerivedFF>& current_field,
                          bool current_field_is_not_ccw);

    IGL_INLINE void solveGaussNewton(integrable_polyvector_fields_parameters &params,
                                     const Eigen::VectorXd &x_initial,
                                     Eigen::VectorXd &x);

    //Compute residuals and Jacobian for Gauss Newton
    IGL_INLINE double RJ(const Eigen::VectorXd &x,
                         const Eigen::VectorXd &x0,
                         const integrable_polyvector_fields_parameters &params,
                         bool doJacs = false);

    IGL_INLINE void RJ_Smoothness(const Eigen::MatrixXd &sol2D,
                                  const double &wSmoothSqrt,
                                  const int startRowInJacobian,
                                  bool doJacs = false,
                                  const int startIndexInVectors = 0);
    IGL_INLINE void RJ_Barrier(const Eigen::MatrixXd &sol2D,
                               const double &s,
                               const double &wBarrierSqrt,
                               const int startRowInJacobian,
                               bool doJacs = false,
                               const int startIndexInVectors = 0);
    IGL_INLINE void RJ_Closeness(const Eigen::MatrixXd &sol2D,
                                 const Eigen::MatrixXd &sol02D,
                                 const double &wCloseUnconstrainedSqrt,
                                 const double &wCloseConstrainedSqrt,
                                 const int startRowInJacobian,
                                 bool doJacs = false,
                                 const int startIndexInVectors = 0);
    IGL_INLINE void RJ_Curl(const Eigen::MatrixXd &sol2D,
                            const double &wCASqrt,
                            const double &wCBSqrt,
                            const int startRowInJacobian,
                            bool doJacs = false,
                            const int startIndexInVectors = 0);
    IGL_INLINE void RJ_QuotCurl(const Eigen::MatrixXd &sol2D,
                                const double &wQuotCurlSqrt,
                                const int startRowInJacobian,
                                bool doJacs = false,
                                const int startIndexInVectors = 0);


  };
};



template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::IntegrableFieldSolverData()
{}


template<typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
precomputeMesh(const Eigen::PlainObjectBase<DerivedV> &_V,
               const Eigen::PlainObjectBase<DerivedF> &_F)
{
  numV = _V.rows();
  numF = _F.rows();
  numVariables = 2*2*numF;
  //Mesh stuff
  igl::edge_topology(_V,_F,E,F2E,E2F);
  numE = E.rows();
  igl::local_basis(_V,_F,B1,B2,FN);
  computeInteriorEdges();
  igl::parallel_transport_angles(_V, _F, FN, E2F, F2E, K);
  EVecNorm.setZero(numE,3);
  for (int k = 0; k<numE; ++k)
    EVecNorm.row(k) = (_V.row(E(k,1))-_V.row(E(k,0))).normalized();
}

template<typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
initializeConstraints(const Eigen::VectorXi& b,
                      const Eigen::PlainObjectBase<DerivedC>& bc,
                      const Eigen::VectorXi& constraint_level)
{
  //save constrained
  Eigen::VectorXi iSorted, constrained_unsorted;
  constrained_unsorted.resize(2*2*b.size());
  is_constrained_face.setZero(numF, 1);
  int ind = 0;
  indInConstrained.setConstant(numF,1,-1);
  for (int k =0; k<b.size(); ++k)
  {
    is_constrained_face[b[k]] = constraint_level[k];
    for (int i=0; i<2;++i)
    {
      int xi = 2*2*b[k] + 2*i +0;
      int yi = 2*2*b[k] + 2*i +1;
      constrained_unsorted[ind++] = xi;
      constrained_unsorted[ind++] = yi;
    }
    indInConstrained[b[k]] = k;
  }
  //sort in descending order (so removing rows will work)
  igl::sort(constrained_unsorted, 1, false, constrained, iSorted);
  constrained_vec3 = bc.template cast<double>();

}

template<typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
makeFieldCCW(Eigen::MatrixXd &sol3D)
{
  //sort ccw
  Eigen::RowVectorXd t;
  Eigen::RowVectorXd all(1,2*2*3);
  Eigen::VectorXi order;
  for (int fi=0; fi<numF; ++fi)
  {
    //take all 4 vectors (including opposites) and pick two that are in ccw order
    all << sol3D.row(fi), -sol3D.row(fi);
	  igl::sort_vectors_ccw(all, FN.row(fi).eval(), order, t);
    //if we are in a constrained face, we need to make sure that the first vector is always the same vector as in the constraints
    if(is_constrained_face[fi])
    {
      const Eigen::RowVector3d &constraint = constrained_vec3.block(indInConstrained[fi], 0, 1, 3);;
      int best_i = -1; double best_score = std::numeric_limits<double>::max();
      for (int i = 0; i<2*2; ++i)
      {
        double score = (t.segment(i*3,3) - constraint).norm();
        if (score<best_score)
        {
          best_score = score;
          best_i = i;
        }
      }
      //do a circshift
      Eigen::RowVectorXd temp = t.segment(0, 3*best_i);
      int s1 = 3*best_i;
      int n2 = 3*best_i;
      int n1 = 3*2*2-n2;
      t.segment(0,n1) = t.segment(s1,n1);
      t.segment(n1, n2) = temp;

    }
    sol3D.row(fi) = t.segment(0, 2*3);
  }
}

template<typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
initializeOriginalVariable(const Eigen::PlainObjectBase<DerivedFF>& original_field)
{
  Eigen::MatrixXd sol2D;
  Eigen::MatrixXd sol3D = original_field.template cast<double>();
  makeFieldCCW(sol3D);
  igl::global2local(B1, B2, sol3D, sol2D);
  xOriginal.setZero(numVariables);
  for (int i =0; i<numF; i++)
    xOriginal.segment(i*2*2, 2*2) = sol2D.row(i);

}

template<typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
computeInteriorEdges()
{
  Eigen::VectorXi isBorderEdge;
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

template<typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
add_Jacobian_to_svector(const int &toplace,
                        const Eigen::MatrixXd &tJac,
                        Eigen::VectorXd &SS_Jac)
{
  int numInnerRows = tJac.rows();
  int numInnerCols = tJac.cols();
  int ind = toplace;
  for (int j=0; j<numInnerRows; ++j)
    for (int k=0; k<numInnerCols; ++k, ++ind)
      SS_Jac(ind) = tJac(j,k);
}

template<typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
add_jac_indices_face(const int numInnerRows,
                     const int numInnerCols,
                     const int startRowInJacobian,
                     const int startIndexInVectors,
                     Eigen::VectorXi &Rows,
                     Eigen::VectorXi &Columns)
{
  for (int fi=0; fi<numF; ++fi)
  {

    int startRow = startRowInJacobian+numInnerRows*fi;
    int startIndex = startIndexInVectors+numInnerRows*numInnerCols*fi;
    face_Jacobian_indices(startRow, startIndex, fi, 2, numInnerRows, numInnerCols, Rows, Columns);
  }
}


template<typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
face_Jacobian_indices(const int &startRow,
                      const int &toplace,
                      const int& fi,
                      const int& half_degree,
                      const int &numInnerRows,
                      const int &numInnerCols,
                      Eigen::VectorXi &rows,
                      Eigen::VectorXi &columns)
{
  int ind = toplace;
  for (int j=0; j<numInnerRows; ++j)
  {
    for (int k=0; k<numInnerCols; ++k, ++ind)
    {
      int iv  = k/2;//which vector (0..half_degree-1) am i at
      int ixy = k%2;//which part (real/imag) am i at
      rows(ind) = startRow+j;
      columns(ind) = 2*half_degree*fi + 2*iv +ixy;
    }
  }
}

template<typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
add_jac_indices_edge(const int numInnerRows,
                     const int numInnerCols,
                     const int startRowInJacobian,
                     const int startIndexInVectors,
                     Eigen::VectorXi &Rows,
                     Eigen::VectorXi &Columns)
{
  for (int ii=0; ii<numInteriorEdges; ++ii)
  {
    // the two faces of the flap
    int a = E2F_int(ii,0);
    int b = E2F_int(ii,1);


    int startRow = startRowInJacobian+numInnerRows*ii;
    int startIndex = startIndexInVectors+numInnerRows*numInnerCols*ii;

    edge_Jacobian_indices(startRow, startIndex, a, b, 2, numInnerRows, numInnerCols, Rows, Columns);

  }
}


template<typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
edge_Jacobian_indices(const int &startRow,
                      const int &toplace,
                      const int& a,
                      const int& b,
                      const int& half_degree,
                      const int &numInnerRows,
                      const int &numInnerCols,
                      Eigen::VectorXi &rows,
                      Eigen::VectorXi &columns)
{
  int ind = toplace;
  for (int j=0; j<numInnerRows; ++j)
  {
    for (int k=0; k<numInnerCols; ++k, ++ind)
    {
      int f   = (k<2*half_degree)?a:b;//which face i am at
      int iv  = (k%(2*half_degree))/2;//which vector (0..half_degree-1) am i at
      int ixy = k%2;//which part (real/imag) am i at
      rows(ind) = startRow+j;
      columns(ind) = 2*half_degree*f + 2*iv +ixy;
    }
  }
}

template<typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
computeJacobianPattern()
{
  num_residuals_smooth = 4*numInteriorEdges;
  num_residuals_close = 4*numF;
  num_residuals_polycurl = 2*numInteriorEdges;
  num_residuals_quotcurl = numInteriorEdges;
  num_residuals_barrier = numF;

  num_residuals = num_residuals_smooth + num_residuals_close + num_residuals_polycurl + num_residuals_quotcurl + num_residuals_barrier;

  residuals.setZero(num_residuals,1);

  //per edge
  numJacElements_smooth = numInteriorEdges*numInnerJacCols_edge*numInnerJacRows_smooth;
  numJacElements_polycurl = numInteriorEdges*numInnerJacCols_edge*numInnerJacRows_polycurl;
  numJacElements_quotcurl = numInteriorEdges*numInnerJacCols_edge*numInnerJacRows_quotcurl;

  //per face
  numJacElements_barrier = numF*numInnerJacCols_face*numInnerJacRows_barrier;
  numJacElements_close = numF*numInnerJacCols_face*numInnerJacRows_close;

  numJacElements = (numJacElements_smooth +numJacElements_polycurl + numJacElements_quotcurl) + (numJacElements_barrier +numJacElements_close);
  //allocate
  II_Jac.setZero(numJacElements);
  JJ_Jac.setZero(numJacElements);
  SS_Jac.setOnes(numJacElements);

  //set stuff (attention: order !)
  int startRowInJacobian = 0;
  int startIndexInVectors = 0;

  //smoothness
  add_jac_indices_edge(numInnerJacRows_smooth,
                       numInnerJacCols_edge,
                       startRowInJacobian,
                       startIndexInVectors,
                       II_Jac,
                       JJ_Jac);
  startRowInJacobian += num_residuals_smooth;
  startIndexInVectors += numJacElements_smooth;

  //closeness
  add_jac_indices_face(numInnerJacRows_close,
                       numInnerJacCols_face,
                       startRowInJacobian,
                       startIndexInVectors,
                       II_Jac,
                       JJ_Jac);
  startRowInJacobian += num_residuals_close;
  startIndexInVectors += numJacElements_close;

  //barrier
  add_jac_indices_face(numInnerJacRows_barrier,
                       numInnerJacCols_face,
                       startRowInJacobian,
                       startIndexInVectors,
                       II_Jac,
                       JJ_Jac);
  startRowInJacobian += num_residuals_barrier;
  startIndexInVectors += numJacElements_barrier;

  //polycurl
  add_jac_indices_edge(numInnerJacRows_polycurl,
                       numInnerJacCols_edge,
                       startRowInJacobian,
                       startIndexInVectors,
                       II_Jac,
                       JJ_Jac);
  startRowInJacobian += num_residuals_polycurl;
  startIndexInVectors += numJacElements_polycurl;

  //quotcurl
  add_jac_indices_edge(numInnerJacRows_quotcurl,
                       numInnerJacCols_edge,
                       startRowInJacobian,
                       startIndexInVectors,
                       II_Jac,
                       JJ_Jac);
  igl::sparse(II_Jac, JJ_Jac, SS_Jac, Jac);
}


template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
computeHessianPattern()
{
  //II_Jac is sorted in ascending order already
  int starti = 0;
  int currI = II_Jac(0);
  for (int ii = 0; ii<II_Jac.rows(); ++ii)
  {
    if(currI != II_Jac(ii))
    {
      starti = ii;
      currI = II_Jac(ii);
    }
    int k1  = II_Jac(ii);
    for (int jj = starti; jj<II_Jac.rows(); ++jj)
    {
      int k2  = II_Jac(jj);
      if (k1 !=k2)
        break;
      indInSS_Hess_1_vec.push_back(ii);
      indInSS_Hess_2_vec.push_back(jj);
      Hess_triplets.push_back(Eigen::Triplet<double> (JJ_Jac(ii),
                                                      JJ_Jac(jj),
                                                      SS_Jac(ii)*SS_Jac(jj)
                                                      )
                              );
    }
  }
  Hess.resize(Jac.cols(),Jac.cols());
  Hess.setFromTriplets(Hess_triplets.begin(), Hess_triplets.end());
  Hess.makeCompressed();
}


template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC>::
computeNewHessValues()
{
  for (int i =0; i<Hess_triplets.size(); ++i)
    Hess_triplets[i] = Eigen::Triplet<double>(Hess_triplets[i].row(),
                                              Hess_triplets[i].col(),
                                              SS_Jac(indInSS_Hess_1_vec[i])*SS_Jac(indInSS_Hess_2_vec[i])
                                              );

  Hess.setFromTriplets(Hess_triplets.begin(), Hess_triplets.end());
}


template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::IntegrableFieldSolver( IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC> &cffsoldata):
data(cffsoldata)
{ };


template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE bool igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
solve(igl::integrable_polyvector_fields_parameters &params,
      Eigen::PlainObjectBase<DerivedFF>& current_field,
      bool current_field_is_not_ccw)
{
  Eigen::MatrixXd sol2D;
  Eigen::MatrixXd sol3D = current_field.template cast<double>();
  if (current_field_is_not_ccw)
    data.makeFieldCCW(sol3D);

  igl::global2local(data.B1, data.B2, sol3D, sol2D);
  Eigen::VectorXd x;
  x.setZero(data.numVariables);
  for (int i =0; i<data.numF; i++)
    x.segment(i*2*2, 2*2) = sol2D.row(i);

  //get x
  solveGaussNewton(params, data.xOriginal, x);
  //get output from x
  for (int i =0; i<data.numF; i++)
    sol2D.row(i) = x.segment(i*2*2, 2*2);
  igl::local2global(data.B1, data.B2, sol2D, sol3D);
  current_field = sol3D.cast<typename DerivedFF::Scalar>();
  return true;
}

template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
solveGaussNewton(integrable_polyvector_fields_parameters &params,
                 const Eigen::VectorXd &x_initial,
                 Eigen::VectorXd &x)
{
  bool converged = false;

  double F;
  Eigen::VectorXd xprev = x;
  Eigen::VectorXd xc = igl::slice(x_initial, data.constrained, 1);
  //  double ESmooth, EClose, ECurl, EQuotCurl, EBarrier;
  for (int innerIter = 0; innerIter<params.numIter; ++innerIter)
  {

    //set constrained entries to those of the initial
    igl::slice_into(xc, data.constrained, 1, xprev);


    //get function, gradients and Hessians
    F = RJ(x, xprev, params, true);

    printf("IntegrableFieldSolver -- Iteration %d\n", innerIter);

    if((data.residuals.array() == std::numeric_limits<double>::infinity()).any())
    {
      std::cerr<<"IntegrableFieldSolver -- residuals: got infinity somewhere"<<std::endl;
      exit(1);
    };
    if((data.residuals.array() != data.residuals.array()).any())
    {
      std::cerr<<"IntegrableFieldSolver -- residuals: got infinity somewhere"<<std::endl;
      exit(1);
    };

    converged = false;

    Eigen::VectorXd rhs = data.Jac.transpose()*data.residuals;

    bool success;
    data.solver.factorize(data.Hess);
    success = data.solver.info() == Eigen::Success;

    if(!success)
      std::cerr<<"IntegrableFieldSolver -- Could not do LU"<<std::endl;

    Eigen::VectorXd direction;

    double error;
    direction = data.solver.solve(rhs);
    error = (data.Hess*direction - rhs).cwiseAbs().maxCoeff();
    if(error> 1e-4)
    {
      std::cerr<<"IntegrableFieldSolver -- Could not solve"<<std::endl;
    }

    // adaptive backtracking
    bool repeat = true;
    int run = 0;
    Eigen::VectorXd cx;
    Eigen::VectorXd tRes;
    double newF;
    while(repeat)
    {
      cx = x - params.gamma*direction;
      newF = RJ(cx, xprev, params);
      if(newF < F)
      {
        repeat = false;
        if(run == 0)
          params.gamma *= 2;
      }
      else
      {
        params.gamma *= 0.5f;
        if(params.gamma<1e-30)
        {
          repeat = false;
          converged = true;
        }
      }
      run++;
    }


    if (!converged)
    {
      xprev = x;
      x = cx;
    }
    else
    {
      std::cerr<<"IntegrableFieldSolver -- Converged"<<std::endl;
      break;
    }
  }


}

template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE double igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
RJ(const Eigen::VectorXd &x,
   const Eigen::VectorXd &x0,
   const integrable_polyvector_fields_parameters &params,
   bool doJacs)
{
  Eigen::MatrixXd sol2D(data.numF,4), sol02D(data.numF,4);
  for (int i =0; i<data.numF; i++)
    sol2D.row(i) = x.segment(i*2*2, 2*2);
  for (int i =0; i<data.numF; i++)
    sol02D.row(i) = x0.segment(i*2*2, 2*2);

  data.SS_Jac.setZero(data.numJacElements);

  //set stuff (attention: order !)
  int startRowInJacobian = 0;
  int startIndexInVectors = 0;

  //smoothness
  RJ_Smoothness(sol2D, sqrt(params.wSmooth), startRowInJacobian, doJacs, startIndexInVectors);
  startRowInJacobian += data.num_residuals_smooth;
  startIndexInVectors += data.numJacElements_smooth;

  //closeness
  RJ_Closeness(sol2D, sol02D, sqrt(params.wCloseUnconstrained), sqrt(params.wCloseConstrained), startRowInJacobian, doJacs, startIndexInVectors);
  startRowInJacobian += data.num_residuals_close;
  startIndexInVectors += data.numJacElements_close;

  //barrier
  RJ_Barrier(sol2D, params.sBarrier, sqrt(params.wBarrier), startRowInJacobian, doJacs, startIndexInVectors);
  startRowInJacobian += data.num_residuals_barrier;
  startIndexInVectors += data.numJacElements_barrier;

  //polycurl
  RJ_Curl(sol2D, params.wCurl, powl(params.wCurl, 2), startRowInJacobian, doJacs, startIndexInVectors);
  startRowInJacobian += data.num_residuals_polycurl;
  startIndexInVectors += data.numJacElements_polycurl;

  //quotcurl
  RJ_QuotCurl(sol2D, sqrt(params.wQuotCurl), startRowInJacobian, doJacs, startIndexInVectors);

  if(doJacs)
  {
    for (int i =0; i<data.numJacElements; ++i)
      data.Jac.coeffRef(data.II_Jac(i), data.JJ_Jac(i)) = data.SS_Jac(i);
    data.computeNewHessValues();
  }

  return data.residuals.transpose()*data.residuals;
}



template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
rj_smoothness_edge(const Eigen::RowVectorXd &vec2D_a,
                   const Eigen::RowVectorXd &vec2D_b,
                   const double &k,
                   const int nA,
                   const int nB,
                   Eigen::VectorXd &residuals,
                   bool do_jac,
                   Eigen::MatrixXd &Jac)
{
  const Eigen::RowVector2d &ua = vec2D_a.segment(0, 2);
  const Eigen::RowVector2d &va = vec2D_a.segment(2, 2);
  const Eigen::RowVector2d &ub = vec2D_b.segment(0, 2);
  const Eigen::RowVector2d &vb = vec2D_b.segment(2, 2);
  const double &xua=ua[0], &yua=ua[1], &xva=va[0], &yva=va[1];
  const double &xub=ub[0], &yub=ub[1], &xvb=vb[0], &yvb=vb[1];

  double xua_2 = xua*xua;
  double xva_2 = xva*xva;
  double yua_2 = yua*yua;
  double yva_2 = yva*yva;
  double xub_2 = xub*xub;
  double xvb_2 = xvb*xvb;
  double yub_2 = yub*yub;
  double yvb_2 = yvb*yvb;

  double sA = sin(nA*k);
  double cA = cos(nA*k);
  double sB = sin(nB*k);
  double cB = cos(nB*k);

  double t1 = xua*yua;
  double t2 = xva*yva;
  double t3 = xub*xvb;
  double t4 = yub*yvb;
  double t5 = xua*xva;
  double t6 = xub*yub;
  double t7 = yua*yva;
  double t8 = xvb*yvb;

  double t9  = xva_2 - yva_2;
  double t10 = xua_2 - yua_2;
  double t11 = xvb_2 - yvb_2;
  double t12 = xub_2 - yub_2;

  double t13 = 2*t1 + 2*t2;

  double t17 = (2*t1*t9 + 2*t2*t10);
  double t19 = (t10*t9 - 4*t5*t7);


  residuals.resize(4, 1);
  residuals <<
  cA*(t10 + t9) - sA*(t13) - t12 - t11,
  sA*(t10 + t9) - 2*t8 - 2*t6 + cA*(t13),
  cB*t19 - (t12)*(t11) - sB*t17 + 4*t3*t4,
  cB*t17 + sB*t19 - 2*t6*(t11) - 2*t8*(t12);


  if (do_jac)
  {
    double t20 = 2*yua*t9 + 4*xua*t2;
    double t21 = 2*xua*t9 - 4*xva*t7;
    double t22 = 2*yva*t10 + 4*t5*yua;
    double t23 = 2*xva*t10 - 4*t1*yva;

    Jac.resize(4,8);
    Jac <<                                                                     2*xua*cA - 2*yua*sA,                                                                     - 2*yua*cA - 2*xua*sA,                                                                     2*xva*cA - 2*yva*sA,                                                                     - 2*yva*cA - 2*xva*sA,                                  -2*xub,                                 2*yub,                                  -2*xvb,                                 2*yvb,
    2*yua*cA + 2*xua*sA,                                                                       2*xua*cA - 2*yua*sA,                                                                     2*yva*cA + 2*xva*sA,                                                                       2*xva*cA - 2*yva*sA,                                  -2*yub,                                -2*xub,                                  -2*yvb,                                -2*xvb,
    cB*(t21) - sB*(t20), - cB*(t20) - sB*(t21), cB*(t23) - sB*(t22), - cB*(t22) - sB*(t23),   4*xvb*t4 - 2*xub*t11, 2*yub*t11 + 4*t3*yvb,   4*xub*t4 - 2*xvb*t12, 2*yvb*t12 + 4*t3*yub,
    cB*(t20) + sB*(t21),   cB*(t21) - sB*(t20), cB*(t22) + sB*(t23),   cB*(t23) - sB*(t22), - 2*yub*t11 - 4*t3*yvb, 4*xvb*t4 - 2*xub*t11, - 2*yvb*t12 - 4*t3*yub, 4*xub*t4 - 2*xvb*t12;
  }
}


template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
RJ_Smoothness(const Eigen::MatrixXd &sol2D,
              const double &wSmoothSqrt,
              const int startRowInJacobian,
              bool doJacs,
              const int startIndexInVectors)
{
  if (wSmoothSqrt ==0)
    return;
  for (int ii=0; ii<data.numInteriorEdges; ++ii)
  {
    // the two faces of the flap
    int a = data.E2F_int(ii,0);
    int b = data.E2F_int(ii,1);

    int k = data.indInteriorToFull[ii];

    Eigen::MatrixXd tJac;
    Eigen::VectorXd tRes;
    rj_smoothness_edge(sol2D.row(a),
                       sol2D.row(b),
                       data.K[k],
                       2*(0+1), //degree of first coefficient
                       2*(1+1), //degree of second coefficient
                       tRes,
                       doJacs,
                       tJac);

    int startRow = startRowInJacobian+data.numInnerJacRows_smooth*ii;
    data.residuals.segment(startRow,data.numInnerJacRows_smooth) = wSmoothSqrt*tRes;

    if(doJacs)
    {
      int startIndex = startIndexInVectors+data.numInnerJacRows_smooth*data.numInnerJacCols_edge*ii;
      data.add_Jacobian_to_svector(startIndex, wSmoothSqrt*tJac,data.SS_Jac);
    }
  }
}


template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
rj_barrier_face(const Eigen::RowVectorXd &vec2D_a,
                const double &s,
                Eigen::VectorXd &residuals,
                bool do_jac,
                Eigen::MatrixXd &Jac)
{

  const Eigen::RowVector2d &ua = vec2D_a.segment(0, 2);
  const Eigen::RowVector2d &va = vec2D_a.segment(2, 2);
  const double &xua=ua[0], &yua=ua[1], &xva=va[0], &yva=va[1];


  double xva_2 = xva*xva;
  double yua_2 = yua*yua;
  double xua_2 = xua*xua;
  double yva_2 = yva*yva;

  double s_2 = s*s;
  double s_3 = s*s_2;

  double t00 = xua*yva;
  double t01 = xva*yua;
  double t05 = t00 - t01;
  double t05_2 = t05*t05;
  double t05_3 = t05*t05_2;

  if (do_jac)
    Jac.setZero(1,4);
  residuals.resize(1, 1);
  if (t05>=s)
    residuals << 0;
  else if (t05<0)
    residuals << std::numeric_limits<double>::infinity();
  else
  {
    residuals << 1/((3*t00 - 3*t01)/s - (3*t05_2)/s_2 + t05_3/s_3) - 1;
    double t03 = (s - t05)*(s - t05);
    double t06 = (3*s_2 - 3*s*t00 + 3*s*t01 + xua_2*yva_2 - 2*xua*t01*yva + xva_2*yua_2);
    double t04 = t06*t06;
    if (do_jac)
      Jac<<
      -(3*s_3*yva*t03)/(t05_2*t04),
      (3*s_3*xva*t03)/(t05_2*t04),
      (3*s_3*yua*t03)/(t05_2*t04),
      -(3*s_3*xua*t03)/(t05_2*t04);
  }
}

template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
RJ_Barrier(const Eigen::MatrixXd &sol2D,
           const double &s,
           const double &wBarrierSqrt,
           const int startRowInJacobian,
           bool doJacs,
           const int startIndexInVectors)
{
  if (wBarrierSqrt ==0)
    return;

  for (int fi=0; fi<data.numF; ++fi)
  {
    Eigen::MatrixXd tJac;
    Eigen::VectorXd tRes;
    rj_barrier_face(sol2D.row(fi),
                    s,
                    tRes,
                    doJacs,
                    tJac);

    int startRow = startRowInJacobian+ data.numInnerJacRows_barrier * fi;
    data.residuals.segment(startRow,data.numInnerJacRows_barrier) = wBarrierSqrt*tRes;

    if(doJacs)
    {
      int startIndex = startIndexInVectors+data.numInnerJacRows_barrier*data.numInnerJacCols_face*fi;
      data.add_Jacobian_to_svector(startIndex, wBarrierSqrt*tJac,data.SS_Jac);
    }
  }
}

template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
RJ_Closeness(const Eigen::MatrixXd &sol2D,
             const Eigen::MatrixXd &sol02D,
             const double &wCloseUnconstrainedSqrt,
             const double &wCloseConstrainedSqrt,
             const int startRowInJacobian,
             bool doJacs,
             const int startIndexInVectors)
{
  if (wCloseUnconstrainedSqrt ==0 && wCloseConstrainedSqrt ==0)
    return;
  for (int fi=0; fi<data.numF; ++fi)
  {
    Eigen::Vector4d weights;
    if (!data.is_constrained_face[fi])
      weights.setConstant(wCloseUnconstrainedSqrt);
    else
    {
      if (data.is_constrained_face[fi]==1)
        //only constrain the first vector
        weights<<wCloseConstrainedSqrt,wCloseConstrainedSqrt,wCloseUnconstrainedSqrt,wCloseUnconstrainedSqrt;
      else
        //either not partial, or partial with 2 constraints
        weights.setConstant(wCloseConstrainedSqrt);
    }

    Eigen::MatrixXd tJac;
    Eigen::VectorXd tRes;
    tJac = Eigen::MatrixXd::Identity(4,4);
    tRes.resize(4, 1); tRes<<(sol2D.row(fi)-sol02D.row(fi)).transpose();
    int startRow = startRowInJacobian+data.numInnerJacRows_close*fi;
    data.residuals.segment(startRow,data.numInnerJacRows_close) = weights.array()*tRes.array();

    if(doJacs)
    {
      int startIndex = startIndexInVectors+data.numInnerJacRows_close*data.numInnerJacCols_face*fi;
      data.add_Jacobian_to_svector(startIndex, weights.asDiagonal()*tJac,data.SS_Jac);
    }

  }
}


template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
rj_polycurl_edge(const Eigen::RowVectorXd &vec2D_a,
                 const Eigen::RowVector2d &ea,
                 const Eigen::RowVectorXd &vec2D_b,
                 const Eigen::RowVector2d &eb,
                 Eigen::VectorXd &residuals,
                 bool do_jac,
                 Eigen::MatrixXd &Jac)
{
  const Eigen::RowVector2d &ua = vec2D_a.segment(0, 2);
  const Eigen::RowVector2d &va = vec2D_a.segment(2, 2);
  const Eigen::RowVector2d &ub = vec2D_b.segment(0, 2);
  const Eigen::RowVector2d &vb = vec2D_b.segment(2, 2);
  const double &xua=ua[0], &yua=ua[1], &xva=va[0], &yva=va[1];
  const double &xub=ub[0], &yub=ub[1], &xvb=vb[0], &yvb=vb[1];
  const double &xea=ea[0], &yea=ea[1];
  const double &xeb=eb[0], &yeb=eb[1];

  const double dua = (xea*xua + yea*yua);
  const double dub = (xeb*xub + yeb*yub);
  const double dva = (xea*xva + yea*yva);
  const double dvb = (xeb*xvb + yeb*yvb);

  const double dua_2 = dua*dua;
  const double dva_2 = dva*dva;
  const double dub_2 = dub*dub;
  const double dvb_2 = dvb*dvb;

  residuals.resize(2, 1);
  residuals << dua_2 - dub_2 + dva_2 - dvb_2,
  dua_2*dva_2 - dub_2*dvb_2 ;



  if (do_jac)
  {

    Jac.resize(2,8);
    Jac << 2*xea*dua,                       2*yea*dua,                       2*xea*dva,                       2*yea*dva,                       -2*xeb*dub,                       -2*yeb*dub,                       -2*xeb*dvb,                       -2*yeb*dvb,
    2*xea*dua*dva_2, 2*yea*dua*dva_2, 2*xea*dua_2*dva, 2*yea*dua_2*dva, -2*xeb*dub*dvb_2, -2*yeb*dub*dvb_2, -2*xeb*dub_2*dvb, -2*yeb*dub_2*dvb;
  }
}

template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
RJ_Curl(const Eigen::MatrixXd &sol2D,
        const double &wCASqrt,
        const double &wCBSqrt,
        const int startRowInJacobian,
        bool doJacs,
        const int startIndexInVectors)
{
  if((wCASqrt==0) &&(wCBSqrt==0))
    return;
  for (int ii=0; ii<data.numInteriorEdges; ++ii)
  {
    // the two faces of the flap
    int a = data.E2F_int(ii,0);
    int b = data.E2F_int(ii,1);

    int k = data.indInteriorToFull[ii];

    // the common edge, a 3 vector
    const Eigen::RowVector3d &ce = data.EVecNorm.row(k);

    // the common edge expressed in local coordinates in the two faces
    // x/y denotes real/imaginary
    double xea = data.B1.row(a).dot(ce);
    double yea = data.B2.row(a).dot(ce);
    Eigen::RowVector2d ea; ea<<xea, yea;
    double xeb = data.B1.row(b).dot(ce);
    double yeb = data.B2.row(b).dot(ce);
    Eigen::RowVector2d eb; eb<<xeb, yeb;


    Eigen::MatrixXd tJac;
    Eigen::VectorXd tRes;
    rj_polycurl_edge(sol2D.row(a),
                     ea,
                     sol2D.row(b),
                     eb,
                     tRes,
                     doJacs,
                     tJac);

    tRes[0] = tRes[0]*wCASqrt;
    tRes[1] = tRes[1]*wCBSqrt;


    int startRow = startRowInJacobian+data.numInnerJacRows_polycurl*ii;
    data.residuals.segment(startRow,data.numInnerJacRows_polycurl) = tRes;

    if(doJacs)
    {
      tJac.row(0) = tJac.row(0)*wCASqrt;
      tJac.row(1) = tJac.row(1)*wCBSqrt;
      int startIndex = startIndexInVectors+data.numInnerJacRows_polycurl*data.numInnerJacCols_edge*ii;
      data.add_Jacobian_to_svector(startIndex, tJac,data.SS_Jac);
    }

  }
}


template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
rj_quotcurl_edge_polyversion(const Eigen::RowVectorXd &vec2D_a,
                             const Eigen::RowVector2d &ea,
                             const Eigen::RowVectorXd &vec2D_b,
                             const Eigen::RowVector2d &eb,
                             Eigen::VectorXd &residuals,
                             bool do_jac,
                             Eigen::MatrixXd &Jac)
{
  const Eigen::RowVector2d &ua = vec2D_a.segment(0, 2);
  const Eigen::RowVector2d &va = vec2D_a.segment(2, 2);
  const Eigen::RowVector2d &ub = vec2D_b.segment(0, 2);
  const Eigen::RowVector2d &vb = vec2D_b.segment(2, 2);
  const double &xua=ua[0], &yua=ua[1], &xva=va[0], &yva=va[1];
  const double &xub=ub[0], &yub=ub[1], &xvb=vb[0], &yvb=vb[1];
  const double &xea=ea[0], &yea=ea[1];
  const double &xeb=eb[0], &yeb=eb[1];

  double dua = (xea*xua + yea*yua);
  double dva = (xea*xva + yea*yva);
  double dub = (xeb*xub + yeb*yub);
  double dvb = (xeb*xvb + yeb*yvb);

  double dua_2 = dua * dua;
  double dva_2 = dva * dva;
  double dub_2 = dub * dub;
  double dvb_2 = dvb * dvb;
  double t00 = (dub_2 - dvb_2);
  double t01 = (dua_2 - dva_2);


  residuals.resize(1, 1);
  residuals << dua*dva*t00 - dub*dvb*t01;

  if (do_jac)
  {
    Jac.resize(1,8);
    Jac <<  xea*dva*t00 - 2*xea*dua*dub*dvb, yea*dva*t00 - 2*yea*dua*dub*dvb, xea*dua*t00 + 2*xea*dub*dva*dvb, yea*dua*t00 + 2*yea*dub*dva*dvb, 2*xeb*dua*dub*dva - xeb*dvb*t01, 2*yeb*dua*dub*dva - yeb*dvb*t01, - xeb*dub*t01 - 2*xeb*dua*dva*dvb, - yeb*dub*t01 - 2*yeb*dua*dva*dvb;
  }
}

template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC>::
RJ_QuotCurl(const Eigen::MatrixXd &sol2D,
            const double &wQuotCurlSqrt,
            const int startRowInJacobian,
            bool doJacs,
            const int startIndexInVectors)
{
  for (int ii=0; ii<data.numInteriorEdges; ++ii)
  {
    // the two faces of the flap
    int a = data.E2F_int(ii,0);
    int b = data.E2F_int(ii,1);

    int k = data.indInteriorToFull[ii];

    // the common edge, a 3 vector
    const Eigen::RowVector3d &ce = data.EVecNorm.row(k);

    // the common edge expressed in local coordinates in the two faces
    // x/y denotes real/imaginary
    double xea = data.B1.row(a).dot(ce);
    double yea = data.B2.row(a).dot(ce);
    Eigen::RowVector2d ea; ea<<xea, yea;
    double xeb = data.B1.row(b).dot(ce);
    double yeb = data.B2.row(b).dot(ce);
    Eigen::RowVector2d eb; eb<<xeb, yeb;


    Eigen::MatrixXd tJac;
    Eigen::VectorXd tRes;
    rj_quotcurl_edge_polyversion(sol2D.row(a),
                                 ea,
                                 sol2D.row(b),
                                 eb,
                                 tRes,
                                 doJacs,
                                 tJac);

    int startRow = startRowInJacobian+ data.numInnerJacRows_quotcurl*ii;
    data.residuals.segment(startRow,data.numInnerJacRows_quotcurl) = wQuotCurlSqrt*tRes;

    if(doJacs)
    {
      int startIndex = startIndexInVectors+data.numInnerJacRows_quotcurl*data.numInnerJacCols_edge*ii;
      data.add_Jacobian_to_svector(startIndex, wQuotCurlSqrt*tJac,data.SS_Jac);
    }
  }
}


template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::integrable_polyvector_fields_precompute(
                                                             const Eigen::PlainObjectBase<DerivedV>& V,
                                                             const Eigen::PlainObjectBase<DerivedF>& F,
                                                             const Eigen::VectorXi& b,
                                                             const Eigen::PlainObjectBase<DerivedC>& bc,
                                                             const Eigen::VectorXi& constraint_level,
                                                             const Eigen::PlainObjectBase<DerivedFF>& original_field,
                                                             igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC> &data)
{
  data.precomputeMesh(V,F);

  data.computeJacobianPattern();
  data.computeHessianPattern();
  data.solver.analyzePattern(data.Hess);

  data.initializeConstraints(b,bc,constraint_level);
  data.initializeOriginalVariable(original_field);
};


template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
IGL_INLINE void igl::integrable_polyvector_fields_solve(igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC> &cffsoldata,
                                                        integrable_polyvector_fields_parameters &params,
                                                        Eigen::PlainObjectBase<DerivedFF>& current_field,
                                                        bool current_field_is_not_ccw)
{
  igl::IntegrableFieldSolver<DerivedV, DerivedF, DerivedFF, DerivedC> cffs(cffsoldata);
  cffs.solve(params, current_field, current_field_is_not_ccw);
};

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template igl::IntegrableFieldSolverData<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >::IntegrableFieldSolverData();
template void igl::integrable_polyvector_fields_solve<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(igl::IntegrableFieldSolverData<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, igl::integrable_polyvector_fields_parameters&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, bool);
template void igl::integrable_polyvector_fields_precompute<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<int, -1, 1, 0, -1, 1> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, igl::IntegrableFieldSolverData<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
#endif
