// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_INTEGRABLE_POLYVECTOR_FIELDS
#define IGL_INTEGRABLE_POLYVECTOR_FIELDS
#include "igl_inline.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace igl {
  // Compute a curl-free frame field from user constraints, optionally starting
  // from a gived frame field (assumed to be interpolating the constraints).
  // Implementation of the paper "Integrable PolyVector Fields", SIGGRAPH 2015.

  // Set of parameters used during solve
  struct integrable_polyvector_fields_parameters;

  // All data necessary for solving. Gets initialized from the original field
  //  and gets updated during each solve.
  template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC> class IntegrableFieldSolverData;


  // Precomputes matrices and necessary initial variables from the mesh and the
  // original field. This function is meant to be called before solve().
  // Inputs:
  //   V                 #V by 3 list of mesh vertex coordinates
  //   F                 #F by 3 list of mesh faces (must be triangles)
  //   b                 #B by 1 list of constrained face indices
  //   bc                #B by 6 list of the representative vectors at the constrained faces
  //   constraint_level  #B by 1 list of "constraint level" flag (level=2 indicates that both directions are constrained,
  //                     level = 1 indicates a partially constrained face, i.e. only the first vector will be constrained)
  //   original_field    #F by 6 list of the representative vectors of the frame field to be used as a starting point
  //                     (up to permutation and sign, stacked horizontally for each face)
  // Returns:
  //   data              an IntegrableFieldSolverData object that holds all intermediate
  //                     data needed by the solve routine, with correctly initialized values.
  //
  template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
  IGL_INLINE void integrable_polyvector_fields_precompute(const Eigen::PlainObjectBase<DerivedV>& V,
                                                    const Eigen::PlainObjectBase<DerivedF>& F,
                                                    const Eigen::VectorXi& b,
                                                    const Eigen::PlainObjectBase<DerivedC>& bc,
                                                    const Eigen::VectorXi& constraint_level,
                                                    const Eigen::PlainObjectBase<DerivedFF>& original_field,
                                                    igl::IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC> &data);


  // Given the current estimate of the field, performes one round of optimization
  // iterations and updates the current estimate. The intermediate data is saved
  // and returned for the next iteration.
  // Inputs:
  //   data                         an IntegrableFieldSolverData object that holds all intermediate
  //                                data needed by the solve routine, with their values at the current time instance.
  //   params                       solver parameters (see below)
  //   current_field                #F by 6 list of the representative vectors of the current frame field to be used as a starting point
  //   current_field_is_not_ccw     boolean, determines whether the representative vectors in the current field are in
  //                                non- ccw order - if true, they will be corrected before being passed to the solver.
  //                                The field returned by this routine is always in ccw order, so this flag typically only
  //                                needs to be set to true during the first call to solve(). If unsure, set to true.
  // Returns:
  //   current_field                updated estimate for the integrable field
  //
template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
  IGL_INLINE void integrable_polyvector_fields_solve(IntegrableFieldSolverData<DerivedV, DerivedF, DerivedFF, DerivedC> &cffsoldata,
    integrable_polyvector_fields_parameters &params,
                                               Eigen::PlainObjectBase<DerivedFF>& current_field,
                                               bool current_field_is_not_ccw);


};


//parameters
struct igl::integrable_polyvector_fields_parameters
{
  // number of optimization iterations
  int numIter;
  //weight for barrier term (ensuring ccw ordering of the vectors per face)
  double wBarrier;
  //the s-parameter of the barrier term (see Schüller et al. 2013, Locally Injective Mappings)
  double sBarrier;
  //weight for the PolyCurl term of the energy
  double wCurl;
  //weight for the PolyQuotient term of the energy
  double wQuotCurl;
  //weight for the smoothness term of the energy
  double wSmooth;
  //weight for the closeness to the original field, for the unconstrained faces (regularization term)
  double wCloseUnconstrained;
  //weight for the closeness to the original field, for the constrained faces
  double wCloseConstrained;
  //the per-iteration reduction factor the smoothness weight
  double redFactor_wsmooth;
  //the step size for the Gauss-Newton optimization
  double gamma;
  //tikhonov regularization term (typically not needed, default value should suffice)
  double tikh_gamma;

  IGL_INLINE integrable_polyvector_fields_parameters();

};

//solver data
template <typename DerivedV, typename DerivedF, typename DerivedFF, typename DerivedC>
class igl::IntegrableFieldSolverData
{
public:
  //Original field
  Eigen::VectorXd xOriginal;

  //Constraints
  Eigen::VectorXi constrained;
  Eigen::VectorXi is_constrained_face;//0 for unconstrained, 1 for constrained once (partial) and 2 for fully constrained
  Eigen::VectorXi indInConstrained;
  Eigen::MatrixXd constrained_vec3;

  //Mesh data
  //edges and normalized edge vectors
  Eigen::MatrixXd EVecNorm;
  Eigen::MatrixXi E, E2F, F2E;
  int numV, numF, numE;
  //interior edge data
  int numInteriorEdges;
  Eigen::Matrix<int,Eigen::Dynamic,2> E2F_int;
  Eigen::VectorXi indInteriorToFull;
  Eigen::VectorXi indFullToInterior;
  //per-edge angles (for parallel transport)
  Eigen::VectorXd K;
  //local bases
  Eigen::PlainObjectBase<DerivedV> B1, B2, FN;

  //Solver Data
  Eigen::VectorXd residuals;
  Eigen::SparseMatrix<double> Jac;
  Eigen::VectorXi II_Jac, JJ_Jac;
  Eigen::VectorXd SS_Jac;
  int numVariables;
  int num_residuals;
  int num_residuals_smooth;
  int num_residuals_close;
  int num_residuals_polycurl;
  int num_residuals_quotcurl;
  int num_residuals_barrier;
  const int numInnerJacCols_edge = 8;
  const int numInnerJacCols_face = 4;
  const int numInnerJacRows_smooth = 4;
  const int numInnerJacRows_polycurl = 2;
  const int numInnerJacRows_quotcurl = 1;
  const int numInnerJacRows_barrier = 1;
  const int numInnerJacRows_close = 4;
  int numJacElements_smooth;
  int numJacElements_polycurl;
  int numJacElements_quotcurl;
  int numJacElements_barrier;
  int numJacElements_close;
  int numJacElements;
  IGL_INLINE void add_jac_indices_face(const int numInnerRows,
                                       const int numInnerCols,
                                       const int startRowInJacobian,
                                       const int startIndexInVectors,
                                       Eigen::VectorXi &Rows,
                                       Eigen::VectorXi &Columns);
  IGL_INLINE void face_Jacobian_indices(const int &startRow,
                                        const int &toplace,
                                        const int& fi,
                                        const int& half_degree,
                                        const int &numInnerRows,
                                        const int &numInnerCols,
                                        Eigen::VectorXi &rows,
                                        Eigen::VectorXi &columns);
  IGL_INLINE void add_Jacobian_to_svector(const int &toplace,
                                          const Eigen::MatrixXd &tJac,
                                          Eigen::VectorXd &SS_Jac);

  IGL_INLINE void add_jac_indices_edge(const int numInnerRows,
                                       const int numInnerCols,
                                       const int startRowInJacobian,
                                       const int startIndexInVectors,
                                       Eigen::VectorXi &Rows,
                                       Eigen::VectorXi &Columns);
  IGL_INLINE void edge_Jacobian_indices(const int &startRow,
                                        const int &toplace,
                                        const int& a,
                                        const int& b,
                                        const int& half_degree,
                                        const int &numInnerRows,
                                        const int &numInnerCols,
                                        Eigen::VectorXi &rows,
                                        Eigen::VectorXi &columns);
  std::vector<int> indInSS_Hess_1_vec;
  std::vector<int> indInSS_Hess_2_vec;
  Eigen::SparseMatrix<double> Hess;
  std::vector<Eigen::Triplet<double> > Hess_triplets;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

  IGL_INLINE void precomputeMesh(const Eigen::PlainObjectBase<DerivedV> &_V,
                                 const Eigen::PlainObjectBase<DerivedF> &_F);
  IGL_INLINE void computeInteriorEdges();
  IGL_INLINE void computeJacobianPattern();
  IGL_INLINE void computeHessianPattern();
  IGL_INLINE void computeNewHessValues();
  IGL_INLINE void initializeOriginalVariable(const Eigen::PlainObjectBase<DerivedFF>& original_field);
  IGL_INLINE void initializeConstraints(const Eigen::VectorXi& b,
                                        const Eigen::PlainObjectBase<DerivedC>& bc,
                                        const Eigen::VectorXi& constraint_level);
  IGL_INLINE void makeFieldCCW(Eigen::MatrixXd &sol3D);

public:
  IGL_INLINE IntegrableFieldSolverData();

  IGL_INLINE IntegrableFieldSolverData(
                                     const Eigen::PlainObjectBase<DerivedV> &_V,
                                     const Eigen::PlainObjectBase<DerivedF> &_F,
                                     const Eigen::VectorXi& b,
                                     const Eigen::VectorXi& constraint_level,
                                     const Eigen::PlainObjectBase<DerivedFF>& original_field);

};

#ifndef IGL_STATIC_LIBRARY
#include "integrable_polyvector_fields.cpp"
#endif


#endif /* defined(IGL_INTEGRABLE_POLYVECTOR_FIELDS) */
