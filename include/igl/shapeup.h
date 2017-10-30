// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SHAPEUP_H
#define IGL_SHAPEUP_H

#include <igl/igl_inline.h>
#include <igl/setdiff.h>
#include <igl/cat.h>
#include <Eigen/Core>
#include <vector>


//This file implements the following algorithm:

//Boaziz et al.
//Shape-Up: Shaping Discrete Geometry with Projections
//Computer Graphics Forum (Proc. SGP) 31(5), 2012

namespace igl
{
    
    struct ShapeupData{

        //input data
        Eigen::MatrixXd P;
        Eigen::VectorXi SD;
        Eigen::MatrixXi S;
        Eigen::VectorXi b;
        double shapeCoeff, closeCoeff;
        std::function<bool(const MatrixXd&, const VectorXi&, const MatrixXi&,  Eigen::MatrixXd&)> local_projection;
        
        //Internally-used matrices
        Eigen::SparseMatrix<double> A, Q, C, E, At, W;
        
         min_quad_with_fixed_data<double> solver_data;
    };
    
    
    
    //This function precomputation the necessary matrices for the ShapeUp process, and prefactorizes them.
    
    //input:
    //  P   #P by 3             point positions
    //  SC  #Set by 1           cardinalities of sets in S
    //  S   #Sets by max(SC)    independent sets where the local projection applies. Values beyond column SC(i) in row S(i,:) are "don't care"
    //  b   #b by 1             boundary (fixed) vertices from P.
    //  w   #Set by 1           weight for each set (used in the global step)
    //  local_projection              function pointer taking (P,SC,S,projP),
    // where the first three parameters are as defined, and "projP" is the output, as a #S by 3*max(SC) function in format xyzxyzxyz, and where it returns the projected points corresponding to each set in S in the same order.
    //  shapecoeff,
    //  closeCoeff,
    //  smoothCoeff             energy coefficients as mentioned in the paper
    
    // Output:
    //  sudata struct ShapeupData     the data necessary to solve the system in shapeup_solve

    template <
    typename DerivedP,
    typename DerivedSC,
    typename DerivedS,
    typename Derivedb,
    typename Derivedw>
    IGL_INLINE void shapeup_precomputation(
                                           const Eigen::PlainObjectBase<DerivedP>& P,
                                           const Eigen::PlainObjectBase<DerivedSC>& SC,
                                           const Eigen::PlainObjectBase<DerivedS>& S,
                                           const Eigen::PlainObjectBase<Derivedb>& b,
                                           const Eigen::PlainObjectBase<Derivedw>& w,
                                           const std::function<bool(const Eigen::PlainObjectBase<DerivedP>&, const Eigen::PlainObjectBase<DerivedSX>&, const Eigen::PlainObjectBase<DerivedS>&,  Eigen::PlainObjectBase<Derivedb&)>& local_projection,
                                           const double shapeCoeff,
                                           const double closeCoeff,
                                           const double smoothCoeff,
                                           struct ShapeupData& sudata);
    
    //This function solve the shapeup project optimization. shapeup_precompute must be called before with the same sudata, or results are unpredictable
    
    //Input:
    //bc                #b by 3 fixed point values corresonding to "b" in sudata
    //NOTE: the input values in P0 don't need to correspond to prescribed values in bc; the iterations will project them automatically (by design).
    //P0                #P by 3 initial solution (point positions)
    //maxIterations     referring to number of local-global pairs.
    //pTolerance        algorithm stops when max(|P_k-P_{k-1}|)<pTolerance.
    //Output:
    //P                 the solution to the problem, corresponding to P0.
    template <
    typename Derivedbc,
    typename DerivedP>
    IGL_INLINE void shapeup_compute(const Eigen::PlainObjectBase<Derivedbc>& bc,
                                   const Eigen::PlainObjectBase<DerivedP>& P0,
                                    const struct ShapeupData& sudata,
                                    const int maxIterations=50,
                                    const double pTolerance=10e-6
                                    const Eigen::PlainObjectBase<DerivedP>& P)
    {
      
    }
}

#ifndef IGL_STATIC_LIBRARY
#include "shapeup.cpp"
#endif

#endif
