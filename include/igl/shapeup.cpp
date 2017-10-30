// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SHAPEUP_H
#define IGL_SHAPEUP_H

#include "shapeup.h"
#include "igl/min_quad_with_fixed.h"
#include <igl/igl_inline.h>
#include <igl/setdiff.h>
#include <igl/cat.h>
#include <Eigen/Core>
#include <vector>



namespace igl
{
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
                                           const Eigen::PlainObjectBase<DerivedS>& E,
                                           const Eigen::PlainObjectBase<Derivedb>& b,
                                           const Eigen::PlainObjectBase<Derivedw>& w,
                                           const std::function<bool(const Eigen::PlainObjectBase<DerivedP>&, const Eigen::PlainObjectBase<DerivedSX>&, const Eigen::PlainObjectBase<DerivedS>&,  Eigen::PlainObjectBase<Derivedb&)>& local_projection,
                                           struct ShapeupData& sudata)
    {
        using namespace std;
        using namespace Eigen;
        sudata.P=P;
        sudata.SC=SC;
        sudata.S=S;
        sudata.E=E;
        sudata.b=b;
        sudata.local_projection=local_projection;
        
        sudata.DShape.conservativeResize(SC.sum(), P.rows());  //Shape matrix (integration);
        sudata.DClose.conservativeResize(h.rows(), P.rows());  //Closeness matrix for positional constraints
        sudata.DSmooth.conservativeResize(EV.rows(), P.rows());  //smoothness matrix
        
        //Building shape matrix
        std::vector<Triplet<double> > DShapeTriplets;
        int currRow=0;
        for (int i=0;i<S.rows();i++){
            double avgCoeff=1.0/(double)SC(i);
            
            for (int j=0;j<SC(i);j++){
                for (int k=0;k<SC(i);k++){
                    if (j==k)
                        DShapeTriplets.push_back(Triplet<double>(currRow+j, S(i,k), (1.0-avgCoeff)));
                    else
                        DShapeTriplets.push_back(Triplet<double>(currRow+j, S(i,k), (-avgCoeff)));
                }
            }
            currRow+=SC(i);
        }
        
        sudata.DShape.setFromTriplets(DShapeTriplets.begin(), DShapeTriplets.end());
        
        //Building closeness matrix
        std::vector<Triplet<double> > DCloseTriplets;
        for (int i=0;i<b.size();i++)
            DCloseTriplets.push_back(Triplet<double>(i,h(i), 1.0));
        
        sudata.DClose.setFromTriplets(DCloseTriplets.begin(), DCloseTriplets.end());
        
        igl::cat(1, sudata.DShape, sudata.DClose, sudata.A);
        //is this allowed? repeating A.
        igl::cat(1, sudata.A, sudata.DSmooth, sudata.A);
        //sudata.At=sudata.A.transpose();  //to save up this expensive computation.
        
        //weight matrix
        vector<Triplet<double> > WTriplets;
        
        //one weight per set in S.
        currRow=0;
        for (int i=0;i<SD.rows();i++){
            for (int j=0;j<SD(i);j++)
                WTriplets.push_back(Triplet<double>(currRow+j,currRow+j,shapeCoeff*w(i)));
            currRow+=SD(i);
        }
        
        for (int i=0;i<b.size();i++)
            WTriplets.push_back(Triplet<double>(SD.sum()+i, SD.sum()+i, closeCoeff));
        
        for (int i=0;i<EV.rows();i++)
            WTriplets.push_back(Triplet<double>(SD.sum()+b.size()+i, SD.sum()+b.size()+i, smoothCoeff));
        
        
        sudata.W.conserativeResize(SD.sum()+b.size()+EV.rows(), SD.sum()+b.size()+EV.rows());
        sudata.W.setFromTriplets(WTriplets.begin(), WTriplets.end());
        
        sudata.At=sudata.A.transpose();
        sudata.Q=sudata.At*sudata.W*sudata.A;
        
        return min_quad_with_fixed_precompute(sudata.Q,VectorXi(),SparseMatrix<double>(),true,solver_data);
    }
    
    
    
    template <
    typename Derivedbc,
    typename DerivedP>
    IGL_INLINE void shapeup_solve(const Eigen::PlainObjectBase<Derivedbc>& bc,
                                    const Eigen::PlainObjectBase<DerivedP>& P0,
                                    const struct ShapeupData& sudata,
                                    Eigen::PlainObjectBase<DerivedP>& P)
    {
        using namespace Eigen;
        using namespace std;
        MatrixXd currP;
        MatrixXd prevP=P0;
        MatrixXd projP;
        MatrixXd b(sudata.A.rows(),3);
        b.block(sudata.Q.rows(), 0, sudata.b.rows(),3)=bc;  //this stays constant throughout the iterations
        
        projP.conservativeResize(sudata.SD.rows(), 3*sudata.SC.maxCoeff());
        for (int i=0;i<maxIterations;i++){
            
            for (int j=0;j<sudata.SC.rows();j++)
                sudata.local_projection(currV, SC,S,projP);
            //constructing the projection part of the right hand side
            int currRow=0;
            for (int i=0;i<sudata.S.rows();i++){
                for (int j=0;j<sudata.SC(i);j++){
                    b.row(currRow++)=projP.block(i, 3*j, 1,3);
                }
            }
   
            //the global solve is independent per dimension
            Eigen::PlainObjectBase<DerivedP> rhs=-sudata.At*sudata.W*b;
            min_quad_with_fixed_solve(sudata.solver_data, rhs,Eigen::PlainObjectBase<DerivedP>(),Eigen::PlainObjectBase<DerivedP>(), currP);

            currV=sudata.solver.solve();
            double currChange=(currP-prevP).lpNorm<Infinity>();
            prevP=currP;
            if (currChange<vTolerance)
                break;
            
        }
    }
}

/*#ifdef IGL_STATIC_LIBRARY
template bool igl::shapeup_solve<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, igl::ARAPData&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template bool igl::arap_precomputation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, igl::ARAPData&);
#endif*/
