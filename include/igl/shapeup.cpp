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
    
    //This function does the precomputation of the necessary operators for the shape-up projection process.
    

    IGL_INLINE void shapeup_precomputation(const Eigen::MatrixXd& V,
                                       const Eigen::VectorXi& D,
                                       const Eigen::MatrixXi& F,
                                       const Eigen::VectorXi& SD,
                                       const Eigen::MatrixXi& S,
                                       const Eigen::VectorXi& h,
                                       const Eigen::VectorXd& w,
                                       const double shapeCoeff,
                                       const double closeCoeff,
                                       struct ShapeupData& sudata)
    {
        using namespace Eigen;
        //The integration solve is separable to x,y,z components
        sudata.V=V; sudata.F=F; sudata.D=D; sudata.SD=SD; sudata.S=S; sudata.h=h; sudata.closeCoeff=closeCoeff; sudata.shapeCoeff=shapeCoeff;
        sudata.Q.conservativeResize(SD.sum(), V.rows());  //Shape matrix (integration);
        sudata.C.conservativeResize(h.rows(), V.rows());        //Closeness matrix for handles
        
        std::vector<Triplet<double> > QTriplets;
        int currRow=0;
        for (int i=0;i<S.rows();i++){
            double avgCoeff=1.0/(double)SD(i);
            
            for (int j=0;j<SD(i);j++){
                for (int k=0;k<SD(i);k++){
                    if (j==k)
                        QTriplets.push_back(Triplet<double>(currRow+j, S(i,k), (1.0-avgCoeff)));
                    else
                        QTriplets.push_back(Triplet<double>(currRow+j, S(i,k), (-avgCoeff)));
                }
            }
            currRow+=SD(i);
        }
        
        sudata.Q.setFromTriplets(QTriplets.begin(), QTriplets.end());
        
        std::vector<Triplet<double> > CTriplets;
        for (int i=0;i<h.size();i++)
            CTriplets.push_back(Triplet<double>(i,h(i), 1.0));
        
        sudata.C.setFromTriplets(CTriplets.begin(), CTriplets.end());
        
        igl::cat(1, sudata.Q, sudata.C, sudata.A);
        sudata.At=sudata.A.transpose();  //to save up this expensive computation.
        
        //weight matrix
        vector<Triplet<double> > WTriplets;
        //std::cout<<"w: "<<w<<std::endl;
        currRow=0;
        for (int i=0;i<SD.rows();i++){
            for (int j=0;j<SD(i);j++)
                WTriplets.push_back(Triplet<double>(currRow+j,currRow+j,shapeCoeff*w(i)));
            currRow+=SD(i);
        }
        for (int i=0;i<h.size();i++)
            WTriplets.push_back(Triplet<double>(SD.sum()+i, SD.sum()+i, closeCoeff));
        sudata.W.resize(SD.sum()+h.size(), SD.sum()+h.size());
        sudata.W.setFromTriplets(WTriplets.begin(), WTriplets.end());
        
        sudata.E=sudata.At*sudata.W*sudata.A;
        sudata.solver.compute(sudata.E);
    }
    
    
    
    IGL_INLINE void shapeup_compute(void (*projection)(int , const hedra::ShapeupData&, const Eigen::MatrixXd& , Eigen::MatrixXd&),
                                    const Eigen::MatrixXd& vh,
                                    const struct ShapeupData& sudata,
                                    Eigen::MatrixXd& currV,
                                    const int maxIterations=50,
                                    const double vTolerance=10e-6)
    {
        using namespace Eigen;
        MatrixXd prevV=currV;
        MatrixXd PV;
        MatrixXd b(sudata.A.rows(),3);
        b.block(sudata.Q.rows(), 0, sudata.h.rows(),3)=vh;  //this stays constant throughout the iterations
        
        //std::cout<<"vh: "<<vh<<std::endl;
        //std::cout<<"V(h(0))"<<currV.row(sudata.h(0))<<std::endl;
        PV.conservativeResize(sudata.SD.rows(), 3*sudata.SD.maxCoeff());
        for (int i=0;i<maxIterations;i++){
            //std::cout<<"A*prevV-b before local projection:"<<(sudata.W*(sudata.A*prevV-b)).squaredNorm()<<std::endl;
            for (int j=0;j<sudata.SD.rows();j++)
                projection(j, sudata, currV, PV);
            //constructing the projection part of the right hand side
            int currRow=0;
            for (int i=0;i<sudata.S.rows();i++){
                for (int j=0;j<sudata.SD(i);j++)
                    b.row(currRow++)=PV.block(i, 3*j, 1,3);
                //currRow+=sudata.SD(i);
            }
            //std::cout<<"A*prevV-b after local projection:"<<(sudata.W*(sudata.A*prevV-b)).squaredNorm()<<std::endl;
            //std::cout<<"A*currV-b:"<<i<<(sudata.A*currV-b)<<std::endl;
            currV=sudata.solver.solve(sudata.At*sudata.W*b);
            //std::cout<<"b: "<<b<<std::endl;
            //std::cout<<"A*cubbV-b after global solve:"<<
            std::cout<<i<<","<<(sudata.W*(sudata.A*currV-b)).squaredNorm()<<std::endl;
            //std::cout<<"b:"<<b.block(b.rows()-1, 0, 1, b.cols())<<std::endl;
            //exit(0);
            double currChange=(currV-prevV).lpNorm<Infinity>();
            //std::cout<<"Iteration: "<<i<<", currChange: "<<currChange<<std::endl;
            prevV=currV;
            if (currChange<vTolerance)
                break;
            
        }
    }
}

#ifndef IGL_STATIC_LIBRARY
#include "shapeup.cpp"
#endif

#endif
