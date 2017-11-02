// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_SHAPEUP_LOCAL_PROJECTIONS_H
#define IGL_SHAPEUP_LOCAL_PROJECTIONS_H

#include <igl/igl_inline.h>
#include <igl/setdiff.h>
#include <igl/cat.h>
#include <Eigen/Core>
#include <vector>


//This file implements several basic lcaol projection functions for the shapeup algorithm in shapeup.h

namespace igl
{
    
    //This projection does nothing but render points into projP. Mostly used for "echoing" the global step
    bool shapeup_identity_projection(const Eigen::MatrixXd& P, const Eigen::VectorXi& SC, const Eigen::MatrixXi& S,  Eigen::MatrixXd& projP){
        projP.conservativeResize(SC.rows(), 3*SC.maxCoeff());
        for (int i=0;i<S.rows();i++){
            Eigen::RowVector3d avgCurrP=Eigen::RowVector3d::Zero();
            for (int j=0;j<SC(i);j++)
                avgCurrP+=P.row(S(i,j))/(double)(SC(i));
            
            for (int j=0;j<SC(i);j++)
                projP.block(i,3*j,1,3)=P.row(S(i,j))-avgCurrP;
        }
    }
    
    
    //the projection assumes that the sets are vertices of polygons in order
    IGL_INLINE void shapeup_regular_face_projection(const Eigen::MatrixXd& P, const Eigen::VectorXi& SC, const Eigen::MatrixXi& S,  Eigen::MatrixXd& projP)
    {
        projP.conservativeResize(SC.rows(), 3*SC.maxCoeff());
        for (int currRow=0;currRow<SC.rows();currRow++){
            //computing average
            int N=SC(currRow);
            const Eigen::RowVectorXi SRow=S.row(currRow);
            Eigen::RowVector3d avgCurrP=Eigen::RowVector3d::Zero();
            Eigen::MatrixXd targetPolygon(N, 3);
            Eigen::MatrixXd sourcePolygon(N, 3);
            for (int j=0;j<N;j++)
                avgCurrP+=P.row(SRow(j))/(double)(N);
            
            for (int j=0;j<N;j++)
                targetPolygon.row(j)=P.row(SRow(j))-avgCurrP;
            
            //creating perfectly regular source polygon
            for (int j=0;j<N;j++)
                sourcePolygon.row(j)<<cos(2*M_PI*(double)j/(double(N))), sin(2*M_PI*(double)j/(double(N))),0.0;
            
            //finding closest similarity transformation between source and target
            Eigen::MatrixXd corrMat=sourcePolygon.transpose()*targetPolygon;
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(corrMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::MatrixXd R=svd.matrixU()*svd.matrixV().transpose();
            //getting scale by edge length change average. TODO: by singular values
            Eigen::VectorXd sourceEdgeLengths(N);
            Eigen::VectorXd targetEdgeLengths(N);
            for (int j=0;j<N;j++){
                sourceEdgeLengths(j)=(sourcePolygon.row((j+1)%N)-sourcePolygon.row(j)).norm();
                targetEdgeLengths(j)=(targetPolygon.row((j+1)%N)-targetPolygon.row(j)).norm();
            }
            double scale=(targetEdgeLengths.cwiseQuotient(sourceEdgeLengths)).mean();
            
            for (int j=0;j<N;j++)
                projP.block(currRow,3*j,1,3)=sourcePolygon.row(j)*R*scale;
        }
        
        
    }
    
}

#ifndef IGL_STATIC_LIBRARY
#include "shapeup_local_projections.cpp"
#endif

#endif
