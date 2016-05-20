// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_READPOLYGONALOFF_H
#define IGL_READPOLYGONALOFF_H
#include "igl_inline.h"
// History:
//  return type changed from void to bool  Alec 18 Sept 2011

#  include <Eigen/Core>
#include <string>
#include <vector>
#include "list_to_matrix.h"
#include <cstdio>

namespace igl
{
    
    
    // reads mesh from an ascii off file of a polygonal mesh
    // Inputs:
    //   str  path to .off file
    // Outputs:
    //   V  eigen double matrix #V by 3
    // D  eigen int vector #F by 1 - face degrees
    //   F  eigen int matrix #F by max(D)
    template <typename DerivedV, typename DerivedD, typename DerivedF>
    IGL_INLINE bool readPolygonalOFF(const std::string str,
                                     Eigen::PlainObjectBase<DerivedV>& V,
                                     Eigen::PlainObjectBase<DerivedD>& D,
                                     Eigen::PlainObjectBase<DerivedF>& F)
    {
        
        using namespace std;
        ifstream FileHandle;
        FileHandle.open(str);
        if (!FileHandle.is_open())
            return false;
        int NumofVertices, NumofFaces, NumofEdges;
        char OFFString[6];
        vector<vector<int> > RawFaces;
        
        FileHandle>>OFFString>>NumofVertices>>NumofFaces>>NumofEdges;
        V.resize(NumofVertices,3);
        RawFaces.resize(NumofFaces);
        D.resize(NumofFaces,1);
        for (int i=0;i<NumofVertices;i++)
            FileHandle>>V(i,0)>>V(i,1)>>V(i,2);
        
        for (int i=0;i<NumofFaces;i++){
            FileHandle>>D(i);
            RawFaces[i].resize(D(i));
            for (int j=0;j<D(i);j++)
                FileHandle>>RawFaces[i][j];
        }
        
        F.resize(NumofFaces,D.maxCoeff());
        F.setConstant(-1);  //to "don't care" vertices
        for (int i=0;i<NumofFaces;i++)
            for (int j=0;j<RawFaces[i].size();j++)
                F(i,j)=RawFaces[i][j];
        
        //Handling non-zero indexed files: Assuming 0 is not an isolated vertex.
        bool FoundZero=false;
        for (int i=0;i<F.rows();i++)
            for (int j=0;j<F.cols();j++)
                if (F(i,j)==0){
                    FoundZero=true;
                    break;
                }
        //int MinIndex=F.minCoeff();
        //cout<<"MinIndex: "<<MinIndex<<endl;
        if (!FoundZero)
            F-=MatrixXi::Constant(F.rows(), F.cols(), 1);
        FileHandle.close();
        return true;
    }
}


#endif


