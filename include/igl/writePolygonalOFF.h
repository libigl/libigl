// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_WRITEPOLYGONALOFF_H
#define IGL_WRITEPOLYGONALOFF_H
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
    
    
    // writes a polygonal mesh into an ascii off file
    // Inputs:
    //   str  path to .off file
    //   V  eigen double matrix #V by 3
    // D  eigen int vector #F by 1 - face degrees
    //   F  eigen int matrix #F by max(D)
    template <typename DerivedV, typename DerivedF>
    IGL_INLINE bool writePolygonalOFF(const std::string str,
                                      const Eigen::PlainObjectBase<DerivedV>& V,
                                      const Eigen::PlainObjectBase<DerivedF>& D,
                                      const Eigen::PlainObjectBase<DerivedF>& F)
    {
        
        using namespace std;
        ofstream FileHandle;
        FileHandle.open(str);
        if (!FileHandle.is_open())
            return false;
        FileHandle<<"OFF"<<endl<<V.rows()<<" "<<F.rows()<<" 0"<<endl;
        FileHandle<<V<<endl;
        MatrixXi FD(D.rows(), D.cols()+F.cols());
        FD<<D, F;
        FileHandle<<FD<<endl;
        FileHandle.close();
        return true;
    }
}


#endif


