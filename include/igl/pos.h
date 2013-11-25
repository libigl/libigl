// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_POS_H
#define IGL_POS_H

#include <Eigen/Core>

#include <vector>

namespace igl 
{
  // Pos - Fake halfedge for fast and easy navigation on triangle meshes with VT and TT adj
template <typename S>
  class Pos
  {
  public:
    // Init the pos by specifying Face,Edge Index and Orientation
    Pos(const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>* F, 
        Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>* FF, 
        Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>* FFi, 
        int fi,
        int ei,
        bool reverse = false
        )
    : F(F), FF(FF), FFi(FFi), fi(fi), ei(ei), reverse(reverse)
    {}

//    // Init the pos by specifying Face,Vertex Index and Orientation
//    Pos(Eigen::MatrixXi& F, 
//        Eigen::MatrixXi& FF, 
//        Eigen::MatrixXi& FFi, 
//        int fi,
//        int vi,
//        bool reverse = false
//        )
//    : F(F), FF(FF), FFi(FFi), fi(fi), reverse(reverse)
//    {
//      ei = -1;
//      for (int i=0;i<3;++i)
//        if (F(fi,i) == vi)
//          ei = i;
//      assert(ei != -1);
//      
//      if (reverse)
//        ei = (ei-1)%3;
//    }
    
    // Change Face
    void flipF()
    {
      int fin = (*FF)(fi,ei);
      int ein = (*FFi)(fi,ei);
      int reversen = !reverse;
      
      fi = fin;
      ei = ein;
      reverse = reversen;
    }
    
    // Change Edge
    void flipE()
    {
      if (!reverse)
        ei = (ei+2)%3; // ei-1
      else
        ei = (ei+1)%3;

      reverse = !reverse;
    }
    
    // Change Vertex
    void flipV()
    {
      reverse = !reverse;
    }
    
    // Get vertex index
    int Vi()
    {
      assert(fi >= 0);
      assert(fi < F->rows());
      assert(ei >= 0);
      assert(ei <= 2);
      
      if (!reverse)
        return (*F)(fi,ei);
      else
        return (*F)(fi,(ei+1)%3);
    }
    
    // Get face index
    int Fi()
    {
      return fi;
    }
    
    bool operator==(Pos& p2)
    {
      return 
      (
       (fi == p2.fi) &&
       (ei == p2.ei) &&
       (reverse == p2.reverse) &&
       (F   == p2.F) &&
       (FF  == p2.FF) &&
       (FFi == p2.FFi)
       );
    }
    
    int fi;
    int ei;
    bool reverse;
    
    const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>*     F;
    Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>*     FF;
    Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>*     FFi;
  };
  
}

#endif
