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
    
    // Change Face
    void flipF()
    {
      if (isBorder())
        return;
      
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
    
    bool isBorder()
    {
      return (*FF)(fi,ei) == -1;
    }
    
    /*!
     * Returns the next edge skipping the border
     *      _________
     *     /\ c | b /\
     *    /  \  |  /  \
     *   / d  \ | / a  \
     *  /______\|/______\
     *          v
     * In this example, if a and d are of-border and the pos is iterating counterclockwise, this method iterate through the faces incident on vertex v,
     * producing the sequence a, b, c, d, a, b, c, ...
     */
    bool NextFE()
    {
      if ( isBorder() ) // we are on a border
      {
        do
        {
					flipF();
					flipE();
        } while (!isBorder());
        flipE();
        return false;
      }
      else
      {
        flipF();
        flipE();
        return true;
      }
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

    // Get edge index
    int Ei()
    {
      return ei;
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
    
  private:
    int fi;
    int ei;
    bool reverse;
    
    const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>*     F;
    const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>*     FF;
    const Eigen::Matrix<S, Eigen::Dynamic, Eigen::Dynamic>*     FFi;
  };
  
}

#endif
