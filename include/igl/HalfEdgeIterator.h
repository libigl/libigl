// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_HALFEDGEITERATOR_H
#define IGL_HALFEDGEITERATOR_H

#include <Eigen/Core>

#include <vector>
#include <igl/igl_inline.h>

namespace igl
{
  // HalfEdgeIterator - Fake halfedge for fast and easy navigation on triangle meshes with vertex_triangle_adjacency and
  // triangle_triangle adjacency
  template <typename DerivedF>
  class HalfEdgeIterator
  {
  public:
    // Init the HalfEdgeIterator by specifying Face,Edge Index and Orientation
    IGL_INLINE HalfEdgeIterator(
        const Eigen::PlainObjectBase<DerivedF>& _F,
        const Eigen::PlainObjectBase<DerivedF>& _FF,
        const Eigen::PlainObjectBase<DerivedF>& _FFi,
        int _fi,
        int _ei,
        bool _reverse = false
        )
    : fi(_fi), ei(_ei), reverse(_reverse), F(_F), FF(_FF), FFi(_FFi)
    {}

    // Change Face
    IGL_INLINE void flipF()
    {
      if (isBorder())
        return;

      int fin = (FF)(fi,ei);
      int ein = (FFi)(fi,ei);
      int reversen = !reverse;

      fi = fin;
      ei = ein;
      reverse = reversen;
    }

    // Change Edge
    IGL_INLINE void flipE()
    {
      if (!reverse)
        ei = (ei+2)%3; // ei-1
      else
        ei = (ei+1)%3;

      reverse = !reverse;
    }

    // Change Vertex
    IGL_INLINE void flipV()
    {
      reverse = !reverse;
    }

    IGL_INLINE bool isBorder()
    {
      return (FF)(fi,ei) == -1;
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
    IGL_INLINE bool NextFE()
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
    IGL_INLINE int Vi()
    {
      assert(fi >= 0);
      assert(fi < F.rows());
      assert(ei >= 0);
      assert(ei <= 2);

      if (!reverse)
        return (*F)(fi,ei);
      else
        return (*F)(fi,(ei+1)%3);
    }

    // Get face index
    IGL_INLINE int Fi()
    {
      return fi;
    }

    // Get edge index
    IGL_INLINE int Ei()
    {
      return ei;
    }


    IGL_INLINE bool operator==(HalfEdgeIterator& p2)
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

    const Eigen::PlainObjectBase<DerivedF>& F;
    const Eigen::PlainObjectBase<DerivedF>& FF;
    const Eigen::PlainObjectBase<DerivedF>& FFi;
  };

}

#endif
