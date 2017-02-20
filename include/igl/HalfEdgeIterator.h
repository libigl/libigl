// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
// Heavily edited and extended for halfedge navigation
// by Christian Schüller <schuellchr@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_HALFEDGEITERATOR_H
#define IGL_HALFEDGEITERATOR_H

#include <Eigen/Core>

#include <vector>
#include <igl/igl_inline.h>

// debug macro
#if defined(IGL_HALFEDGE_ITERATOR_DEBUG_0)
  #define HE_ITER_DEBUG(action) std::cout << action << " (" << Fi() << "/" << Vi0() << "/" << Vi1() << "/" << state.boundary << ")" << std::endl;
#elif defined(IGL_HALFEDGE_ITERATOR_DEBUG_1)
  #define HE_ITER_DEBUG(action) std::cout << action << " (" << state.fi << "/" << state.ei << "/" << state.reverse << "/" << state.boundary << ")" << std::endl;
#else
  #define HE_ITER_DEBUG(action)
#endif

namespace igl
{
  // HalfEdgeIterator - Fake halfedge for fast and easy navigation on manifold triangle meshes with triangle_triangle adjacency
  // Example:
  //   Eigen::MatrixXi TT,TTi;
  //   igl::triangle_triangle_adjacency(Triangles,TT,TTi);
  //   igl::HalfEdgeIterator<MatI> heIter(T,TT,TTi,0,0);

  template <typename DerivedF>
  class HalfEdgeIterator
  {
  public:
    struct State
    {
      int fi;
      int ei;
      bool reverse;
      bool boundary;
    };

    // Init the HalfEdgeIterator by specifying Face,Edge Index and Orientation
    IGL_INLINE HalfEdgeIterator(
        const Eigen::PlainObjectBase<DerivedF>& _F,
        const Eigen::PlainObjectBase<DerivedF>& _FF,
        const Eigen::PlainObjectBase<DerivedF>& _FFi,
        int _fi,
        int _ei,
        bool _reverse = false
        )
    : F(_F), FF(_FF), FFi(_FFi)
    {
      state.fi = _fi;
      state.ei = _ei;
      state.reverse = _reverse;
      state.boundary = false;
      HE_ITER_DEBUG("Constructor");
    }

    // Init the HalfEdgeIterator by another
    IGL_INLINE HalfEdgeIterator(
      const Eigen::PlainObjectBase<DerivedF>& _F,
      const Eigen::PlainObjectBase<DerivedF>& _FF,
      const Eigen::PlainObjectBase<DerivedF>& _FFi,
      const HalfEdgeIterator& other)
      : F(_F),FF(_FF),FFi(_FFi)
    {
      state.fi = other.fi;
      state.ei = other.ei;
      state.reverse = other.reverse;
      state.boundary = false;
      HE_ITER_DEBUG("Constructor");
    }

    // Set current face and edge index
    IGL_INLINE bool init(int faceIndex,int edgeIndex,bool reverse = false)
    {
      if(state.fi >= 0 && state.fi < F.rows() && state.ei >= 0 && state.ei <= 2)
      {
        state.fi = faceIndex;
        state.ei = edgeIndex;
        state.reverse = reverse;
        state.boundary = false;

        HE_ITER_DEBUG("Init");
        return true;
      }

      HE_ITER_DEBUG("Init failed");
      return false;
    }

    IGL_INLINE State getState() const
    {
      return state;
    }

    IGL_INLINE void setState(const State& state)
    {
      this->state = state;
      HE_ITER_DEBUG("Set state");
    }

    // Change Vertex
    IGL_INLINE void flipV()
    {
      state.reverse = !state.reverse;
      HE_ITER_DEBUG("Flip Vertex");
    }

    // Change Edge
    IGL_INLINE void flipE()
    {
      state.ei = Eif();
      state.reverse = !state.reverse;
      state.boundary = false;
      HE_ITER_DEBUG("Flip edge");
    }

    // Change to other Halfedge
    // Like flipF() but also works for boundary edges
    IGL_INLINE void flipHE()
    {
      state.boundary = !flipF() && !state.boundary;
      HE_ITER_DEBUG("Flip halfedge");
    }

    // Change Face
    IGL_INLINE bool flipF()
    {
      int fin = Fif();

      // check if not boundary face
      if(fin != -1)
      {
        state.ei = FFi(state.fi,state.ei);
        state.fi = fin;
        state.reverse = !state.reverse;
        HE_ITER_DEBUG("Flip face");
        return true;
      }

      HE_ITER_DEBUG("Flip face failed - boundary");
      return false;
    }

    // Return if vertex is on boundary
    IGL_INLINE bool isBoundaryV() const
    {
      HalfEdgeIterator<DerivedF> iter(F,FF,FFi,0,0);
      HalfEdgeIterator<DerivedF> end(F,FF,FFi,0,0);
      iter.setState(state);
      
      if(state.reverse != state.boundary)
      {
        iter.flipHE();
      }

      bool isBoundary = false;
      end = iter;
      do
      {
        isBoundary = isBoundary || iter.isBoundaryE();
        iter.iterHE();
      } while(iter != end && !isBoundary);

      return isBoundary;
    }

    // Return if edge is on boundary
    IGL_INLINE bool isBoundaryE() const
    {
      return FF(state.fi,state.ei) == -1;
    }

    // Todo
    /*IGL_INLINE int isBoundaryF() const
    {
    }*/

    // Return if halfedge is on boundary
    IGL_INLINE bool isBoundaryHE() const
    {
      return state.boundary;
    }

    // Todo
    /*IGL_INLINE std::vector<int> neighbourV() const
    {
    }*/

    // Todo
    /*IGL_INLINE std::vector<int> bool neighbourF() const
    {
    }*/

    // Move to next halfedge such that Vi0 becomes Vi1
    // Can also be used to travel along boundary halfedges
    IGL_INLINE void nextHE()
    {
      if(state.boundary)
      {
        if(state.reverse)
        {
          flipV();
          do
          {
            flipF();
            flipE();
          } while(!isBoundaryE());
          flipHE();
        }
        else
        {
          do
          {
            flipF();
            flipE();
          } while(!isBoundaryE());
          flipHE();
          flipV();
        }
      }
      else
      {
        state.ei = (state.ei+1)%3;
      }
      HE_ITER_DEBUG("Next halfedge");
    }

    // Returns the next halfedge around the current vertex v, including boundaries
    IGL_INLINE void iterHE()
    {
      if(state.reverse != state.boundary)
      {
        nextHE();
        flipV();
        flipHE();
      }
      else
      {
        flipHE();
        nextHE();
        flipV();
      }
      HE_ITER_DEBUG("Iterate halfedge");
    }

    /*!
     * Returns the next edge around the current vertex v, skipping the boundary
     *      _________
     *     /\ c | b /\
     *    /  \  |  /  \
     *   / d  \ | / a  \
     *  /______\|/______\
     *          v
     * In this example, if a and d are of-boundary and the pos is iterating counterclockwise,
     * this method iterate through the faces incident on vertex v,
     * producing the sequence a, b, c, d, a, b, c, ...
     */
    IGL_INLINE bool nextFE()
    {
      if ( isBoundaryE() ) // we are on a boundary
      {
        do
        {
          flipF();
          flipE();
        } while (!isBoundaryE());
        flipE();
        HE_ITER_DEBUG("Next face edge - border");
        return false;
      }
      else
      {
        flipF();
        flipE();
        HE_ITER_DEBUG("Next face edge");
        return true;
      }
    }

    // Get inner triangle vertex index
    IGL_INLINE int Vii() const
    {
      return !state.reverse ? state.ei : (state.ei+1)%3;
    }

    // Get vertex index
    IGL_INLINE int Vi() const
    {
      assert(state.fi >= 0 && state.fi < F.rows() && state.ei >= 0 && state.ei <= 2);
      return F(state.fi,Vii());
    }

    // Get inner triangle flipped vertex index
    IGL_INLINE int Viif() const
    {
      return !state.reverse ? (state.ei+1)%3 : state.ei;
    }

    // Get flipped vertex index
    IGL_INLINE int Vif() const
    {
      assert(state.fi >= 0 && state.fi < F.rows() && state.ei >= 0 && state.ei <= 2);
      return F(state.fi,Viif());
    }

    // Get inner triangle vertex index at halfedge start
    IGL_INLINE int Vii0() const
    {
      return !state.boundary ? state.ei : (state.ei+1)%3;
    }

    // Get vertex index at halfedge start
    IGL_INLINE int Vi0() const
    {
      assert(state.fi >= 0 && state.fi < F.rows() && state.ei >= 0 && state.ei <= 2);
      return F(state.fi,Vii0());
    }

    // Get inner triangle vertex index at halfedge end
    IGL_INLINE int Vii1() const
    {
      return !state.boundary ? (state.ei+1)%3 : state.ei;
    }

    // Get vertex index at halfedge end
    IGL_INLINE int Vi1() const
    {
      assert(state.fi >= 0 && state.fi < F.rows() && state.ei >= 0 && state.ei <= 2);
      return F(state.fi,Vii1());
    }

    // Get edge index
    IGL_INLINE int Ei() const
    {
      return state.ei;
    }

    // Get flipped edge index
    IGL_INLINE int Eif() const
    {
      return !state.reverse ? (state.ei+2)%3 : (state.ei+1)%3;;
    }

    // Get flipped halfedge index
    IGL_INLINE int HEi() const
    {
      return FFi(state.fi,state.ei);
    }

    // Get face index
    IGL_INLINE int Fi() const
    {
      return state.fi;
    }

    // Get flipped face index
    IGL_INLINE int Fif() const
    {
      return FF(state.fi,state.ei);
    }

    IGL_INLINE HalfEdgeIterator& operator=(const HalfEdgeIterator& p2)
    {
      assert((F == p2.F) && (FF  == p2.FF) && (FFi == p2.FFi));

      state.fi = p2.state.fi;
      state.ei = p2.state.ei;
      state.reverse = p2.state.reverse;
      state.boundary = p2.state.boundary;

      HE_ITER_DEBUG("Assigment");

      return *this;
    }

    IGL_INLINE bool operator==(HalfEdgeIterator& p2) const
    {
      return
      (
       (state.fi == p2.state.fi) &&
       (state.ei == p2.state.ei) &&
       (state.reverse == p2.state.reverse) &&
       (state.boundary == p2.state.boundary) &&
       (F   == p2.F) &&
       (FF  == p2.FF) &&
       (FFi == p2.FFi)
       );
    }

    IGL_INLINE bool operator!=(HalfEdgeIterator& p2) const
    {
      return !(*this == p2);
    }

  private:
    State state;

    const Eigen::PlainObjectBase<DerivedF>& F;
    const Eigen::PlainObjectBase<DerivedF>& FF;
    const Eigen::PlainObjectBase<DerivedF>& FFi;
  };

}

#endif
