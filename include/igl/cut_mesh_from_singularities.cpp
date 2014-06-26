// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>, Olga Diamanti <olga.diam@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "cut_mesh_from_singularities.h"

#include <vector>
#include <deque>


namespace igl {
  template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedM,
  typename DerivedS,
  typename DerivedO
  >
  class MeshCutter
  {
  protected:
    const Eigen::PlainObjectBase<DerivedV> &V;
    const Eigen::PlainObjectBase<DerivedF> &F;
    const Eigen::PlainObjectBase<DerivedS> &Handle_Singular;
    const Eigen::PlainObjectBase<DerivedS> &Handle_SingularDegree;
    const Eigen::PlainObjectBase<DerivedM> &Handle_MMatch;

    Eigen::VectorXi F_visited;
    Eigen::PlainObjectBase<DerivedF> TT;
    Eigen::PlainObjectBase<DerivedF> TTi;

  protected:

    inline bool IsRotSeam(const int f0,const int edge)
    {
      unsigned char MM = Handle_MMatch(f0,edge);
      return (MM!=0);
    }

    inline void FloodFill(const int start, Eigen::PlainObjectBase<DerivedO> &Handle_Seams)
    {
      std::deque<int> d;
      ///clean the visited flag
      F_visited(start) = true;
      d.push_back(start);

      while (!d.empty())
      {
        int f = d.at(0); d.pop_front();
        for (int s = 0; s<3; s++)
        {
          int g = TT(f,s); // f->FFp(s);
          int j = TTi(f,s); // f->FFi(s);

          if (j == -1)
          {
            g = f;
            j = s;
          }

          if ((!(IsRotSeam(f,s))) && (!(IsRotSeam(g,j)))  && (!F_visited(g)) )
          {
            Handle_Seams(f,s)=false;
            Handle_Seams(g,j)=false;
            F_visited(g) = true;
            d.push_back(g);
          }
        }
      }
    }

    inline void Retract(Eigen::PlainObjectBase<DerivedO> &Handle_Seams)
    {
      std::vector<int> e(V.rows(),0); // number of edges per vert

      for (unsigned f=0; f<F.rows(); f++)
      {
        for (int s = 0; s<3; s++)
        {
          if (Handle_Seams(f,s))
            if (TT(f,s)<=f)
            {
              e[ F(f,s) ] ++;
              e[ F(f,(s+1)%3) ] ++;
            }
        }
      }

      bool over=true;
      int guard = 0;
      do
      {
        over = true;
        for (int f = 0; f<F.rows(); f++) //if (!f->IsD())
        {
          for (int s = 0; s<3; s++)
          {
            if (Handle_Seams(f,s))
              if (!(IsRotSeam(f,s))) // never retract rot seams
              {
                if (e[ F(f,s) ] == 1) {
                  // dissolve seam
                  Handle_Seams(f,s)=false;
                  if (TT(f,s) != -1)
                    Handle_Seams(TT(f,s),TTi(f,s))=false;

                  e[ F(f,s)] --;
                  e[ F(f,(s+1)%3) ] --;
                  over = false;
                }
              }
          }
        }

        if (guard++>10000)
          over = true;

      } while (!over);
    }

  public:

    inline MeshCutter(const Eigen::PlainObjectBase<DerivedV> &V_,
               const Eigen::PlainObjectBase<DerivedF> &F_,
               const Eigen::PlainObjectBase<DerivedM> &Handle_MMatch_,
               const Eigen::PlainObjectBase<DerivedS> &Handle_Singular_,
               const Eigen::PlainObjectBase<DerivedS> &Handle_SingularDegree_):
    V(V_),
    F(F_),
    Handle_MMatch(Handle_MMatch_),
    Handle_Singular(Handle_Singular_),
    Handle_SingularDegree(Handle_SingularDegree_)
    {
      triangle_triangle_adjacency(V,F,TT,TTi);
    };

    inline void cut(Eigen::PlainObjectBase<DerivedO> &Handle_Seams)
    {
      F_visited.setConstant(F.rows(),0);
      Handle_Seams.setConstant(F.rows(),3,1);

      int index=0;
      for (unsigned f = 0; f<F.rows(); f++)
      {
        if (!F_visited(f))
        {
          index++;
          FloodFill(f, Handle_Seams);
        }
      }

      Retract(Handle_Seams);

      for (unsigned int f=0;f<F.rows();f++)
        for (int j=0;j<3;j++)
          if (IsRotSeam(f,j))
            Handle_Seams(f,j)=true;

    }




  };
}
template <typename DerivedV,
  typename DerivedF,
  typename DerivedM,
  typename DerivedS,
  typename DerivedO>
IGL_INLINE void igl::cut_mesh_from_singularities(const Eigen::PlainObjectBase<DerivedV> &V,
                                                 const Eigen::PlainObjectBase<DerivedF> &F,
                                                 const Eigen::PlainObjectBase<DerivedM> &Handle_MMatch,
                                                 const Eigen::PlainObjectBase<DerivedS> &isSingularity,
                                                 const Eigen::PlainObjectBase<DerivedS> &singularityIndex,
                                                 Eigen::PlainObjectBase<DerivedO> &Handle_Seams)
{
  igl::MeshCutter< DerivedV, DerivedF, DerivedM, DerivedS, DerivedO> mc(V, F, Handle_MMatch, isSingularity, singularityIndex);
  mc.cut(Handle_Seams);

}
#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif
