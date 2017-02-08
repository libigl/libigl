// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Christian Schüller <schuellchr@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include <igl/cut_mesh_simple.h>

#include <igl/triangle_triangle_adjacency.h>
#include <igl/HalfEdgeIterator.h>

#include <map>

namespace igl
{
  template <typename DerivedS,typename DerivedI>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<typename DerivedI::Scalar>& cut)
  {
    std::vector<std::vector<DerivedI::Scalar>> cutVertices;
    cut_mesh(V.derived(),F.derived(),cut,cutVertices);
  }

  template <typename DerivedS,typename DerivedI>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<typename DerivedI::Scalar>& cut,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices)
  {
    std::vector<int> cutVerticesLink;
    cut_mesh(V.derived(),F.derived(),cut,cutVertices,cutVerticesLink);
  }

  template <typename DerivedS,typename DerivedI>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<typename DerivedI::Scalar>& cut,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices,
                std::vector<int>& cutVerticesLink)
  {
    std::vector<std::vector<DerivedI::Scalar>> cuts;
    cuts.push_back(cut);
    std::vector<std::vector<int>> cutVerticesLinkTemp;
    cut_mesh<DerivedS,DerivedI>(V.derived(),F.derived(),cuts,cutVertices,cutVerticesLinkTemp);
    cutVerticesLink = cutVerticesLinkTemp[0];
  }

  template <typename DerivedS,typename DerivedI>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<typename DerivedI::Scalar>& cut,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices,
                std::vector<int>& cutVerticesLink,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutHalfedges,
                std::vector<int>& cutHalfedgesLink)
  {
    std::vector<std::vector<DerivedI::Scalar>> cuts;
    cuts.push_back(cut);
    std::vector<std::vector<int>> cutVerticesLinkTemp;
    std::vector<std::vector<int>> cutHalfedgesLinkTemp;
    cut_mesh(V,F,cuts,cutVertices,cutVerticesLinkTemp,cutHalfedges,cutHalfedgesLinkTemp);
    cutVerticesLink = cutVerticesLinkTemp[0];
    cutHalfedgesLink = cutHalfedgesLinkTemp[0];
  }

  template <typename DerivedS,typename DerivedI>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                 Eigen::MatrixBase<DerivedI>& F,
                 const std::vector<std::vector<typename DerivedI::Scalar>>& cuts,
                 std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices,
                 std::vector<std::vector<int>>& cutVerticesLink)
  {
    std::vector<std::vector<typename DerivedI::Scalar>> cutHalfedges;
    std::vector<std::vector<int>> cutHalfedgesLink;
    cut_mesh(V,F,cuts,cutVertices,cutVerticesLink,cutHalfedges,cutHalfedgesLink);
  }

  template <typename DerivedS,typename DerivedI>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<std::vector<typename DerivedI::Scalar>>& cuts,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices,
                std::vector<std::vector<int>>& cutVerticesLink,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutHalfedges,
                std::vector<std::vector<int>>& cutHalfedgesLink)
  {
    cutVertices.clear();
    cutVerticesLink.clear();
    cutHalfedges.clear();
    cutHalfedgesLink.clear();

    DerivedI TT,TTi;
    igl::triangle_triangle_adjacency(F.derived(),TT,TTi);
    igl::HalfEdgeIterator<DerivedI> heIter(F.derived(),TT,TTi,0,0);

    DerivedI TTHelp = TT;

    std::map<int,int> cutHalfEdgesLinks;
    std::map<std::pair<int,int>,int> cutHalfEdgesPairLinks;
    std::vector<DerivedI::Scalar> cutVerexIds;
    std::vector<igl::HalfEdgeIterator<DerivedI>::State> cutHalfIter;
    for(auto c : cuts)
    {
      if(c.size() == 0)
      {
        std::cout << "Empty cut vertex sequence provided!" << std::endl;
        continue;
      }

      // find outgoing halfedge from the cut starting vertex
      int f=-1,e=-1;
      for(int t=0;t<F.rows();t++)
      {
        for(int i=0;i<3;i++)
        {
          if(F(t,i) == c[0])
          {
            e = i;
            f = t;
            t = F.rows(); // break outer loop
            break;
          }
        }
      }

      cutVerticesLink.push_back(std::vector<int>());
      cutHalfedgesLink.push_back(std::vector<int>());

      // l-function to remember cut vertices
      auto rememberVertex = [&](){
        auto itv = cutHalfEdgesLinks.insert({heIter.Vi(),-1});
        if(itv.second)
        {
          cutHalfIter.push_back(heIter.getState());
          cutVerexIds.push_back(heIter.Vi());
          itv.first->second = cutHalfIter.size()-1;
        }
        cutVerticesLink.back().push_back(itv.first->second);
      };

      // l-function to remember halfedge pairs 
      auto rememberHalfedge = [&](){
        int vA = heIter.Vi0();
        int vB = heIter.Vi1();
        if(vA > vB) std::swap(vA,vB);
        auto ith = cutHalfEdgesPairLinks.insert({{vA,vB},-1});
        if(ith.second)
        {
          cutHalfedges.push_back(std::vector<DerivedI::Scalar>(4));
          cutHalfedges.back()[0] = heIter.Fi();
          cutHalfedges.back()[1] = heIter.Ei();
          cutHalfedges.back()[2] = heIter.Fif();
          cutHalfedges.back()[3] = heIter.HEi();

          ith.first->second = cutHalfedges.size()-1;
        }
        cutHalfedgesLink.back().push_back(ith.first->second);
      };

      // travel along the cut and detach faces
      heIter.init(f,e);
      if(heIter.Vi() != c[0])
        heIter.flipV();
      for(int i=1;i<c.size();i++)
      {
        int startId = heIter.Vif();
        while(heIter.Vif() != c[i])
        {
          heIter.iterHE();
          assert(startId != heIter.Vif() && "Provided vertex id cut sequence is not valid!");
        }

        // make edge a boundary
        TTHelp(heIter.Fi(),heIter.Ei()) = -1;
        if(heIter.Fif() != -1)
        {
          TTHelp(heIter.Fif(),heIter.HEi()) = -1;
        }

        rememberVertex();
        rememberHalfedge();

        heIter.flipV();
        heIter.flipHE();
      }

      // remember last vertex
      rememberVertex();
    }

    // properly detach all cut vertices
    for(int h=0;h<cutHalfIter.size();h++)
    {
      std::vector<std::vector<std::pair<int,int>>> sharedFIds(1);
      heIter.setState(cutHalfIter[h]);
      do
      {
        sharedFIds[sharedFIds.size()-1].push_back(std::make_pair(heIter.Fi(),heIter.Vii()));

        if(heIter.isBoundaryHE())
        {
          // skip boundary halfedges
          heIter.iterHE();
        }
        if(TTHelp(heIter.Fi(),heIter.Ei()) == -1)
        {
          // add a new vertex for all independent components
          sharedFIds.push_back(std::vector<std::pair<int,int>>());
        }

        heIter.iterHE();

      } while(heIter.Fi() != cutHalfIter[h].fi);

      // add new vertices
      int offset = V.rows();
      V.derived().conservativeResize(offset+sharedFIds.size()-2,3);
      cutVertices.push_back(std::vector<typename DerivedI::Scalar>());
      cutVertices.back().push_back(cutVerexIds[h]);
      for(int i=1;i<sharedFIds.size()-1;i++)
      {
        V.row(offset+i-1) = V.row(cutVerexIds[h]);// +RowVec3S::Random()*0.01;
        cutVertices.back().push_back(offset+i-1);
      }

      // relink new vertices to triangles
      for(int i=1;i<sharedFIds.size()-1;i++)
      {
        for(auto& f : sharedFIds[i])
        {
          F(f.first,f.second) = offset+i-1;
        }
      }
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif
