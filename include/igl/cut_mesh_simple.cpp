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

  template <typename DerivedS,typename DerivedI,int Option>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<typename DerivedI::Scalar>& cut,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices,
                std::vector<int> cutVerticesLink,
                Eigen::Matrix<typename DerivedI::Scalar,Eigen::Dynamic,4,Option>& cutHalfedges,
                std::vector<int> cutHalfedgesLink)
  {

  }

  template <typename DerivedS,typename DerivedI>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                 Eigen::MatrixBase<DerivedI>& F,
                 const std::vector<std::vector<typename DerivedI::Scalar>>& cuts,
                 std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices,
                 std::vector<std::vector<int>>& cutVerticesLink)
  {
    cutVertices.clear();
    cutVerticesLink.clear();

    DerivedI TT,TTi;
    igl::triangle_triangle_adjacency(F.derived(),TT,TTi);
    igl::HalfEdgeIterator<DerivedI> heIter(F.derived(),TT,TTi,0,0);

    DerivedI TTHelp = TT;

    std::map<int,int> cutHalfEdgesLinks;
    std::vector<DerivedI::Scalar> cutVerexIds;
    std::vector<igl::HalfEdgeIterator<DerivedI>::State> cutHalfEdges;
    for(auto c : cuts)
    {
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

        // remember vertex of halfedge
        auto it = cutHalfEdgesLinks.insert({heIter.Vi(),-1});
        if(it.second)
        {
          cutHalfEdges.push_back(heIter.getState());
          cutVerexIds.push_back(heIter.Vi());
          it.first->second = cutHalfEdges.size()-1;
        }
        cutVerticesLink.back().push_back(it.first->second);

        heIter.flipV();
        heIter.flipHE();
      }

      // remember last vertex
      auto it = cutHalfEdgesLinks.insert({heIter.Vi(),-1});
      if(it.second)
      {
        cutHalfEdges.push_back(heIter.getState());
        cutVerexIds.push_back(heIter.Vi());
        it.first->second = cutHalfEdges.size()-1;
      }
      cutVerticesLink.back().push_back(it.first->second);
    }

    // properly detach all cut vertices
    for(int h=0;h<cutHalfEdges.size();h++)
    {
      vector<vector<pair<int,int>>> sharedFIds(1);
      heIter.setState(cutHalfEdges[h]);
      do
      {
        sharedFIds[sharedFIds.size()-1].push_back(make_pair(heIter.Fi(),heIter.Vii()));

        if(heIter.isBoundaryHE())
        {
          // skip boundary halfedges
          heIter.iterHE();
        }
        if(TTHelp(heIter.Fi(),heIter.Ei()) == -1)
        {
          // add a new vertex for all independent components
          sharedFIds.push_back(vector<pair<int,int>>());
        }

        heIter.iterHE();

      } while(heIter.Fi() != cutHalfEdges[h].fi);

      // add new vertices
      int offset = V.rows();
      V.derived().conservativeResize(offset+sharedFIds.size()-2,3);
      cutVertices.push_back(vector<typename DerivedI::Scalar>());
      cutVertices.back().push_back(cutVerexIds[h]);
      for(int i=1;i<sharedFIds.size()-1;i++)
      {
        V.row(offset+i-1) = V.row(cutVerexIds[h]);// +RowVec3S::Random()*0.05;
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

  template <typename DerivedS,typename DerivedI,int Option>
  void cut_mesh(Eigen::MatrixBase<DerivedS>& V,
                Eigen::MatrixBase<DerivedI>& F,
                const std::vector<std::vector<typename DerivedI::Scalar>>& cuts,
                std::vector<std::vector<typename DerivedI::Scalar>>& cutVertices,
                std::vector<std::vector<int>> cutVerticesLink,
                std::vector<Eigen::Matrix<typename DerivedI::Scalar,Eigen::Dynamic,4,Option>>& cutHalfedges,
                std::vector<std::vector<int>>& cutHalfedgesLink)
  {

  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
#endif
