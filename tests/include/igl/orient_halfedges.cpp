#include <test_common.h>
#include <igl/orient_halfedges.h>
#include <igl/is_border_vertex.h>
#include <igl/edges.h>
#include <igl/remove_unreferenced.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/EPS.h>

#include <vector>


TEST_CASE("orient_halfedges: sanity checks", "[igl]")
{
  const auto meshes = test_common::manifold_meshes();

  for(const auto& mesh : meshes) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(test_common::data_path(mesh), V, F);
    Eigen::MatrixXi I;
    igl::remove_unreferenced(Eigen::MatrixXd(V), Eigen::MatrixXi(F), V, F,
     I);

    Eigen::MatrixXi TT, TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);
    // Fix mis-match convention
    {
      Eigen::PermutationMatrix<3,3> perm(3);
      perm.indices() = Eigen::Vector3i(1,2,0);
      TT = (TT*perm).eval();
      TTi = (TTi*perm).eval();
      for(int i=0;i<TTi.rows();i++) {
        for(int j=0;j<TTi.cols();j++) {
          TTi(i,j)=TTi(i,j)==-1?-1:(TTi(i,j)+3-1)%3;
        }
      }
    }

    Eigen::MatrixXi E, oE;
    igl::orient_halfedges(F, E, oE);

    const int m = E.maxCoeff()+1;
    std::vector<bool> b(m);
    for(int i=0; i<F.rows(); ++i) {
      for(int j=0; j<3; ++j) {
        b[E(i,j)] = TT(i,j)<0;
      }
    }
    int nb = 0;
    for(int i=0; i<b.size(); ++i) {
      if(b[i]) {
        ++nb;
      }
    }


    //Perform a variety of sanity checks.

    //Correct number of edges
    Eigen::MatrixXi sanityEdges;
    igl::edges(F, sanityEdges);
    REQUIRE(m == sanityEdges.rows());

    //All border halfedges edges have orientation 1 only, all others
    // orientations 1 and -1, so oE must sum to the number of border
    // edges.
    REQUIRE(nb == oE.array().sum());

    //Every border halfedge has orientation 1. Every other edge appeats
    // with orientation 1 and -1.
    std::vector<int> appeared1(m,0), appearedm1(m,0);
    for(int i=0; i<F.rows(); ++i) {
      for(int j=0; j<3; ++j) {
        if(oE(i,j)==1) {
          ++appeared1[E(i,j)];
        } else if(oE(i,j)==-1) {
          ++appearedm1[E(i,j)];
        } else {
          REQUIRE(false); //Only 1 and -1 should occur.
        }
      }
    }
    for(int i=0; i<m; ++i) {
      REQUIRE(appeared1[i]==1);
      if(b[i]) {
        REQUIRE(appearedm1[i]==0);
      } else {
        REQUIRE(appearedm1[i]==1);
      }
    }

    //Two opposite halfedges always map to the same edge
    for(int i=0; i<F.rows(); ++i) {
      for(int j=0; j<3; ++j) {
        if(TT(i,j)>=0) {
          REQUIRE(E(i,j) == E(TT(i,j),TTi(i,j)));
        }
      }
    }
  }
}

