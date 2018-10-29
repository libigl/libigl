#include <test_common.h>
#include <igl/cut_mesh.h>
#include <igl/cut_to_disk.h>
#include <igl/edges.h>

#include <array>
#include <iostream>
#include <vector>
#include <set>

namespace cut_to_disk_helper {
  template<typename DerivedV, typename DerivedF>
  void assert_is_disk(
      const Eigen::PlainObjectBase<DerivedV>& V,
      const Eigen::PlainObjectBase<DerivedF>& F,
      const std::vector<std::vector<int>>& cuts) {
    using namespace igl;
    std::set<std::array<int, 2>> cut_edges;
    for (const auto& cut : cuts) {
      const size_t cut_len = cut.size();
      for (size_t i=0; i<cut_len-1; i++) {
        std::array<int, 2> e{cut[i], cut[i+1]};
        if (e[0] > e[1]) {
          std::swap(e[0], e[1]);
        }
        cut_edges.insert(e);
      }
    }

    const size_t num_faces = F.rows();
    Eigen::MatrixXi cut_mask(num_faces, 3);
    cut_mask.setZero();
    for (size_t i=0; i<num_faces; i++) {
      std::array<int, 2> e0{F(i, 0), F(i, 1)};
      std::array<int, 2> e1{F(i, 1), F(i, 2)};
      std::array<int, 2> e2{F(i, 2), F(i, 0)};
      if (e0[0] > e0[1]) std::swap(e0[0], e0[1]);
      if (e1[0] > e1[1]) std::swap(e1[0], e1[1]);
      if (e2[0] > e2[1]) std::swap(e2[0], e2[1]);

      if (cut_edges.find(e0) != cut_edges.end()) {
        cut_mask(i, 0) = 1;
      }
      if (cut_edges.find(e1) != cut_edges.end()) {
        cut_mask(i, 1) = 1;
      }
      if (cut_edges.find(e2) != cut_edges.end()) {
        cut_mask(i, 2) = 1;
      }
    }

    Eigen::MatrixXd V2;
    Eigen::MatrixXi F2;
    igl::cut_mesh(V, F, cut_mask, V2, F2);

    Eigen::MatrixXi E2;
    edges(F2, E2);
    const auto euler = V2.rows() - E2.rows() + F2.rows();
    CHECK ((1 == euler || 2 == euler));
  }
}

TEST_CASE("cut_to_disk: simple_tet", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXi F(4, 3);
  F << 0, 2, 1,
       0, 3, 2,
       1, 2, 3,
       0, 1, 3;
  std::vector<std::vector<int>> cuts;
  cut_to_disk(F, cuts);
  REQUIRE (cuts.size() == 0);
}

TEST_CASE("cut_to_disk: two_disconnected_tet", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXi F(8, 3);
  F << 0, 2, 1,
       0, 3, 2,
       1, 2, 3,
       0, 1, 3,
       4, 6, 5,
       4, 7, 6,
       5, 6, 7,
       4, 5, 7;
  std::vector<std::vector<int>> cuts;
  cut_to_disk(F, cuts);
  REQUIRE (cuts.size() == 0);
}

TEST_CASE("cut_to_disk: simple_square", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXi F(2, 3);
  F << 0, 1, 2,
       2, 1, 3;
  std::vector<std::vector<int>> cuts;
  cut_to_disk(F, cuts);
  REQUIRE (cuts.size() == 1);
  REQUIRE (cuts[0].size() == 5);
  REQUIRE (cuts[0][4] == cuts[0][0]);
}

TEST_CASE("cut_to_disk: torus", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  test_common::load_mesh("TinyTorus.obj", V, F);

  std::vector<std::vector<int>> cuts;
  cut_to_disk(F, cuts);
  REQUIRE (cuts.size() == 2);

  cut_to_disk_helper::assert_is_disk(V, F, cuts);
}

TEST_CASE("cut_to_disk: cube", "[igl]")
{
  using namespace igl;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  test_common::load_mesh("cube.obj", V, F);

  std::vector<std::vector<int>> cuts;
  cut_to_disk(F, cuts);
  REQUIRE (cuts.size() == 0);

  cut_to_disk_helper::assert_is_disk(V, F, cuts);
}
