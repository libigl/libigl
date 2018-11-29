#include <test_common.h>

#include <igl/copyleft/cgal/outer_facet.h>

namespace OuterFacetHelper {

/**
 * Check if the outer facet is indeed valid.
 * Assumption: mesh is closed.
 */
template<typename DerivedV, typename DerivedF>
void assert_outer_facet_is_correct(
        const Eigen::PlainObjectBase<DerivedV>& V,
        const Eigen::PlainObjectBase<DerivedF>& F,
        size_t fid, bool flipped) {
    // Todo.
}

TEST_CASE("OuterFacet: Simple", "[igl/copyleft/cgal]")
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    test_common::load_mesh("cube.obj", V, F);

    const size_t num_faces = F.rows();

    Eigen::VectorXi I(num_faces);
    I.setLinSpaced(num_faces, 0, num_faces-1);

    size_t fid = num_faces + 1;
    bool flipped;
    igl::copyleft::cgal::outer_facet(V, F, I, fid, flipped);

    REQUIRE (num_faces > fid);
    REQUIRE (!flipped);
}

TEST_CASE("OuterFacet: DuplicatedOppositeFaces", "[igl/copyleft/cgal]")
    Eigen::MatrixXd V;
    Eigen::MatrixXi F1;
    test_common::load_mesh("cube.obj", V, F1);

    Eigen::MatrixXi F2 = F1;
    F2.col(0).swap(F2.col(1));

    Eigen::MatrixXi F(F1.rows()*2, F1.cols());
    F << F1, F2;

    Eigen::VectorXi I(F.rows());
    I.setLinSpaced(F.rows(), 0, F.rows()-1);

    size_t fid = F.rows() + 1;
    bool flipped;
    igl::copyleft::cgal::outer_facet(V, F, I, fid, flipped);

    REQUIRE (F.rows() > fid);
    REQUIRE (!flipped);
}

TEST_CASE("OuterFacet: FullyDegnerated", "[igl/copyleft/cgal]")
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    test_common::load_mesh("degenerated.obj", V, F);

    Eigen::VectorXi I(F.rows());
    I.setLinSpaced(F.rows(), 0, F.rows()-1);

    size_t fid = F.rows() + 1;
    bool flipped;
    igl::copyleft::cgal::outer_facet(V, F, I, fid, flipped);

    REQUIRE (F.rows() > fid);
    REQUIRE (!flipped);
}

TEST_CASE("OuterFacet: InvertedNormal", "[igl/copyleft/cgal]")
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    test_common::load_mesh("cube.obj", V, F);
    F.col(0).swap(F.col(1));

    Eigen::VectorXi I(F.rows());
    I.setLinSpaced(F.rows(), 0, F.rows()-1);

    size_t fid = F.rows() + 1;
    bool flipped;
    igl::copyleft::cgal::outer_facet(V, F, I, fid, flipped);

    REQUIRE (F.rows() > fid);
    REQUIRE (flipped);
}

TEST_CASE("OuterFacet: SliverTet", "[igl/copyleft/cgal]")
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    test_common::load_mesh("sliver_tet.ply", V, F);

    Eigen::VectorXi I(F.rows());
    I.setLinSpaced(F.rows(), 0, F.rows()-1);

    size_t fid = F.rows() + 1;
    bool flipped;
    igl::copyleft::cgal::outer_facet(V, F, I, fid, flipped);

    REQUIRE (F.rows() > fid);
    REQUIRE (!flipped);
}

}
