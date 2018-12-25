#include <test_common.h>
#include <iostream>
#include <Eigen/Dense>

#include <igl/copyleft/cgal/peel_outer_hull_layers.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/RemeshSelfIntersectionsParam.h>
#include <igl/per_face_normals.h>
#include <igl/remove_unreferenced.h>
#include <igl/writeOBJ.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

TEST_CASE("copyleft_cgal_peel_outer_hull_layers: TwoCubes", "[igl/copyleft/cgal]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    test_common::load_mesh("two-boxes-bad-self-union.ply", V, F);
    REQUIRE (V.rows() == 486);
    REQUIRE (F.rows() == 708);

    typedef CGAL::Exact_predicates_exact_constructions_kernel K;
    typedef K::FT Scalar;
    typedef Eigen::Matrix<Scalar,
            Eigen::Dynamic,
            Eigen::Dynamic> MatrixXe;

    MatrixXe Vs;
    Eigen::MatrixXi Fs, IF;
    Eigen::VectorXi J, IM;
    igl::copyleft::cgal::RemeshSelfIntersectionsParam param;
    igl::copyleft::cgal::remesh_self_intersections(V, F, param, Vs, Fs, IF, J, IM);

    std::for_each(Fs.data(),Fs.data()+Fs.size(),
            [&IM](int & a){ a=IM(a); });
    MatrixXe Vt;
    Eigen::MatrixXi Ft;
    igl::remove_unreferenced(Vs,Fs,Vt,Ft,IM);
    const size_t num_faces = Ft.rows();

    Eigen::VectorXi I, flipped;
    size_t num_peels = igl::copyleft::cgal::peel_outer_hull_layers(Vt, Ft, I, flipped);

    Eigen::MatrixXd vertices(Vt.rows(), Vt.cols());
    std::transform(Vt.data(), Vt.data() + Vt.rows() * Vt.cols(),
            vertices.data(), [](Scalar v) { return CGAL::to_double(v); });
    igl::writeOBJ("debug.obj", vertices, Ft);

    REQUIRE (I.rows() == num_faces);
    REQUIRE (I.minCoeff() == 0);
    REQUIRE (I.maxCoeff() == 1);
}

TEST_CASE("PeelOuterHullLayers: CubeWithFold", "[igl/copyleft/cgal]")
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> V;
    Eigen::MatrixXi F;
    test_common::load_mesh("cube_with_fold.ply", V, F);

    typedef CGAL::Exact_predicates_exact_constructions_kernel K;
    typedef K::FT Scalar;
    typedef Eigen::Matrix<Scalar,
            Eigen::Dynamic,
            Eigen::Dynamic> MatrixXe;

    MatrixXe Vs;
    Eigen::MatrixXi Fs, IF;
    Eigen::VectorXi J, IM;
    igl::copyleft::cgal::RemeshSelfIntersectionsParam param;
    igl::copyleft::cgal::remesh_self_intersections(V, F, param, Vs, Fs, IF, J, IM);

    std::for_each(Fs.data(),Fs.data()+Fs.size(),
            [&IM](int & a){ a=IM(a); });
    MatrixXe Vt;
    Eigen::MatrixXi Ft;
    igl::remove_unreferenced(Vs,Fs,Vt,Ft,IM);

    Eigen::VectorXi I, flipped;
    size_t num_peels = igl::copyleft::cgal::peel_outer_hull_layers(Vt, Ft, I, flipped);
}
