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

TEST(copyleft_cgal_peel_outer_hull_layers, TwoCubes) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    test_common::load_mesh("two-boxes-bad-self-union.ply", V, F);
    ASSERT_EQ(486, V.rows());
    ASSERT_EQ(708, F.rows());

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

    ASSERT_EQ(num_faces, I.rows());
    ASSERT_EQ(0, I.minCoeff());
    ASSERT_EQ(1, I.maxCoeff());
}

TEST(PeelOuterHullLayers, CubeWithFold) {
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
