#include <test_common.h>

#include <catch2/catch.hpp>

#include <igl/readMSH.h>

#include <set>


TEST_CASE("readMSH","[igl]")
{
    Eigen::MatrixXd X;
    Eigen::MatrixXi Tri;
    Eigen::MatrixXi Tet;
    Eigen::VectorXi TriTag;
    Eigen::VectorXi TetTag;

    std::vector<std::string> XFields;
    std::vector<std::string> EFields;

    std::vector<Eigen::MatrixXd> XF;
    std::vector<Eigen::MatrixXd> TriF;
    std::vector<Eigen::MatrixXd> TetF;

    REQUIRE(igl::readMSH(test_common::data_path("sphere_lowres_TMS_1-0001_Magstim_70mm_Fig8_nii_scalar.msh"), 
        X, Tri, Tet, TriTag, TetTag, XFields, XF, EFields, TriF, TetF));

    REQUIRE(X.cols() == 3);
    REQUIRE(X.rows() == (398+4506));

    REQUIRE(Tri.cols() == 3);
    REQUIRE(Tri.rows() == 8988);
    REQUIRE(TriTag.rows() == 8988);

    // determine all tags
    std::set<int> tri_tags_unique;
    for(size_t i=0; i<TriTag.rows(); ++i) tri_tags_unique.insert(TriTag(i));
    REQUIRE(tri_tags_unique.size()==6);

    // make sure we have tags 1001-1006
    for(int i=1;i<6;++i)
        REQUIRE(tri_tags_unique.find(i+1000)!=std::end(tri_tags_unique));

    REQUIRE(Tet.cols() == 4);
    REQUIRE(Tet.rows() == 25937);
    REQUIRE(TetTag.rows() == 25937);
    // determine all tags
    std::set<int> tet_tags_unique;
    for(size_t i=0; i<TetTag.rows(); ++i) tet_tags_unique.insert(TetTag(i));
    REQUIRE(tet_tags_unique.size()==6);

    // make sure we have tags 1-6
    for(int i=1;i<6;++i)
        REQUIRE(tet_tags_unique.find(i)!=std::end(tet_tags_unique));

    REQUIRE(XFields.size()==0);
    REQUIRE(EFields.size()==1);

    REQUIRE(EFields[0]=="normE");

    //make sure field sizes are correct
    REQUIRE(XF.size()==0);
    REQUIRE(TriF.size()==1);
    REQUIRE(TetF.size()==1);

    // normE , scalar field
    REQUIRE(TriF[0].cols()==1);
    REQUIRE(TriF[0].rows()==8988);
    REQUIRE(TetF[0].cols()==1);
    REQUIRE(TetF[0].rows()==25937);
}

