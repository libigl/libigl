#include <test_common.h>
#include <catch2/catch.hpp>

#include <igl/readMSH.h>
#include <igl/writeMSH.h>

#include <set>

TEST_CASE("writeMSH","[igl]")
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

    // load source data
    REQUIRE(igl::readMSH(test_common::data_path("sphere_lowres_TMS_1-0001_Magstim_70mm_Fig8_nii_scalar.msh"), 
        X, Tri, Tet, TriTag, TetTag, XFields, XF, EFields, TriF, TetF));

    // save check data
    REQUIRE(igl::writeMSH("test_binary_sphere_2.msh",
        X, Tri, Tet, TriTag, TetTag, XFields, XF, EFields, TriF, TetF));

    // load again
    Eigen::MatrixXd _X;
    Eigen::MatrixXi _Tri;
    Eigen::MatrixXi _Tet;
    Eigen::VectorXi _TriTag;
    Eigen::VectorXi _TetTag;

    std::vector<std::string> _XFields;
    std::vector<std::string> _EFields;

    std::vector<Eigen::MatrixXd> _XF;
    std::vector<Eigen::MatrixXd> _TriF;
    std::vector<Eigen::MatrixXd> _TetF;

    REQUIRE(igl::readMSH("test_binary_sphere_2.msh", 
        _X, _Tri, _Tet, _TriTag, _TetTag, _XFields, _XF, _EFields, _TriF, _TetF));

    REQUIRE(_X.size()   == X.size());
    REQUIRE(_Tri.size() == Tri.size());
    REQUIRE(_Tet.size() == Tet.size());
    REQUIRE(_TriTag.size() == TriTag.size());
    REQUIRE(_TetTag.size() == TetTag.size());
    REQUIRE(_XFields.size() == XFields.size());
    REQUIRE(_XF.size() == XF.size());
    REQUIRE(_EFields.size() == EFields.size());
    REQUIRE(_TriF.size() == TriF.size());
    REQUIRE(_TetF.size() == TetF.size());

    REQUIRE(_X.cols() == 3);
    REQUIRE(_X.rows() == (398+4506));

    REQUIRE(_Tri.cols() == 3);
    REQUIRE(_Tri.rows() == 8988);
    REQUIRE(_TriTag.rows() == 8988);

    REQUIRE(_Tet.cols() == 4);
    REQUIRE(_Tet.rows() == 25937);
    REQUIRE(_TetTag.rows() == 25937);

    //make sure field sizes are correct
    REQUIRE(_XF.size()==0);
    REQUIRE(_TriF.size()==1);
    REQUIRE(_TetF.size()==1);

    // normE , scalar field
    REQUIRE(_TriF[0].cols()==1);
    REQUIRE(_TriF[0].rows()==8988);

    REQUIRE(_TetF[0].cols()==1);
    REQUIRE(_TetF[0].rows()==25937);

    // check the contents too
    for(size_t i=0;i<X.rows();++i)
        for(size_t j=0;j<X.cols();++j)
            REQUIRE(_X(i,j) == X(i,j));

    for(size_t i=0;i<Tri.rows();++i)
        for(size_t j=0;j<Tri.cols();++j)
            REQUIRE(_Tri(i,j) == Tri(i,j));

    for(size_t i=0;i<Tet.rows();++i)
        for(size_t j=0;j<Tet.cols();++j)
            REQUIRE(_Tet(i,j) == Tet(i,j));

    for(size_t i=0;i<XFields.size();++i)
    {
        REQUIRE(XFields[i]==_XFields[i]);
        REQUIRE(XF[i].rows()==_XF[i].rows());
        REQUIRE(XF[i].cols()==_XF[i].cols());
        for(size_t j=0;j<XF[i].rows();++j)
            for(size_t k=0;k<XF[i].cols();++k)
                REQUIRE(XF[i](j,k)==_XF[i](j,k));
    }

    for(size_t i=0;i<EFields.size();++i)
    {
        REQUIRE(EFields[i]==_EFields[i]);
        REQUIRE(TriF[i].rows()==_TriF[i].rows());
        REQUIRE(TriF[i].cols()==_TriF[i].cols());

        for(size_t j=0;j<TriF[i].rows();++j)
            for(size_t k=0;k<TriF[i].cols();++k)
                REQUIRE(TriF[i](j,k)==_TriF[i](j,k));
        
        REQUIRE(TetF[i].rows()==_TetF[i].rows());
        REQUIRE(TetF[i].cols()==_TetF[i].cols());
        for(size_t j=0;j<TetF[i].rows();++j)
            for(size_t k=0;k<TetF[i].cols();++k)
                REQUIRE(TetF[i](j,k)==_TetF[i](j,k));
    }
}
