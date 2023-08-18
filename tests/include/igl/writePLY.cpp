#include <test_common.h>
#include <igl/readPLY.h>
#include <igl/writePLY.h>
#include <igl/edges.h>
#include <fstream>
#include <string>
#include <vector>

TEST_CASE("writePLY: bunny.ply", "[igl]")
{
    std::ifstream f(test_common::data_path("bunny.ply"));
    REQUIRE (f.good());
    f.close();

    Eigen::MatrixXd V1,N1,UV1,VD1,FD1,ED1;
    std::vector<std::string> Vheader1,Fheader1,Eheader1,comments1;

    Eigen::MatrixXi F1,E1;

    // load test data first
    REQUIRE (igl::readPLY(test_common::data_path("bunny.ply"), V1, F1, E1, N1, UV1, VD1,Vheader1, FD1,Fheader1, ED1,Eheader1,comments1));

    // add more data
    Vheader1.push_back("dummy_data");
    Eigen::VectorXd dummy_data(V1.rows());
    for(size_t i=0;i<V1.rows();++i)
        dummy_data(i)=(double)i;
    Eigen::MatrixXd VD2(V1.rows(),VD1.cols()+1);
    VD2<<VD1,dummy_data;

    Fheader1.push_back("face_data");
    Eigen::VectorXd face_data(F1.rows());
    for(size_t i=0;i<F1.rows();++i)
        face_data(i)=(double)i;



    // there is no face data in the input file
    Eigen::MatrixXd FD2(F1.rows(),FD1.cols()+1);
    FD2<<face_data;

    //input file have no edge data
    REQUIRE (E1.rows() == 0);
    REQUIRE (E1.cols() == 0);
    REQUIRE (ED1.rows() == 0);
    REQUIRE (Eheader1.empty());    

    Eigen::MatrixXi E2;
    std::vector<std::string> Eheader2;

    //generate edges
    igl::edges(F1,E2);

    //generate edge data
    Eheader2.push_back("edge_data");
    Eigen::VectorXd edge_data(E2.rows());
    for(size_t i=0;i<E2.rows();++i)
        edge_data(i)=(double)i;
    
    // there is no edge data in the input file
    Eigen::MatrixXd ED2(E2.rows(),1);
    ED2<<edge_data;

    // test that saving preserves all the data, including new data column
    REQUIRE (igl::writePLY("writePLY_test_bunny.ply", V1, F1, E2, N1, UV1, VD2, Vheader1, FD2,Fheader1, ED2, Eheader2, comments1, igl::FileEncoding::Binary));

    Eigen::MatrixXd V,N,UV,VD,FD,ED;
    Eigen::MatrixXi F,E;
    std::vector<std::string> Vheader,Fheader,Eheader,comments;

    // test that saving preserves all the data
    REQUIRE (igl::readPLY("writePLY_test_bunny.ply", V, F, E, N, UV, VD,Vheader, FD,Fheader, ED,Eheader, comments));

    REQUIRE (V.rows() == 35947);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 69451);
    REQUIRE (F.cols() == 3);

    // generated edges
    REQUIRE (E.rows() == E2.rows());
    REQUIRE (E.cols() == E2.cols());

    // no normals or texture coordinates
    REQUIRE (N.rows() == 0);
    REQUIRE (N.cols() == 0);
    REQUIRE (UV.rows() == 0);
    REQUIRE (UV.cols() == 0);

    // this bunny have additional data
    REQUIRE (VD.rows() == 35947);
    REQUIRE (VD.cols() == 3);

    // the dummy column contents check
    for(size_t i=0;i<V.rows();++i)
        REQUIRE (VD(i,2) == (double)i);

    REQUIRE (Vheader.size() == 3);
    REQUIRE (Vheader[0] == "confidence" );
    REQUIRE (Vheader[1] == "intensity" );
    REQUIRE (Vheader[2] == "dummy_data" );

    // Face datashould have only one column
    REQUIRE (Fheader.size() == 1);
    REQUIRE (Fheader[0] == "face_data" );
    REQUIRE (FD.rows() == F.rows());
    REQUIRE (FD.cols() == 1);

    // the dummy column contents check
    for(size_t i=0;i<F.rows();++i)
        REQUIRE (FD(i,0) == (double)i);

    // Edge data should have only one column
    REQUIRE (Eheader.size() == 1);
    REQUIRE (Eheader[0] == "edge_data" );
    REQUIRE (ED.rows() == E.rows());
    REQUIRE (ED.cols() == 1);

    // the dummy column contents check
    for(size_t i=0;i<E.rows();++i)
        REQUIRE (ED(i,0) == (double)i);


    // there are comments
    REQUIRE (comments.size() == 2);
}



TEST_CASE("writePLY: bunny.ply float", "[igl]")
{
    std::ifstream f(test_common::data_path("bunny.ply"));
    REQUIRE (f.good());
    f.close();

    Eigen::MatrixXf V1,N1,UV1,VD1,FD1,ED1;
    std::vector<std::string> Vheader1,Fheader1,Eheader1,comments1;

    Eigen::MatrixXi F1,E1;

    // load test data first
    REQUIRE (igl::readPLY(test_common::data_path("bunny.ply"), V1, F1, E1, N1, UV1, VD1,Vheader1, FD1,Fheader1, ED1,Eheader1,comments1));

    // add more data
    Vheader1.push_back("dummy_data");
    Eigen::VectorXf dummy_data(V1.rows());
    for(size_t i=0;i<V1.rows();++i)
        dummy_data(i)=(double)i;
    Eigen::MatrixXf VD2(V1.rows(),VD1.cols()+1);
    VD2<<VD1,dummy_data;


    // test that saving preserves all the data, including new data column
    REQUIRE (igl::writePLY("writePLY_test_bunny_float.ply", V1, F1, E1, N1, UV1, VD2, Vheader1, FD1,Fheader1, ED1, Eheader1, comments1, igl::FileEncoding::Binary));

    Eigen::MatrixXf V,N,UV,VD,FD,ED;
    Eigen::MatrixXi F,E;
    std::vector<std::string> Vheader,Fheader,Eheader,comments;

    // test that saving preserves all the data
    REQUIRE (igl::readPLY("writePLY_test_bunny_float.ply", V, F, E, N, UV, VD,Vheader, FD,Fheader, ED,Eheader, comments));

    REQUIRE (V.rows() == 35947);
    REQUIRE (V.cols() == 3);
    REQUIRE (F.rows() == 69451);
    REQUIRE (F.cols() == 3);
    // no edge data
    REQUIRE (E.rows() == 0);
    REQUIRE (E.cols() == 0);

    // no normals or texture coordinates
    REQUIRE (N.rows() == 0);
    REQUIRE (N.cols() == 0);
    REQUIRE (UV.rows() == 0);
    REQUIRE (UV.cols() == 0);

    // this bunny have additonal data
    REQUIRE (VD.rows() == 35947);
    REQUIRE (VD.cols() == 3);

    // the dummy column contents check
    for(size_t i=0;i<V.rows();++i)
        REQUIRE (VD(i,2) == (double)i);

    REQUIRE (Vheader.size() == 3);
    REQUIRE (Vheader[0] == "confidence" );
    REQUIRE (Vheader[1] == "intensity" );
    REQUIRE (Vheader[2] == "dummy_data" );

    // no Face data or edge data
    REQUIRE (FD.rows() == 0);
    REQUIRE (FD.cols() == 0);
    REQUIRE (Fheader.size() == 0);

    REQUIRE (ED.rows() == 0);
    REQUIRE (ED.cols() == 0);
    REQUIRE (Eheader.size() == 0);

    // there are comments
    REQUIRE (comments.size() == 2);
}
