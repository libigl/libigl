#include <test_common.h>

#include <catch2/catch.hpp>
#include <igl/MshLoader.h>


TEST_CASE("MshLoader","[igl]")
{
    igl::MshLoader msh_loader(test_common::data_path("sphere_lowres_TMS_1-0001_Magstim_70mm_Fig8_nii_scalar.msh"));

    //
    // size of nodes - 3xNodes (x,y,z)
    //
    REQUIRE( msh_loader.get_nodes().size()    == (398+4506)*3 );

    // 
    // size of elements: number of triangles x 3 
    //                   number of tetrahedra x 4
    // 
    REQUIRE( msh_loader.get_elements().size()         == (8988*3 + 25937*4) );
    // all these should have the same size
    REQUIRE( msh_loader.get_elements_lengths().size() == (8988 + 25937) );
    REQUIRE( msh_loader.get_elements_ids().size()     == (8988 + 25937) );
    REQUIRE( msh_loader.get_elements_types().size()   == (8988 + 25937) );

    // the test file should have identity map
    REQUIRE( msh_loader.is_element_map_identity() );

    // there should be 12 tags: 6 surface structures and 6 volume structures 
    // index 1st structure tag: Volumes
    msh_loader.index_structures(1); 

    REQUIRE( msh_loader.is_element_map_identity() );
    REQUIRE( msh_loader.get_structures().size() == 6+6 );

    // verify that all triangle elements have length of 3 and 
    // tetrahedral elements have length 4

    for(auto t =msh_loader.get_structures().begin(); 
             t!=msh_loader.get_structures().end(); ++t )
    {
        REQUIRE(( t->el_type == igl::MshLoader::ELEMENT_TRI || 
                  t->el_type == igl::MshLoader::ELEMENT_TET ));

        auto element_range = msh_loader.get_structure_index().equal_range(*t);

        for(auto e=element_range.first; e != element_range.second; ++e)
        { 
            if( t->el_type == igl::MshLoader::ELEMENT_TRI )
            {
                REQUIRE( msh_loader.get_elements_lengths()[e->second] == 3 );
            } else {
                // t->first.el_type==MshLoader::ELEMENT_TET
                REQUIRE( msh_loader.get_elements_lengths()[e->second] == 4 );
            }
        }
    }
}

TEST_CASE("MshLoader_Non_Sequential", "[igl]")
{
    // minimal .msh file with non-sequential tags
    std::string tmp_path = test_common::data_path("_tmp_nonseq_tag.msh");
    {
        std::ofstream out(tmp_path);
        out << "$MeshFormat\n"
            << "2.2 0 8\n"
            << "$EndMeshFormat\n"
            << "$Nodes\n"
            << "3\n"
            << "100 0.0 0.0 0.0\n"
            << "200 1.0 0.0 0.0\n"
            << "300 0.0 1.0 0.0\n"
            << "$EndNodes\n"
            << "$Elements\n"
            << "1\n"
            << "50 2 0 100 200 300\n"
            << "$EndElements\n";
    }

    igl::MshLoader msh_loader(tmp_path);

    // 3 nodes stored
    REQUIRE(msh_loader.get_nodes().size() == 9);

    // node tag 100 -> dense index 0: position (0, 0, 0)
    REQUIRE(msh_loader.get_nodes()[0] == 0.0);  // x
    REQUIRE(msh_loader.get_nodes()[1] == 0.0);  // y
    REQUIRE(msh_loader.get_nodes()[2] == 0.0);  // z

    // node tag 200 -> dense index 1: position (1, 0, 0)
    REQUIRE(msh_loader.get_nodes()[3] == 1.0);  // x
    REQUIRE(msh_loader.get_nodes()[4] == 0.0);  // y
    REQUIRE(msh_loader.get_nodes()[5] == 0.0);  // z

    // node tag 300 -> dense index 2: position (0, 1, 0)
    REQUIRE(msh_loader.get_nodes()[6] == 0.0);  // x
    REQUIRE(msh_loader.get_nodes()[7] == 1.0);  // y
    REQUIRE(msh_loader.get_nodes()[8] == 0.0);  // z

    REQUIRE(msh_loader.get_elements().size() == 3);
    REQUIRE(msh_loader.get_elements_lengths()[0] == 3);
    REQUIRE(msh_loader.get_elements_types()[0] == igl::MshLoader::ELEMENT_TRI);

    REQUIRE(msh_loader.get_elements()[0] == 0);
    REQUIRE(msh_loader.get_elements()[1] == 1);
    REQUIRE(msh_loader.get_elements()[2] == 2);

    REQUIRE(!msh_loader.is_element_map_identity());

    // cleanup
    std::remove(tmp_path.c_str());
}