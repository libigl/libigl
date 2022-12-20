#include <test_common.h>

#include <catch2/catch.hpp>
#include <igl/MshLoader.h>
#include <igl/MshSaver.h>

TEST_CASE("MshSaver","[igl]")
{
    igl::MshLoader msh_loader1(test_common::data_path("sphere_lowres_TMS_1-0001_Magstim_70mm_Fig8_nii_scalar.msh"));
    REQUIRE( msh_loader1.is_element_map_identity());
    
    igl::MshSaver msh_saver("test_binary_sphere.msh",true);
    msh_saver.save_mesh( msh_loader1.get_nodes(), 
        msh_loader1.get_elements(),
        msh_loader1.get_elements_lengths(), 
        msh_loader1.get_elements_types(),
        msh_loader1.get_elements_tags()[1]);

    for(size_t i=0;i<msh_loader1.get_node_fields_names().size();++i)
    {
        if(msh_loader1.get_node_fields_components()[i]==3)
            msh_saver.save_vector_field(msh_loader1.get_node_fields_names()[i],msh_loader1.get_node_fields()[i]);
        else
            msh_saver.save_scalar_field(msh_loader1.get_node_fields_names()[i],msh_loader1.get_node_fields()[i]);
    }

    for(size_t i=0;i<msh_loader1.get_element_fields_names().size();++i)
    {
        if(msh_loader1.get_element_fields_components()[i]==3)
            msh_saver.save_elem_vector_field(msh_loader1.get_element_fields_names()[i],msh_loader1.get_element_fields()[i]);
        else
            msh_saver.save_elem_scalar_field(msh_loader1.get_element_fields_names()[i],msh_loader1.get_element_fields()[i]);
    }

    igl::MshLoader msh_loader("test_binary_sphere.msh");

    for(size_t i=0;i<msh_loader1.get_elements().size();++i)
        REQUIRE(msh_loader.get_elements()[i] == msh_loader1.get_elements()[i]);

    for(size_t i=0;i<msh_loader1.get_elements_lengths().size();++i)
        REQUIRE(msh_loader.get_elements_lengths()[i] == msh_loader1.get_elements_lengths()[i]);

    for(size_t i=0;i<msh_loader1.get_elements_types().size();++i)
        REQUIRE(msh_loader.get_elements_types()[i] == msh_loader1.get_elements_types()[i]);

    for(size_t j=0;j<2;++j)
        for(size_t i=0;i<msh_loader1.get_elements_tags()[j].size();++i)
            REQUIRE(msh_loader.get_elements_tags()[j][i] == msh_loader1.get_elements_tags()[j][i]);

    REQUIRE(msh_loader.get_node_fields_names().size() == msh_loader1.get_node_fields_names().size());

    for(size_t i=0;i<msh_loader1.get_node_fields_names().size();++i)
    {
        REQUIRE(msh_loader.get_node_fields_names()[i] == msh_loader1.get_node_fields_names()[i]);
        REQUIRE(msh_loader.get_node_fields_components()[i] == msh_loader.get_node_fields_components()[i]);

        for(size_t j=0;j<msh_loader1.get_node_fields()[i].size();++j)
            REQUIRE(msh_loader1.get_node_fields()[i][j] == msh_loader.get_node_fields()[i][j]);
    }

    REQUIRE(msh_loader.get_element_fields_names().size() == msh_loader1.get_element_fields_names().size());
    for(size_t i=0;i<msh_loader1.get_element_fields_names().size();++i)
    {
        REQUIRE(msh_loader.get_element_fields_names()[i] == msh_loader1.get_element_fields_names()[i]);
        REQUIRE(msh_loader.get_element_fields_components()[i] == msh_loader1.get_element_fields_components()[i]);

        for(size_t j=0;j<msh_loader1.get_element_fields()[i].size();++j)
            REQUIRE(msh_loader1.get_element_fields()[i][j] == msh_loader.get_element_fields()[i][j]);
    }

}

