#include <test_common.h>

#include <catch2/catch.hpp>

#include <igl/stb/write_image.h>
#include <igl/stb/read_image.h>

namespace {
    // magick -size 17x29 gradient:white-black -depth 8 test.png
    // xxd -i test.png
    #include "test_png.h"
    // convert test.png test.bmp
    #include "test_bmp.h"
    //convert -quality 100% test.png test.jpg
    #include "test_jpg.h"
    // convert test.png test.tga
    #include "test_tga.h"
    // convert test.png test.pgm
    #include "test_pgm.h"

    #include "test_reference.h"

    static bool dump_to_file(const unsigned char *buf,int size,const char* fname)
    {
        //generate file:
        FILE *out;
        if(out=fopen(fname,"wb"))
        {
            bool r=(fwrite(buf,1,size,out)==size);
            r=(fclose(out)==0)&&r;
            return r;
        } else {
            return false;
        }
    }
}

TEST_CASE("read_image","[igl/stb]")
{
    // convert -verbose test.png -depth 8 -colorspace Gray -format pgm test.pgm
    // xxd -i test.pgm # skip header
    REQUIRE(dump_to_file(test_png,test_png_len,"test.png"));
    REQUIRE(dump_to_file(test_bmp,test_bmp_len,"test.bmp"));
    REQUIRE(dump_to_file(test_tga,test_tga_len,"test.tga"));
    REQUIRE(dump_to_file(test_jpg,test_jpg_len,"test.jpg"));
    REQUIRE(dump_to_file(test_pgm,test_pgm_len,"test.pgm"));

    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
    Eigen::Map<const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor > > ref_image(ref_image_raw, 29, 17) ;

    for(auto i:{"test.png","test.bmp","test.tga","test.jpg", "test.pgm"}) 
    {
        WHEN(i) {
            REQUIRE(igl::stb::read_image(i,R,G,B,A));

            //check the size
            REQUIRE(R.rows()==17);
            REQUIRE(R.cols()==29);

            REQUIRE(G.rows()==17);
            REQUIRE(G.cols()==29);

            REQUIRE(B.rows()==17);
            REQUIRE(B.cols()==29);

            REQUIRE(A.rows()==17);
            REQUIRE(A.cols()==29);

            // check the contents, it is transposed and upside down 
            for(int i=0;i<17;++i)
                for(int j=0;j<29;++j)
                {
                    REQUIRE( (unsigned int)R(i,28-j) == (unsigned int)ref_image(j,i) );
                    REQUIRE( (unsigned int)G(i,28-j) == (unsigned int)ref_image(j,i) );
                    REQUIRE( (unsigned int)A(i,j) == 255 );
                }
        }
    }
}


TEST_CASE("write_image","[igl/stb]")
{
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(17,29),G(17,29),B(17,29),A(17,29),A2(17,29);

    for(int i=0;i<17;++i)
        for(int j=0;j<29;++j)
        {
            R(i,j)=(i+j*17)%255;
            G(i,j)=(i+j*17+1)%255;
            B(i,j)=(i+j*17+2)%255;
            A(i,j)=i+j;
            A2(i,j)=255;
        }

    for(auto i:{"check.png","check.bmp","check.tga","check.jpg"}) 
    {
        WHEN(i) {

            if(i=="check.png" || i=="check.tga")
                REQUIRE(igl::stb::write_image(R, G, B, A, i, 100));
            else // Alpha channel somehow affects other colors in BMP, even though it's not saved
                REQUIRE(igl::stb::write_image(R, G, B, A2, i, 100));

            // read it back
            Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> _R,_G,_B,_A;
            REQUIRE(igl::stb::read_image(i, _R,_G,_B,_A));

            //only PNG and TGA files support alpha channel
            if(i=="check.png" || i=="check.tga")
                REQUIRE((A.array()==_A.array()).all());

            if(i=="check.jpg")
            {
                // somehow jpeg loses information even with 100% quality
                REQUIRE((R.cast<float>()-_R.cast<float>()).norm()<R.size()/20.);
                REQUIRE((G.cast<float>()-_G.cast<float>()).norm()<G.size()/20.);
                REQUIRE((B.cast<float>()-_B.cast<float>()).norm()<B.size()/20.);
                
            } else {
                REQUIRE((R.array()==_R.array()).all());
                REQUIRE((G.array()==_G.array()).all());
                REQUIRE((B.array()==_B.array()).all());
            }
        }
    }
}