#include <test_common.h>
#include <igl/opengl/glfw/map_texture.h>

TEST_CASE("map_texture: identity","[igl/glfw]") 
{
  // 2 triangle quad
  Eigen::MatrixXd V(4,3);
  V<<0,0,0,
     1,0,0,
     1,1,0,
     0,1,0;
  Eigen::MatrixXi F(2,3);
  F<<0,1,2,
     0,2,3;
  // Random RGBA texture
  const int in_w = 16;
  const int in_h = 16;
  const int in_nc = 4;
  using ArrayXuc = Eigen::Array<unsigned char,Eigen::Dynamic,1>;
  ArrayXuc in_rgba = ArrayXuc::Random(in_h*in_w*in_nc,1);
  std::vector<unsigned char> out_rgba;
  int out_w, out_h, out_nc;
  igl::opengl::glfw::map_texture(
    V,F,V,
    in_rgba.data(),
    in_w,in_h,in_nc,
    out_rgba,
    out_w,out_h,out_nc);
  REQUIRE(out_w == in_w);
  REQUIRE(out_h == in_h);
  REQUIRE(out_nc == in_nc);
  {
    // Map in_rgb
    Eigen::Map<ArrayXuc> out_rgba_map(out_rgba.data(),out_w*out_h*out_nc,1);
    REQUIRE(out_rgba_map.isApprox(in_rgba,0));
  }
}

TEST_CASE("map_texture: transpose","[igl/glfw]") 
{
  // 2 triangle quad
  Eigen::MatrixXd V(4,3);
  V<<0,0,0,
     1,0,0,
     1,1,0,
     0,1,0;
  Eigen::MatrixXi F(2,3);
  F<<0,1,2,
     0,2,3;
  // Random RGBA texture
  const int in_w = 16;
  const int in_h = 16;
  const int in_nc = 4;
  using ArrayXuc = Eigen::Array<unsigned char,Eigen::Dynamic,1>;
  ArrayXuc in_rgba = ArrayXuc::Random(in_h*in_w*in_nc,1);
  std::vector<unsigned char> out_rgba;
  int out_w, out_h, out_nc;
  Eigen::MatrixXd U(4,3);
  U<< 
    1,1,0,
    1,0,0,
    0,0,0,
    0,1,0;
  igl::opengl::glfw::map_texture(
    V,F,U,
    in_rgba.data(),
    in_w,in_h,in_nc,
    out_rgba,
    out_w,out_h,out_nc);
  REQUIRE(out_w == in_w);
  REQUIRE(out_h == in_h);
  REQUIRE(out_nc == in_nc);
  {
    // Treat each 4 unsigned chars as a single int32_t and then transpose
    using FourChars = std::int32_t;
    Eigen::Map<Eigen::Array<FourChars,Eigen::Dynamic,Eigen::Dynamic>>
      pixel_map(reinterpret_cast<FourChars*>(out_rgba.data()), out_w,out_h);
    pixel_map.transposeInPlace();
    // Array of chars
    Eigen::Map<ArrayXuc> out_rgba_map(out_rgba.data(),out_w*out_h*out_nc,1);
    REQUIRE(out_rgba_map.isApprox(in_rgba,0));
  }
}
