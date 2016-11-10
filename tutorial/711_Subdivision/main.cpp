#include <igl/read_triangle_mesh.h>
#include <igl/loop.h>
#include <igl/upsample.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/viewer/Viewer.h>
#include <Eigen/Core>
#include <iostream>

#include "tutorial_shared_path.h"

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace igl;
  Eigen::MatrixXi F,UF,LF,FF;
  Eigen::MatrixXd V,UV,LV,FV;
  bool show_swept_volume = false;
  read_triangle_mesh(
      TUTORIAL_SHARED_PATH "/decimated-knight.off",V,F);
  cout<<R"(Usage:
1  Original mesh
2  In-plane upsampled mesh
3  Loop subdivided mesh
4  False barycentric subdivision
)";
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V,F);
  viewer.data.set_face_based(true);
  igl::upsample(V,F,UV,UF);
  igl::loop(V,F,LV,LF);
  igl::false_barycentric_subdivision(V,F,FV,FF);

  viewer.callback_key_down =
    [&](igl::viewer::Viewer & viewer, unsigned char key, int mod)->bool
    {
      switch(key)
      {
        default:
          return false;
        case '1':
        {
          viewer.data.clear();
          viewer.data.set_mesh(V,F);
          break;
        }
        case '2':
        {
          viewer.data.clear();
          viewer.data.set_mesh(UV,UF);
          break;
        }
        case '3':
        {
          viewer.data.clear();
          viewer.data.set_mesh(LV,LF);
          break;
        }
        case '4':
        {
          viewer.data.clear();
          viewer.data.set_mesh(FV,FF);
          break;
        }
      }
      viewer.data.set_face_based(true);
      return true;
    };
  viewer.launch();
}
