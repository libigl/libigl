#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/cgal/CSGTree.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <Eigen/Core>

#include "tutorial_shared_path.h"

int main(int argc, char * argv[])
{
  using namespace Eigen;
  using namespace igl::copyleft::cgal;
  using namespace std;
  using namespace igl;
  cout<<R"(
[,]  Toggle between boolean sub-tree operations
)";

  MatrixXi FA,FB,FC,FD,FE;
  MatrixXd VA,VB,VC,VD,VE;
  // Read in inputs as double precision floating point meshes
  read_triangle_mesh(TUTORIAL_SHARED_PATH "/cube.obj"     ,VA,FA);
  read_triangle_mesh(TUTORIAL_SHARED_PATH "/sphere.obj"   ,VB,FB);
  read_triangle_mesh(TUTORIAL_SHARED_PATH "/xcylinder.obj",VC,FC);
  read_triangle_mesh(TUTORIAL_SHARED_PATH "/ycylinder.obj",VD,FD);
  read_triangle_mesh(TUTORIAL_SHARED_PATH "/zcylinder.obj",VE,FE);
  igl::opengl::glfw::Viewer viewer;

  int num_views = 5+4;
  int view_id = num_views-1;
  const auto & update = [&]()
  {
    viewer.data().clear();
    // CSGTree templated on type of F
    VectorXd I;
    const auto & set_mesh = 
      [&](const MatrixXd & V, const MatrixXi & F, const int i)
    {
      viewer.data().set_mesh(V,F);
      I = VectorXd::Constant(F.rows(),1,i);
    };
    switch(view_id)
    {
      case 0:
        set_mesh(VA,FA,5);
        break;
      case 1:
        set_mesh(VB,FB,4);
        break;
      case 2:
        set_mesh(VC,FC,3);
        break;
      case 3:
        set_mesh(VD,FD,2);
        break;
      case 4:
        set_mesh(VE,FE,1);
        break;
      default:
      {
        CSGTree M;
        Matrix<MatrixXi::Index,Dynamic,1> J;
        switch(view_id)
        {
          case 5:
            // Compute result of (A ∩ B)
            M = {{VA,FA},{VB,FB},"i"};
            J = M.J().array()+0;
            break;
          case 6:
            // Compute result of (C ∪ D)
            M = {{VC,FC},{VD,FD},"u"};
            J = M.J().array()+FA.rows()+FB.rows();
            break;
          case 7:
            // Compute result of (C ∪ D) ∪ E
            M = {{{VC,FC},{VD,FD},"u"},{VE,FE},"u"};
            J = M.J().array()+FA.rows()+FB.rows();
            break;
          case 8:
            // Compute result of (A ∩ B) \ ((C ∪ D) ∪ E)
            M = {{{VA,FA},{VB,FB},"i"},{{{VC,FC},{VD,FD},"u"},{VE,FE},"u"},"m"};
            J = M.J().array()+0;
            break;
          default:
            assert(false && "unknown view id");
        }
        viewer.data().set_mesh(M.cast_V<MatrixXd>(),M.F());
        I.resize(M.F().rows(),1);
        // Compute colors based on original facets
        for(int f = 0;f<M.F().rows();f++)
        {
          const int j = J(f);
          I(f) = 
            (int)(j<FA.rows())+
            (int)(j<FA.rows()+FB.rows())+
            (int)(j<FA.rows()+FB.rows()+FC.rows())+
            (int)(j<FA.rows()+FB.rows()+FC.rows()+FD.rows())+
            (int)(j<FA.rows()+FB.rows()+FC.rows()+FD.rows()+FE.rows());
        }
      }
    }

    MatrixXd C;
    jet(I,1,5,C);
    viewer.data().set_colors(C);
  };
  update();

  viewer.callback_key_down = 
    [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)->bool
    {
      switch(key)
      {
        case ']':
          view_id = (view_id+1)%num_views;
          break;
        case '[':
          view_id = (view_id+num_views-1)%num_views;
          break;
        default:
          return false;
      }
      update();
      return true;
    };
  viewer.launch();
}
