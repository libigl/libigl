#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/fast_find_self_intersections.h>
#include <igl/unique.h>
#include <igl/remove_unreferenced.h>
#include <igl/get_seconds.h>
#include <igl/AABB.h>
#include <igl/writePLY.h>
#include <igl/intersection_blocking_collapse_edge_callbacks.h>
#include <igl/qslim_optimal_collapse_edge_callbacks.h>
#include <igl/per_vertex_point_to_plane_quadrics.h>
#include <igl/STR.h>
#include <igl/connect_boundary_to_infinity.h>
#include <igl/decimate.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/point_simplex_squared_distance.h>
#include <igl/edge_flaps.h>
#include <igl/decimate_callback_types.h>
#include <igl/find.h>

int main(int argc, char *argv[])
{
  IGL_TICTOC_LAMBDA;
  Eigen::MatrixXd V,V0;
  Eigen::MatrixXi F,F0;
  igl::read_triangle_mesh( 
      argc<=1 ? TUTORIAL_SHARED_PATH "/octopus-low.mesh" :
      argv[1],V,F);
  V0 = V;F0 = F;
  ///////////////////////////////////////////////////////////
  /// Before collapsing starts
  ///////////////////////////////////////////////////////////
  tictoc();
  igl::AABB<Eigen::MatrixXd, 3> * tree = new igl::AABB<Eigen::MatrixXd, 3>();
  tree->init(V,F);
  printf("tree->init: %g\n",tictoc());
  tree->validate();

  const int target_m = F.rows() * 0.1 + 1;
  Eigen::MatrixXd dV[2];
  Eigen::MatrixXd dC[2];
  Eigen::MatrixXi dF[2];
  Eigen::RowVector3d gray(0.9,0.9,0.9);
  for(auto pass : {0,1})
  {
    Eigen::MatrixXd VO;
    Eigen::MatrixXi FO;
    igl::connect_boundary_to_infinity(V,F,VO,FO);
    Eigen::VectorXi EMAP;
    Eigen::MatrixXi E,EF,EI;
    igl::edge_flaps(FO,E,EMAP,EF,EI);

    igl::decimate_cost_and_placement_callback cost_and_placement;
    igl::decimate_pre_collapse_callback  pre_collapse;
    igl::decimate_post_collapse_callback post_collapse;

    // Quadrics per vertex
    typedef std::tuple<Eigen::MatrixXd,Eigen::RowVectorXd,double> Quadric;
    std::vector<Quadric> quadrics;
    igl::per_vertex_point_to_plane_quadrics(VO,FO,EMAP,EF,EI,quadrics);
    // State variables keeping track of edge we just collapsed
    int v1 = -1;
    int v2 = -1;
    // Callbacks for computing and updating metric
    igl::qslim_optimal_collapse_edge_callbacks(
      E,quadrics,v1,v2, cost_and_placement, pre_collapse,post_collapse);

    if(pass == 1)
    {
      igl::intersection_blocking_collapse_edge_callbacks(
        pre_collapse, post_collapse, // These will get copied as needed
        tree,
        pre_collapse, post_collapse);
    }

    int m = F.rows();
    const int orig_m = m;
    Eigen::MatrixXd U;
    Eigen::MatrixXi G;
    Eigen::VectorXi J,I;
    tictoc();
    const bool ret = igl::decimate(
      VO, FO,
      cost_and_placement,
      igl::max_faces_stopping_condition(m,orig_m,target_m),
      pre_collapse,
      post_collapse,
      E, EMAP, EF, EI,
      U, G, J, I);
    G = G(igl::find((J.array()<orig_m).eval()), Eigen::all).eval();
    {
      Eigen::VectorXi _;
      igl::remove_unreferenced(Eigen::MatrixXd(U),Eigen::MatrixXi(G),U,G,_);
    }
    printf("qslim-%22s in %g secs\n",pass?"-intersection-blocking":"",tictoc());
    igl::writePLY(STR("out-"<<pass<<".ply"),U,G);
    dV[pass] = U;
    dF[pass] = G;
    {
      Eigen::VectorXi BI;
      {
        Eigen::MatrixXd EV;
        Eigen::MatrixXi IF,EE;
        Eigen::VectorXi EI;
        igl::fast_find_self_intersections(dV[pass],dF[pass],true,false,IF,EV,EE,EI);
        igl::unique(IF,BI);
      }
      printf("  # self-intersections: %d\n",(int)BI.size());
      dC[pass] = gray.replicate(dF[pass].rows(),1);
      dC[pass](BI,Eigen::all) = 
        Eigen::RowVector3d(0.95,0.15,0.15).replicate(BI.size(),1);
    }
  }
  

  tree->validate();
  assert(tree == tree->root());
  tree = tree->root();

  igl::opengl::glfw::Viewer vr;
  vr.callback_key_pressed = [&](decltype(vr) &,unsigned int key, int mod)
  {
    switch(key)
    {
      case ',':
      case '.':
        vr.data_list[vr.selected_data_index].is_visible = false;
        vr.selected_data_index += (key==','?-1:1) + vr.data_list.size();
        vr.selected_data_index %= vr.data_list.size();
        vr.data_list[vr.selected_data_index].is_visible = true;
        switch(vr.selected_data_index)
        {
          case 0:
            printf("Original mesh\n");
            break;
          case 1:
            printf("Qslim mesh\n");
            break;
          case 2:
            printf("Qslim mesh with intersection blocking\n");
            break;
        }
        return true;
      default:
        return false;
    }
  };

  vr.data().set_mesh(V0,F0);
  vr.data().show_lines = false;
  vr.data().set_face_based(true);
  vr.data().is_visible = false;
  vr.data().set_colors(gray);
  vr.append_mesh();
  vr.data().set_mesh(dV[0],dF[0]);
  vr.data().show_lines = true;
  vr.data().set_face_based(true);
  vr.data().is_visible = false;
  vr.data().set_colors(dC[0]);
  vr.append_mesh();
  vr.data().set_mesh(dV[1],dF[1]);
  vr.data().show_lines = true;
  vr.data().set_face_based(true);
  vr.data().set_colors(dC[1]);
  vr.launch();


}
