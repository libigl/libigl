#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_unreferenced.h>
#include <igl/get_seconds.h>
#include <igl/AABB.h>
#include <igl/barycenter.h>
#include <igl/writePLY.h>
#include <igl/intersection_blocking_collapse_edge_callbacks.h>
#include <igl/qslim_optimal_collapse_edge_callbacks.h>
#include <igl/per_vertex_point_to_plane_quadrics.h>
#include <igl/STR.h>
#include <igl/connect_boundary_to_infinity.h>
#include <igl/box_faces.h>
#include <igl/quad_edges.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/colon.h>
#include <igl/writeDMAT.h>
#include <igl/decimate.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/is_self_intersecting.h>
#include <igl/point_simplex_squared_distance.h>
#include <igl/edge_flaps.h>
#include <igl/decimate_callback_types.h>
#include <igl/decimate_trivial_callbacks.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/collapse_edge.h>
#include <igl/matlab_format.h>
#include <igl/randperm.h>
#include <igl/avg_edge_length.h>
#include <igl/find.h>
#include <limits>
#include <deque>

const int MAX_RUNS = 10;
using AABB = igl::AABB<Eigen::MatrixXd,3>;

template <typename DerivedV, int DIM>
void validate(
    const igl::AABB<DerivedV,DIM> * root,
    const std::vector<igl::AABB<DerivedV,DIM> > & leaves)
{
  // Check that all leaves are in the tree
  for(int i = 0;i<leaves.size();i++)
  {
    const auto * leaf = &leaves[i];
    assert(leaf->m_primitive == i);
    assert(leaf->root() == root);
  }
  return root->validate();
}

void vis(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const AABB & tree)
{
  igl::opengl::glfw::Viewer vr;
  vr.data().set_mesh(V,F);
  vr.data().set_face_based(true);
  int depth = 0;
  Eigen::MatrixXd TV;
  Eigen::MatrixXi TQ;
  Eigen::VectorXi TD;
  igl::box_faces(tree,TV,TQ,TD);
  const auto update_edges = [&]()
  {
    Eigen::MatrixXi TQd = TQ(igl::find((TD.array()==depth).eval()),Eigen::all);
    Eigen::MatrixXi TE;
    igl::quad_edges(TQd,TE);
    //Eigen::MatrixXi TF(TQd.rows()*2,3);
    //TF<< 
    //  TQd.col(0),TQd.col(1),TQd.col(2),
    //  TQd.col(0),TQd.col(2),TQd.col(3);
    //vr.append_mesh();
    //vr.data().set_mesh(TV,TF);
    vr.data().set_edges(TV,TE,Eigen::RowVector3d(1,1,1));
    vr.data().line_width = 2;
  };
  update_edges();
  vr.callback_key_pressed = [&](decltype(vr) &,unsigned int key, int mod)
  {
    switch(key)
    {
      case ',':
      case '.':
        depth = std::max(depth+(key==','?-1:1),0);
        update_edges();
        return true;
      default:
        return false;
    }
  };
  vr.launch();
}



///
/// @param[in] orig_pre_collapse  Original pre-collapse callback
/// @param[in] orig_post_collapse Original post-collapse callback
///

int main(int argc, char *argv[])
{
  IGL_TICTOC_LAMBDA;
#if true
  Eigen::MatrixXd V,V0;
  Eigen::MatrixXi F,F0;
  const bool use_test = argc<=1;
  igl::read_triangle_mesh(use_test?"/Users/alecjacobson/Downloads/simple-intersection-collapse-2.ply":argv[1],V,F);
  V0 = V;F0 = F;
  ///////////////////////////////////////////////////////////
  /// Before collapsing starts
  ///////////////////////////////////////////////////////////
  tictoc();
  igl::AABB<Eigen::MatrixXd, 3> * tree = new igl::AABB<Eigen::MatrixXd, 3>();
  tree->init(V,F);
  printf("tree->init: %g\n",tictoc());
  tree->validate();

  // Dummy
  const int target_m = F.rows() * 0.1 + 1;
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
    //cost_and_placement = igl::shortest_edge_and_midpoint;
    //igl::decimate_trivial_callbacks(pre_collapse,post_collapse);

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

    printf("pass %d in %g sec\n",pass,tictoc());
    igl::writePLY(STR("out-"<<pass<<".ply"),U,G);
    if(pass == 1)
    {
      V = U;
      F = G;
    }
  }
  

  vis(V,F,*tree);

  tree->validate();
  assert(tree == tree->root());
  tree = tree->root();
  delete tree;
  exit(1);

  //Eigen::MatrixXd K = Eigen::MatrixXd::Constant(F.rows(),3,0.9);
  //for(const auto f : new_one_ring)
  //{
  //  K.row(f) = Eigen::RowVector3d(0.3,0.9,0.3);
  //}
  //Eigen::MatrixXd K0 = Eigen::MatrixXd::Constant(F.rows(),3,0.9);
  //for(const auto f : old_one_ring)
  //{
  //  K0.row(f) = Eigen::RowVector3d(0.9,0.3,0.3);
  //}


  igl::opengl::glfw::Viewer vr;
  //vr.data().set_mesh(V0,F0);
  //vr.data().set_face_based(true);
  //vr.data().set_colors(K0);
  //vr.data().show_faces = false;
  //vr.data().line_color = Eigen::RowVector4f(1,1,1,1);
  //vr.data().double_sided = true;
  //vr.append_mesh();
  vr.data().set_mesh(V,F);
  //vr.data().set_colors(K);
  vr.data().set_face_based(true);
  vr.data().show_lines = false;
  vr.data().double_sided = true;

  //Eigen::MatrixXd TV;
  //Eigen::MatrixXi TQ,TE;
  //box_faces(big_box,0,TV,TQ);
  //quad_edges(TQ,TE);
  //vr.append_mesh();
  //vr.data().set_edges(TV,TE,Eigen::RowVector3d(1,1,1));
  //vr.data().line_width = 2;

  
  vr.launch();

#else
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(argv[1],V,F);

  tictoc();
  igl::AABB<Eigen::MatrixXd, 3> tree;
  tree.init(V,F);
  printf("static: %g\n",tictoc());

  Eigen::MatrixXd BC;
  igl::barycenter(V,F,BC);
  Eigen::VectorXi I;
  Eigen::VectorXd sqrD;
  Eigen::MatrixXd C;
  tictoc();
  printf("static: %g\n",tictoc());
  tree.squared_distance(V,F,BC,sqrD,I,C);
  Eigen::VectorXi J = igl::colon<int>(0,F.rows()-1);
  assert(I.isApprox(J,0));
  printf("tree.squared_distance(): %g secs\n",tictoc());
  printf("  surface_area: %g\n",tree.internal_surface_area());
  printf("  is_root(): %d\n",tree.is_root());
  printf("  size: %d\n",tree.size());
  printf("  height: %d/%d\n",tree.height(),(int)F.rows());


  igl::AABB<Eigen::MatrixXd, 3> * dynamic = nullptr;
  std::vector<igl::AABB<Eigen::MatrixXd, 3> > leaves;
  for(auto rotation_amount : {1, 2, 3})
  {
    // if vector is storing objects, must clear first.
    leaves.clear();
    // tree is now invalid, but deleting is safe.
    delete dynamic;
    {
      printf("\n--------------------------------\n\n");
      printf("rotation_amount: %d\n",rotation_amount);
      tictoc();
      // The root starts as the first one which will be self-inserted
      leaves.resize(F.rows());
      dynamic = leaves.data();
      for(int i = 0;i<F.rows();i++)
      {
        auto * leaf = &leaves[i];
        // Use the idiotic .init()
        leaf->init(V,F,Eigen::MatrixXi(),(Eigen::VectorXi(1)<<i).finished());
        dynamic = dynamic->insert(leaf)->root();

        if(rotation_amount==1)
        {
          const bool ret = leaf->rotate();
        }

        if(rotation_amount>=2)
        {
          leaf->rotate_lineage();
        }
      }
      if(rotation_amount==3)
      {
        std::deque<igl::AABB<Eigen::MatrixXd, 3> *> bfs;
        bfs.push_back(dynamic);
        std::vector<igl::AABB<Eigen::MatrixXd, 3> *> all_nodes;
        while(!bfs.empty())
        {
          auto * node = bfs.back();
          bfs.pop_back();
          if(node->m_left)
          {
            bfs.push_back(node->m_left);
          }
          if(node->m_right)
          {
            bfs.push_back(node->m_right);
          }
          all_nodes.push_back(node);
        }
        while(!all_nodes.empty())
        {
          auto * node = all_nodes.back();
          all_nodes.pop_back();
          assert(node);
          const bool ret = node->rotate();
        }
      }
      printf("dynamic %g\n",tictoc());
      tictoc();
      dynamic->squared_distance(V,F,BC,sqrD,I,C);
      assert(I.isApprox(J,0));
      printf("tree.squared_distance(): %g secs\n",tictoc());
      printf("  surface_area: %g\n",dynamic->internal_surface_area());
      printf("  is_root(): %d\n",dynamic->is_root());
      printf("  size: %d\n",dynamic->size());
      printf("  height: %d/%d\n",dynamic->height(),(int)F.rows());
      //print(dynamic);
      validate(dynamic, leaves);
    }
    printf("********\n");
  }

  vis(V,F,*dynamic);

#endif
}
