#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/get_seconds.h>
#include <igl/AABB.h>
#include <igl/barycenter.h>
#include <igl/writePLY.h>
#include <igl/intersection_blocking_collapse_edge_callbacks.h>
#include <igl/STR.h>
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
  if(use_test)
  {
    igl::decimate_pre_collapse_callback  pre_collapse;
    igl::decimate_post_collapse_callback post_collapse;
    intersection_blocking_collapse_edge_callbacks(
      tree,
      pre_collapse,
      post_collapse);
    
    Eigen::VectorXi EMAP;
    Eigen::MatrixXi E,EF,EI;
    igl::edge_flaps(F,E,EMAP,EF,EI);
    igl::min_heap< std::tuple<double,int,int> > Q;
    Eigen::VectorXi EQ;
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(E.rows(),3);
    // Try to do the collapses.
    std::vector<int> edges_to_collapse;
    for(const auto edge : {Eigen::RowVector2i(3,6),Eigen::RowVector2i(6,18)})
    {
      int e = -1;
      const bool found = 
        ( (E.array().col(0)==edge[0] && E.array().col(1)==edge[1])||
          (E.array().col(0)==edge[1] && E.array().col(1)==edge[0])).maxCoeff(&e);
      assert(found);
      edges_to_collapse.push_back(e);
    }
    for(auto & e : edges_to_collapse)
    {
      printf("E(e=%d,:) = %d,%d\n",e,E(e,0),E(e,1));
      bool collapsed = true;
      int e1,e2,f1,f2;
      // pre_collapse only has access to C not p
      C.row(e) = 0.5*(V.row(E(e,0))+V.row(E(e,1)));
      if(pre_collapse(V,F,E,EMAP,EF,EI,Q,EQ,C,e))
      {
        collapsed = igl::collapse_edge(
          e,C.row(e).eval(),
          /*Nsv,Nsf,Ndv,Ndf,*/
          V,F,E,EMAP,EF,EI,e1,e2,f1,f2);
      }else
      {
        // Aborted by pre collapse callback
        collapsed = false;
      }
      printf("!!!!!!!!!!!!!!!!!collapsed? %s\n",collapsed?"✅":"❌");
      post_collapse(V,F,E,EMAP,EF,EI,Q,EQ,C,e,e1,e2,f1,f2,collapsed);
    }
    printf("tree: %p\n",tree);
  }else
  {
    const int target_m = F.rows() * 0.1 + 1;
    for(auto pass : {0,1})
    {
      igl::decimate_pre_collapse_callback  pre_collapse;
      igl::decimate_post_collapse_callback post_collapse;
      igl::decimate_trivial_callbacks(pre_collapse,post_collapse);
      if(pass == 1)
      {
        intersection_blocking_collapse_edge_callbacks(
          pre_collapse, post_collapse, // These will get copied as needed
          tree,
          pre_collapse, post_collapse);
      }
      Eigen::MatrixXd VO = V;
      Eigen::MatrixXi FO = F;
#     warning "Should connect_boundary_to_infinity"
      Eigen::VectorXi EMAP;
      Eigen::MatrixXi E,EF,EI;
      igl::edge_flaps(FO,E,EMAP,EF,EI);
      int m = F.rows();
      const int orig_m = m;
      Eigen::MatrixXd U;
      Eigen::MatrixXi G;
      Eigen::VectorXi J,I;
      tictoc();
      const bool ret = igl::decimate(
        VO, FO,
        igl::shortest_edge_and_midpoint,
        igl::max_faces_stopping_condition(m,orig_m,target_m),
        pre_collapse,
        post_collapse,
        E, EMAP, EF, EI,
        U, G, J, I);
      printf("pass %d in %g sec\n",pass,tictoc());
      igl::writePLY(STR("out-"<<pass<<".ply"),U,G);
      if(pass == 1)
      {
        V = U;
        F = G;
      }
    }
  }

  vis(V,F,*tree);

  tree->validate();
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
  igl::read_triangle_mesh(argc>1?argv[1]:TUTORIAL_SHARED_PATH "/armadillo.obj",V,F);
  // make into soup
  V = V(Eigen::Map<Eigen::VectorXi>(F.data(),F.size()), Eigen::all).eval();
  F = Eigen::Map<Eigen::MatrixXi>(igl::colon<int>(0,V.rows()-1).data(),V.rows()/3,3).eval();

  //F = F.topRows(4).eval();
  //
  Eigen::MatrixXd BC;
  igl::barycenter(V,F,BC);
  Eigen::VectorXd sqrD;
  Eigen::VectorXi I;
  Eigen::VectorXi J = igl::colon<int>(0,F.rows()-1);
  Eigen::MatrixXd C;

  IGL_TICTOC_LAMBDA;
  tictoc();
  igl::AABB<Eigen::MatrixXd, 3> * tree = new igl::AABB<Eigen::MatrixXd, 3>();
  tree->init(V,F);
  printf("tree.init(): %g secs\n",tictoc());
  tictoc();
  tree->squared_distance(V,F,BC,sqrD,I,C);
  assert(I.isApprox(J,0));
  printf("tree.squared_distance(): %g secs\n",tictoc());
  printf("  surface_area: %g\n",tree->internal_surface_area());
  printf("  is_root(): %d\n",tree->is_root());
  printf("  size: %d\n",size(tree));
  printf("  height: %d/%d\n",height(tree),(int)F.rows());
  //print(&tree);
  tree->validate(tree);

  // Gather list of pointers to leaves
  std::vector<igl::AABB<Eigen::MatrixXd,3>*> leaves = tree->gather_leaves(F.rows());
  for(auto * leaf : leaves) { assert(leaf); }
  //print(tree);
  printf("--------------------------------\n");

  // detach and insert each leaf
  const double h = igl::avg_edge_length(V,F);
  const double pad = h;
  {
    // Gather list of pointers to leaves
    std::vector<igl::AABB<Eigen::MatrixXd,3>*> leaves = tree->gather_leaves(F.rows());
    for(auto * leaf : leaves) { assert(leaf); }
    tree = tree->pad(leaves,pad,2);
  }

  tictoc();
  tree->squared_distance(V,F,BC,sqrD,I,C);
  assert(I.isApprox(J,0));
  printf("tree.squared_distance(): %g secs\n",tictoc());
  printf("  surface_area: %g\n",tree->internal_surface_area());
  printf("  is_root(): %d\n",tree->is_root());
  printf("  size: %d\n",size(tree));
  printf("  height: %d/%d\n",height(tree),(int)F.rows());
  //print(tree);


  Eigen::VectorXi RV;
  Eigen::VectorXi RF;
  // Perturb a small subset of the triangles
  {
    {
      igl::randperm(F.rows(),RF);
      RF = RF.topRows(std::min(12,(int)F.rows())).eval();
      RV.resize(RF.size()*3);
      RV << RF, RF.array()+F.rows(), RF.array()+2*F.rows();
    }
    Eigen::MatrixXd TF = 0.1*h*Eigen::MatrixXd::Random(RF.size(),3);
    Eigen::MatrixXd TV(RV.rows(),3);
    TV<<TF,TF,TF;
    V(RV,Eigen::all) += TV;
    igl::barycenter(V,F,BC);
  }
  const int qi = RF(0);


  {
    tictoc();
    for(int i = 0;i<RF.size();i++)
    {
      tree = leaves[RF(i)]->update_primitive(V,F,pad)->root();
    }
    printf("        tree.refit                :          %g secs\n",tictoc());
    tictoc();
    for(int r = 0;r<MAX_RUNS;r++)
    {
      tree->squared_distance(V,F,Eigen::MatrixXd(BC.row(qi)),sqrD,I,C);
    }
    assert(I(0) == qi);
    printf("%d,%g ← tree.squared_distance(0,…):          %g secs\n",I(0),sqrD(0),tictoc()/MAX_RUNS);
  }
  {
    tictoc();
    for(int r = 0;r<MAX_RUNS;r++)
    {
      igl::point_mesh_squared_distance(Eigen::MatrixXd(BC.row(qi)),V,F,sqrD,I,C);
    }
    assert(I(0) == qi);
    printf("%d,%g ← point_mesh_squared_distance(0,…):    %g secs\n",I(0),sqrD(0),tictoc()/MAX_RUNS);
  }

  {
    tictoc();
    double min_dist;
    int min_i;
    for(int r = 0;r<MAX_RUNS;r++)
    {
      min_dist = std::numeric_limits<double>::infinity();
      min_i = -1;
      Eigen::RowVector3d q = Eigen::RowVector3d(BC.row(qi));
      for(int i = 0;i<F.rows();i++)
      {
        Eigen::RowVector3d c;
        double d;
        igl::point_simplex_squared_distance<3>(q,V,F,i,d,c);
        if(d < min_dist)
        {
          min_dist = d;
          min_i = i;
        }
      }
    }
    printf("%d,%g ← point_simplex_squared_distance(0,…): %g secs\n",min_i,min_dist,tictoc()/MAX_RUNS);
    assert(min_i == qi);
  }

  
  vis(V,F,*tree);
  delete tree;


  //igl::AABB<Eigen::MatrixXd, 3> * dynamic = nullptr;
  //std::vector<igl::AABB<Eigen::MatrixXd, 3> > leaves;
  //for(auto rotation_amount : {1, 2, 3})
  //{
  //  // if vector is storing objects, must clear first.
  //  leaves.clear();
  //  // tree is now invalid, but deleting is safe.
  //  delete dynamic;
  //  {
  //    printf("\n--------------------------------\n\n");
  //    printf("rotation_amount: %d\n",rotation_amount);
  //    tictoc();
  //    // The root starts as the first one which will be self-inserted
  //    leaves.resize(F.rows());
  //    dynamic = leaves.data();
  //    for(int i = 0;i<F.rows();i++)
  //    {
  //      auto * leaf = &leaves[i];
  //      // Use the idiotic .init()
  //      leaf->init(V,F,Eigen::MatrixXi(),(Eigen::VectorXi(1)<<i).finished());
  //      dynamic = dynamic->insert(leaf)->root();

  //      if(rotation_amount==1)
  //      {
  //        const bool ret = leaf->rotate();
  //      }

  //      if(rotation_amount>=2)
  //      {
  //        leaf->rotate_lineage();
  //      }
  //    }
  //    if(rotation_amount==3)
  //    {
  //      std::deque<igl::AABB<Eigen::MatrixXd, 3> *> bfs;
  //      bfs.push_back(dynamic);
  //      std::vector<igl::AABB<Eigen::MatrixXd, 3> *> all_nodes;
  //      while(!bfs.empty())
  //      {
  //        auto * node = bfs.back();
  //        bfs.pop_back();
  //        if(node->m_left)
  //        {
  //          bfs.push_back(node->m_left);
  //        }
  //        if(node->m_right)
  //        {
  //          bfs.push_back(node->m_right);
  //        }
  //        all_nodes.push_back(node);
  //      }
  //      while(!all_nodes.empty())
  //      {
  //        auto * node = all_nodes.back();
  //        all_nodes.pop_back();
  //        assert(node);
  //        const bool ret = node->rotate();
  //      }
  //    }
  //    printf("dynamic %g\n",tictoc());
  //    tictoc();
  //    dynamic->squared_distance(V,F,BC,sqrD,I,C);
  //    assert(I.isApprox(J,0));
  //    printf("tree.squared_distance(): %g secs\n",tictoc());
  //    printf("  surface_area: %g\n",dynamic->internal_surface_area());
  //    printf("  is_root(): %d\n",dynamic->is_root());
  //    printf("  size: %d\n",size(dynamic));
  //    printf("  height: %d/%d\n",height(dynamic),(int)F.rows());
  //    //print(dynamic);
  //    validate(dynamic, leaves);
  //  }
  //  printf("********\n");
  //}

  //vis(V,F,*dynamic);

#endif
}
