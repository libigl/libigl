#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/unique_simplices.h>
#include <igl/get_seconds.h>
#include <igl/AABB.h>
#include <igl/barycenter.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/colon.h>
#include <igl/point_simplex_squared_distance.h>
#include <igl/matlab_format.h>
#include <limits>
#include <igl/randperm.h>
#include <igl/avg_edge_length.h>
#include <igl/find.h>
#include <deque>

const int MAX_RUNS = 10;
using AABB = igl::AABB<Eigen::MatrixXd,3>;

template <typename Scalar, int Dim>
void pad_box(const Scalar pad, Eigen::AlignedBox<Scalar,Dim> & box)
{
  box.min().array() -= pad;
  box.max().array() += pad;
}

template <typename DerivedV, int DIM>
void validate(const igl::AABB<DerivedV,DIM> * tree, int depth = 0)
{
  if(tree->is_leaf())
  {
    assert(tree->m_primitive >= 0 || tree->is_root());
  }
  if(tree->m_left)
  {
    assert(tree->m_box.contains(tree->m_left->m_box));
    assert(tree->m_left->m_parent == tree);
    validate(tree->m_left,depth+1);
  }
  if(tree->m_right)
  {
    assert(tree->m_box.contains(tree->m_right->m_box));
    assert(tree->m_right->m_parent == tree);
    validate(tree->m_right,depth+1);
  }
}

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
  return validate(root);
}

template <typename DerivedV, int DIM>
void print(const igl::AABB<DerivedV,DIM> * tree, int depth = 0)
{
  const auto indent = std::string(depth*2,' ');
  printf("%s%p",indent.c_str(),tree);
  if(tree->is_leaf())
  {
    printf(" [%d]",tree->m_primitive);
  }
  printf("\n");
  if(tree->m_left)
  {
    assert(tree->m_box.contains(tree->m_left->m_box));
    assert(tree->m_left->m_parent == tree);
    print(tree->m_left,depth+1);
  }
  if(tree->m_right)
  {
    assert(tree->m_box.contains(tree->m_right->m_box));
    assert(tree->m_right->m_parent == tree);
    print(tree->m_right,depth+1);
  }
}


template <typename DerivedV, int DIM>
int size(const igl::AABB<DerivedV,DIM> * tree)
{
  return 1 + (tree->m_left?size(tree->m_left):0) + (tree->m_right?size(tree->m_right):0);
}
template <typename DerivedV, int DIM>
int height(const igl::AABB<DerivedV,DIM> * tree)
{
  return 1 + std::max((tree->m_left?height(tree->m_left):0),(tree->m_right?height(tree->m_right):0));
}


template <typename DerivedV, int DIM>
void box_faces(
  const igl::AABB<DerivedV,DIM> & tree,
  Eigen::MatrixXd & P,
  Eigen::MatrixXi & Q,
  Eigen::VectorXi & D)
{
  static_assert(DIM == 3,"Assumes 3D");
  const int num_nodes = size(&tree);
  P.resize(8*num_nodes,DIM);
  Q.resize(6*num_nodes,4);
  D.resize(6*num_nodes);
  int d = 0;
  int p = 0;
  int q = 0;
  std::vector<std::pair<const igl::AABB<DerivedV,DIM> *,int> > stack;
  stack.push_back({&tree,0});
  while(!stack.empty())
  {
    const auto pair = stack.back();
    const auto * node = pair.first;
    const int depth = pair.second;
    D(d++) = depth;
    D(d++) = depth;
    D(d++) = depth;
    D(d++) = depth;
    D(d++) = depth;
    D(d++) = depth;
    stack.pop_back();
    const auto & box = node->m_box;
    auto min_corner = box.min();
    auto max_corner = box.max();
    // shrink by 3%
    min_corner = min_corner + 0.03*(max_corner-min_corner);
    max_corner = max_corner - 0.03*(max_corner-min_corner);

    Q.row(q++) << p+0,p+1,p+2,p+3;
    Q.row(q++) << p+0,p+1,p+5,p+4;
    Q.row(q++) << p+1,p+2,p+6,p+5;
    Q.row(q++) << p+2,p+3,p+7,p+6;
    Q.row(q++) << p+3,p+0,p+4,p+7;
    Q.row(q++) << p+4,p+5,p+6,p+7;
    P.row(p++) = min_corner;
    P.row(p++) = Eigen::RowVector3d(max_corner[0],min_corner[1],min_corner[2]);
    P.row(p++) = Eigen::RowVector3d(max_corner[0],max_corner[1],min_corner[2]);
    P.row(p++) = Eigen::RowVector3d(min_corner[0],max_corner[1],min_corner[2]);
    P.row(p++) = Eigen::RowVector3d(min_corner[0],min_corner[1],max_corner[2]);
    P.row(p++) = Eigen::RowVector3d(max_corner[0],min_corner[1],max_corner[2]);
    P.row(p++) = max_corner;
    P.row(p++) = Eigen::RowVector3d(min_corner[0],max_corner[1],max_corner[2]);
    if(node->m_left)
    {
      stack.push_back({node->m_left,depth+1});
    }
    if(node->m_right)
    {
      stack.push_back({node->m_right,depth+1});
    }
  }
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
  box_faces(tree,TV,TQ,TD);
  const auto update_edges = [&]()
  {
    //std::cout<<igl::matlab_format(TV,"TV")<<std::endl;
    //std::cout<<igl::matlab_format_index(TQ,"TQ")<<std::endl;
    //std::cout<<igl::matlab_format_index(TD,"TD")<<std::endl;

    Eigen::MatrixXi TQd = TQ(igl::find((TD.array()==depth).eval()),Eigen::all);

    Eigen::MatrixXi TE(4*TQd.rows(),2);
    TE <<
      TQd.col(0), TQd.col(1),
      TQd.col(1), TQd.col(2),
      TQd.col(2), TQd.col(3),
      TQd.col(3), TQd.col(0);
    igl::unique_simplices(Eigen::MatrixXi(TE),TE);

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


int main(int argc, char *argv[])
{
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
  validate(tree);

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

}
