#include <test_common.h>
#include <igl/intersection_blocking_collapse_edge_callbacks.h>
#include <igl/edge_flaps.h>
#include <igl/decimate_callback_types.h>
#include <igl/AABB.h>
#include <igl/collapse_edge.h>
#include <igl/min_heap.h>
#include <tuple>

TEST_CASE("intersection_blocking_collapse_edge_callbacks: simple", "[igl]")
{
  // A mesh where collapsing edge e1 will create a self-intersection and
  // collapsed e2 will not
  Eigen::MatrixXd V(21,3);
  V<< 
    0,20,0,
    0,30,0,
    10,10,0,
    10,20,-30,
    10,30,0,
    20,0,0,
    20,10,30,
    20,20,0,
    30,0,0,
    30,10,0,
    10,17,-15,
    7,23,15,
    13,23,15,
    10,30,30,
    20,30,0,
    0,10,3,
    20,5,15,
    25,0,0,
    25,5,15,
    25,10,15,
    30,5,0;
  Eigen::MatrixXi F(22,3);
  F<< 
    5,17,16,
    8,20,18,
    8,18,17,
    9,19,20,
    6,16,18,
    6,18,19,
    18,16,17,
    19,18,20,
    2,5,16,
    7,6,19,
    2,16,6,
    7,19,9,
    2,3,0,
    0,3,1,
    3,4,1,
    2,6,3,
    6,7,3,
    3,7,4,
    10,11,12,
    12,11,13,
    4,7,14,
    2,0,15;
  igl::AABB<Eigen::MatrixXd, 3> * tree = new igl::AABB<Eigen::MatrixXd, 3>();
  tree->init(V,F);


  igl::decimate_pre_collapse_callback  pre_collapse;
  igl::decimate_post_collapse_callback post_collapse;
  igl::intersection_blocking_collapse_edge_callbacks(
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
    REQUIRE(found);
    edges_to_collapse.push_back(e);
  }
  std::vector<bool> should_collapse = {false,true};
  for(int i=0;i<edges_to_collapse.size();i++)
  {
    const auto e = edges_to_collapse[i];
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
    REQUIRE(collapsed == should_collapse[i]);
    post_collapse(V,F,E,EMAP,EF,EI,Q,EQ,C,e,e1,e2,f1,f2,collapsed);
  }

  {
    auto * root = tree->root();
    const auto leaves = tree->gather_leaves();
    for(int f = 0;f<F.rows();f++)
    {
      if((F.row(f).array() == IGL_COLLAPSE_EDGE_NULL).all())
      {
        continue;
      }
      REQUIRE(f < leaves.size());
      auto * leaf = leaves[f];
      REQUIRE(root == leaf->root());
      // check containment all the way up to
      auto * node = leaf;
      while(node)
      {
        REQUIRE(node->m_box.contains(V.row(F(f,0)).transpose()));
        REQUIRE(node->m_box.contains(V.row(F(f,1)).transpose()));
        REQUIRE(node->m_box.contains(V.row(F(f,2)).transpose()));
        auto * parent = node->m_parent;
        if(parent)
        {
          REQUIRE(parent->m_box.contains(node->m_box));
        }
        node = parent;
      }
    }
  }

  tree = tree->root();
  delete tree;
}
