#include <igl/AABB.h>
#include <igl/box_faces.h>
#include <igl/colon.h>
#include <igl/find.h>
#include <igl/get_seconds.h>
#include <igl/per_face_normals.h>
#include <igl/quad_edges.h>
#include <igl/read_triangle_mesh.h>
#include <igl/unproject_ray.h>
#include <igl/opengl/glfw/Viewer.h>
#include <limits>
#include <deque>
#include <stdio.h>



int main(int argc, char *argv[])
{
  Eigen::MatrixXd V,N;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    argc>1?argv[1]:TUTORIAL_SHARED_PATH "/decimated-knight.off",V,F);
  // Make mesh into disconnected soup
  V = V(Eigen::Map<Eigen::VectorXi>(F.data(),F.size()), Eigen::all).eval();
  F = Eigen::Map<Eigen::MatrixXi>(igl::colon<int>(0,V.rows()-1).data(),V.rows()/3,3).eval();
  // Cache normals
  igl::per_face_normals(V,F,N);

  // Build an initial tree. Store as pointer so we can update root.
  auto * tree = new igl::AABB<Eigen::MatrixXd, 3>();
  tree->init(V,F);
  // Gather leaf pointers (these are stable)
  const auto leaves = tree->gather_leaves(F.rows());

  // Set up visualization and interaction
  igl::opengl::glfw::Viewer vr;
  Eigen::MatrixXd C = Eigen::RowVector3d(0.9,0.9,0.9).replicate(F.rows(),1);
  vr.data().set_mesh(V,F);
  vr.data().set_face_based(true);
  vr.data().set_colors(C);

  // draw the boxes of the tree at a given depth
  int depth = 0;
  Eigen::MatrixXd TV;
  Eigen::MatrixXi TQ;
  Eigen::VectorXi TD;
  const auto update_edges = [&]()
  {
    Eigen::MatrixXi TQd = TQ(igl::find((TD.array()==depth).eval()),Eigen::all);
    Eigen::MatrixXi TE;
    igl::quad_edges(TQd,TE);
    vr.data().set_edges(TV,TE,Eigen::RowVector3d(1,1,1));
    vr.data().line_width = 2;
  };
  const auto update_tree_vis = [&]()
  {
    igl::box_faces(*tree,TV,TQ,TD);
    update_edges();
  };

  // Update the tree for each triangle in hits which might have moved.
  const auto update_tree = [&](const std::vector<igl::Hit> & hits)
  {
    for(const auto hit : hits)
    {
      const int fid = hit.id;
      Eigen::AlignedBox<double,3> box;
      box.extend(V.row(F(fid,0)).transpose());
      box.extend(V.row(F(fid,1)).transpose());
      box.extend(V.row(F(fid,2)).transpose());
      // This is O(height) which is typically O(log n)
      tree = leaves[fid]->update(box,0)->root();
    }
    // update the visualization. This is O(n) 🙃
    update_tree_vis();
  };

  update_tree_vis();

  std::vector<igl::Hit> hits;
  Eigen::RowVector3d dir;
  Eigen::RowVector3d red(1,0.2,0.2);
  vr.callback_pre_draw = [&](decltype(vr) &)->bool
  {
    // While mouse is down, pull front (back) faces towards (away from) cursor.
    if(hits.size())
    {
      for(const auto hit : hits)
      {
        for(int c : {0,1,2})
        {
          const int fid = hit.id;
          V.row(F(fid,c)) += 
            ((N.row(fid).dot(dir))>0?1:-1) *
            0.001*dir.cast<double>().transpose();
          C.row(fid) = red;
        }
      }
      // Things have moved, so update tree ~O(#hits log n)
      update_tree(hits);
      // update viewer. This is O(n) 🙃
      vr.data().set_vertices(V);
      vr.data().set_colors(C);
    }
    return false;
  };
  vr.callback_mouse_down = [&](decltype(vr) &, int button, int modifier)->bool
  {
    // Shoot a ray through cursor into mesh accelerated by current AABB tree
    const Eigen::Vector2f pos(
      vr.current_mouse_x,vr.core().viewport(3)-vr.current_mouse_y);
    const auto model = vr.core().view;
    const auto proj = vr.core().proj;
    const auto viewport = vr.core().viewport;
    Eigen::RowVector3d src;
    {
      Eigen::Vector3f src_f,dir_f;
      igl::unproject_ray(pos,model,proj,viewport,src_f,dir_f);
      src = src_f.cast<double>().transpose();
      dir = dir_f.cast<double>().transpose();
    }
    if(tree->intersect_ray(V,F,src,dir,hits))
    {
      vr.core().is_animating = true;
      return true;
    }
    return false;
  };
  vr.callback_mouse_up = [&](decltype(vr) &, int button, int modifier)->bool
  {
    hits.clear();
    vr.core().is_animating = false;
    return false;
  };
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
  printf("  ,/. to change tree-vis depth\n");
  printf("  Click and hold to push triangles\n");
  vr.launch();
  // delete the tree (and all subtrees/leaves)
  delete tree;
}
