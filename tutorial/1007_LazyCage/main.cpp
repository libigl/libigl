// 1007_LazyCage
//
// A "lazy cage" is a coarse mesh that encloses an input surface, obtained the
// lazy way: offset the input outward by some amount σ, then decimate the dense
// offset down to a target face count. The catch is choosing σ. Too small and
// there is no room to decimate to the target without the coarsening mesh poking
// back through the input (or through itself); large enough and it succeeds but
// the cage is looser than necessary.
//
// igl::lazy_cage bisects on σ to find the *smallest* offset for which the dense
// offset can be decimated to `num_faces` while igl::block_self_intersections and
// igl::block_intersections_with_input(input) refuse every offending collapse.
// The result is a tight, self-intersection-free cage that provably encloses the
// input.
//
// The isosurfacing grid resolution (--grid) trades cage fidelity for speed.
//
// Usage: 1007_LazyCage [mesh] [--faces N] [--grid G] [--max-sigma frac]
//                      [--iters K] [--qslim] [--metric sigma|volume|area]
//                      [--sparse]
// With no --grid, the isosurfacing resolution is chosen automatically from the
// target face count (igl::lazy_cage_default_grid_size). --sparse switches the
// dense signed-distance grid for a Lipschitz-octree + sparse marching cubes
// extraction (useful for very dense grids / large inputs).
// Press ' ' to toggle the cage wireframe.
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/lazy_cage.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/centroid.h>
#include <igl/doublearea.h>
#include <igl/get_seconds.h>
#include <igl/predicates/find_self_intersections.h>
#include <igl/predicates/find_intersections.h>
#include <Eigen/Core>
#include <cstdlib>
#include <string>

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  std::string mesh_path = TUTORIAL_SHARED_PATH "/decimated-knight.off";
  int num_faces = 100;
  int grid_size = 0;   // 0 => automatic from num_faces
  double max_sigma_frac = 0.1;
  int num_iters = 12;
  bool use_qslim = false;
  igl::LazyCageMetric metric = igl::LAZY_CAGE_METRIC_SIGMA;
  igl::LazyCageGridMode grid_mode = igl::LAZY_CAGE_GRID_DENSE;
  for(int i = 1;i<argc;i++)
  {
    const std::string a = argv[i];
    if(a == "--faces" && i+1<argc) { num_faces = std::atoi(argv[++i]); }
    else if(a == "--grid" && i+1<argc) { grid_size = std::atoi(argv[++i]); }
    else if(a == "--max-sigma" && i+1<argc) { max_sigma_frac = std::atof(argv[++i]); }
    else if(a == "--iters" && i+1<argc) { num_iters = std::atoi(argv[++i]); }
    else if(a == "--qslim") { use_qslim = true; }
    else if(a == "--sparse") { grid_mode = igl::LAZY_CAGE_GRID_SPARSE; }
    else if(a == "--metric" && i+1<argc)
    {
      const std::string m = argv[++i];
      if(m == "volume") { metric = igl::LAZY_CAGE_METRIC_VOLUME; }
      else if(m == "area") { metric = igl::LAZY_CAGE_METRIC_SURFACE_AREA; }
      else { metric = igl::LAZY_CAGE_METRIC_SIGMA; }
    }
    else if(a == "--batch") { /* handled below */ }
    else if(!a.empty() && a[0] != '-') { mesh_path = a; }
  }
  igl::read_triangle_mesh(mesh_path, V,F);
  const double diag = igl::bounding_box_diagonal(V);
  const int grid_used =
    grid_size>0 ? grid_size : igl::lazy_cage_default_grid_size(num_faces);
  const char * metric_name =
    metric==igl::LAZY_CAGE_METRIC_VOLUME ? "volume" :
    metric==igl::LAZY_CAGE_METRIC_SURFACE_AREA ? "surface-area" : "sigma";
  printf("input: %ld vertices, %ld faces\n",(long)V.rows(),(long)F.rows());
  printf("target cage faces: %d, grid: %d%s, decimation: %s, metric: %s, %s\n",
    num_faces,grid_used,grid_size>0?"":" (auto)",
    use_qslim?"qslim":"shortest-edge/midpoint",metric_name,
    grid_mode==igl::LAZY_CAGE_GRID_SPARSE?"sparse":"dense");

  Eigen::MatrixXd CV;
  Eigen::MatrixXi CF;
  double sigma = 0;
  const double t = igl::get_seconds();
  const bool ok = igl::lazy_cage(
    V,F,num_faces,grid_size,max_sigma_frac*diag,num_iters,use_qslim,metric,
    grid_mode,CV,CF,sigma);
  double vol = 0; { Eigen::RowVector3d c; igl::centroid(CV,CF,c,vol); }
  double area = 0; { Eigen::VectorXd dblA; igl::doublearea(CV,CF,dblA); area = 0.5*dblA.sum(); }
  printf("lazy_cage: %s  sigma=%g (%.4f*diag)  cage faces=%ld  "
         "volume=%g  area=%g  (%.2fs)\n",
    ok?"reached target":"FAILED (target not reached in range)",
    sigma,sigma/diag,(long)CF.rows(),std::abs(vol),area,igl::get_seconds()-t);

  // Verify the cage: self-intersection free, disjoint from the input.
  {
    Eigen::MatrixXi IF; Eigen::Array<bool,Eigen::Dynamic,1> CP;
    igl::predicates::find_self_intersections(CV,CF,false,IF,CP);
    printf("  cage self-intersections: %ld triangle pairs\n",(long)IF.rows());
    Eigen::MatrixXi IFio; Eigen::Array<bool,Eigen::Dynamic,1> CPio;
    igl::predicates::find_intersections(CV,CF,V,F,false,IFio,CPio);
    printf("  cage-vs-input intersections: %ld triangle pairs\n",(long)IFio.rows());
  }

  bool batch = std::getenv("LIBIGL_NO_VIEWER") != nullptr;
  for(int i = 1;i<argc;i++)
  {
    if(std::string(argv[i]) == "--batch") { batch = true; }
  }
  if(batch) { return 0; }

  igl::opengl::glfw::Viewer vr;
  vr.data().set_mesh(V,F);
  vr.data().set_face_based(true);
  vr.data().show_lines = false;
  vr.append_mesh();
  vr.data().set_mesh(CV,CF);
  vr.data().show_lines = true;
  vr.data().set_face_based(true);
  // Draw the cage as a translucent shell around the input.
  vr.data().set_colors(Eigen::RowVector3d(0.2,0.6,1.0));
  bool cage_visible = true;
  const int cage_id = vr.selected_data_index;
  vr.callback_key_pressed =
    [&](decltype(vr)&, unsigned int key, int)->bool
  {
    if(key==' ')
    {
      cage_visible = !cage_visible;
      vr.data_list[cage_id].is_visible = cage_visible;
      return true;
    }
    return false;
  };
  vr.launch();
  return 0;
}
