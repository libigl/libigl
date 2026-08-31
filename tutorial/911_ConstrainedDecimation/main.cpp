// 911_ConstrainedDecimation
//
// Decimate (qslim) a mesh while optionally enforcing several geometric
// constraints, each expressed as an igl::decimate callback "decorator" that can
// be freely cascaded:
//
//   1. block self-intersections            (igl::block_self_intersections)
//   2. stay inside an offset "clearance"   (igl::block_intersections_with_input)
//      shape (a static obstacle mesh)
//   3. stay strictly outside the input      (a *custom* decorator defined below
//      surface, a "barrier on distance"     using a signed-distance barrier)
//
// With all three enabled the decimated armadillo remains free of
// self-intersections and stays inside the thin shell between the original
// surface and its outward offset.
//
// Press ' ' to cycle through the results.
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/qslim.h>
#include <igl/decimate_callback_types.h>
#include <igl/block_self_intersections.h>
#include <igl/block_intersections_with_input.h>
#include <igl/AABB.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_edge_normals.h>
#include <igl/signed_distance.h>
#include <igl/voxel_grid.h>
#include <igl/marching_cubes.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/placeholders.h>
#include <igl/predicates/find_self_intersections.h>
#include <igl/unique.h>
#include <igl/get_seconds.h>
#include <igl/STR.h>
#include <Eigen/Core>
#include <cstdlib>
#include <memory>
#include <tuple>
#include <vector>

// A custom decorator: reject any collapse whose merged vertex would be pushed
// (more than `tol`) to the *inside* of the static surface (VB,FB). This is a
// simple distance barrier keeping the decimated surface on/outside the input —
// "strict outside decimation".
//
// It demonstrates that callers can author their own decorators against the same
// interface used by the shipped ones.
static igl::decimate_pre_post_collapse_callbacks_decorator block_inside_surface(
  const Eigen::MatrixXd & VB,
  const Eigen::MatrixXi & FB,
  const double tol)
{
  using Tree = igl::AABB<Eigen::MatrixXd,3>;
  auto tree = std::make_shared<Tree>();
  tree->init(VB,FB);
  // Angle-weighted pseudonormal data for robust signed distance.
  auto V = std::make_shared<Eigen::MatrixXd>(VB);
  auto F = std::make_shared<Eigen::MatrixXi>(FB);
  auto FN = std::make_shared<Eigen::MatrixXd>();
  auto VN = std::make_shared<Eigen::MatrixXd>();
  auto EN = std::make_shared<Eigen::MatrixXd>();
  auto EMAP = std::make_shared<Eigen::VectorXi>();
  {
    Eigen::MatrixXi E;
    igl::per_face_normals(*V,*F,*FN);
    igl::per_vertex_normals(
      *V,*F,igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,*FN,*VN);
    igl::per_edge_normals(
      *V,*F,igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,*FN,*EN,E,*EMAP);
  }
  return [=](
    const Eigen::MatrixXd & /*V*/,
    const Eigen::MatrixXi & /*F*/,
    const int /*orig_m*/,
    igl::decimate_pre_collapse_callback  & pre_collapse,
    igl::decimate_post_collapse_callback & /*post_collapse*/)
  {
    igl::decimate_pre_collapse_callback inner_pre = pre_collapse;
    pre_collapse = [=](
      const Eigen::MatrixXd & Vd,
      const Eigen::MatrixXi & Fd,
      const Eigen::MatrixXi & E,
      const Eigen::VectorXi & EMAP_d,
      const Eigen::MatrixXi & EF,
      const Eigen::MatrixXi & EI,
      const igl::min_heap<std::tuple<double,int,int>> & Q,
      const Eigen::VectorXi & EQ,
      const Eigen::MatrixXd & C,
      const int e) -> bool
    {
      if(!inner_pre(Vd,Fd,E,EMAP_d,EF,EI,Q,EQ,C,e)) { return false; }
      const Eigen::RowVector3d q = C.row(e);
      const double s = igl::signed_distance_pseudonormal(
        *tree,*V,*F,*FN,*VN,*EN,*EMAP,q);
      // Positive is outside the closed surface; reject going inside.
      return s >= -tol;
    };
  };
}

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(
    argc>1 ? argv[1] : TUTORIAL_SHARED_PATH "/armadillo.obj", V,F);
  printf("input: %ld vertices, %ld faces\n",(long)V.rows(),(long)F.rows());

  const double diag = igl::bounding_box_diagonal(V);
  const double offset = 0.02*diag;

  // Build an outward offset surface (OV,OF) via a background grid signed
  // distance + marching cubes. This is used as a static "clearance" obstacle:
  // the decimated mesh is not allowed to poke outside of it.
  Eigen::MatrixXd OV;
  Eigen::MatrixXi OF;
  {
    double t = igl::get_seconds();
    Eigen::MatrixXd GV;
    Eigen::RowVector3i side;
    igl::voxel_grid(V,offset*2.0,64,1,GV,side);
    Eigen::VectorXd S;
    {
      Eigen::VectorXi I; Eigen::MatrixXd C,N;
      igl::signed_distance(
        GV,V,F,igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
        -std::numeric_limits<double>::infinity(),
         std::numeric_limits<double>::infinity(),
        S,I,C,N);
    }
    igl::marching_cubes(S,GV,side(0),side(1),side(2),offset,OV,OF);
    printf("offset surface: %ld faces (%.2fs)\n",(long)OF.rows(),
      igl::get_seconds()-t);
  }

  const int target_m = std::max<int>(1000,int(F.rows()*0.05));

  // Assemble the shippable decorators.
  igl::decimate_pre_post_collapse_callbacks_decorator dec_self =
    igl::block_self_intersections();
  igl::decimate_pre_post_collapse_callbacks_decorator dec_offset =
    igl::block_intersections_with_input(OV,OF);
  igl::decimate_pre_post_collapse_callbacks_decorator dec_inside =
    block_inside_surface(V,F,1e-4*diag);

  struct Result { std::string name; Eigen::MatrixXd U; Eigen::MatrixXi G; };
  std::vector<Result> results;

  const auto run = [&](
    const std::string & name,
    const std::vector<igl::decimate_pre_post_collapse_callbacks_decorator> & decs)
  {
    double t = igl::get_seconds();
    Eigen::MatrixXd U; Eigen::MatrixXi G; Eigen::VectorXi J,I;
    igl::qslim(V,F,target_m,decs,U,G,J,I);
    const double secs = igl::get_seconds()-t;

    // Report self-intersections.
    int num_si = 0;
    {
      Eigen::MatrixXi IF;
      Eigen::Array<bool,Eigen::Dynamic,1> CP;
      igl::predicates::find_self_intersections(U,G,false,IF,CP);
      Eigen::VectorXi u; igl::unique(IF,u); num_si = u.size();
    }
    // Report shell containment: signed distance to input in [~0, offset].
    int num_inside=0, num_outside=0;
    {
      Eigen::VectorXd Sd; Eigen::VectorXi Id; Eigen::MatrixXd Cd,Nd;
      igl::signed_distance(
        U,V,F,igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
        -std::numeric_limits<double>::infinity(),
         std::numeric_limits<double>::infinity(),Sd,Id,Cd,Nd);
      num_inside  = (Sd.array() < -1e-3*diag).count();
      num_outside = (Sd.array() > offset+1e-3*diag).count();
    }
    printf("%-28s %6ld F  %5.2fs  self-int:%5d  inside-input:%5d  "
           "outside-offset:%5d\n",
      name.c_str(),(long)G.rows(),secs,num_si,num_inside,num_outside);
    results.push_back({name,U,G});
  };

  run("unconstrained",              {});
  run("+self",                      {dec_self});
  run("+self +inside(input)",       {dec_self,dec_inside});
  run("+self +inside +offset(all)", {dec_self,dec_inside,dec_offset});

  bool batch = std::getenv("LIBIGL_NO_VIEWER") != nullptr;
  for(int i = 1;i<argc;i++)
  {
    if(std::string(argv[i]) == "--batch") { batch = true; }
  }
  if(batch) { return 0; }

  // Viewer: cycle results with the space bar.
  igl::opengl::glfw::Viewer vr;
  int shown = 0;
  const auto show = [&](int i)
  {
    shown = (i%(int)results.size()+(int)results.size())%(int)results.size();
    vr.data().clear();
    vr.data().set_mesh(results[shown].U,results[shown].G);
    vr.data().set_face_based(true);
    printf("showing [%d/%d]: %s\n",shown+1,(int)results.size(),
      results[shown].name.c_str());
  };
  vr.callback_key_pressed =
    [&](decltype(vr)&, unsigned int key, int)->bool
  {
    if(key==' '){ show(shown+1); return true; }
    return false;
  };
  show(0);
  vr.launch();
  return 0;
}
