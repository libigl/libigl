// 1006_ConstrainedDecimation
//
// "Progressive hulls" [Sander et al. 2000] decimate a closed mesh while placing
// every new vertex *outside* the current surface, so the coarse output strictly
// *contains* the input. That makes a practically useful question well posed:
//
//   can we build a coarse enclosing hull that stays within a given clearance —
//   i.e., inside an outward offset of the input — and is free of
//   self-intersections?
//
// This tutorial answers yes by cascading igl::decimate callback "decorators"
// onto the progressive-hulls decimation:
//
//   1. igl::block_self_intersections        keep the hull free of self-intersections
//   2. igl::block_intersections_with_input   keep the hull inside an outward
//                                            offset (a static "clearance" mesh)
//   3. a *custom* decorator (defined below)  cap how far the hull may bulge out
//                                            (a signed-distance barrier)
//
// Unconstrained, the progressive hull bulges well past a tight offset (and can
// self-intersect); with the clearance decorator on, no collapse is allowed to
// push the hull outside the offset shell.
//
// Press ' ' to cycle through the results.
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/progressive_hulls.h>
#include <igl/decimate_callback_types.h>
#include <igl/block_self_intersections.h>
#include <igl/block_intersections_with_input.h>
#include <igl/AABB.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_edge_normals.h>
#include <igl/signed_distance.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/predicates/find_self_intersections.h>
#include <igl/unique.h>
#include <igl/get_seconds.h>
#include <Eigen/Core>
#include <cstdlib>
#include <limits>
#include <memory>
#include <tuple>
#include <vector>

// A custom decorator demonstrating the same interface as the shipped ones:
// reject any collapse whose merged vertex would be placed farther than
// `max_dist` *outside* the surface (VB,FB). For progressive hulls (which only
// ever place vertices outside) this is a pure "barrier on distance" that caps
// how far the enclosing hull may bulge away from the input.
static igl::decimate_pre_post_collapse_callbacks_decorator cap_distance_outside(
  const Eigen::MatrixXd & VB,
  const Eigen::MatrixXi & FB,
  const double max_dist)
{
  using Tree = igl::AABB<Eigen::MatrixXd,3>;
  auto tree = std::make_shared<Tree>();
  tree->init(VB,FB);
  auto V  = std::make_shared<Eigen::MatrixXd>(VB);
  auto F  = std::make_shared<Eigen::MatrixXi>(FB);
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
      // Positive is outside; reject placements more than max_dist outside.
      return s <= max_dist;
    };
  };
}

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // Args: first non-flag arg is the mesh path; --offset <frac> and
  // --target <frac> are fractions of the bbox diagonal / input face count.
  std::string mesh_path = TUTORIAL_SHARED_PATH "/armadillo.obj";
  double offset_frac = 0.02;
  double target_frac = 0.02;
  for(int i = 1;i<argc;i++)
  {
    const std::string a = argv[i];
    if(a == "--offset" && i+1<argc) { offset_frac = std::atof(argv[++i]); }
    else if(a == "--target" && i+1<argc) { target_frac = std::atof(argv[++i]); }
    else if(a == "--batch") { /* handled below */ }
    else if(!a.empty() && a[0] != '-') { mesh_path = a; }
  }
  igl::read_triangle_mesh(mesh_path, V,F);
  printf("input: %ld vertices, %ld faces\n",(long)V.rows(),(long)F.rows());

  const double diag = igl::bounding_box_diagonal(V);
  const double offset = offset_frac*diag;

  // Build the outward "clearance" offset by moving each vertex along its
  // (outward) normal by `offset`. Exact-resolution and fast for any offset.
  // Mild self-overlap at sharp convexities is harmless: it is only used as a
  // static obstacle for winding-number distance / intersection queries.
  Eigen::MatrixXd OV;
  Eigen::MatrixXi OF = F;
  {
    Eigen::MatrixXd N;
    igl::per_vertex_normals(V,F,N);
    OV = V + offset*N;
    printf("offset: %g (%.4f*diag), outward clearance shape %ld faces\n",
      offset,offset_frac,(long)OF.rows());
  }

  const int target_m = std::max<int>(100,int(F.rows()*target_frac));
  printf("target faces: %d (%.3f*input)\n",target_m,target_frac);

  const double tol = 1e-3*diag;

  // Shipped decorators + one custom decorator.
  igl::decimate_pre_post_collapse_callbacks_decorator dec_self =
    igl::block_self_intersections();
  igl::decimate_pre_post_collapse_callbacks_decorator dec_clearance =
    igl::block_intersections_with_input(OV,OF);
  igl::decimate_pre_post_collapse_callbacks_decorator dec_cap =
    cap_distance_outside(V,F,offset);

  struct Result {
    std::string name; Eigen::MatrixXd U; Eigen::MatrixXi G;
    int num_si=0, not_contained=0, outside_offset=0;
  };
  std::vector<Result> results;

  const auto run = [&](
    const std::string & name,
    const std::vector<igl::decimate_pre_post_collapse_callbacks_decorator> & decs)
  {
    double t = igl::get_seconds();
    Eigen::MatrixXd U; Eigen::MatrixXi G; Eigen::VectorXi J;
    igl::copyleft::progressive_hulls(V,F,target_m,decs,U,G,J);
    const double secs = igl::get_seconds()-t;

    // Self-intersections of the hull.
    int num_si = 0;
    {
      Eigen::MatrixXi IF;
      Eigen::Array<bool,Eigen::Dynamic,1> CP;
      igl::predicates::find_self_intersections(U,G,false,IF,CP);
      Eigen::VectorXi u; igl::unique(IF,u); num_si = u.size();
    }
    // Containment: input vertices should be inside/on the hull (signed distance
    // to the hull <= 0). Count those that poke outside the hull.
    int not_contained = 0;
    {
      Eigen::VectorXd S; Eigen::VectorXi I; Eigen::MatrixXd C,N;
      igl::signed_distance(
        V,U,G,igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
        -std::numeric_limits<double>::infinity(),
         std::numeric_limits<double>::infinity(),S,I,C,N);
      not_contained = (S.array() > tol).count();
    }
    // Clearance: hull vertices should stay inside the outward offset (signed
    // distance to the offset mesh <= 0). Count those that poke outside it.
    int outside_offset = 0;
    {
      Eigen::VectorXd S; Eigen::VectorXi I; Eigen::MatrixXd C,N;
      igl::signed_distance(
        U,OV,OF,igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
        -std::numeric_limits<double>::infinity(),
         std::numeric_limits<double>::infinity(),S,I,C,N);
      outside_offset = (S.array() > tol).count();
    }
    printf("%-26s %6ld F  %6.2fs  self-int:%5d  input-outside-hull:%5d  "
           "hull-outside-offset:%5d\n",
      name.c_str(),(long)G.rows(),secs,num_si,not_contained,outside_offset);
    results.push_back({name,U,G,num_si,not_contained,outside_offset});
  };

  // "unconstrained" contains the input but may self-intersect and bulge past
  // the offset. "+self" removes self-intersections. "+self +clearance" also
  // keeps the hull within the outward offset (growing the face count as
  // collapses that would breach the offset are refused). "+self +custom-cap"
  // uses the custom signed-distance barrier instead of the offset mesh.
  run("unconstrained",            {});
  run("+self",                    {dec_self});
  run("+self +clearance",         {dec_self,dec_clearance});
  run("+self +custom-cap",        {dec_self,dec_cap});

  bool batch = std::getenv("LIBIGL_NO_VIEWER") != nullptr;
  for(int i = 1;i<argc;i++)
  {
    if(std::string(argv[i]) == "--batch") { batch = true; }
  }
  if(batch) { return 0; }

  // Viewer: cycle results with the space bar; the input is shown for reference.
  igl::opengl::glfw::Viewer vr;
  int shown = 0;
  const auto show = [&](int i)
  {
    shown = (i%(int)results.size()+(int)results.size())%(int)results.size();
    vr.data().clear();
    vr.data().set_mesh(results[shown].U,results[shown].G);
    vr.data().set_face_based(true);
    vr.data().show_lines = true;
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
