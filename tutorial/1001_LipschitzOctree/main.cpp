#include <igl/read_triangle_mesh.h>
#include <igl/lipschitz_octree.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/AABB.h>
#include <igl/fast_winding_number.h>
#include <igl/signed_distance.h>
#include <igl/grid.h>
#include <igl/colormap.h>
#include <igl/matlab_format.h>
#include <igl/marching_cubes.h>
#include <igl/unique_sparse_voxel_corners.h>
#include <igl/get_seconds.h>


// Edges of a cube box with minimum corner at origin and side length h.
void box_edges(
  const Eigen::RowVector3d & origin,
  const double h,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & E)
{
  V.resize(8,3);
  V<< 0,0,0,  0,0,1,  0,1,0,  0,1,1,  1,0,0,  1,0,1,  1,1,0,  1,1,1;
  V *= h;
  V.rowwise() += origin;
  E.resize(12,2);
  E<< 0,1,  0,2,  0,4,  1,3,  1,5,  2,3,  2,6,  3,7,  4,5,  4,6,  5,7,  6,7;
}

// Edges of cubes placed at all corners of octree with origin and cell side length h.
void all_box_edges(
  const Eigen::RowVector3d & origin,
  const double h,
  const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & ijk,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & E)
{
  V.resize(ijk.rows()*8,3);
  E.resize(ijk.rows()*12,2);
  for(int c = 0;c<ijk.rows();c++)
  {
    Eigen::MatrixXd Vc;
    Eigen::MatrixXi Ec;
    box_edges(
      origin + h * Eigen::RowVector3d(ijk(c,1),ijk(c,0),ijk(c,2)), h, Vc,Ec);
    V.middleRows(c*8,8) = Vc;
    E.middleRows(c*12,12) = Ec.array() + c*8;
  }
}

int main(int argc, char * argv[])
{
  IGL_TICTOC_LAMBDA;
  // Read mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(!igl::read_triangle_mesh
     (argc>1?argv[1]: TUTORIAL_SHARED_PATH "/decimated-knight.off",V,F)) {
    std::cout << "Failed to load mesh." << std::endl;
  }

  // Prepare unsigned and signed distance function handles
  igl::AABB<Eigen::MatrixXd,3> aabb;
  aabb.init(V,F);
  const double bbd = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
  double offset_factor = 0.05;
  igl::FastWindingNumberBVH fwn_bvh;
  igl::fast_winding_number(V.cast<float>().eval(), F, 2, fwn_bvh);
  const std::function<double(const Eigen::RowVector3d &)>
    sdf = [&]( const Eigen::RowVector3d & p) -> double
  {
    double w = fast_winding_number(fwn_bvh, 2, p.template cast<float>().eval());         
    double s = 1.-2.*std::abs(w);  
    int i;
    Eigen::RowVector3d c;
    return s*std::sqrt(aabb.squared_distance(V,F,p,i,c)) - offset_factor*bbd;
  };
  // This can't be auto when using static library or compiler gets confused.
  const std::function<double(const Eigen::RowVector3d &)>
    udf = [&]( const Eigen::RowVector3d & p) -> double
  {
    return std::abs(sdf(p));
  };

  // Centered bounding cube (root of octree) with padding.
  double h0 = (V.colwise().maxCoeff() - V.colwise().minCoeff()).maxCoeff();
  int max_depth = 5;
  // Pad bounding box by max_depth
  //{
  //  double leaf_h = h0 / (1 << max_depth);
  //  h0 += leaf_h*2;
  //}
  // Pad to root
  h0 *= 2;
  Eigen::RowVector3d origin = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff()).array() - 0.5*h0;


  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  const int main_index = viewer.selected_data_index;
  Eigen::MatrixXd rV;
  Eigen::MatrixXi rE;
  box_edges(origin,h0,rV,rE);
  viewer.data().set_edges(rV,rE,Eigen::RowVector3d(1,1,0));
  viewer.data().set_face_based(true);
  viewer.data().show_lines = false;
  viewer.append_mesh();
  const int aux_index = viewer.selected_data_index;

  const auto update = [&]()
  {
    if(max_depth<=8)
    {
      // For comparison, how long would it take to compute on the dense grid
      tictoc();
      Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> GV;
      Eigen::RowVector3i res(1<<max_depth,1<<max_depth,1<<max_depth);
      igl::grid(res,GV);
      GV *= h0;
      GV.rowwise() += origin;
      Eigen::VectorXd GS(GV.rows());
      igl::parallel_for( GV.rows(), [&](const int u)
      {
        // evaluate the function at the corner
        GS(u) = sdf(GV.row(u));
      },1000);
      printf("                    n³ grid: %0.7f seconds\n",tictoc());
    }else
    {
      printf("                    n³ grid: [omitted]\n");
    }

    tictoc();
    // Identify the octree cells that could contain the surface
    Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> ijk;
    igl::lipschitz_octree(origin,h0,max_depth,udf,ijk);
    printf("           lipschitz_octree: %0.7f seconds\n",tictoc());

    // Gather the corners of those leaf cells
    const double h = h0 / (1 << max_depth);
    Eigen::Matrix<int,Eigen::Dynamic,8,Eigen::RowMajor> J;
    Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> unique_ijk;
    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> unique_corner_positions;
    igl::unique_sparse_voxel_corners(origin,h0,max_depth,ijk,unique_ijk,J,unique_corner_positions);
    printf("unique_sparse_voxel_corners: %0.7f seconds\n",tictoc());
    /// Evaluate the signed distance function at the corners
    Eigen::VectorXd S(unique_corner_positions.rows());
    //for(int u = 0;u<unique_corner_positions.rows();u++)
    igl::parallel_for(
      unique_corner_positions.rows(),
      [&](const int u)
      {
        // evaluate the function at the corner
        S(u) = sdf(unique_corner_positions.row(u));
      },1000);
      printf("                        sdf: %0.7f seconds\n",tictoc());
    // Run marching cubes on the sparse set of leaf cells
    Eigen::Matrix<double,Eigen::Dynamic,3> mV;
    Eigen::Matrix<int,Eigen::Dynamic,3> mF;
    igl::marching_cubes( S,unique_corner_positions,J, 0, mV,mF);
    printf("             marching_cubes: %0.7f seconds\n",tictoc());

    // Visualize
    viewer.data(aux_index).clear();
    viewer.data(aux_index).set_mesh(mV,mF);
    viewer.data(aux_index).set_face_based(true);
    viewer.data().show_lines = max_depth <= 7;
    viewer.data(aux_index).set_colors(Eigen::RowVector3d(0.941176,0.305882,0.282353));
    if(max_depth <= 7)
    {
      // If tree is shallow enough, visualize the edges of the octree cells
      Eigen::MatrixXd gV;
      Eigen::MatrixXi gE;
      all_box_edges(origin,h,ijk,gV,gE);
      viewer.data(aux_index).set_edges(gV,gE,Eigen::RowVector3d(1,1,1));
      Eigen::MatrixXd point_colors;
      igl::colormap(
          igl::ColorMapType::COLOR_MAP_TYPE_TURBO,
          S,-S.array().abs().maxCoeff(),S.array().abs().maxCoeff(),
          point_colors);
      viewer.data(aux_index).add_points(unique_corner_positions,point_colors);
      viewer.data(aux_index).point_size = 8;
    }
    printf("                     |ijk| : %ld\n",ijk.rows());
  };

  viewer.callback_key_pressed =
  [&](igl::opengl::glfw::Viewer & viewer, unsigned int key, int mod)->bool
  {
    switch(key) {
      default:
        return false;
      case '=':
      case '+':
      case '-':
      case '_':
        max_depth = std::max(0,max_depth + ((key == '-'||key=='_') ? -1 : 1));
        printf("max_depth = %d\n",max_depth);
        update();
        break;
      case '}':
      case '{':
        offset_factor =
          std::min(std::max(-1.0,offset_factor+(key=='{'?-0.01:0.01)),1.0);
        printf("offset_factor = %f\n",offset_factor);
        update();
        break;
    }
    return true;
  };
  update();
  // Trivial batched function
  const std::function<Eigen::VectorXd(const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> &)>
    udf_batch = [&](const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & P) -> Eigen::VectorXd
  {
    Eigen::VectorXd U(P.rows());
    for(int i = 0;i<P.rows();i++)
    {
      U(i) = udf(P.row(i));
    }
    return U;
  };
  {
    // The batched version should produce the same size output
    Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> ijk_batch;
    igl::lipschitz_octree<true>(origin,h0,max_depth,udf_batch,ijk_batch);
    printf("               |ijk_batch| : %ld\n",ijk_batch.rows());
  }


  std::cout << R"(Usage:
  -,+  Decrease,Increase the depth of the octree and recompute
  {,}  Decrease,Increase the offset factor for the signed distance function
)";
  viewer.launch();
  return 0;
}
