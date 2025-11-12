#include <igl/read_triangle_mesh.h>
#include <igl/box_simplices.h>
#include <igl/eytzinger_aabb.h>
#include <igl/lipschitz_octree.h>
#include <igl/unique_sparse_voxel_corners.h>
#include <igl/eytzinger_aabb_sdf.h>
#include <igl/unique.h>
#include <igl/readDMAT.h>
#include <igl/grid.h>
#include <igl/edges.h>
#include <igl/marching_cubes.h>
#include <igl/parallel_for.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/matlab_format.h>
#include <igl/variable_radius_offset.h>
#include <igl/get_seconds.h>
#include <igl/SphereMeshWedge.h>

#include <limits>

#include <Eigen/Geometry>

// What kind of offset is being computed?
enum class OffsetType
{
  EDGE_OFFSET = 0,
  TRIANGLE_OFFSET = 1,
  VARIABLE_RADIUS_OFFSET = 2
} offset = OffsetType::VARIABLE_RADIUS_OFFSET;

///////////////////////////////////////////////////////////////////////
/// Some helper functions for triangle and line segment distance
///////////////////////////////////////////////////////////////////////
template <typename T> inline T sign(const T & x){ return (x>T(0)) - (x<T(0)); }

inline double dot2(const Eigen::RowVector3d & v)
{
  return v.dot(v);
}

inline double clamp(const double & x, const double & a, const double & b)
{
  return std::max(a, std::min(b, x));
}

// https://iquilezles.org/articles/triangledistance/
template <typename T1, typename T2, typename T3, typename Tp>
typename Tp::Scalar udTriangle(
    const T1 & v1, const T2 & v2, const T3 & v3, const Tp & p)
{
  using vec3 = Eigen::Matrix<typename Tp::Scalar, 3, 1>;
  // prepare data    
  vec3 v21 = v2 - v1; vec3 p1 = p - v1;
  vec3 v32 = v3 - v2; vec3 p2 = p - v2;
  vec3 v13 = v1 - v3; vec3 p3 = p - v3;
  vec3 nor =  v21.cross( v13 );

  return sqrt( // inside/outside test    
    (sign(v21.cross(nor).dot(p1)) + 
     sign(v32.cross(nor).dot(p2)) + 
     sign(v13.cross(nor).dot(p3))<2.0) 
     ?
     // 3 edges    
     std::min( std::min( 
     dot2(v21*clamp(v21.dot(p1)/dot2(v21),0.0,1.0)-p1), 
     dot2(v32*clamp(v32.dot(p2)/dot2(v32),0.0,1.0)-p2) ), 
     dot2(v13*clamp(v13.dot(p3)/dot2(v13),0.0,1.0)-p3) )
     :
     // 1 face    
     nor.dot(p1)*nor.dot(p1)/dot2(nor) );
}
    
template <typename Tp, typename Ta, typename Tb>
typename Tp::Scalar udLineSegment(const Tp & p, const Ta & a, const Tb & b)
{
  using vec3 = Eigen::Matrix<typename Tp::Scalar, 3, 1>;
  vec3 pa = p - a;
  vec3 ba = b - a;
  const auto h = clamp( pa.dot(ba)/ba.dot(ba), 0.0, 1.0 );
  return ( pa - ba*h ).norm();
}

int main(int argc, char * argv[])
{
  IGL_TICTOC_LAMBDA;
  // Octree depth
  int max_depth = 7;

  // Read mesh
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> V;
  Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> F;
  if(!igl::read_triangle_mesh
     (argc>1?argv[1]: TUTORIAL_SHARED_PATH "/decimated-knight.off",V,F)) {
    std::cout << "Failed to load mesh." << std::endl;
  }
  // Normalize
  V.rowwise() -= 0.5* (V.colwise().maxCoeff() + V.colwise().minCoeff());
  const double scale = 0.8/V.array().abs().maxCoeff();
  V *= scale;

  // Per-vertex radius for VARIABLE_RADIUS_OFFSET
  // if second arg exists read R from .dmat
  Eigen::VectorXd R;
  if(argc>2)
  {
    igl::readDMAT( argv[2],R);
    R *= scale;
  }else
  {
    R = V.col(1);
    R.array() -= R.minCoeff();
    R /= R.maxCoeff();
    R.array() += 1.0;
    R = (R.array().pow(4)).eval();
    R *= 0.01;
  }

  ////////////////////////////////////////////////////////////////////////
  /// Box up simplices (and, for VARIABLE_RADIUS_OFFSET, precompute data)
  ////////////////////////////////////////////////////////////////////////
  tictoc();
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> PB1,PB2;
  Eigen::Matrix<int,Eigen::Dynamic,2,Eigen::RowMajor> E;
  igl::edges(F,E);

  double isovalue = 0.02;
  std::function<double(const Eigen::RowVector3d &,const int i)> primitive;
  std::vector<igl::SphereMeshWedge<double>> data;
  switch(offset)
  {
    case OffsetType::EDGE_OFFSET:
    {
      igl::box_simplices(V,E,PB1,PB2);
      primitive = [&](const Eigen::RowVector3d & p, const int i)
      {
        return udLineSegment(p, V.row(E(i,0)), V.row(E(i,1)));
      };
      break;
    }
    case OffsetType::TRIANGLE_OFFSET:
    {
      igl::box_simplices(V,F,PB1,PB2);
      primitive = [&](const Eigen::RowVector3d & p, const int i)
      {
        return udTriangle( V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)), p);
      };
      break;
    }
    case OffsetType::VARIABLE_RADIUS_OFFSET:
    {
      isovalue = 0;
      igl::variable_radius_offset(V,F,R,PB1,PB2,data,primitive);
      break;
    }
  }
  printf("%-20s: %g secs\n","PB",tictoc());

  ////////////////////////////////////////////////////////////////////////
  /// Put boxes into eytzinger AABB tree
  ////////////////////////////////////////////////////////////////////////
  tictoc();
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> B1,B2;
  Eigen::VectorXi leaf;
  igl::eytzinger_aabb(PB1,PB2,B1,B2,leaf);
  printf("%-20s: %g secs\n","eytzinger_aabb",tictoc());
  //{
  //  Eigen::VectorXi U;
  //  igl::unique(leaf,U);
  //  printf("%d → %d\n",PB1.rows(),U.size()-2);
  //}

  igl::opengl::glfw::Viewer viewer;
  Eigen::Matrix<double,Eigen::Dynamic,3> mV,oV;
  Eigen::Matrix<int,Eigen::Dynamic,3> mF,oF;
  viewer.data().set_mesh(V,F);
  viewer.data().set_face_based(true);
  viewer.data().show_lines = false;

  if(offset == OffsetType::VARIABLE_RADIUS_OFFSET)
  {
    viewer.data().set_colors(R);
  }
  viewer.append_mesh();
  const int m_index = viewer.selected_data_index;
  viewer.append_mesh();
  const int o_index = viewer.selected_data_index;

  viewer.data(m_index).set_face_based(true);
  viewer.data(m_index).show_lines = true;
  viewer.data(o_index).set_face_based(true);
  viewer.data(o_index).show_lines = false;


  const auto update = [&]()
  {
    ////////////////////////////////////////////////////////////////////////
    /// If the grid isn't that big then compute a dense M.C. mesh
    /// Call batched eytzinger_aabb_sdf directly on grid nodes.
    ////////////////////////////////////////////////////////////////////////
    tictoc();
    const int ns = (1 << max_depth) + 1;
    if(ns <= 128+1)
    {
      Eigen::RowVector3i res(ns,ns,ns);
      Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> GV;
      igl::grid(res,GV);
      GV *= 2;
      GV.array() -= 1;
      Eigen::VectorXd S;
      igl::eytzinger_aabb_sdf(GV,primitive,B1,B2,leaf,S);
      printf("%-20s: %g secs\n","eytzinger_aabb_sdf",tictoc());

      tictoc();
      igl::marching_cubes(S,GV,res(0),res(1),res(2),isovalue,mV,mF);
      printf("%-20s: %g secs\n","marching_cubes",tictoc());
    }else
    {
      mV.resize(0,3); 
      mF.resize(0,3);
    }

    ////////////////////////////////////////////////////////////////////////
    /// Prepare an unsigned distance and signed distance function handle.
    ////////////////////////////////////////////////////////////////////////
    tictoc();
    const std::function<double(const Eigen::RowVector3d &)>
      sdf = [&](const Eigen::RowVector3d & p) -> double
    {
      const std::function<double(const int)> primitive_p = [&](const int j)
      {
        return primitive(p,j);
      };
      double f;
      igl::eytzinger_aabb_sdf(p,primitive_p,B1,B2,leaf,f);
      return f - isovalue;
    };
    const std::function<double(const Eigen::RowVector3d &)>
      udf = [&](const Eigen::RowVector3d & p) -> double
    {
      return std::abs(sdf(p));
    };
    ////////////////////////////////////////////////////////////////////////
    /// Use udf to build sparse voxel octree around the zero level set
    ////////////////////////////////////////////////////////////////////////
    Eigen::RowVector3d origin(-1,-1,-1);
    const double h0 = 2;
    Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> ijk;
    igl::lipschitz_octree( origin,h0,max_depth,udf,ijk);
    printf("%-20s: %g secs\n","lipschitz_octree",tictoc());

    ////////////////////////////////////////////////////////////////////////
    /// Gather the corners of those leaf cells, compute sdf there and run 
    /// (sparse) marching cubes
    ////////////////////////////////////////////////////////////////////////
    {
      tictoc();
      // Gather the corners of those leaf cells
      const double h = h0 / (1 << max_depth);
      Eigen::Matrix<int,Eigen::Dynamic,8,Eigen::RowMajor> J;
      Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> unique_ijk;
      Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> unique_corner_positions;
      igl::unique_sparse_voxel_corners(origin,h0,max_depth,ijk,unique_ijk,J,unique_corner_positions);
      //printf("unique_sparse_voxel_corners: %0.7f seconds\n",tictoc());
      printf("%-20s: %g secs\n","unique_sparse_vo...",tictoc());
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
        //printf("                        sdf: %0.7f seconds\n",tictoc());
        printf("%-20s: %g secs\n","sdf",tictoc());
      // Run marching cubes on the sparse set of leaf cells
      igl::marching_cubes( S,unique_corner_positions,J, 0, oV,oF);
      //printf("             marching_cubes: %0.7f seconds\n",tictoc());
      printf("%-20s: %g secs\n","marching_cubes",tictoc());
    }

    viewer.data(m_index).clear();
    viewer.data(m_index).set_mesh(mV,mF);
    viewer.data(m_index).set_colors(Eigen::RowVector3d(0.8,0.2,0.2));

    viewer.data(o_index).clear();
    viewer.data(o_index).set_mesh(oV,oF);
    viewer.data(o_index).set_colors(Eigen::RowVector3d(0.2,0.8,0.2));
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
    }
    return true;
  };
  update();

  viewer.launch();
  return 0;
}
