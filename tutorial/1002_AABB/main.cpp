#include <igl/read_triangle_mesh.h>
#include <igl/box_simplices.h>
#include <igl/eytzinger_aabb.h>
#include <igl/lipschitz_octree.h>
#include <igl/unique_sparse_voxel_corners.h>
#include <igl/eytzinger_aabb_sdf.h>
#include <igl/unique.h>
#include <igl/grid.h>
#include <igl/edges.h>
#include <igl/marching_cubes.h>
#include <igl/parallel_for.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/matlab_format.h>
#include <igl/get_seconds.h>

#include <Eigen/QR>

#include <limits>

#include <Eigen/Geometry>
//template <typename T>
//inline T sign(const T & x)
//{
//  return (x > T(0)) - (x < T(0));
//}
#define sign(x) ((x>0) - (x<0))


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
    
template <typename Deriveda, typename Derivedba>
inline double sdRoundCone(
  const Eigen::RowVector3d & p, 
  const Eigen::MatrixBase<Deriveda> & a, 
  const double & r1, 
  const double & r2,
  const Eigen::MatrixBase<Derivedba> & ba,
  const double & l2,
  const double & rr,
  const double & a2,
  const double & il2)
{
  // sampling dependant computations
  Eigen::RowVector3d pa = p - a;
  double y = pa.dot(ba);
  double z = y - l2;
  double x2 = ( pa*l2 - ba*y ).squaredNorm();
  double y2 = y*y*l2;
  double z2 = z*z*l2;

  // single square root!
  double k = sign(rr)*rr*rr*x2;
  if( sign(z)*a2*z2>k ) return  sqrt(x2 + z2)        *il2 - r2;
  if( sign(y)*a2*y2<k ) return  sqrt(x2 + y2)        *il2 - r1;
                        return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
}

template <typename Deriveda, typename Derivedb>
inline double sdRoundCone(
  const Eigen::RowVector3d & p, 
  const Eigen::MatrixBase<Deriveda> & a, 
  const Eigen::MatrixBase<Derivedb> & b, 
  const double & r1, 
  const double & r2)
{
  // sampling independent computations (only depend on shape)
  Eigen::RowVector3d ba = b - a;
  double l2 = ba.dot(ba);
  double rr = r1 - r2;
  double a2 = l2 - rr*rr;
  double il2 = 1.0/l2;
  return sdRoundCone(p,a,r1,r2,ba,l2,rr,a2,il2);
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


inline double sdPlane(
  const Eigen::RowVector3d & P,
  const Eigen::RowVector3d & A,
  const Eigen::RowVector3d & B,
  const Eigen::RowVector3d & C)
{
  return (P-A).dot( ((B-A).cross(C-A)).normalized() );
}

// This has the correct zero-level set and 1-Lipschitz continuity but 
// incorrect value away from the zero-level set. When combined with the
// cones below I _think_ ∇f is still correct when f=0.
//
// If it was important to be correct everywhere then a straightforward way
// to do it would be to compute distances to triangles and determine
// insidnesss based on each plane.
inline double sdSkewedExtrudedTriangle(
  const Eigen::RowVector3d & P,
  const Eigen::Matrix<double,3,3,Eigen::RowMajor> & V,
  const Eigen::Matrix<double,3,3,Eigen::RowMajor> & T)
{
  return  std::max(sdPlane(P,V.row(0),V.row(1),V.row(2)),
          std::max(sdPlane(P,T.row(2),T.row(1),T.row(0)),
          std::max(sdPlane(P,V.row(1),V.row(0),T.row(0)),
          std::max(sdPlane(P,V.row(2),V.row(1),T.row(1)),
                   sdPlane(P,V.row(0),V.row(2),T.row(2))))));
};

class sdRoundTriangle
{
  public: 
  // Fields
  enum 
  {
    BIG_VERTEX = 0,
    BIG_EDGE = 1,
    NO_TRIANGLE = 2,
    FULL = 3
  } flavor;
  Eigen::Matrix<double,3,3,Eigen::RowMajor> V;
  Eigen::Matrix<double,3,1> r;
  Eigen::Matrix<double,3,3,Eigen::RowMajor> EV;
  Eigen::Matrix<double,3,1> l,l2,rr,a2,il2;
  int max_i;
  Eigen::Matrix<double,5,4,Eigen::RowMajor> planes;

  sdRoundTriangle(){}

  sdRoundTriangle(
    const Eigen::RowVector3d & V0,
    const Eigen::RowVector3d & V1,
    const Eigen::RowVector3d & V2,
    const double r0,
    const double r1,
    const double r2)
  {
    // Internal copy
    V.row(0) = V0;
    V.row(1) = V1;
    V.row(2) = V2;
    r(0) = r0;
    r(1) = r1;
    r(2) = r2;

    flavor = FULL;
    // By default use full

    EV.row(0) = V.row(2) - V.row(1);
    EV.row(1) = V.row(0) - V.row(2);
    EV.row(2) = V.row(1) - V.row(0);
    l = EV.rowwise().norm();
    l2 = l.array().square();
    rr << r(1) - r(2), r(2) - r(0), r(0) - r(1);
    a2 = l2.array() - rr.array().square();
    il2 = 1.0/l2.array();
  
    /////////////////////////////////////////////
    /// BIG_VERTEX ?
    /////////////////////////////////////////////
    {
      r.maxCoeff(&max_i);
      int j = (max_i+1)%3;
      int k = (max_i+2)%3;
      if((l(k) + r(j) < r(max_i)) && (l(j) + r(k) < r(max_i)))
      {
        flavor = BIG_VERTEX;
      }
    }

    /////////////////////////////////////////////
    /// BIG_EDGE ?
    /////////////////////////////////////////////
    if(flavor == FULL)
    {
      // Case where one edge's roundCone containes the others
      for(int e = 0;e<3;e++)
      {
        const int i = (e+1)%3;
        const int j = (e+2)%3;
        const int k = (e+3)%3;
        const double s = ::sdRoundCone(V.row(k),V.row(i),V.row(j),r(i),r(j));
        if(-s > r(k))
        {
          flavor = BIG_EDGE;
          max_i = i;
          break;
        }
      }
    }

    if(flavor == FULL && !compute_planes())
    {
      flavor = NO_TRIANGLE;
    }
  }

  bool compute_planes()
  {
    // Non-degenerate case
    const Eigen::RowVector3d & a = V.row(0);
    const Eigen::RowVector3d & b = V.row(1);
    const Eigen::RowVector3d & c = V.row(2);
    const double & ra = r(0);
    const double & rb = r(1);
    const double & rc = r(2);
    Eigen::Matrix<double,2,3,Eigen::RowMajor> A;
    A<<
      b-a,
      c-a;
    const Eigen::Vector2d d(rb-ra,rc-ra);
    const Eigen::RowVector3d N = (A.row(0).cross(A.row(1))).normalized();
    //const Eigen::CompleteOrthogonalDecomposition<decltype(A)> cod(A);
    const Eigen::RowVector3d n0 = A.completeOrthogonalDecomposition().solve(d);
    const double qA = N.squaredNorm();
#warning "qB is 0 by construction"
    const double qB = 2 * N.dot(n0);
    const double qC = n0.squaredNorm() - 1;
    const double qD = qB*qB - 4*qA*qC;
    double t_sol_1;
    double t_sol_2;
    Eigen::RowVector3d                        n1;
    Eigen::RowVector3d                        n2;

    if(qD<0) { return false; }

    t_sol_1 = (-qB + std::sqrt(qD)) / (2*qA);
    n1 = -(t_sol_1 * N + n0);
    const Eigen::Matrix<double,3,3,Eigen::RowMajor> T  = V + r * n1;

    const auto plane_equation = [](
        const Eigen::RowVector3d & a,
        const Eigen::RowVector3d & b,
        const Eigen::RowVector3d & c)->Eigen::RowVector4d
    {
      Eigen::RowVector3d n = (b-a).cross(c-a).normalized();
      n.normalize();
      double d = -n.dot(a);
      return Eigen::RowVector4d(n(0),n(1),n(2),d);
    };

    planes.row(0) = plane_equation(V.row(0),V.row(1),V.row(2));
    planes.row(1) = plane_equation(T.row(2),T.row(1),T.row(0));
    planes.row(2) = plane_equation(V.row(1),V.row(0),T.row(0));
    planes.row(3) = plane_equation(V.row(2),V.row(1),T.row(1));
    planes.row(4) = plane_equation(V.row(0),V.row(2),T.row(2));
    //std::cout<<igl::matlab_format(Eigen::MatrixXd(planes),"planes")<<std::endl;
    //exit(1);

    //t_sol_2 = (-qB - std::sqrt(qD)) / (2*qA);
    //n2 = -(t_sol_2 * N + n0);
    //B  = V + r * n2;
    return true;
  }

  double sdRoundCone(const Eigen::RowVector3d & p, const int i, const int j) const
  {
    const int e = (j+1)%3;
    return ::sdRoundCone(
      p,
      V.row(i),
      r(i),
      r(j),
      EV.row(e),
       l2(e),
       rr(e),
       a2(e),
      il2(e)
      );
  }

  double operator()(const Eigen::RowVector3d & p) const
  {
    if(flavor == BIG_VERTEX)
    {
      // Case 0: Vertex i
      return (p - V.row(max_i)).norm() - r(max_i);
    }

    if(flavor == BIG_EDGE)
    {
      const int i = max_i;
      const int j = (i+1)%3;
      // Case 1: Edge e
      return sdRoundCone(p,i,j);
    }

    double s = std::numeric_limits<double>::infinity();
    if(flavor == FULL)
    {
      // This is possibly the bottleneck and could be turned into precomputed
      // plane equations.
      
      // signed distance to triangle plane (this is immediately recomputed later in
      // sdSkewedExtrudedTriangle...)
      const auto plane_sdf = [](
        const Eigen::RowVector3d & p,
        const Eigen::RowVector4d & plane)
      {
        return plane.head<3>().dot(p) + plane(3);
      };
      double d0 = plane_sdf(p, planes.row(0));
      double planes_s = -std::abs(d0);
      // Reflect if necessary so that q is always on negative side of plane
      Eigen::RowVector3d q = p - (d0 - planes_s) * planes.row(0).head<3>();
      // Other planes (for negative side slab, by symmetry)
      for(int i = 1;i<planes.rows();i++)
      {
        planes_s = std::max(planes_s,plane_sdf(q, planes.row(i)));
      }
      s = std::min(s,planes_s);

      //s = std::min(s,sdSkewedExtrudedTriangle(q,V,T));
      //s = std::min(s,sdSkewedExtrudedTriangle(p,B,V));
    }
    assert(flavor == FULL || flavor == NO_TRIANGLE);

    for(int e = 0;e<3;e++)
    {
      const int i = (e+1)%3;
      const int j = (e+2)%3;
      s = std::min(s,sdRoundCone(p,i,j));
    }
    return s;
  }
};

int main(int argc, char * argv[])
{
  IGL_TICTOC_LAMBDA;
  // Read mesh
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> V;
  Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> F;
  if(!igl::read_triangle_mesh
     (argc>1?argv[1]: TUTORIAL_SHARED_PATH "/decimated-knight.off",V,F)) {
    std::cout << "Failed to load mesh." << std::endl;
  }

  V.rowwise() -= 0.5* (V.colwise().maxCoeff() + V.colwise().minCoeff());
  V /= V.array().abs().maxCoeff();
  V *= 0.8;


  const int ns = 512+1;
  Eigen::RowVector3i res(ns,ns,ns);

  //Eigen::VectorXd R = (0.1+0.9*(0.5+(Eigen::VectorXd::Random(V.rows())*0.5).array()))*((2.0/(res(0)-1))*2);
  Eigen::VectorXd R = V.col(1);
  R.array() -= R.minCoeff();
  R /= R.maxCoeff();
  R.array() += 1.0;
  R = (R.array().pow(4)).eval();
  R *= 0.01;


  tictoc();
  Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> PB1,PB2;
  Eigen::Matrix<int,Eigen::Dynamic,2,Eigen::RowMajor> E;
  igl::edges(F,E);

  enum class OffsetType
  {
    EDGE_OFFSET = 0,
    TRIANGLE_OFFSET = 1,
    VARIABLE_RADIUS_OFFSET = 2
  } offset = OffsetType::VARIABLE_RADIUS_OFFSET;

  double isovalue = 0.16;(2.0/(res(0)-1))*2;
  switch(offset)
  {
    case OffsetType::EDGE_OFFSET:
    {
      igl::box_simplices(V,E,PB1,PB2);
      break;
    }
    case OffsetType::TRIANGLE_OFFSET:
    {
      igl::box_simplices(V,F,PB1,PB2);
      break;
    }
    case OffsetType::VARIABLE_RADIUS_OFFSET:
    {
      isovalue = 0;
      PB1.setConstant(F.rows(),3,std::numeric_limits<double>::infinity());
      PB2.setConstant(F.rows(),3,-std::numeric_limits<double>::infinity());
      for(int f = 0;f<F.rows();f++)
      {
        for(int c = 0;c<F.cols();c++)
        {
          for(int d = 0;d<V.cols();d++)
          {
            PB1(f,d) = std::min(PB1(f,d),V(F(f,c),d)-R(F(f,c)));
            PB2(f,d) = std::max(PB2(f,d),V(F(f,c),d)+R(F(f,c)));
          }
        }
      }
      break;
    }
  }
  printf("%-20s: %g secs\n","PB",tictoc());






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


  tictoc();
  std::vector<sdRoundTriangle> data(F.rows());
  for(int f = 0;f<F.rows();f++)
  {
    data[f] = sdRoundTriangle(
      V.row(F(f,0)),V.row(F(f,1)),V.row(F(f,2)),
      R(F(f,0)),R(F(f,1)),R(F(f,2)));
  }

  std::function<double(const Eigen::RowVector3d &,const int i)> primitive;
  switch(offset)
  {
    case OffsetType::EDGE_OFFSET:
    {
      primitive = [&](const Eigen::RowVector3d & p, const int i)
      {
        return udLineSegment(p, V.row(E(i,0)), V.row(E(i,1)));
      };
      break;
    }
    case OffsetType::TRIANGLE_OFFSET:
    {
      primitive = [&](const Eigen::RowVector3d & p, const int i)
      {
        return udTriangle( V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)), p);
      };
      break;
    }
    case OffsetType::VARIABLE_RADIUS_OFFSET:
    {
      primitive = [&](const Eigen::RowVector3d & p, const int i)
      {
        return data[i](p);
      };
      break;
    }
  }


  Eigen::Matrix<double,Eigen::Dynamic,3> mV;
  Eigen::Matrix<int,Eigen::Dynamic,3> mF;
  if(ns <= 64+1)
  {
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
  }

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
    return f;
  };
  const std::function<double(const Eigen::RowVector3d &)>
    udf = [&](const Eigen::RowVector3d & p) -> double
  {
    return std::abs(sdf(p));
    //// This made performance worse.
    //const std::function<double(const int)> primitive_p = [&](const int j)
    //{
    //  const double d = udTriangle( V.row(F(j,0)), V.row(F(j,1)), V.row(F(j,2)), p);
    //  const Eigen::RowVector3d r(R(F(j,0)),R(F(j,1)),R(F(j,2)));
    //  const double min_r = r.minCoeff();
    //  const double max_r = r.maxCoeff();
    //  if(d > max_r)
    //  {
    //    return d - max_r;
    //  }else if(d < min_r)
    //  {
    //    return d - min_r;
    //  }
    //  return 0.0;
    //};
    //double f;
    //igl::eytzinger_aabb_sdf(p,primitive_p,B1,B2,leaf,f);
    //return f;
  };
  Eigen::RowVector3d origin(-1,-1,-1);
  const double h0 = 2;
  const int max_depth = floor(log2(res(0)));
  Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> ijk;
  igl::lipschitz_octree( origin,h0,max_depth,udf,ijk);
  printf("%-20s: %g secs\n","lipschitz_octree",tictoc());

  Eigen::Matrix<double,Eigen::Dynamic,3> oV;
  Eigen::Matrix<int,Eigen::Dynamic,3> oF;
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


  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V,F);
  viewer.data().set_face_based(true);
  viewer.data().show_lines = false;

  viewer.append_mesh();
  viewer.data().set_mesh(mV,mF);
  viewer.data().set_face_based(true);
  viewer.data().show_lines = true;
  viewer.data().show_faces = true;
  viewer.data().set_colors(Eigen::RowVector3d(0.8,0.2,0.2));

  viewer.append_mesh();
  viewer.data().set_mesh(oV,oF);
  viewer.data().set_face_based(true);
  viewer.data().show_lines = true;
  viewer.data().show_faces = true;
  viewer.data().set_colors(Eigen::RowVector3d(0.2,0.8,0.2));

  viewer.launch();
  return 0;
}
