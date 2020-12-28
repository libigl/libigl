// Something bad is happening when I try to make this a function. Maybe
// something is not inlining? It ends up 1.25× slower :-(
//
// Even if I make it a lambda with no arguments (all capture by reference [&])
// and call it immediately I get a 1.25× slow-down. 
// 
// Maybe keeping it out of a function allows the compiler to optimize with the
// loop? But then I guess that measn this function is not getting inlined? Or
// that it's not getting optimized after inlining?
//
// Minimal benchmarking code:
#ifdef false
#include <igl/grid.h>
#include <igl/marching_cubes.h>
#include <igl/get_seconds.h>
#include <Eigen/Core>
#include <cstdio>
int main(int argc, char * argv[])
{
  const auto & tictoc = []()
  {
    static double t_start = igl::get_seconds();
    double diff = igl::get_seconds()-t_start;
    t_start += diff;
    return diff;
  };
  tictoc();
  Eigen::MatrixXd GV;
  const int s = 256;
  igl::grid(Eigen::RowVector3i(s,s,s),GV);
  GV.array() *= 2.0;
  GV.array() -= 1.0;
  const auto f = [&](const double x, const double y, const double z)->double
  {
    const double R = sqrt(x*x+y*y+z*z);
    const double s = atan2(sqrt(x*x+y*y),z);
    const double p = atan2(y,x);
    return pow(sin(s),2.)*(pow(cos(12.*s),3.)*0.1+pow(sin(6.*p),2)*0.2)+(R-0.5);
  };
  Eigen::VectorXd S(GV.rows());
  for(int i = 0;i<GV.rows();i++) { S(i) = f(GV(i,0),GV(i,1),GV(i,2)); }
  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF;
  igl::marching_cubes(S,GV,s,s,s,0,SV,SF);
  tictoc();
  const int MAX_RUNS = 10;
  for(int r = 0;r<MAX_RUNS;r++) { igl::marching_cubes(S,GV,s,s,s,0,SV,SF); }
  printf("mc: %g secs\n",tictoc()/MAX_RUNS);
}
#endif

//march_cube(GV,cS,cI,isovalue,V,n,F,m,E2V);
//
//namespace igl{
//inline void march_cube(
//  const DerivedGV & GV,
//  const Eigen::Matrix<Scalar,8,1> & cS,
//  const Eigen::Matrix<Index,8,1> & cI,
//  const typename DerivedS::Scalar & isovalue,
//  Eigen::PlainObjectBase<DerivedV> &V,
//  Index & n,
//  Eigen::PlainObjectBase<DerivedF> &F,
//  Index & m,
//  std::unordered_map<int64_t,int> & E2V)
//{

// These consts get stored reasonably
#include "marching_cubes_tables.h"

  // Seems this is also successfully inlined
  const auto ij2vertex =
    [&E2V,&V,&n,&GV]
      (const Index & i, const Index & j, const Scalar & t)->Index
  {
    // Seems this is successfully inlined.
    const auto ij2key = [](int32_t i,int32_t j)
    {
      if(i>j){ std::swap(i,j); }
      std::int64_t ret = 0;
      ret |= i;
      ret |= static_cast<std::int64_t>(j) << 32;
      return ret;
    };
    const auto key = ij2key(i,j);
    const auto it = E2V.find(key);
    int v = -1;
    if(it == E2V.end())
    {
      // new vertex
      if(n==V.rows()){ V.conservativeResize(V.rows()*2+1,V.cols()); }
      V.row(n) = GV.row(i) + t*(GV.row(j) - GV.row(i));
      v = n;
      E2V[key] = v;
      n++;
    }else
    {
      v = it->second;
    }
    return v;
  };

    int c_flags = 0;
    for(int c = 0; c < 8; c++)
    {
      if(cS(c) > isovalue){ c_flags |= 1<<c; }
    }
    //Find which edges are intersected by the surface
    int e_flags = aiCubeEdgeFlags[c_flags];
    //If the cube is entirely inside or outside of the surface, then there will be no intersections
    if(e_flags == 0) { return; }
    //Find the point of intersection of the surface with each edge
    //Then find the normal to the surface at those points
    Eigen::Matrix<Index,12,1> edge_vertices;
    for(int e = 0; e < 12; e++)
    {
#ifndef NDEBUG
      edge_vertices[e] = -1;
#endif
      //if there is an intersection on this edge
      if(e_flags & (1<<e))
      {
        // find crossing point assuming linear interpolation along edges
        const Scalar & a = cS(a2eConnection[e][0]);
        const Scalar & b = cS(a2eConnection[e][1]);
        Scalar t;
        {
          const Scalar delta = b-a;
          if(delta == 0) { t = 0.5; }
          t = (isovalue - a)/delta;
        };
        // record global index into local table
        edge_vertices[e] = 
          ij2vertex(cI(a2eConnection[e][0]),cI(a2eConnection[e][1]),t);
        assert(edge_vertices[e] >= 0);
        assert(edge_vertices[e] < n);
      }
    }
    // Insert the triangles that were found.  There can be up to five per cube
    for(int f = 0; f < 5; f++)
    {
      if(a2fConnectionTable[c_flags][3*f] < 0) break;
      if(m==F.rows()){ F.conservativeResize(F.rows()*2+1,F.cols()); }
      assert(edge_vertices[a2fConnectionTable[c_flags][3*f+0]]>=0);
      assert(edge_vertices[a2fConnectionTable[c_flags][3*f+1]]>=0);
      assert(edge_vertices[a2fConnectionTable[c_flags][3*f+2]]>=0);
      F.row(m) <<
        edge_vertices[a2fConnectionTable[c_flags][3*f+0]],
        edge_vertices[a2fConnectionTable[c_flags][3*f+1]],
        edge_vertices[a2fConnectionTable[c_flags][3*f+2]];
      m++;
    }
//};
