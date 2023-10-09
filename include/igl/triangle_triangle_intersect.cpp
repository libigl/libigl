#include "triangle_triangle_intersect.h"
#include "PI.h"
#include "tri_tri_intersect.h"
#include "ray_triangle_intersect.h"
#include "barycentric_coordinates.h"
#include "matlab_format.h"
#include <Eigen/Geometry>
#include <iostream>
#include <iomanip>
#include <stdio.h>
// std::signbit
#include <cmath>

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedE,
  typename DerivedEMAP,
  typename DerivedEF,
  typename DerivedEI,
  typename Derivedp>
IGL_INLINE bool igl::triangle_triangle_intersect(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<DerivedE> & E,
  const Eigen::MatrixBase<DerivedEMAP> & EMAP,
  const Eigen::MatrixBase<DerivedEI> & EF,
  const Eigen::MatrixBase<DerivedEF> & EI,
  const int f,
  const int c,
  const Eigen::MatrixBase<Derivedp> & p,
  const int g)
{
  // I'm leaving this debug printing stuff in for a bit until I trust this
  // better. 
  constexpr bool stinker = false;
  //const bool stinker = (f==1492 && g==1554);
  if(stinker) { printf("👀\n"); }
  bool found_intersection = false;
  // So edge opposite F(f,c) is the outer edge.
  const int o = EMAP(f + c*F.rows());
  // Do they share the edge opposite c?
  if((EF(o,0) == f && EF(o,1) == g) || (EF(o,1) == f && EF(o,0) == g))
  {
    if(stinker) { printf("⚠️ shares an edge\n"); }
    // Only intersects if the dihedral angle is zero (precondition: no zero
    // area triangles before or after collapse)
    const auto vg10 = (V.row(F(g,1))-V.row(F(g,0))).template head<3>();
    const auto vg20 = (V.row(F(g,2))-V.row(F(g,0))).template head<3>();
    const auto ng = vg10.cross(vg20);
    const int fo = EF(o,0) == f ? EI(o,0) : EI(o,1);
    const auto vf1p = (V.row(F(f,(fo+1)%3))-p).template head<3>();
    const auto vf2p = (V.row(F(f,(fo+2)%3))-p).template head<3>();
    const auto nf = vf1p.cross(vf2p);
    const auto o_vec_un = (V.row(E(o,1))-V.row(E(o,0))).template head<3>();
    const auto o_vec = o_vec_un.stableNormalized();

    const auto dihedral_angle = igl::PI - std::atan2(o_vec.dot(ng.cross(nf)),ng.dot(nf));
    if(dihedral_angle > 1e-8)
    {
      return false;
    }
    // Triangles really really might intersect.
    found_intersection = true;
  }else
  {
    if(stinker) { printf("does not share an edge\n"); }
    // Do they share a vertex?
    int sf,sg;
    bool found_shared_vertex = false;
    for(sf = 0;sf<3;sf++)
    {
      if(sf == c){ continue;}
      for(sg = 0;sg<3;sg++)
      {
        if(F(f,sf) == F(g,sg))
        {
          found_shared_vertex = true;
          break;
        }
      }
      if(found_shared_vertex) { break;} 
    }
    if(found_shared_vertex)
    {
      if(stinker) { printf("⚠️ shares an vertex\n"); }
      // If they share a vertex and intersect, then an opposite edge must
      // stab through the other triangle.

      // intersect_triangle1 needs non-const inputs.
      Eigen::RowVector3d g0 = V.row(F(g,0)).template cast<double>();
      Eigen::RowVector3d g1 = V.row(F(g,1)).template cast<double>();
      Eigen::RowVector3d g2 = V.row(F(g,2)).template cast<double>();
      Eigen::RowVector3d fs;
      if(((sf+1)%3)==c)
      {
        fs = p;
      }else
      {
        fs = V.row(F(f,(sf+1)%3));
      }
      Eigen::RowVector3d fd;
      if( ((sf+2)%3)==c )
      {
        fd = p.template cast<double>();
      }else
      {
        fd = V.row(F(f,(sf+2)%3)).template cast<double>();
      }
      Eigen::RowVector3d fdir = fd - fs;
      double t,u,v;

      if(stinker)
      {
        std::cout<<"T = ["<<g0<<";" <<g1<<";"<<g2<<"];"<<std::endl;
        std::cout<<"src = [" <<fs<<"];"<<std::endl;
        std::cout<<"dir = [" <<fdir<<"];"<<std::endl;
      }
      // p = (1-u-v)*a + u*b + v*c
      const auto bary = [](
        const Eigen::RowVector3d & p,
        const Eigen::RowVector3d & a,
        const Eigen::RowVector3d & b,
        const Eigen::RowVector3d & c,
        double & u,
        double & v)
      {
        const auto v0 = (b-a).eval();
        const auto v1 = (c-a).eval();
        const auto v2 = (p-a).eval();
        const double d00 = v0.dot(v0);
        const double d01 = v0.dot(v1);
        const double d11 = v1.dot(v1);
        const double d20 = v2.dot(v0);
        const double d21 = v2.dot(v1);
        const double denom = d00 * d11 - d01 * d01;
        u = (d11 * d20 - d01 * d21) / denom;
        v = (d00 * d21 - d01 * d20) / denom;
        // Equivalent:
        //Eigen::RowVector3d l;
        //igl::barycentric_coordinates(p,a,b,c,l);
        //u = l(1); v = l(2);
      };

      // Does the segment (A,B) intersect the triangle (0,0),(1,0),(0,1)?
      const auto intersect_unit = [](
        const Eigen::RowVector2d & A,
        const Eigen::RowVector2d & B) -> bool
      {
        // Check if P is inside (0,0),(1,0),(0,1) triangle
        const auto inside_unit = []( const Eigen::RowVector2d & P) -> bool
        {
          return P(0) >= 0 && P(1) >= 0 && P(0) + P(1) <= 1;
        };
        if(inside_unit(A) || inside_unit(B)) { return true; }

        const auto open_interval_contains_zero = [](
          const double a, const double b) -> bool
        {
          // handle case where either is 0.0 or -0.0
          if(a==0 || b==0) { return false; }
          return std::signbit(a) != std::signbit(b);
        };

        // Now check if the segment intersects any of the edges.
        // Does A-B intesect X-axis?
        if(open_interval_contains_zero(A(1),B(1)))
        {
          assert((A(1) - B(1)) != 0);
          // A and B are on opposite sides of the X-axis
          const double t = A(1) / (A(1) - B(1));
          const double x = A(0) + t * (B(0) - A(0));
          if(x >= 0 && x <= 1)
          {
            return true;
          }
        }
        // Does A-B intesect Y-axis?
        if(open_interval_contains_zero(A(0),B(0)))
        {
          assert((A(0) - B(0)) != 0);
          // A and B are on opposite sides of the Y-axis
          const double t = A(0) / (A(0) - B(0));
          const double y = A(1) + t * (B(1) - A(1));
          if(y >= 0 && y <= 1)
          {
            return true;
          }
        }
        // Does A-B intersect the line x+y=1?
        if(open_interval_contains_zero(A(0) + A(1),B(0) + B(1)))
        {
          assert((A(0) + A(1) - (B(0) + B(1))) != 0);
          const double t = (A(0) + A(1)) / ((A(0) + A(1)) - (B(0) + B(1)));
          const double x = A(0) + t * (B(0) - A(0));
          if(x >= 0 && x <= 1)
          {
            return true;
          }
        }
        return false;
      };


      //if(intersect_triangle1(
      //      fs.data(),fdir.data(),
      //      g0.data(),g1.data(),g2.data(),
      //      &t,&u,&v))
      bool coplanar = false;
      double epsilon = 1e-14;
      if(ray_triangle_intersect(
        fs,fdir,
        g0,g1,g2,
        epsilon,
        t,u,v,coplanar))
      {
        found_intersection = t > 0 && t<1+epsilon;
      }else if(coplanar)
      {
        // deal with coplanar
        Eigen::RowVector2d s2,d2;
        bary(fs,g0,g1,g2,s2(0),s2(1));
        bary(fd,g0,g1,g2,d2(0),d2(1));
        found_intersection = intersect_unit(s2,d2);
      }

      if(!found_intersection)
      {
        Eigen::RowVector3d fv[3];
        fv[0] = V.row(F(f,0)).template cast<double>();
        fv[1] = V.row(F(f,1)).template cast<double>();
        fv[2] = V.row(F(f,2)).template cast<double>();
        fv[c] = p.template cast<double>();
        Eigen::RowVector3d gs = V.row(F(g,(sg+1)%3)).template cast<double>();
        Eigen::RowVector3d gd = V.row(F(g,(sg+2)%3)).template cast<double>();
        Eigen::RowVector3d gdir = gd - gs;
        if(ray_triangle_intersect(
              gs,gdir,
              fv[0],fv[1],fv[2],
              epsilon,
              t,u,v,coplanar))
        {
          found_intersection = t > 0 && t<1+epsilon;
        }else if(coplanar)
        {
          // deal with coplanar
          //assert(false);
          Eigen::RowVector2d s2,d2;
          bary(gs,fv[0],fv[1],fv[2],s2(0),s2(1));
          bary(gd,fv[0],fv[1],fv[2],d2(0),d2(1));
          found_intersection = intersect_unit(s2,d2);
        }
      }
    }else
    {
      bool coplanar;
      Eigen::RowVector3d i1,i2;
      found_intersection = igl::tri_tri_intersection_test_3d(
                V.row(F(g,0)).template cast<double>(), 
                V.row(F(g,1)).template cast<double>(), 
                V.row(F(g,2)).template cast<double>(),
                            p.template cast<double>(),
          V.row(F(f,(c+1)%3)).template cast<double>(),
          V.row(F(f,(c+2)%3)).template cast<double>(),
          coplanar,
          i1,i2);
      if(stinker) { printf("tri_tri_intersection_test_3d says %s\n",found_intersection?"☠️":"✅"); }
    }
  }
  if(stinker) { printf("%s\n",found_intersection?"☠️":"✅"); }
  return found_intersection;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template bool igl::triangle_triangle_intersect<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, int, int, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, int);
#endif
