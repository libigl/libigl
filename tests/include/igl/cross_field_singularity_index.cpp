#include <test_common.h>
#include <igl/cross_field_singularity_index.h>
#include <igl/find_cross_field_singularities.h>
#include <igl/per_face_normals.h>
#include <igl/local_basis.h>
#include <igl/euler_characteristic.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_border_vertex.h>
#include <igl/PI.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cstdlib>

namespace
{
  // Longitude/latitude sphere: pole 0, nlat-1 rings of nlon vertices, pole last
  void uv_sphere(int nlat, int nlon, Eigen::MatrixXd & V, Eigen::MatrixXi & F)
  {
    V.resize(2 + (nlat-1)*nlon, 3);
    V.row(0) << 0,0,1;
    auto ring = [&](int i, int j){ return 1 + (i-1)*nlon + ((j % nlon) + nlon) % nlon; };
    for(int i = 1;i<nlat;i++)
    {
      const double th = igl::PI*i/nlat;
      for(int j = 0;j<nlon;j++)
      {
        const double ph = 2*igl::PI*j/nlon;
        V.row(ring(i,j)) << std::sin(th)*std::cos(ph), std::sin(th)*std::sin(ph), std::cos(th);
      }
    }
    const int last = V.rows()-1;
    V.row(last) << 0,0,-1;
    F.resize(2*nlon + 2*nlon*(nlat-2), 3);
    int f = 0;
    for(int j = 0;j<nlon;j++){ F.row(f++) << 0, ring(1,j), ring(1,j+1); }
    for(int i = 1;i<nlat-1;i++)
    {
      for(int j = 0;j<nlon;j++)
      {
        const int a = ring(i,j), b = ring(i,j+1), c = ring(i+1,j+1), d = ring(i+1,j);
        F.row(f++) << a,d,c;
        F.row(f++) << a,c,b;
      }
    }
    for(int j = 0;j<nlon;j++){ F.row(f++) << last, ring(nlat-1,j+1), ring(nlat-1,j); }
    assert(f == F.rows());
  }

  void torus(int nu, int nv, Eigen::MatrixXd & V, Eigen::MatrixXi & F)
  {
    const double R = 2.0, r = 0.7;
    V.resize(nu*nv,3);
    auto id = [&](int i, int j){ return ((i%nu+nu)%nu)*nv + ((j%nv+nv)%nv); };
    for(int i = 0;i<nu;i++)
    {
      const double u = 2*igl::PI*i/nu;
      for(int j = 0;j<nv;j++)
      {
        const double v = 2*igl::PI*j/nv;
        V.row(id(i,j)) << (R + r*std::cos(v))*std::cos(u), (R + r*std::cos(v))*std::sin(u), r*std::sin(v);
      }
    }
    F.resize(2*nu*nv,3);
    int f = 0;
    for(int i = 0;i<nu;i++)
    {
      for(int j = 0;j<nv;j++)
      {
        F.row(f++) << id(i,j), id(i+1,j), id(i+1,j+1);
        F.row(f++) << id(i,j), id(i+1,j+1), id(i,j+1);
      }
    }
    // Orient outward
    for(int k = 0;k<F.rows();k++)
    {
      const Eigen::RowVector3d c = (V.row(F(k,0))+V.row(F(k,1))+V.row(F(k,2)))/3.0;
      const Eigen::RowVector3d e1 = V.row(F(k,1))-V.row(F(k,0));
      const Eigen::RowVector3d e2 = V.row(F(k,2))-V.row(F(k,0));
      const Eigen::RowVector3d n = e1.cross(e2);
      Eigen::RowVector3d axis(c(0),c(1),0); axis.normalize();
      const Eigen::RowVector3d out = c - R*axis;
      if(n.dot(out) < 0){ std::swap(F(k,1),F(k,2)); }
    }
  }

  // Geometry and field of libigl PR #2546 / issue #2526: a bipyramid whose faces
  // carry the "latitude" direction, so each pole is a full-turn singularity.
  void bipyramid(int n, Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::MatrixXd & PD1, Eigen::MatrixXd & PD2, double h = 1.0)
  {
    V.resize(n+2,3); F.resize(2*n,3); PD1.resize(2*n,3); PD2.resize(2*n,3);
    V.row(0) << 0,0,h;
    V.row(1) << 0,0,-h;
    for(int i = 0;i<n;i++)
    {
      V.row(i+2) << std::cos(2*igl::PI*i/n), std::sin(2*igl::PI*i/n), 0;
      F.row(i) << 0, i+2, (i+1)%n+2;
      F.row(i+n) << 1, (i+1)%n+2, i+2;
    }
    for(int i = 0;i<n;i++)
    {
      PD1.row(i) = V.row(F(i,2)) - V.row(F(i,1));
      PD1.row(i+n) = PD1.row(i);
      PD2.row(i) = V.row(0) - (V.row(F(i,1)) + V.row(F(i,2)))/2;
      PD2.row(i+n) = V.row(1) - (V.row(F(i+n,1)) + V.row(F(i+n,2)))/2;
    }
    PD1.rowwise().normalize();
    PD2.rowwise().normalize();
  }

  void icosahedron(Eigen::MatrixXd & V, Eigen::MatrixXi & F)
  {
    const double p = (1.0+std::sqrt(5.0))/2.0;
    V.resize(12,3);
    V << -1, p, 0,  1, p, 0,  -1,-p, 0,  1,-p, 0,
          0,-1, p,  0, 1, p,   0,-1,-p,  0, 1,-p,
          p, 0,-1,  p, 0, 1,  -p, 0,-1, -p, 0, 1;
    V.rowwise().normalize();
    F.resize(20,3);
    F << 0,11,5, 0,5,1, 0,1,7, 0,7,10, 0,10,11,
         1,5,9, 5,11,4, 11,10,2, 10,7,6, 7,1,8,
         3,9,4, 3,4,2, 3,2,6, 3,6,8, 3,8,9,
         4,9,5, 2,4,11, 6,2,10, 8,6,7, 9,8,1;
    for(int k = 0;k<F.rows();k++)
    {
      const Eigen::RowVector3d c = (V.row(F(k,0))+V.row(F(k,1))+V.row(F(k,2)))/3.0;
      const Eigen::RowVector3d e1 = V.row(F(k,1))-V.row(F(k,0));
      const Eigen::RowVector3d e2 = V.row(F(k,2))-V.row(F(k,0));
      const Eigen::RowVector3d n = e1.cross(e2);
      if(n.dot(c) < 0){ std::swap(F(k,1),F(k,2)); }
    }
  }

  // Planar polar grid in z=0: center vertex 0, `rings` rings of m vertices, CCW faces (normal +z)
  void polar_disk(int m, int rings, Eigen::MatrixXd & V, Eigen::MatrixXi & F)
  {
    V.resize(1 + rings*m, 3);
    V.row(0) << 0,0,0;
    auto id = [&](int i, int j){ return 1 + (i-1)*m + ((j%m)+m)%m; };
    for(int i = 1;i<=rings;i++)
    {
      for(int j = 0;j<m;j++)
      {
        // slight twist per ring so no edge is radial
        const double a = 2*igl::PI*j/m + 0.05*i;
        V.row(id(i,j)) << i*std::cos(a)/rings, i*std::sin(a)/rings, 0;
      }
    }
    F.resize(m + 2*m*(rings-1), 3);
    int f = 0;
    for(int j = 0;j<m;j++){ F.row(f++) << 0, id(1,j), id(1,j+1); }
    for(int i = 1;i<rings;i++)
    {
      for(int j = 0;j<m;j++)
      {
        F.row(f++) << id(i,j), id(i+1,j), id(i+1,j+1);
        F.row(f++) << id(i,j), id(i+1,j+1), id(i,j+1);
      }
    }
    assert(f == F.rows());
  }

  // Deterministic pseudo-random tangent field, PD2 = n × PD1
  void random_tangent_field(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, Eigen::MatrixXd & PD1, Eigen::MatrixXd & PD2)
  {
    Eigen::MatrixXd N;
    igl::per_face_normals(V,F,N);
    std::srand(2526);
    Eigen::MatrixXd R = Eigen::MatrixXd::Random(F.rows(),3);
    PD1.resize(F.rows(),3);
    PD2.resize(F.rows(),3);
    for(int f = 0;f<F.rows();f++)
    {
      const Eigen::RowVector3d n = N.row(f);
      Eigen::RowVector3d d = R.row(f) - R.row(f).dot(n)*n;
      if(d.norm() < 1e-3){ d = Eigen::RowVector3d(n(1),-n(0),0); }
      d.normalize();
      PD1.row(f) = d;
      PD2.row(f) = n.cross(d);
    }
  }

  // Replace each face's representative by another of the N symmetric directions
  Eigen::MatrixXd resymmetrize(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, const Eigen::MatrixXd & PD1, int N)
  {
    Eigen::MatrixXd Nrm;
    igl::per_face_normals(V,F,Nrm);
    Eigen::MatrixXd out(PD1.rows(),3);
    for(int f = 0;f<F.rows();f++)
    {
      const double a = 2*igl::PI*((f*7)%N)/N;
      const Eigen::RowVector3d n = Nrm.row(f);
      const Eigen::RowVector3d d = PD1.row(f);
      out.row(f) = std::cos(a)*d + std::sin(a)*n.cross(d);
    }
    return out;
  }

  int mod4(int x){ return ((x%4)+4)%4; }
}

TEST_CASE("cross_field_singularity_index: bipyramid poles (issue #2526)", "[igl]")
{
  // Resolved regime: the field turns by the apex angle (< 45°) across each
  // radial edge, so each pole must be reported as a full turn (+4 quarter turns)
  for(int n = 6;n<=16;n++)
  {
    INFO("n = " << n);
    Eigen::MatrixXd V,PD1,PD2; Eigen::MatrixXi F;
    bipyramid(n,V,F,PD1,PD2);
    Eigen::VectorXi I;
    igl::cross_field_singularity_index(V,F,PD1,I);
    REQUIRE(I.size() == n+2);
    REQUIRE(I(0) == 4);
    REQUIRE(I(1) == 4);
    for(int v = 2;v<n+2;v++){ REQUIRE(I(v) == 0); }
    REQUIRE(I.sum() == 4*2);
    // The Z/4 residue reported by the existing function is consistent (and 0 here,
    // which is why these poles are invisible to it)
    Eigen::VectorXi S,SI;
    igl::find_cross_field_singularities(V,F,PD1,PD2,S,SI,false);
    for(int v = 0;v<n+2;v++){ REQUIRE(SI(v) == mod4(I(v))); }
  }
  // Under-resolved regime (apex angle > 45°): the field aliases and the per-vertex
  // values change, but Poincaré–Hopf still holds exactly
  for(int n = 3;n<=5;n++)
  {
    INFO("n = " << n);
    Eigen::MatrixXd V,PD1,PD2; Eigen::MatrixXi F;
    bipyramid(n,V,F,PD1,PD2);
    Eigen::VectorXi I;
    igl::cross_field_singularity_index(V,F,PD1,I);
    REQUIRE(I.sum() == 4*2);
  }
}

TEST_CASE("cross_field_singularity_index: residuals on the wrap boundary", "[igl]")
{
  // Hexagonal bipyramid with apex height 2^(-1/4): every radial edge has an apex
  // angle of exactly 45°, so each residual sits on the ±π/4 wrap boundary and
  // floating point decides its sign edge by edge. Whatever the per-vertex outcome,
  // the indices must be integers and Poincaré–Hopf must hold exactly, because
  // each residual enters its two endpoint stars with opposite signs.
  Eigen::MatrixXd V,PD1,PD2; Eigen::MatrixXi F;
  bipyramid(6,V,F,PD1,PD2,std::pow(2.0,-0.25));
  Eigen::VectorXi I;
  igl::cross_field_singularity_index(V,F,PD1,I);
  REQUIRE(I.sum() == 4*2);
  // Perturbing the field by a tiny rotation to either side of the boundary
  // resolves the ambiguity consistently: the poles are full turns on one side
  Eigen::MatrixXd Nrm;
  igl::per_face_normals(V,F,Nrm);
  for(const double eps : {-1e-9, 1e-9})
  {
    // rotate every face's direction by eps·(face parity) so radial residuals move off ±π/4
    Eigen::MatrixXd P(F.rows(),3);
    for(int f = 0;f<F.rows();f++)
    {
      const double a = eps*((f%6)%2 == 0 ? 1.0 : -1.0);
      const Eigen::RowVector3d n = Nrm.row(f), d = PD1.row(f);
      P.row(f) = std::cos(a)*d + std::sin(a)*n.cross(d);
    }
    Eigen::VectorXi J;
    igl::cross_field_singularity_index(V,F,P,J);
    REQUIRE(J.sum() == 4*2);
  }
}

TEST_CASE("cross_field_singularity_index: icosahedron local basis field", "[igl]")
{
  // Case from libigl-python-bindings test_find_cross_field_singularities, where
  // the Poincaré–Hopf assertion had to be commented out: the existing function
  // returns residues summing to 16 because it stores -1 as 3.
  Eigen::MatrixXd V; Eigen::MatrixXi F;
  icosahedron(V,F);
  Eigen::MatrixXd B1,B2,B3;
  igl::local_basis(V,F,B1,B2,B3);
  Eigen::VectorXi I;
  igl::cross_field_singularity_index(V,F,B1,I);
  REQUIRE(I.sum() == 4*2);
  REQUIRE(I.cwiseAbs().maxCoeff() == 1);
  REQUIRE(I.cwiseAbs().sum() == 12);
  Eigen::VectorXi S,SI;
  igl::find_cross_field_singularities(V,F,B1,B2,S,SI,false);
  for(int v = 0;v<V.rows();v++){ REQUIRE(SI(v) == mod4(I(v))); }
  REQUIRE(SI.sum() == 16);
}

TEST_CASE("cross_field_singularity_index: Poincaré–Hopf on random fields", "[igl]")
{
  struct Case { std::string name; Eigen::MatrixXd V; Eigen::MatrixXi F; };
  std::vector<Case> cases;
  { Case c; c.name = "uv_sphere"; uv_sphere(14,28,c.V,c.F); cases.push_back(c); }
  { Case c; c.name = "torus"; torus(20,12,c.V,c.F); cases.push_back(c); }
  for(const auto & name : test_common::closed_manifold_meshes())
  {
    Case c; c.name = name;
    igl::read_triangle_mesh(test_common::data_path(name),c.V,c.F);
    cases.push_back(c);
  }
  for(const auto & c : cases)
  {
    INFO("mesh = " << c.name);
    REQUIRE(igl::is_edge_manifold(c.F));
    const std::vector<bool> border = igl::is_border_vertex(c.F);
    REQUIRE(std::none_of(border.begin(),border.end(),[](bool b){ return b; }));
    const int chi = igl::euler_characteristic(c.F);
    Eigen::MatrixXd PD1,PD2;
    random_tangent_field(c.V,c.F,PD1,PD2);
    for(const int N : {1,2,4,6})
    {
      INFO("N = " << N);
      Eigen::VectorXi I;
      igl::cross_field_singularity_index(c.V,c.F,PD1,N,I);
      REQUIRE(I.sum() == N*chi);
      // Independent of which of the N symmetric directions represents each face
      Eigen::VectorXi J;
      igl::cross_field_singularity_index(c.V,c.F,resymmetrize(c.V,c.F,PD1,N),N,J);
      REQUIRE(I == J);
    }
    // The existing function's output is exactly the Z/4 residue of the true index
    Eigen::VectorXi I,S,SI;
    igl::cross_field_singularity_index(c.V,c.F,PD1,I);
    igl::find_cross_field_singularities(c.V,c.F,PD1,PD2,S,SI,false);
    for(int v = 0;v<c.V.rows();v++)
    {
      REQUIRE(SI(v) == mod4(I(v)));
      REQUIRE(S(v) == (mod4(I(v)) != 0));
    }
    // Default overload is N = 4
    Eigen::VectorXi I4;
    igl::cross_field_singularity_index(c.V,c.F,PD1,4,I4);
    REQUIRE(I == I4);
  }
}

TEST_CASE("cross_field_singularity_index: prescribed winding on a planar disk", "[igl]")
{
  // Field angle (k/N)·arg(z) has an isolated singularity of index k/N at the
  // origin and is smooth elsewhere; the mesh has valence 24 at the center so
  // every |k| <= 5 is resolved. Checks sign, magnitude, N = 1, 2, 4, and that
  // boundary vertices report 0.
  Eigen::MatrixXd V; Eigen::MatrixXi F;
  polar_disk(24,4,V,F);
  const std::vector<bool> border = igl::is_border_vertex(F);
  for(const int N : {1,2,4})
  {
    for(int k = -5;k<=5;k++)
    {
      INFO("N = " << N << ", k = " << k);
      Eigen::MatrixXd PD1(F.rows(),3), PD2(F.rows(),3);
      for(int f = 0;f<F.rows();f++)
      {
        const Eigen::RowVector3d c = (V.row(F(f,0))+V.row(F(f,1))+V.row(F(f,2)))/3.0;
        const double a = double(k)/N*std::atan2(c(1),c(0));
        PD1.row(f) << std::cos(a), std::sin(a), 0;
        PD2.row(f) << -std::sin(a), std::cos(a), 0;
      }
      Eigen::VectorXi I;
      igl::cross_field_singularity_index(V,F,PD1,N,I);
      REQUIRE(I(0) == k);
      for(int v = 1;v<V.rows();v++){ REQUIRE(I(v) == 0); }
      if(N == 4)
      {
        Eigen::VectorXi S,SI;
        igl::find_cross_field_singularities(V,F,PD1,PD2,S,SI,false);
        REQUIRE(SI(0) == mod4(k));
      }
    }
  }
  // Boundary vertices are always 0, even for a field that is singular there
  REQUIRE(std::count(border.begin(),border.end(),true) == 24);
}

TEST_CASE("cross_field_singularity_index: accepts unnormalized and non-tangent input", "[igl]")
{
  Eigen::MatrixXd V,PD1,PD2; Eigen::MatrixXi F;
  bipyramid(10,V,F,PD1,PD2);
  Eigen::MatrixXd Nrm;
  igl::per_face_normals(V,F,Nrm);
  Eigen::MatrixXd skew = 3.0*PD1 + 0.5*Nrm;
  Eigen::VectorXi I,J;
  igl::cross_field_singularity_index(V,F,PD1,I);
  igl::cross_field_singularity_index(V,F,skew,J);
  REQUIRE(I == J);
  // Matrix (rather than vector) output type
  Eigen::MatrixXi M;
  igl::cross_field_singularity_index(V,F,PD1,M);
  REQUIRE(M.rows() == V.rows());
  REQUIRE(M.cols() == 1);
  REQUIRE(M.col(0) == I);
}
