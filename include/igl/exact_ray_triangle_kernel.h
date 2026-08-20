// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_EXACT_RAY_TRIANGLE_KERNEL_H
#define IGL_EXACT_RAY_TRIANGLE_KERNEL_H
//
// Internal (detail) header shared by exact_ray_triangle_intersect and
// exact_ray_triangle_hit_compare.
//
// The exact ray-triangle predicates only need the sign of a handful of 3x3
// determinants whose columns are either the raw floating point direction vector
// `dir` or an exact difference of two input points (e.g. B-A). Because the
// direction is a free vector, these determinants cannot be reduced to Shewchuk's
// point-based orient3d without first forming source+dir (which would round and
// perturb the exact ray). Instead we evaluate the exact sign using standard
// adaptive expansion arithmetic (Shewchuk, "Adaptive Precision Floating-Point
// Arithmetic and Fast Robust Geometric Predicates", 1997). A cheap floating
// point filter handles the common non-degenerate case; the expansion path is
// only taken when the filter cannot certify the sign.
//
// Everything here operates on plain `double`: the error-free transforms
// (two_product via fma, two_sum, ...) are only exact under IEEE-754 double
// semantics. Callers static_assert that inputs are float or double (both promote
// to double exactly), preserving exactness with respect to the given
// coordinates. Only signs are ever produced; no exact construction is exposed.
//
#include <vector>
#include <cmath>
#include <cfloat>

namespace igl
{
  namespace exact_ray_triangle_kernel
  {
    // A nonoverlapping expansion sorted in order of increasing magnitude
    // (Shewchuk).
    using Expansion = std::vector<double>;

    inline void two_sum(double a, double b, double &x, double &y)
    {
      x = a + b;
      const double bv = x - a;
      const double av = x - bv;
      const double br = b - bv;
      const double ar = a - av;
      y = ar + br;
    }

    inline void two_diff(double a, double b, double &x, double &y)
    {
      x = a - b;
      const double bv = a - x;
      const double av = x + bv;
      const double br = bv - b;
      const double ar = a - av;
      y = ar + br;
    }

    inline void two_product(double a, double b, double &x, double &y)
    {
      x = a * b;
      y = std::fma(a, b, -x);
    }

    // h = e + f, where e and f are each nonoverlapping and increasing in
    // magnitude.
    inline Expansion fast_expansion_sum(const Expansion & e, const Expansion & f)
    {
      if(e.empty()){ return f; }
      if(f.empty()){ return e; }
      Expansion h;
      h.reserve(e.size() + f.size());
      size_t ei = 0, fi = 0;
      double Q, hh;
      {
        const double en = e[ei], fn = f[fi];
        if((fn > en) == (fn > -en)){ Q = en; ei++; }else{ Q = fn; fi++; }
      }
      while(ei < e.size() && fi < f.size())
      {
        const double en = e[ei], fn = f[fi];
        double addend;
        if((fn > en) == (fn > -en)){ addend = en; ei++; }else{ addend = fn; fi++; }
        two_sum(Q, addend, Q, hh);
        if(hh != 0.0){ h.push_back(hh); }
      }
      while(ei < e.size()){ two_sum(Q, e[ei], Q, hh); ei++; if(hh != 0.0){ h.push_back(hh); } }
      while(fi < f.size()){ two_sum(Q, f[fi], Q, hh); fi++; if(hh != 0.0){ h.push_back(hh); } }
      if(Q != 0.0 || h.empty()){ h.push_back(Q); }
      return h;
    }

    // h = e * b for a single double b.
    inline Expansion scale_expansion(const Expansion & e, double b)
    {
      Expansion h;
      if(e.empty() || b == 0.0){ return h; }
      h.reserve(e.size() * 2);
      double Q, hh, p1, p0, sum;
      two_product(e[0], b, Q, hh);
      if(hh != 0.0){ h.push_back(hh); }
      for(size_t i = 1; i < e.size(); i++)
      {
        two_product(e[i], b, p1, p0);
        two_sum(Q, p0, sum, hh);
        if(hh != 0.0){ h.push_back(hh); }
        // Fast_Two_Sum(p1, sum): |p1| >= |sum|
        Q = p1 + sum;
        const double bv = Q - p1;
        hh = sum - bv;
        if(hh != 0.0){ h.push_back(hh); }
      }
      if(Q != 0.0 || h.empty()){ h.push_back(Q); }
      return h;
    }

    inline Expansion expansion_product(const Expansion & e, const Expansion & f)
    {
      Expansion h;
      for(size_t i = 0; i < f.size(); i++)
      {
        const Expansion t = scale_expansion(e, f[i]);
        h = fast_expansion_sum(h, t);
      }
      return h;
    }

    inline Expansion expansion_negate(const Expansion & e)
    {
      Expansion h(e);
      for(double & x : h){ x = -x; }
      return h;
    }

    inline Expansion expansion_diff(const Expansion & a, const Expansion & b)
    {
      return fast_expansion_sum(a, expansion_negate(b));
    }

    inline int expansion_sign(const Expansion & e)
    {
      if(e.empty()){ return 0; }
      const double b = e.back();
      return b > 0.0 ? 1 : (b < 0.0 ? -1 : 0);
    }

    // A column of the 3x3 matrix. Either the raw vector p (diff==false) or the
    // exact difference p-q (diff==true).
    struct Column
    {
      double p[3];
      double q[3];
      bool diff;
    };

    inline Expansion column_entry(const Column & c, int r)
    {
      Expansion e;
      if(!c.diff)
      {
        if(c.p[r] != 0.0){ e.push_back(c.p[r]); }
        return e;
      }
      double x, y;
      two_diff(c.p[r], c.q[r], x, y);
      if(y != 0.0){ e.push_back(y); }
      if(x != 0.0){ e.push_back(x); }
      return e;
    }

    // Exact determinant expansion of the matrix whose columns are c0,c1,c2 (rows
    // are the three coordinates).
    inline Expansion det3_expansion(
      const Column & c0, const Column & c1, const Column & c2)
    {
      Expansion E[3][3];
      for(int r = 0; r < 3; r++)
      {
        E[r][0] = column_entry(c0, r);
        E[r][1] = column_entry(c1, r);
        E[r][2] = column_entry(c2, r);
      }
      const Expansion m0 = expansion_diff(
        expansion_product(E[1][1], E[2][2]), expansion_product(E[1][2], E[2][1]));
      const Expansion m1 = expansion_diff(
        expansion_product(E[1][0], E[2][2]), expansion_product(E[1][2], E[2][0]));
      const Expansion m2 = expansion_diff(
        expansion_product(E[1][0], E[2][1]), expansion_product(E[1][1], E[2][0]));
      return fast_expansion_sum(
        expansion_diff(expansion_product(E[0][0], m0), expansion_product(E[0][1], m1)),
        expansion_product(E[0][2], m2));
    }

    // Sign of det([c0 c1 c2]) with a floating point filter and exact fallback.
    inline int det3_sign(const Column & c0, const Column & c1, const Column & c2)
    {
      double M[3][3];
      const Column * col[3] = {&c0, &c1, &c2};
      for(int c = 0; c < 3; c++)
      {
        for(int r = 0; r < 3; r++)
        {
          M[r][c] = col[c]->diff ? (col[c]->p[r] - col[c]->q[r]) : col[c]->p[r];
        }
      }
      const double det =
          M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
        - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
        + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
      const double a[3][3] = {
        {std::fabs(M[0][0]), std::fabs(M[0][1]), std::fabs(M[0][2])},
        {std::fabs(M[1][0]), std::fabs(M[1][1]), std::fabs(M[1][2])},
        {std::fabs(M[2][0]), std::fabs(M[2][1]), std::fabs(M[2][2])}};
      const double perm =
          a[0][0]*(a[1][1]*a[2][2] + a[1][2]*a[2][1])
        + a[0][1]*(a[1][0]*a[2][2] + a[1][2]*a[2][0])
        + a[0][2]*(a[1][0]*a[2][1] + a[1][1]*a[2][0]);
      // Generously conservative dynamic filter. When it cannot certify the sign
      // we fall back to the always-correct expansion evaluation, so the exact
      // answer is never compromised by the choice of constant.
      const double errbound = 256.0 * DBL_EPSILON * perm;
      if(det > errbound){ return 1; }
      if(det < -errbound){ return -1; }
      return expansion_sign(det3_expansion(c0, c1, c2));
    }

    inline Column raw(const double * v)
    {
      Column c;
      c.p[0] = v[0]; c.p[1] = v[1]; c.p[2] = v[2];
      c.diff = false;
      return c;
    }

    inline Column difference(const double * p, const double * q)
    {
      Column c;
      c.p[0] = p[0]; c.p[1] = p[1]; c.p[2] = p[2];
      c.q[0] = q[0]; c.q[1] = q[1]; c.q[2] = q[2];
      c.diff = true;
      return c;
    }

    // Exact boolean: does the ray {O + t*d : t >= 0} transversally intersect the
    // triangle (A,B,C)?  Returns false for the coplanar (parallel) degeneracy.
    inline bool intersects(
      const double * O, const double * d,
      const double * A, const double * B, const double * C)
    {
      const Column cd = raw(d);
      const Column eAB = difference(B, A);
      const Column eAC = difference(C, A);
      // Denominator D' = det(d, B-A, C-A).  Zero => ray parallel to triangle plane.
      const int Dp = det3_sign(cd, eAB, eAC);
      if(Dp == 0){ return false; }
      // Edge tests: g_i = det(d, Vi-O, Vj-O); for an interior/boundary crossing
      // each shares the sign of D' (or is zero).  Any strictly opposite sign is a
      // miss.
      const int g0 = det3_sign(cd, difference(A, O), difference(B, O));
      if(g0 == -Dp){ return false; }
      const int g1 = det3_sign(cd, difference(B, O), difference(C, O));
      if(g1 == -Dp){ return false; }
      const int g2 = det3_sign(cd, difference(C, O), difference(A, O));
      if(g2 == -Dp){ return false; }
      // t = -a/D' with a = det(B-A, C-A, O-A) = n.(O-A).  t >= 0 <=> a == 0 or
      // sign(a) == -sign(D').
      const int a = det3_sign(eAB, eAC, difference(O, A));
      if(!(a == 0 || a == -Dp)){ return false; }
      return true;
    }

    // Fill num/den expansions with t = num/den for a single triangle, normalized
    // so that den > 0.  (t = -a/D'.)
    inline void hit_fraction(
      const double * O, const double * d,
      const double * A, const double * B, const double * C,
      Expansion & num, Expansion & den)
    {
      const Column cd = raw(d);
      const Column eAB = difference(B, A);
      const Column eAC = difference(C, A);
      const Expansion a  = det3_expansion(eAB, eAC, difference(O, A));
      const Expansion Dp = det3_expansion(cd, eAB, eAC);
      num = expansion_negate(a);
      den = Dp;
      if(expansion_sign(den) < 0)
      {
        num = expansion_negate(num);
        den = expansion_negate(den);
      }
    }

    // Exact comparison of the hit distances t1,t2 of two triangles along a common
    // ray. Returns -1 if t1<t2, 0 if equal, +1 if t1>t2.
    inline int compare(
      const double * O, const double * d,
      const double * A1, const double * B1, const double * C1,
      const double * A2, const double * B2, const double * C2)
    {
      Expansion n1, e1, n2, e2;
      hit_fraction(O, d, A1, B1, C1, n1, e1);
      hit_fraction(O, d, A2, B2, C2, n2, e2);
      // t1 < t2  <=>  n1*e2 < n2*e1  (denominators e1,e2 > 0)
      const Expansion lhs = expansion_product(n1, e2);
      const Expansion rhs = expansion_product(n2, e1);
      return expansion_sign(expansion_diff(lhs, rhs));
    }
  }
}
#endif
