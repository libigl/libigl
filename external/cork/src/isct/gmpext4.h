// +-------------------------------------------------------------------------
// | gmpext4.h
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#pragma once

#ifdef _WIN32
#pragma warning(disable: 4800)
#pragma warning(disable: 4244)
#include <mpirxx.h>
#pragma warning(default: 4244)
#pragma warning(default: 4800)
#else
#include <gmpxx.h>
#endif

#include <iostream>


namespace GMPExt4 {

// types for k-vectors in R4:
//      Ext4_k

struct GmpExt4_1 {
    mpz_class e0;
    mpz_class e1;
    mpz_class e2;
    mpz_class e3;
};

struct GmpExt4_2 {
    mpz_class e01;
    mpz_class e02;
    mpz_class e03;
    mpz_class e12;
    mpz_class e13;
    mpz_class e23;
};

struct GmpExt4_3 {
    mpz_class e012;
    mpz_class e013;
    mpz_class e023;
    mpz_class e123;
};


// ********************************
// Output Routines
inline std::ostream& operator<<(std::ostream &out, const GmpExt4_1 &ext)
{
    return out << '[' << ext.e0
               << ',' << ext.e1
               << ',' << ext.e2
               << ',' << ext.e3 << ']';
}
inline std::ostream& operator<<(std::ostream &out, const GmpExt4_2 &ext)
{
    return out << '[' << ext.e01
               << ',' << ext.e02
               << ',' << ext.e03
               << ',' << ext.e12
               << ',' << ext.e13
               << ',' << ext.e23 << ']';
}
inline std::ostream& operator<<(std::ostream &out, const GmpExt4_3 &ext)
{
    return out << '[' << ext.e012
               << ',' << ext.e013
               << ',' << ext.e023
               << ',' << ext.e123 << ']';
}


// ********************************
// A neg takes a k-vector and returns its negation
// neg(X,Y) is safe for X=Y
inline
void neg(GmpExt4_1 &out, const GmpExt4_1 &in)
{
    out.e0 = -in.e0;
    out.e1 = -in.e1;
    out.e2 = -in.e2;
    out.e3 = -in.e3;
}
inline
void neg(GmpExt4_2 &out, const GmpExt4_2 &in)
{
    out.e01 = -in.e01;
    out.e02 = -in.e02;
    out.e03 = -in.e03;
    out.e12 = -in.e12;
    out.e13 = -in.e13;
    out.e23 = -in.e23;
}
inline
void neg(GmpExt4_3 &out, const GmpExt4_3 &in)
{
    out.e012 = -in.e012;
    out.e013 = -in.e013;
    out.e023 = -in.e023;
    out.e123 = -in.e123;
}


// ********************************
// A dual operation takes a k-vector and returns a (4-k)-vector
// A reverse dual operation inverts the dual operation
// dual(X,Y) is not safe for X=Y (same with revdual)
inline
void dual(GmpExt4_1 &out, const GmpExt4_3 &in)
{
    out.e0 =  in.e123;
    out.e1 = -in.e023;
    out.e2 =  in.e013;
    out.e3 = -in.e012;
}
inline
void dual(GmpExt4_2 &out, const GmpExt4_2 &in)
{
    out.e01 =  in.e23;
    out.e02 = -in.e13;
    out.e03 =  in.e12;
    out.e12 =  in.e03;
    out.e13 = -in.e02;
    out.e23 =  in.e01;
}
inline
void dual(GmpExt4_3 &out, const GmpExt4_1 &in)
{
    out.e012 =  in.e3;
    out.e013 = -in.e2;
    out.e023 =  in.e1;
    out.e123 = -in.e0;
}
inline
void revdual(GmpExt4_1 &out, const GmpExt4_3 &in)
{
    out.e0 = -in.e123;
    out.e1 =  in.e023;
    out.e2 = -in.e013;
    out.e3 =  in.e012;
}
inline
void revdual(GmpExt4_2 &out, const GmpExt4_2 &in)
{
    out.e01 =  in.e23;
    out.e02 = -in.e13;
    out.e03 =  in.e12;
    out.e12 =  in.e03;
    out.e13 = -in.e02;
    out.e23 =  in.e01;
}
inline
void revdual(GmpExt4_3 &out, const GmpExt4_1 &in)
{
    out.e012 = -in.e3;
    out.e013 =  in.e2;
    out.e023 = -in.e1;
    out.e123 =  in.e0;
}


// ********************************
// A join takes a j-vector and a k-vector and returns a (j+k)-vector
inline
void join(GmpExt4_2 &out, const GmpExt4_1 &lhs, const GmpExt4_1 &rhs)
{
    out.e01 = (lhs.e0 * rhs.e1) - (rhs.e0 * lhs.e1);
    out.e02 = (lhs.e0 * rhs.e2) - (rhs.e0 * lhs.e2);
    out.e03 = (lhs.e0 * rhs.e3) - (rhs.e0 * lhs.e3);
    out.e12 = (lhs.e1 * rhs.e2) - (rhs.e1 * lhs.e2);
    out.e13 = (lhs.e1 * rhs.e3) - (rhs.e1 * lhs.e3);
    out.e23 = (lhs.e2 * rhs.e3) - (rhs.e2 * lhs.e3);
}
inline
void join(GmpExt4_3 &out, const GmpExt4_2 &lhs, const GmpExt4_1 &rhs)
{
    out.e012 = (lhs.e01 * rhs.e2) - (lhs.e02 * rhs.e1) + (lhs.e12 *rhs.e0);
    out.e013 = (lhs.e01 * rhs.e3) - (lhs.e03 * rhs.e1) + (lhs.e13 *rhs.e0);
    out.e023 = (lhs.e02 * rhs.e3) - (lhs.e03 * rhs.e2) + (lhs.e23 *rhs.e0);
    out.e123 = (lhs.e12 * rhs.e3) - (lhs.e13 * rhs.e2) + (lhs.e23 *rhs.e1);
}
inline
void join(GmpExt4_3 &out, const GmpExt4_1 &lhs, const GmpExt4_2 &rhs)
{
    join(out, rhs, lhs);
    // no negation since swapping the arguments requires two
    // swaps of 1-vectors
}


// ********************************
// A meet takes a j-vector and a k-vector and returns a (j+k-4)-vector
inline
void meet(GmpExt4_2 &out, const GmpExt4_3 &lhs, const GmpExt4_3 &rhs)
{
    GmpExt4_2 out_dual;
    GmpExt4_1 lhs_dual;
    GmpExt4_1 rhs_dual;
    dual(lhs_dual, lhs);
    dual(rhs_dual, rhs);
    join(out_dual, lhs_dual, rhs_dual);
    revdual(out, out_dual);
}
inline
void meet(GmpExt4_1 &out, const GmpExt4_2 &lhs, const GmpExt4_3 &rhs)
{
    GmpExt4_3 out_dual;
    GmpExt4_2 lhs_dual;
    GmpExt4_1 rhs_dual;
    dual(lhs_dual, lhs);
    dual(rhs_dual, rhs);
    join(out_dual, lhs_dual, rhs_dual);
    revdual(out, out_dual);
}
inline
void meet(GmpExt4_1 &out, const GmpExt4_3 &lhs, const GmpExt4_2 &rhs)
{
    GmpExt4_3 out_dual;
    GmpExt4_1 lhs_dual;
    GmpExt4_2 rhs_dual;
    dual(lhs_dual, lhs);
    dual(rhs_dual, rhs);
    join(out_dual, lhs_dual, rhs_dual);
    revdual(out, out_dual);
}


// ********************************
// An inner product takes two k-vectors and produces a single number
inline
mpz_class inner(const GmpExt4_1 &lhs, const GmpExt4_1 &rhs)
{
    return lhs.e0 * rhs.e0 +
           lhs.e1 * rhs.e1 +
           lhs.e2 * rhs.e2 +
           lhs.e3 * rhs.e3;
}
inline
mpz_class inner(const GmpExt4_2 &lhs, const GmpExt4_2 &rhs)
{
    return lhs.e01 * rhs.e01 +
           lhs.e02 * rhs.e02 +
           lhs.e03 * rhs.e03 +
           lhs.e12 * rhs.e12 +
           lhs.e13 * rhs.e13 +
           lhs.e23 * rhs.e23;
}
inline
mpz_class inner(const GmpExt4_3 &lhs, const GmpExt4_3 &rhs)
{
    return lhs.e012 * rhs.e012 +
           lhs.e013 * rhs.e013 +
           lhs.e023 * rhs.e023 +
           lhs.e123 * rhs.e123;
}






} // end namespace GMPExt4
