// +-------------------------------------------------------------------------
// | absext4.h
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

/*
 *
 *  AbsExt4
 *
 *      Support for performing exterior calculus in R4
 *      
 *      Used for performing filter computations.
 *      All numbers are positive! (overflow/underflow are concerns)
 *          and all numbers are kept positive
 *
 */

#include "ext4.h"

#include <cmath>
 
namespace AbsExt4 {

using namespace Ext4;
using std::fabs;

struct AbsExt4_1 {
    union {
        double v[4];
        struct {
            double e0;
            double e1;
            double e2;
            double e3;
        };
    };
};

struct AbsExt4_2 {
    union {
        double v[6];
        struct {
            double e01;
            double e02;
            double e03;
            double e12;
            double e13;
            double e23;
        };
    };
};

struct AbsExt4_3 {
    union {
        double v[4];
        struct {
            double e012;
            double e013;
            double e023;
            double e123;
        };
    };
};

// ************************
// An abs takes a k-vector and returns a version with
//  absolute values of all the coordinates taken
inline void abs(AbsExt4_1 &out, const Ext4_1 &in) {
    for(int i=0; i<4; i++)
        out.v[i] = fabs(in.v[i]);
}
inline void abs(AbsExt4_2 &out, const Ext4_2 &in) {
    for(int i=0; i<6; i++)
        out.v[i] = fabs(in.v[i]);
}
inline void abs(AbsExt4_3 &out, const Ext4_3 &in) {
    for(int i=0; i<4; i++)
        out.v[i] = fabs(in.v[i]);
}


// ************************
// A neg takes a k-vector and returns its negation
// neg(X,Y) is safe for X=Y
inline void neg(AbsExt4_1 &out, const AbsExt4_1 &in) {
    out = in;
}
inline void neg(AbsExt4_2 &out, const AbsExt4_2 &in) {
    out = in;
}
inline void neg(AbsExt4_3 &out, const AbsExt4_3 &in) {
    out = in;
}


// ************************
// A dual operation takes a k-vector and returns a (4-k)-vector
// A reverse dual operation inverts the dual operation
// dual(X,Y) is not safe for X=Y (same with revdual)
inline void dual(AbsExt4_1 &out, const AbsExt4_3 &in) {
    out.e0 = in.e123;
    out.e1 = in.e023;
    out.e2 = in.e013;
    out.e3 = in.e012;
}
inline void dual(AbsExt4_2 &out, const AbsExt4_2 &in) {
    out.e01 = in.e23;
    out.e02 = in.e13;
    out.e03 = in.e12;
    out.e12 = in.e03;
    out.e13 = in.e02;
    out.e23 = in.e01;
}
inline void dual(AbsExt4_3 &out, const AbsExt4_1 &in) {
    out.e012 = in.e3;
    out.e013 = in.e2;
    out.e023 = in.e1;
    out.e123 = in.e0;
}
inline void revdual(AbsExt4_1 &out, const AbsExt4_3 &in) {
    out.e0 = in.e123;
    out.e1 = in.e023;
    out.e2 = in.e013;
    out.e3 = in.e012;
}
inline void revdual(AbsExt4_2 &out, const AbsExt4_2 &in) {
    out.e01 = in.e23;
    out.e02 = in.e13;
    out.e03 = in.e12;
    out.e12 = in.e03;
    out.e13 = in.e02;
    out.e23 = in.e01;
}
inline void revdual(AbsExt4_3 &out, const AbsExt4_1 &in) {
    out.e012 = in.e3;
    out.e013 = in.e2;
    out.e023 = in.e1;
    out.e123 = in.e0;
}


// ************************
// A join takes a j-vector and a k-vector and returns a (j+k)-vector
inline void join(AbsExt4_2 &out, const AbsExt4_1 &lhs, const AbsExt4_1 &rhs) {
    out.e01 = (lhs.e0 * rhs.e1) + (rhs.e0 * lhs.e1);
    out.e02 = (lhs.e0 * rhs.e2) + (rhs.e0 * lhs.e2);
    out.e03 = (lhs.e0 * rhs.e3) + (rhs.e0 * lhs.e3);
    out.e12 = (lhs.e1 * rhs.e2) + (rhs.e1 * lhs.e2);
    out.e13 = (lhs.e1 * rhs.e3) + (rhs.e1 * lhs.e3);
    out.e23 = (lhs.e2 * rhs.e3) + (rhs.e2 * lhs.e3);
}
inline void join(AbsExt4_3 &out, const AbsExt4_2 &lhs, const AbsExt4_1 &rhs) {
    out.e012 = (lhs.e01 * rhs.e2) + (lhs.e02 * rhs.e1) + (lhs.e12 *rhs.e0);
    out.e013 = (lhs.e01 * rhs.e3) + (lhs.e03 * rhs.e1) + (lhs.e13 *rhs.e0);
    out.e023 = (lhs.e02 * rhs.e3) + (lhs.e03 * rhs.e2) + (lhs.e23 *rhs.e0);
    out.e123 = (lhs.e12 * rhs.e3) + (lhs.e13 * rhs.e2) + (lhs.e23 *rhs.e1);
}
inline void join(AbsExt4_3 &out, const AbsExt4_1 &lhs, const AbsExt4_2 &rhs) {
    join(out, rhs, lhs);
}


// ************************
// A meet takes a j-vector and a k-vector and returns a (j+k-4)-vector
inline void meet(AbsExt4_2 &out, const AbsExt4_3 &lhs, const AbsExt4_3 &rhs) {
    AbsExt4_2 out_dual;
    AbsExt4_1 lhs_dual;
    AbsExt4_1 rhs_dual;
    dual(lhs_dual, lhs);
    dual(rhs_dual, rhs);
    join(out_dual, lhs_dual, rhs_dual);
    revdual(out, out_dual);
}
inline void meet(AbsExt4_1 &out, const AbsExt4_2 &lhs, const AbsExt4_3 &rhs) {
    AbsExt4_3 out_dual;
    AbsExt4_2 lhs_dual;
    AbsExt4_1 rhs_dual;
    dual(lhs_dual, lhs);
    dual(rhs_dual, rhs);
    join(out_dual, lhs_dual, rhs_dual);
    revdual(out, out_dual);
}
inline void meet(AbsExt4_1 &out, const AbsExt4_3 &lhs, const AbsExt4_2 &rhs) {
    AbsExt4_3 out_dual;
    AbsExt4_1 lhs_dual;
    AbsExt4_2 rhs_dual;
    dual(lhs_dual, lhs);
    dual(rhs_dual, rhs);
    join(out_dual, lhs_dual, rhs_dual);
    revdual(out, out_dual);
}


// ************************
// An inner product takes two k-vectors and produces a single number
inline double inner(const AbsExt4_1 &lhs, const AbsExt4_1 &rhs) {
    double acc = 0.0;
    for(int i=0; i<4; i++)
        acc += lhs.v[i] * rhs.v[i];
    return acc;
}
inline double inner(const AbsExt4_2 &lhs, const AbsExt4_2 &rhs) {
    double acc = 0.0;
    for(int i=0; i<6; i++)
        acc += lhs.v[i] * rhs.v[i];
    return acc;
}
inline double inner(const AbsExt4_3 &lhs, const AbsExt4_3 &rhs) {
    double acc = 0.0;
    for(int i=0; i<4; i++)
        acc += lhs.v[i] * rhs.v[i];
    return acc;
}




} // end namespace AbsExt4


