// +-------------------------------------------------------------------------
// | bbox.h
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

#include "vec.h"

#include <cfloat>

// NOTE on usage of BBoxes
//  all BBoxes are initialized so that
//      convex(BBox(), bb) == bb
//  for any bb


// **************************************************************************
// *  BBox2 stores 2-dimensional axis aligned bounding boxes
// **************************************************************************
template<class N>
class BBox2 {
public: // data
    Vec2<N> minp, maxp;
public: // constructors
    BBox2(const Vec2<N> &minpp, const Vec2<N> &maxpp) :
        minp(minpp), maxp(maxpp) {}
    inline BBox2(const BBox2<N> &bb) : minp(bb.minp), maxp(bb.maxp) {}
    inline BBox2(); // specialized implementations for float/double only
};

template<> inline
BBox2<float>::BBox2() :
    minp( FLT_MAX, FLT_MAX),
    maxp(-FLT_MAX,-FLT_MAX)
{}
template<> inline
BBox2<double>::BBox2() :
    minp( DBL_MAX, DBL_MAX),
    maxp(-DBL_MAX,-DBL_MAX)
{}

template<class N>
inline bool isEmpty(const BBox2<N> &bb) {
    return bb.maxp[0] < bb.minp[0] ||
           bb.maxp[1] < bb.minp[1];
}
template<class N>
inline bool isIn(const Vec2<N> &p, const BBox2<N> &bb) {
    return bb.minp[0] <= p[0] && p[0] <= bb.maxp[0] &&
           bb.minp[1] <= p[1] && p[1] <= bb.maxp[1];
}
template<class N>
inline bool hasIsct(const BBox2<N> &lhs, const BBox2<N> &rhs) {
    return lhs.minp[0] <= rhs.maxp[0] && lhs.maxp[0] >= rhs.minp[0] &&
           lhs.minp[1] <= rhs.maxp[1] && lhs.maxp[1] >= rhs.minp[1];
}
template<class N>
inline BBox2<N> convex(const BBox2<N> &lhs, const BBox2<N> &rhs) {
    return BBox2<N>(min(lhs.minp, rhs.minp), max(lhs.maxp, rhs.maxp));
}
template<class N>
inline BBox2<N> isct(const BBox2<N> &lhs, const BBox2<N> &rhs) {
    return BBox2<N>(max(lhs.minp, rhs.minp), min(lhs.maxp, rhs.maxp));
}

template<class N>
inline Vec2<N> dim(const BBox2<N> &bb) {
    return bb.maxp - bb.minp;
}
template<class N>
inline N perimeter(const BBox2<N> &bb) {
    Vec2<N> d = dim(bb);
    return 2*(d[0] + d[1]);
}

template<class N>
inline std::ostream& operator<<(std::ostream &out, const BBox2<N> &bb) {
    return out << "[min" << bb.minp << ";max" << bb.maxp << ']';
}

// **************************************************************************
// *  BBox3 stores 3-dimensional axis aligned bounding boxes
// **************************************************************************
template<class N>
class BBox3 {
public: // data
    Vec3<N> minp, maxp;
public: // constructors
    BBox3(const Vec3<N> &minpp, const Vec3<N> &maxpp) :
        minp(minpp), maxp(maxpp) {}
    inline BBox3(const BBox3<N> &bb) : minp(bb.minp), maxp(bb.maxp) {}
    inline BBox3(); // specialized implementations for float/double only
};

template<> inline
BBox3<float>::BBox3() :
    minp( FLT_MAX, FLT_MAX, FLT_MAX),
    maxp(-FLT_MAX,-FLT_MAX,-FLT_MAX)
{}
template<> inline
BBox3<double>::BBox3() :
    minp( DBL_MAX, DBL_MAX, DBL_MAX),
    maxp(-DBL_MAX,-DBL_MAX,-DBL_MAX)
{}

template<class N>
inline bool isEmpty(const BBox3<N> &bb) {
    return bb.maxp[0] < bb.minp[0] ||
           bb.maxp[1] < bb.minp[1] ||
           bb.maxp[2] < bb.minp[2];
}
template<class N>
inline bool isIn(const Vec3<N> &p, const BBox3<N> &bb) {
    return bb.minp[0] <= p[0] && p[0] <= bb.maxp[0] &&
           bb.minp[1] <= p[1] && p[1] <= bb.maxp[1] &&
           bb.minp[2] <= p[2] && p[2] <= bb.maxp[2];
}
template<class N>
inline bool hasIsct(const BBox3<N> &lhs, const BBox3<N> &rhs) {
    return lhs.minp[0] <= rhs.maxp[0] && lhs.maxp[0] >= rhs.minp[0] &&
           lhs.minp[1] <= rhs.maxp[1] && lhs.maxp[1] >= rhs.minp[1] &&
           lhs.minp[2] <= rhs.maxp[2] && lhs.maxp[2] >= rhs.minp[2];
}
template<class N>
inline BBox3<N> convex(const BBox3<N> &lhs, const BBox3<N> &rhs) {
    return BBox3<N>(min(lhs.minp, rhs.minp), max(lhs.maxp, rhs.maxp));
}
template<class N>
inline BBox3<N> isct(const BBox3<N> &lhs, const BBox3<N> &rhs) {
    return BBox3<N>(max(lhs.minp, rhs.minp), min(lhs.maxp, rhs.maxp));
}

template<class N>
inline Vec3<N> dim(const BBox3<N> &bb) {
    return bb.maxp - bb.minp;
}
template<class N>
inline N surfaceArea(const BBox3<N> &bb) {
    Vec3<N> d = dim(bb);
    return 2*(d[1]*d[2] + d[0]*d[2] + d[0]*d[1]);
}

template<class N>
inline std::ostream& operator<<(std::ostream &out, const BBox3<N> &bb) {
    return out << "[min" << bb.minp << ";max" << bb.maxp << ']';
}


// **************************************************************************
// *  BBox4 stores 4-dimensional axis aligned bounding boxes
// **************************************************************************
template<class N>
class BBox4 {
public: // data
    Vec4<N> minp, maxp;
public: // constructors
    BBox4(const Vec4<N> &minpp, const Vec4<N> &maxpp) :
        minp(minpp), maxp(maxpp) {}
    inline BBox4(const BBox4<N> &bb) : minp(bb.minp), maxp(bb.maxp) {}
    inline BBox4();  // specialized implementations for float/double only
};

template<> inline
BBox4<float>::BBox4() :
    minp( FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX),
    maxp(-FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX)
{}
template<> inline
BBox4<double>::BBox4() :
    minp( DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX),
    maxp(-DBL_MAX,-DBL_MAX,-DBL_MAX,-DBL_MAX)
{}

template<class N>
inline bool isEmpty(const BBox4<N> &bb) {
    return bb.maxp[0] < bb.minp[0] ||
           bb.maxp[1] < bb.minp[1] ||
           bb.maxp[2] < bb.minp[2] ||
           bb.maxp[3] < bb.minp[3];
}
template<class N>
inline bool isIn(const Vec4<N> &p, const BBox4<N> &bb) {
    return bb.minp[0] <= p[0] && p[0] <= bb.maxp[0] &&
           bb.minp[1] <= p[1] && p[1] <= bb.maxp[1] &&
           bb.minp[2] <= p[2] && p[2] <= bb.maxp[2] &&
           bb.minp[3] <= p[3] && p[3] <= bb.maxp[3];
}
template<class N>
inline bool hasIsct(const BBox4<N> &lhs, const BBox4<N> &rhs) {
    return lhs.minp[0] <= rhs.maxp[0] && lhs.maxp[0] >= rhs.minp[0] &&
           lhs.minp[1] <= rhs.maxp[1] && lhs.maxp[1] >= rhs.minp[1] &&
           lhs.minp[2] <= rhs.maxp[2] && lhs.maxp[2] >= rhs.minp[2] &&
           lhs.minp[3] <= rhs.maxp[3] && lhs.maxp[3] >= rhs.minp[3];
}
template<class N>
inline BBox4<N> convex(const BBox4<N> &lhs, const BBox4<N> &rhs) {
    return BBox4<N>(min(lhs.minp, rhs.minp), max(lhs.maxp, rhs.maxp));
}
template<class N>
inline BBox4<N> isct(const BBox4<N> &lhs, const BBox4<N> &rhs) {
    return BBox4<N>(max(lhs.minp, rhs.minp), min(lhs.maxp, rhs.maxp));
}

template<class N>
inline Vec4<N> dim(const BBox4<N> &bb) {
    return bb.maxp - bb.minp;
}

template<class N>
inline std::ostream& operator<<(std::ostream &out, const BBox4<N> &bb) {
    return out << "[min" << bb.minp << ";max" << bb.maxp << ']';
}



// define some common useful versions of the boxes
typedef BBox2<float> BBox2f;
typedef BBox3<float> BBox3f;
typedef BBox4<float> BBox4f;

typedef BBox2<double> BBox2d;
typedef BBox3<double> BBox3d;
typedef BBox4<double> BBox4d;


