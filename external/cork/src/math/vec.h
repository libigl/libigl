// +-------------------------------------------------------------------------
// | vec.h
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

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "prelude.h"

// **************************************************************************
// *  Vec2 stores 2-dimensional vectors using whatever base type you choose
// **************************************************************************
template<class N>
class Vec2 {
public: // names
    union {
        N v[2];
        struct {
            N x;
            N y;
        };
        struct {
            N s;
            N t;
        };
    };
// +---------------------------------
public: // constructors
    inline Vec2() {}
    inline Vec2(const N &n0, const N &n1) : x(n0), y(n1) {}
    // allows implicit conversion based on generics...
    template<class T>
    inline Vec2(const Vec2<T> &cp) : x(cp.x), y(cp.y) {}
    // build from short array in memory
    inline Vec2(const N *ns) : x(ns[0]), y(ns[1]) {}
// +---------------------------------
public: // indexing
    inline N& operator[](uint i)       { return v[i]; }
    inline N  operator[](uint i) const { return v[i]; }
// +---------------------------------
public: // destructive arithmetic
    inline Vec2<N>& operator+=(const Vec2<N> &rhs);
    inline Vec2<N>& operator-=(const Vec2<N> &rhs);
    inline Vec2<N>& operator*=(const N &rhs);
    inline Vec2<N>& operator/=(const N &rhs);
};
// +---------------------------------
// comparison operators
template<class N>
inline bool operator==(const Vec2<N> &lhs, const Vec2<N> &rhs);
template<class N>
inline bool operator!=(const Vec2<N> &lhs, const Vec2<N> &rhs);
// +---------------------------------
// non-destructive arithmetic
template<class N>
inline Vec2<N> operator+(const Vec2<N> &lhs, const Vec2<N> &rhs);
template<class N>
inline Vec2<N> operator-(const Vec2<N> &lhs, const Vec2<N> &rhs);
template<class N>
inline Vec2<N> operator*(const Vec2<N> &lhs, const N &rhs);
template<class N>
inline Vec2<N> operator*(const N &lhs, const Vec2<N> &rhs);
template<class N>
inline Vec2<N> operator/(const Vec2<N> &lhs, const N &rhs);
template<class N>
inline Vec2<N> operator-(const Vec2<N> &vec);
template<class N>
inline Vec2<N> abs(const Vec2<N> &vec);
// +---------------------------------
// component-wise max/min
template<class N>
inline Vec2<N> max(const Vec2<N> &lhs, const Vec2<N> &rhs);
template<class N>
inline Vec2<N> min(const Vec2<N> &lhs, const Vec2<N> &rhs);
// +---------------------------------
// dimension logic for projections
template<class N>
inline uint maxDim(const Vec2<N> &vec);
template<class N>
inline uint minDim(const Vec2<N> &vec);
template<class N>
inline N proj(uint dim, const Vec2<N> &vec);
// +---------------------------------
// more sophisticated vector arithmetic
template<class N>
inline N dot(const Vec2<N> &lhs, const Vec2<N> &rhs);
// determinant of 2 vec2s
template<class N>
inline N det(const Vec2<N> &lhs, const Vec2<N> &rhs);
// +---------------------------------
// collapsing measurements
template<class N>
inline N len2(const Vec2<N> &vec);
template<class N>
inline N len(const Vec2<N> &vec);
template<class N>
inline N max(const Vec2<N> &vec);
template<class N>
inline N min(const Vec2<N> &vec);
// +---------------------------------
// normalization functions
template<class N>
inline void normalize(Vec2<N> &vec);
template<class N>
inline Vec2<N> normalized(const Vec2<N> &vec);
// +---------------------------------
// output/input stream operators
template<class N>
inline std::ostream& operator<<(std::ostream &out, const Vec2<N> &vec) {
    return out << '[' << vec.x << ',' << vec.y << ']';
}

// **************************************************************************
// *  Vec3 stores 3-dimensional vectors using whatever base type you choose
// **************************************************************************
template<class N>
class Vec3 {
public: // names
    union {
        N v[3];
        struct {
            N x;
            N y;
            N z;
        };
        struct {
            N r;
            N g;
            N b;
        };
    };
// +---------------------------------
public: // constructors
    inline Vec3() {}
    inline Vec3(const N &n0, const N &n1, const N &n2) : x(n0), y(n1), z(n2) {}
    // allows implicit conversion based on generics...
    template<class T>
    inline Vec3(const Vec3<T> &cp) : x(cp.x), y(cp.y), z(cp.z) {}
    // build from short array in memory
    inline Vec3(const N *ns) : x(ns[0]), y(ns[1]), z(ns[2]) {}
// +---------------------------------
public: // indexing
    inline N& operator[](uint i)       { return v[i]; }
    inline N  operator[](uint i) const { return v[i]; }
// +---------------------------------
public: // destructive arithmetic
    inline Vec3<N>& operator+=(const Vec3<N> &rhs);
    inline Vec3<N>& operator-=(const Vec3<N> &rhs);
    inline Vec3<N>& operator*=(const N &rhs);
    inline Vec3<N>& operator/=(const N &rhs);
};
// +---------------------------------
// comparison operators
template<class N>
inline bool operator==(const Vec3<N> &lhs, const Vec3<N> &rhs);
template<class N>
inline bool operator!=(const Vec3<N> &lhs, const Vec3<N> &rhs);
// +---------------------------------
// non-destructive arithmetic
template<class N>
inline Vec3<N> operator+(const Vec3<N> &lhs, const Vec3<N> &rhs);
template<class N>
inline Vec3<N> operator-(const Vec3<N> &lhs, const Vec3<N> &rhs);
template<class N>
inline Vec3<N> operator*(const Vec3<N> &lhs, const N &rhs);
template<class N>
inline Vec3<N> operator*(const N &lhs, const Vec3<N> &rhs);
template<class N>
inline Vec3<N> operator/(const Vec3<N> &lhs, const N &rhs);
template<class N>
inline Vec3<N> operator-(const Vec3<N> &vec);
template<class N>
inline Vec3<N> abs(const Vec3<N> &vec);
// +---------------------------------
// component-wise max/min
template<class N>
inline Vec3<N> max(const Vec3<N> &lhs, const Vec3<N> &rhs);
template<class N>
inline Vec3<N> min(const Vec3<N> &lhs, const Vec3<N> &rhs);
// +---------------------------------
// dimension logic for projections
template<class N>
inline uint maxDim(const Vec3<N> &vec);
template<class N>
inline uint minDim(const Vec3<N> &vec);
template<class N>
inline Vec2<N> proj(uint dim, const Vec3<N> &vec);
// +---------------------------------
// more sophisticated vector arithmetic
template<class N>
inline N dot(const Vec3<N> &lhs, const Vec3<N> &rhs);
template<class N>
inline Vec3<N> cross(const Vec3<N> &lhs, const Vec3<N> &rhs);
// determinant of 3 Vec3s
template<class N>
inline N det(const Vec3<N> &v0, const Vec3<N> &v1, const Vec3<N> &v2);
// +---------------------------------
// collapsing measurements
template<class N>
inline N len2(const Vec3<N> &vec);
template<class N>
inline N len(const Vec3<N> &vec);
template<class N>
inline N max(const Vec3<N> &vec);
template<class N>
inline N min(const Vec3<N> &vec);
// +---------------------------------
// normalization functions
template<class N>
inline void normalize(Vec3<N> &vec);
template<class N>
inline Vec3<N> normalized(const Vec3<N> &vec);
// +---------------------------------
// output/input stream operators
template<class N>
inline std::ostream& operator<<(std::ostream &out, const Vec3<N> &vec) {
    return out << '[' << vec.x << ',' << vec.y << ',' << vec.z << ']';
}

// **************************************************************************
// *  Vec4 stores 4-dimensional vectors using whatever base type you choose
// **************************************************************************
template<class N>
class Vec4 {
public: // names
    union {
        N v[4];
        struct {
            N x;
            N y;
            N z;
            N w;
        };
        struct {
            N r;
            N g;
            N b;
            N a;
        };
    };
// +---------------------------------
public: // constructors
    inline Vec4() {}
    inline Vec4(const N &n0, const N &n1, const N &n2, const N &n3) :
        x(n0), y(n1), z(n2), w(n3) {}
    // allows implicit conversion based on generics...
    template<class T>
    inline Vec4(const Vec4<T> &cp) : x(cp.x), y(cp.y), z(cp.z), w(cp.w) {}
    // build from short array in memory
    inline Vec4(const N *ns) : x(ns[0]), y(ns[1]), z(ns[2]), w(ns[3]) {}
// +---------------------------------
public: // indexing
    inline N& operator[](uint i)       { return v[i]; }
    inline N  operator[](uint i) const { return v[i]; }
// +---------------------------------
public: // destructive arithmetic
    inline Vec4<N>& operator+=(const Vec4<N> &rhs);
    inline Vec4<N>& operator-=(const Vec4<N> &rhs);
    inline Vec4<N>& operator*=(const N &rhs);
    inline Vec4<N>& operator/=(const N &rhs);
};
// +---------------------------------
// comparison operators
template<class N>
inline bool operator==(const Vec4<N> &lhs, const Vec4<N> &rhs);
template<class N>
inline bool operator!=(const Vec4<N> &lhs, const Vec4<N> &rhs);
// +---------------------------------
// non-destructive arithmetic
template<class N>
inline Vec4<N> operator+(const Vec4<N> &lhs, const Vec4<N> &rhs);
template<class N>
inline Vec4<N> operator-(const Vec4<N> &lhs, const Vec4<N> &rhs);
template<class N>
inline Vec4<N> operator*(const Vec4<N> &lhs, const N &rhs);
template<class N>
inline Vec4<N> operator*(const N &lhs, const Vec4<N> &rhs);
template<class N>
inline Vec4<N> operator/(const Vec4<N> &lhs, const N &rhs);
template<class N>
inline Vec4<N> operator-(const Vec4<N> &vec);
template<class N>
inline Vec4<N> abs(const Vec4<N> &vec);
// +---------------------------------
// component-wise max/min
template<class N>
inline Vec4<N> max(const Vec4<N> &lhs, const Vec4<N> &rhs);
template<class N>
inline Vec4<N> min(const Vec4<N> &lhs, const Vec4<N> &rhs);
// +---------------------------------
// dimension logic for projections
template<class N>
inline uint maxDim(const Vec4<N> &vec);
template<class N>
inline uint minDim(const Vec4<N> &vec);
template<class N>
inline Vec3<N> proj(uint dim, const Vec4<N> &vec);
// +---------------------------------
// more sophisticated vector arithmetic
template<class N>
inline N dot(const Vec4<N> &lhs, const Vec4<N> &rhs);
// determinant of 4 Vec4s
template<class N>
inline N det(const Vec4<N> &v0, const Vec4<N> &v1,
             const Vec4<N> &v2, const Vec4<N> &v3);
// +---------------------------------
// collapsing measurements
template<class N>
inline N len2(const Vec4<N> &vec);
template<class N>
inline N len(const Vec4<N> &vec);
template<class N>
inline N max(const Vec4<N> &vec);
template<class N>
inline N min(const Vec4<N> &vec);
// +---------------------------------
// normalization functions
template<class N>
inline void normalize(Vec4<N> &vec);
template<class N>
inline Vec4<N> normalized(const Vec4<N> &vec);
// +---------------------------------
// output/input stream operators
template<class N>
inline std::ostream& operator<<(std::ostream &out, const Vec4<N> &vec) {
    return out << '[' << vec.x << ',' << vec.y << ','
                      << vec.z << ',' << vec.w << ']';
}

// **************************************************************************
// *  Aliases for Common Specializations
// **************************************************************************
typedef Vec2<int> Vec2i;

typedef Vec2<float> Vec2f;
typedef Vec3<float> Vec3f;
typedef Vec4<float> Vec4f;

typedef Vec2<double> Vec2d;
typedef Vec3<double> Vec3d;
typedef Vec4<double> Vec4d;

// **************************************************************************
// *  Conversions Between Dimensions
// **************************************************************************
// +---------------------------------
// Cartesian -> Homogeneous
template<class N>
inline Vec3<N> toHom(const Vec2<N> &v2) {
    return Vec3<N>(v2.x, v2.y, 1);
}
template<class N>
inline Vec4<N> toHom(const Vec3<N> &v3) {
    return Vec4<N>(v3.x, v3.y, v3.z, 1);
}
/*
// Homogeneous -> Cartesian (BEWARE when last coordinate == 0)
template<class N>
inline Vec2<N> cartProj(const Vec3<N> &v3) {
    return Vec2<N>(v3.x, v3.y) / v3.z;
}
template<class N>
inline Vec3<N> cartProj(const Vec4<N> &v4) {
    return Vec3<N>(v4.x, v4.y, v4.z) / v4.w;
}*/


// **************************************************************************
// **************************************************************************
// ***************************  IMPLEMENTATION  *****************************
// **************************************************************************
// **************************************************************************

// **************************************************************************
// *  Vec2 IMPLEMENTATION
// **************************************************************************
// +---------------------------------
// destructive arithmetic
template<class N>
inline Vec2<N>& Vec2<N>::operator+=(const Vec2<N> &rhs) {
    x += rhs.x;
    y += rhs.y;
    return *this;
}
template<class N>
inline Vec2<N>& Vec2<N>::operator-=(const Vec2<N> &rhs) {
    x -= rhs.x;
    y -= rhs.y;
    return *this;
}
template<class N>
inline Vec2<N>& Vec2<N>::operator*=(const N &rhs) {
    x *= rhs;
    y *= rhs;
    return *this;
}
template<class N>
inline Vec2<N>& Vec2<N>::operator/=(const N &rhs) {
    x /= rhs;
    y /= rhs;
    return *this;
}
// +---------------------------------
// comparison operators
template<class N>
inline bool operator==(const Vec2<N> &lhs, const Vec2<N> &rhs) {
    return lhs.x == rhs.x && lhs.y == rhs.y;
}
template<class N>
inline bool operator!=(const Vec2<N> &lhs, const Vec2<N> &rhs) {
    return lhs.x != rhs.x || lhs.y != rhs.y;
}
// +---------------------------------
// non-destructive arithmetic
template<class N>
inline Vec2<N> operator+(const Vec2<N> &lhs, const Vec2<N> &rhs) {
    Vec2<N> res(lhs);
    return res += rhs;
}
template<class N>
inline Vec2<N> operator-(const Vec2<N> &lhs, const Vec2<N> &rhs) {
    Vec2<N> res(lhs);
    return res -= rhs;
}
template<class N>
inline Vec2<N> operator*(const Vec2<N> &lhs, const N &rhs) {
    Vec2<N> res(lhs);
    return res *= rhs;
}
template<class N>
inline Vec2<N> operator*(const N &lhs, const Vec2<N> &rhs) {
    Vec2<N> res(rhs);
    return res *= lhs;
}
template<class N>
inline Vec2<N> operator/(const Vec2<N> &lhs, const N &rhs) {
    Vec2<N> res(lhs);
    return res /= rhs;
}
template<class N>
inline Vec2<N> operator-(const Vec2<N> &vec) {
    return Vec2<N>(-vec.x, -vec.y);
}
template<>
inline Vec2<double> abs(const Vec2<double> &vec) {
    return Vec2<double>(fabs(vec.x), fabs(vec.y));
}
template<>
inline Vec2<float> abs(const Vec2<float> &vec) {
    return Vec2<float>(fabs(vec.x), fabs(vec.y));
}
// +---------------------------------
// component-wise max/min
template<class N>
inline Vec2<N> max(const Vec2<N> &lhs, const Vec2<N> &rhs) {
    return Vec2<N>(std::max(lhs.x, rhs.x),
                   std::max(lhs.y, rhs.y));
}
template<class N>
inline Vec2<N> min(const Vec2<N> &lhs, const Vec2<N> &rhs) {
    return Vec2<N>(std::min(lhs.x, rhs.x),
                   std::min(lhs.y, rhs.y));
}
// +---------------------------------
// dimension logic for projections
template<class N>
inline uint maxDim(const Vec2<N> &vec) {
    return (vec.x >= vec.y)? 0 : 1;
}
template<class N>
inline uint minDim(const Vec2<N> &vec) {
    return (vec.x <= vec.y)? 0 : 1;
}
template<class N>
inline N proj(uint dim, const Vec2<N> &vec) {
    return (dim == 0)? vec.y : vec.x;
}
// +---------------------------------
// more sophisticated vector arithmetic
template<class N>
inline N dot(const Vec2<N> &lhs, const Vec2<N> &rhs) {
    return lhs.x*rhs.x + lhs.y*rhs.y;
}
// determinant of 2 vec2s
template<class N>
inline N det(const Vec2<N> &lhs, const Vec2<N> &rhs) {
    return lhs.x*rhs.y - lhs.y*rhs.x;
}
// +---------------------------------
// collapsing measurements
template<class N>
inline N len2(const Vec2<N> &vec) {
    return vec.x*vec.x + vec.y*vec.y;
}
template<class N>
inline N len(const Vec2<N> &vec) {
    return std::sqrt(len2(vec));
}
template<class N>
inline N max(const Vec2<N> &vec) {
    return std::max(vec.x, vec.y);
}
template<class N>
inline N min(const Vec2<N> &vec) {
    return std::min(vec.x, vec.y);
}
// +---------------------------------
// normalization functions
template<class N>
inline void normalize(Vec2<N> &vec) {
    vec /= len(vec);
}
template<class N>
inline Vec2<N> normalized(const Vec2<N> &vec) {
    return vec / len(vec);
}

// **************************************************************************
// *  Vec3 IMPLEMENTATION
// **************************************************************************
// +---------------------------------
// destructive arithmetic
template<class N>
inline Vec3<N>& Vec3<N>::operator+=(const Vec3<N> &rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
}
template<class N>
inline Vec3<N>& Vec3<N>::operator-=(const Vec3<N> &rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
}
template<class N>
inline Vec3<N>& Vec3<N>::operator*=(const N &rhs) {
    x *= rhs;
    y *= rhs;
    z *= rhs;
    return *this;
}
template<class N>
inline Vec3<N>& Vec3<N>::operator/=(const N &rhs) {
    x /= rhs;
    y /= rhs;
    z /= rhs;
    return *this;
}
// +---------------------------------
// comparison operators
template<class N>
inline bool operator==(const Vec3<N> &lhs, const Vec3<N> &rhs) {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}
template<class N>
inline bool operator!=(const Vec3<N> &lhs, const Vec3<N> &rhs) {
    return lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z;
}
// +---------------------------------
// non-destructive arithmetic
template<class N>
inline Vec3<N> operator+(const Vec3<N> &lhs, const Vec3<N> &rhs) {
    Vec3<N> res(lhs);
    return res += rhs;
}
template<class N>
inline Vec3<N> operator-(const Vec3<N> &lhs, const Vec3<N> &rhs) {
    Vec3<N> res(lhs);
    return res -= rhs;
}
template<class N>
inline Vec3<N> operator*(const Vec3<N> &lhs, const N &rhs) {
    Vec3<N> res(lhs);
    return res *= rhs;
}
template<class N>
inline Vec3<N> operator*(const N &lhs, const Vec3<N> &rhs) {
    Vec3<N> res(rhs);
    return res *= lhs;
}
template<class N>
inline Vec3<N> operator/(const Vec3<N> &lhs, const N &rhs) {
    Vec3<N> res(lhs);
    return res /= rhs;
}
template<class N>
inline Vec3<N> operator-(const Vec3<N> &vec) {
    return Vec3<N>(-vec.x, -vec.y, -vec.z);
}
template<>
inline Vec3<double> abs(const Vec3<double> &vec) {
    return Vec3<double>(fabs(vec.x), fabs(vec.y), fabs(vec.z));
}
template<>
inline Vec3<float> abs(const Vec3<float> &vec) {
    return Vec3<float>(fabs(vec.x), fabs(vec.y), fabs(vec.z));
}
// +---------------------------------
// component-wise max/min
template<class N>
inline Vec3<N> max(const Vec3<N> &lhs, const Vec3<N> &rhs) {
    return Vec3<N>(std::max(lhs.x, rhs.x),
                   std::max(lhs.y, rhs.y),
                   std::max(lhs.z, rhs.z));
}
template<class N>
inline Vec3<N> min(const Vec3<N> &lhs, const Vec3<N> &rhs) {
    return Vec3<N>(std::min(lhs.x, rhs.x),
                   std::min(lhs.y, rhs.y),
                   std::min(lhs.z, rhs.z));
}
// +---------------------------------
// dimension logic for projections
template<class N>
inline uint maxDim(const Vec3<N> &vec) {
    return (vec.x >= vec.y)? ((vec.x >= vec.z)? 0 : 2) :
                             ((vec.y >= vec.z)? 1 : 2);
}
template<class N>
inline uint minDim(const Vec3<N> &vec) {
    return (vec.x <= vec.y)? ((vec.x <= vec.z)? 0 : 2) :
                             ((vec.y <= vec.z)? 1 : 2);
}
template<class N>
inline Vec2<N> proj(uint dim, const Vec3<N> &vec) {
    return Vec2<N>( vec.v[(dim+1)%3], vec.v[(dim+2)%3] );
}
// +---------------------------------
// more sophisticated vector arithmetic
template<class N>
inline N dot(const Vec3<N> &lhs, const Vec3<N> &rhs) {
    return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
}
template<class N>
inline Vec3<N> cross(const Vec3<N> &lhs, const Vec3<N> &rhs) {
    return Vec3<N>(lhs.y*rhs.z - lhs.z*rhs.y,
                   lhs.z*rhs.x - lhs.x*rhs.z,
                   lhs.x*rhs.y - lhs.y*rhs.x);
}
// determinant of 3 Vec3s
template<class N>
inline N det(const Vec3<N> &v0, const Vec3<N> &v1, const Vec3<N> &v2) {
    N xy = v0.x*v1.y - v0.y*v1.x;
    N xz = v0.x*v1.z - v0.z*v1.x;
    N yz = v0.y*v1.z - v0.z*v1.y;
    return xy*v2.z - xz*v2.y + yz*v2.x;
}
// +---------------------------------
// collapsing measurements
template<class N>
inline N len2(const Vec3<N> &vec) {
    return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
}
template<class N>
inline N len(const Vec3<N> &vec) {
    return std::sqrt(len2(vec));
}
template<class N>
inline N max(const Vec3<N> &vec) {
    return std::max(vec.x, std::max(vec.y, vec.z));
}
template<class N>
inline N min(const Vec3<N> &vec) {
    return std::min(vec.x, std::min(vec.y, vec.z));
}
// +---------------------------------
// normalization functions
template<class N>
inline void normalize(Vec3<N> &vec) {
    vec /= len(vec);
}
template<class N>
inline Vec3<N> normalized(const Vec3<N> &vec) {
    return vec / len(vec);
}

// **************************************************************************
// *  Vec4 IMPLEMENTATION
// **************************************************************************
// destructive arithmetic
template<class N>
inline Vec4<N>& Vec4<N>::operator+=(const Vec4<N> &rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    w += rhs.w;
    return *this;
}
template<class N>
inline Vec4<N>& Vec4<N>::operator-=(const Vec4<N> &rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    w -= rhs.w;
    return *this;
}
template<class N>
inline Vec4<N>& Vec4<N>::operator*=(const N &rhs) {
    x *= rhs;
    y *= rhs;
    z *= rhs;
    w *= rhs;
    return *this;
}
template<class N>
inline Vec4<N>& Vec4<N>::operator/=(const N &rhs) {
    x /= rhs;
    y /= rhs;
    z /= rhs;
    w /= rhs;
    return *this;
}
// +---------------------------------
// comparison operators
template<class N>
inline bool operator==(const Vec4<N> &lhs, const Vec4<N> &rhs) {
    return lhs.x == rhs.x && lhs.y == rhs.y &&
           lhs.z == rhs.z && lhs.w == rhs.w;
}
template<class N>
inline bool operator!=(const Vec4<N> &lhs, const Vec4<N> &rhs) {
    return lhs.x != rhs.x || lhs.y != rhs.y ||
           lhs.z != rhs.z || lhs.w != rhs.w;
}
// +---------------------------------
// non-destructive arithmetic
template<class N>
inline Vec4<N> operator+(const Vec4<N> &lhs, const Vec4<N> &rhs) {
    Vec4<N> res(lhs);
    return res += rhs;
}
template<class N>
inline Vec4<N> operator-(const Vec4<N> &lhs, const Vec4<N> &rhs) {
    Vec4<N> res(lhs);
    return res -= rhs;
}
template<class N>
inline Vec4<N> operator*(const Vec4<N> &lhs, const N &rhs) {
    Vec4<N> res(lhs);
    return res *= rhs;
}
template<class N>
inline Vec4<N> operator*(const N &lhs, const Vec4<N> &rhs) {
    Vec4<N> res(rhs);
    return res *= lhs;
}
template<class N>
inline Vec4<N> operator/(const Vec4<N> &lhs, const N &rhs) {
    Vec4<N> res(lhs);
    return res /= rhs;
}
template<class N>
inline Vec4<N> operator-(const Vec4<N> &vec) {
    return Vec4<N>(-vec.x, -vec.y, -vec.z, -vec.w);
}
template<>
inline Vec4<double> abs(const Vec4<double> &vec) {
    return Vec4<double>(fabs(vec.x), fabs(vec.y), fabs(vec.z), fabs(vec.w));
}
template<>
inline Vec4<float> abs(const Vec4<float> &vec) {
    return Vec4<float>(fabs(vec.x), fabs(vec.y), fabs(vec.z), fabs(vec.w));
}
// +---------------------------------
// component-wise max/min
template<class N>
inline Vec4<N> max(const Vec4<N> &lhs, const Vec4<N> &rhs) {
    return Vec4<N>(std::max(lhs.x, rhs.x),
                   std::max(lhs.y, rhs.y),
                   std::max(lhs.z, rhs.z),
                   std::max(lhs.w, rhs.w));
}
template<class N>
inline Vec4<N> min(const Vec4<N> &lhs, const Vec4<N> &rhs) {
    return Vec4<N>(std::min(lhs.x, rhs.x),
                   std::min(lhs.y, rhs.y),
                   std::min(lhs.z, rhs.z),
                   std::min(lhs.w, rhs.w));
}
// +---------------------------------
// dimension logic for projections
template<class N>
inline uint maxDim(const Vec4<N> &vec) {
    return (vec.x >= vec.y)?
                ((vec.x >= vec.z)? ((vec.x >= vec.w)? 0 : 3) :
                                   ((vec.z >= vec.w)? 2 : 3)) :
                ((vec.y >= vec.z)? ((vec.y >= vec.w)? 1 : 3) :
                                   ((vec.z >= vec.w)? 2 : 3));
}
template<class N>
inline uint minDim(const Vec4<N> &vec) {
    return (vec.x <= vec.y)?
                ((vec.x <= vec.z)? ((vec.x <= vec.w)? 0 : 3) :
                                   ((vec.z <= vec.w)? 2 : 3)) :
                ((vec.y <= vec.z)? ((vec.y <= vec.w)? 1 : 3) :
                                   ((vec.z <= vec.w)? 2 : 3));
}
template<class N>
inline Vec3<N> proj(uint dim, const Vec4<N> &vec) {
    return Vec3<N>( vec.v[(dim+1)%4],
                    vec.v[(dim+2)%4],
                    vec.v[(dim+3)%4] );
}
// +---------------------------------
// more sophisticated vector arithmetic
template<class N>
inline N dot(const Vec4<N> &lhs, const Vec4<N> &rhs) {
    return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z + lhs.w*rhs.w;
}
// determinant of 4 Vec4s
template<class N>
inline N det(const Vec4<N> &v0, const Vec4<N> &v1,
             const Vec4<N> &v2, const Vec4<N> &v3) {
    N xy01 = v0.x*v1.y - v0.y*v1.x;     N yz01 = v0.y*v1.z - v0.z*v1.y;
    N xz01 = v0.x*v1.z - v0.z*v1.x;     N yw01 = v0.y*v1.w - v0.w*v1.y;
    N xw01 = v0.x*v1.w - v0.w*v1.x;     N zw01 = v0.z*v1.w - v0.w*v1.z;
    
    N xy23 = v2.x*v3.y - v2.y*v3.x;     N yz23 = v2.y*v3.z - v2.z*v3.y;
    N xz23 = v2.x*v3.z - v2.z*v3.x;     N yw23 = v2.y*v3.w - v2.w*v3.y;
    N xw23 = v2.x*v3.w - v2.w*v3.x;     N zw23 = v2.z*v3.w - v2.w*v3.z;
    
    return xy01*zw23 - xz01*yw23 + xw01*yz23
         + yz01*xw23 - yw01*xz23 + zw01*xy23;
}
// +---------------------------------
// collapsing measurements
template<class N>
inline N len2(const Vec4<N> &vec) {
    return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z + vec.w*vec.w;
}
template<class N>
inline N len(const Vec4<N> &vec) {
    return std::sqrt(len2(vec));
}
template<class N>
inline N max(const Vec4<N> &vec) {
    return std::max(std::max(vec.x, vec.y), std::max(vec.z, vec.w));
}
template<class N>
inline N min(const Vec4<N> &vec) {
    return std::min(std::min(vec.x, vec.y), std::min(vec.z, vec.w));
}
// +---------------------------------
// normalization functions
template<class N>
inline void normalize(Vec4<N> &vec) {
    vec /= len(vec);
}
template<class N>
inline Vec4<N> normalized(const Vec4<N> &vec) {
    return vec / len(vec);
}

