// +-------------------------------------------------------------------------
// | ray.h
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

template <class N>
struct Ray3
{
    Vec3<N> p; // point of origin
    Vec3<N> r; // ray direction
    
    Ray3<N>() {}
    Ray3<N>(const Vec3<N> &point, const Vec3<N> &dir) :
        p(point), r(dir)
    {}
    template<class T>
    Ray3<N>(const Ray3<T> &cp) : p(cp.p), r(cp.r) {}
};


template<class N>
inline std::ostream& operator<<(std::ostream &out, const Ray3<N> &ray) {
    return out << '[' << ray.p << ';' << ray.r << ']';
}


typedef Ray3<float> Ray3f;
typedef Ray3<double> Ray3d;

