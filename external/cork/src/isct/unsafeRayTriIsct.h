// +-------------------------------------------------------------------------
// | unsafeRayTriIsct.h
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

#include "ray.h"
#include "vec.h"

struct isct_ray_triangle_result
{
    double t;
    Vec3d bary; 
};

inline
bool isct_ray_triangle(
    Ray3d               ray,
    Vec3d               va,
    Vec3d               vb,
    Vec3d               vc,
    double              *t,
    Vec3d               *bary
) {
    ENSURE(t);
    ENSURE(bary);
    
    // re-center the problem at the base point of the ray
    va -= ray.p;
    vb -= ray.p;
    vc -= ray.p;
    
    // Then compute volumes of tetrahedra spanning
    //  * the base point / ray direction line segment
    //  * an edge of the triangle
    // Keeping orientations in mind...
    double volAB  =   det(va, vb, ray.r);
    double volBC  =   det(vb, vc, ray.r);
    double volCA  = - det(va, vc, ray.r);
    // then also compute the volume of tet with the entire triangle as a face...
    double volABC =   det(va, vb, vc);
    // if any of the signs of the edge tests
    // disagree with the sign of the whole triangle, then
    // the ray does not pass through the triangle
    if(volAB * volABC < 0 ||
       volBC * volABC < 0 ||
       volCA * volABC < 0)      return false;
    
    // otherwise, compute the t - value for the ray to intersect
    // if this is negative, then the client can detect that the
    // ray would have to travel backwards to hit the triangle in question.
    double edgeSum = volAB + volBC + volCA;
    if(edgeSum == 0)            return false;
    *t = volABC / (volAB + volBC + volCA);
    if(*t <= 0)                 return false;
    
    *bary = Vec3d(volBC/edgeSum, volCA/edgeSum, volAB/edgeSum);
                                return true;
}



