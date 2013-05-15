// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#ifndef __EMBREE_TRIANGLE_LIGHT_H__
#define __EMBREE_TRIANGLE_LIGHT_H__

#include "lights/light.h"
#include "shapes/triangle.h"

namespace embree
{
  /*! Implements a triangle shaped area light. */
  class TriangleLight : public AreaLight
  {
    /*! Construction from triangle vertices and radiance. */
    TriangleLight (const Vec3f& v0, const Vec3f& v1, const Vec3f& v2, const Col3f& L)
      : v0(v0), v1(v1), v2(v2), e1(v0-v1), e2(v2-v0), Ng(cross(e1,e2)), L(L), tri(new Triangle(v0,v1,v2)) { }

  public:

    /*! Construction from parameter container. */
    TriangleLight(const Parms& parms) {
      v0 = parms.getVec3f("v0");
      v1 = parms.getVec3f("v1");
      v2 = parms.getVec3f("v2");
      e1 = v0-v1;
      e2 = v2-v0;
      Ng = cross(e1,e2);
      L = parms.getCol3f("L");
      tri = new Triangle(v0,v1,v2);
    }

    Ref<Light> transform(const AffineSpace3f& xfm) const {
      return new TriangleLight(xfmPoint(xfm,v0),xfmPoint(xfm,v1),xfmPoint(xfm,v2),L);
    }

    /*! Returns the shape of the triangle light. */
    Ref<Shape> shape() { return tri.cast<Shape>(); }

    /*! Intersects the triangle with a ray of origin O and directio D. */
    __forceinline void intersect(const Vec3f& O, const Vec3f& D, float& t, float& u, float& v, float& w) const
    {
      Vec3f C = v0 - O;
      Vec3f R = cross(D,C);
      float det = dot(Ng,D);
      float rcp_det = rcp(det);
      t = dot(Ng,C) * rcp_det;
      u = dot(R,e2) * rcp_det;
      v = dot(R,e1) * rcp_det;
      w = float(one)-u-v;
    }

    Col3f Le(const DifferentialGeometry& dg, const Vec3f& wo) const {
      return L;
    }

    Col3f eval(const DifferentialGeometry& dg, const Vec3f& wi) const {
      float t,u,v,w; intersect(dg.P,wi,t,u,v,w);
      if (0.0f <= t && min(u,v,w) >= 0) return L;
      return zero;
    }

    Col3f sample(const DifferentialGeometry& dg, Sample3f& wi, float& tMax, const Vec2f& s) const
    {
      Vec3f d = uniformSampleTriangle(s.x,s.y,tri->v0,tri->v1,tri->v2)-dg.P;
      tMax = length(d);
      float dDotNg = dot(d,Ng);
      if (dDotNg >= 0) return zero;
      wi = Sample3f(d*rcp(tMax),2.0f*tMax*tMax*tMax*rcp(fabs(dDotNg)));
      return L;
    }

    float pdf(const DifferentialGeometry& dg, const Vec3f& wi) const {
      float t,u,v,w; intersect(dg.P,wi,t,u,v,w);
      if (t < 0.0f || min(u,v,w) < 0) return zero;
      return 2.0f*t*t*rcp(abs(dot(wi,Ng)));
    }

  public:
    Vec3f v0;                //!< First vertex of the triangle
    Vec3f v1;                //!< Second vertex of the triangle
    Vec3f v2;                //!< Third vertex of the triangle

  private:
    Vec3f e1;               //!< First edge of the tringle (v0-v1)
    Vec3f e2;               //!< Second edge of the triangle (v2-v0)
    Vec3f Ng;               //!< Normal of triangle light
    Col3f L;                //!< Radiance (W/(m^2*sr))
    Ref<Triangle> tri;      //!< Triangle shape
  };
}

#endif
