// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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

#pragma once

#include "bezier1i.h"
#include "common/ray.h"
#include "geometry/filter.h"

namespace embree
{
  namespace isa
  {
    /* code block for correct cone intersection */
#if 0
/* subdivide 3 levels at once */ 
    const BezierCurve3D curve2D(v0,v1,v2,v3,0.0f,1.0f,0);
    const avx4f a = curve2D.eval(coeff0[0],coeff0[1],coeff0[2],coeff0[3]);
    const avx4f b = curve2D.eval(coeff1[0],coeff1[1],coeff1[2],coeff1[3]); // FIXME: can be calculated from p0 by shifting
    const avx3f a3(a.x,a.y,a.z);
    const avx3f b3(b.x,b.y,b.z);
    
    const avxf  rl0 = 1.0f/length(b3-a3); // FIXME: multiply equation with this
    const avx3f p0 = a3, d0 = (b3-a3)*rl0;
    const avxf  r0 = a.w, dr = (b.w-a.w)*rl0;
    const float rl1 = 1.0f/length(ray.dir); // FIXME: normalization not required
    const avx3f p1 = ray.org, d1 = ray.dir*rl1;
    
    const avx3f dp = p1-p0;
    const avxf dpdp = dot(dp,dp);
    const avxf d1d1 = dot(d1,d1);
    const avxf d0d1 = dot(d0,d1);
    const avxf d0dp = dot(d0,dp);
    const avxf d1dp = dot(d1,dp);
    const avxf R = r0 + d0dp*dr;
    const avxf A = d1d1 - sqr(d0d1) * (1.0f+dr*dr);
    const avxf B = 2.0f * (d1dp - d0d1*(d0dp + R*dr));
    const avxf C = dpdp - (sqr(d0dp) + sqr(R));
    const avxf D = B*B - 4.0f*A*C;
    avxb valid = D >= 0.0f;
    if (none(valid)) return;
    
    const avxf Q = sqrt(D);
    //const avxf t0 = (-B-Q)*rcp2A;
    //const avxf t1 = (-B+Q)*rcp2A;
    const avxf t0 = (-B-Q)/(2.0f*A);
    const avxf u0 = d0dp+t0*d0d1;
    const avxf t = t0*rl1;
    const avxf u = u0*rl0;
    valid &= (ray.tnear < t) & (t < ray.tfar) & (0.0f <= u) & (u <= 1.0f);
#endif
    
    /*! Intersector for a single ray with a bezier curve. */
    struct Bezier1iIntersector1Slow
    {
      typedef Bezier1i Primitive;
      
      struct Precalculations {
        __forceinline Precalculations (const Ray& ray) {}
      };
      
      __forceinline bool intersect_box(const BBox3fa& box, const Ray& ray) {
        return max(box.lower.x,box.lower.y) <= 0.0f && 0.0f <= min(box.upper.x,box.upper.y);
      }
      
      static __forceinline void intersect(const Precalculations& pre, Ray& ray, const Bezier1i& curve_in, const void* geom)
      {
        /* load bezier curve control points */
        STAT3(normal.trav_prims,1,1,1);
        const Vec3fa v0 = curve_in.p[0];
        const Vec3fa v1 = curve_in.p[1];
        const Vec3fa v2 = curve_in.p[2];
        const Vec3fa v3 = curve_in.p[3];
        
        /* transform control points into ray space */
        LinearSpace3fa ray_space = rcp(frame(ray.dir)); // FIXME: calculate once per ray
        Vec3fa w0 = xfmVector(ray_space,v0-ray.org); w0.w = v0.w;
        Vec3fa w1 = xfmVector(ray_space,v1-ray.org); w1.w = v1.w;
        Vec3fa w2 = xfmVector(ray_space,v2-ray.org); w2.w = v2.w;
        Vec3fa w3 = xfmVector(ray_space,v3-ray.org); w3.w = v3.w;
        
        /* hit information */
        float ray_u = 0.0f;
        float ray_tfar = ray.tfar;
        bool hit = false;
        
        /* push first curve onto stack */
        BezierCurve3D stack[32];
        new (&stack[0]) BezierCurve3D(w0,w1,w2,w3,0.0f,1.0f,4);
        size_t sptr = 1;
        
        while (true) 
        {
        pop:
          if (sptr == 0) break;
          BezierCurve3D curve = stack[--sptr];
          
          while (curve.depth)
          {
            BezierCurve3D curve0,curve1;
            curve.subdivide(curve0,curve1);
            
            BBox3fa bounds0 = curve0.bounds();
            BBox3fa bounds1 = curve1.bounds();
            bool hit0 = intersect_box(bounds0,ray);
            bool hit1 = intersect_box(bounds1,ray);
            
            if (!hit0 && !hit1) goto pop;
            else if (likely(hit0 != hit1)) {
              if (hit0) { curve = curve0; continue; } 
              else      { curve = curve1; continue; }
            } else {
              curve = curve0;
              stack[sptr++] = curve1;
            }
          }
          
          /* approximative intersection with cone */
          const Vec3fa v = curve.v3-curve.v0;
          const Vec3fa w = -curve.v0;
          const float d0 = w.x*v.x + w.y*v.y;
          const float d1 = v.x*v.x + v.y*v.y;
          const float u = clamp(d0/d1,0.0f,1.0f);
          const Vec3fa p = curve.v0 + u*v;
          const float d2 = p.x*p.x + p.y*p.y; 
          const float r2 = p.w*p.w;
          if (unlikely(d2 > r2)) continue;
          const float t = p.z;
          if (unlikely(t < ray.tnear || t > ray_tfar)) continue;
          ray_u = curve.t0+u*(curve.t1-curve.t0);
          ray_tfar = t;
          hit = true;
        }
        
        /* compute final hit data */
        if (likely(hit)) 
        {
          BezierCurve3D curve(v0,v1,v2,v3,0.0f,1.0f,0);
          Vec3fa P,T; curve.eval(ray.u,P,T);
          ray.u = ray_u;
          ray.v = 1.0f;
          ray.tfar = ray_tfar;
          ray.Ng = T;
          ray.geomID = curve_in.geomID;
          ray.primID = curve_in.primID;
        }
      }
      
      static __forceinline void intersect(const Precalculations& pre, Ray& ray, const Bezier1i* curves, size_t num, void* geom)
      {
        for (size_t i=0; i<num; i++)
          intersect(pre,ray,curves[i],geom);
      }
      
      static __forceinline bool occluded(const Precalculations& pre, Ray& ray, const Bezier1i& curve_in, const void* geom) {
        return false;
      }
      
      static __forceinline bool occluded(const Precalculations& pre, Ray& ray, const Bezier1i* curves, size_t num, void* geom) 
      {
        for (size_t i=0; i<num; i++) 
          if (occluded(pre,ray,curves[i],geom))
            return true;
        
        return false;
      }
    };
  }
}
