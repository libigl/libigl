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

#define USE_RECURSIVE_INTERSECTION 0

namespace embree
{
  namespace isa
  {
    /*! Intersector for a single ray with a bezier curve. */
    struct Bezier1iIntersector1
    {
      typedef Bezier1i Primitive;
      
      struct Mailbox {
        avxi ids;      
        unsigned int index;
        
        __forceinline Mailbox() {
          ids = -1;
          index = 0;
        };
      };
      
      struct Precalculations 
      {
        __forceinline Precalculations (const Ray& ray)
          : ray_space(rcp(frame(ray.dir))) {}
        
        LinearSpace3fa ray_space;
        Mailbox mbox;
      };
      
      
      static __forceinline bool intersectCylinder(const Ray& ray,const Vec3fa &v0,const Vec3fa &v1,const Vec3fa &v2,const Vec3fa &v3)
      {
        const Vec3fa cyl_dir = v3 - v0;
        const float dist_v1 = length((v1 - v0) - (dot(v1-v0,cyl_dir)*cyl_dir));
        const float dist_v2 = length((v2 - v0) - (dot(v2-v0,cyl_dir)*cyl_dir));
        
        const float w0 = v0.w;
        const float w1 = v1.w + dist_v1;
        const float w2 = v2.w + dist_v2;
        const float w3 = v3.w;
        
        const float radius = max(max(w0,w1),max(w2,w3));
        const Vec3fa cyl_org = v0;
        const Vec3fa delta_org = ray.org - cyl_org;
        const Vec3fa g0 = ray.dir - (dot(ray.dir,cyl_dir) * cyl_dir);
        const Vec3fa g1 = delta_org - (dot(delta_org,cyl_dir) * cyl_dir);
        
        const float A = dot(g0,g0);
        const float B = 2.0f * dot(g0,g1);
        const float C = dot(g1,g1) - radius*radius;
        
        const float D = B*B - 4*A*C;
        if (D <= 0.0f) return false;
        
        return true;
      }
      
      static __forceinline bool intersectBoxes(const Ray& ray,const Vec3fa &v0,const Vec3fa &v1,const Vec3fa &v2,const Vec3fa &v3)
      {
        BezierCurve3D curve3D(v0,v1,v2,v3,0.0f,1.0f,0);
        const avx4f p0 = curve3D.eval(coeff0[0],coeff0[1],coeff0[2],coeff0[3]);
        const avx4f p1 = curve3D.eval(coeff1[0],coeff1[1],coeff1[2],coeff1[3]);
        const avx3f d0(p0.x,p0.y,p0.z);
        const avx3f d1(p1.x,p1.y,p1.z);
        const avx3f d_min = min(d0,d1) - avx3f(p0.w);
        const avx3f d_max = max(d0,d1) + avx3f(p1.w);
        
        const Vec3fa ray_rdir = rcp_safe(ray.dir);
        const avx3f t_min = (d_min - avx3f(ray.org)) * avx3f(ray_rdir);
        const avx3f t_max = (d_max - avx3f(ray.org)) * avx3f(ray_rdir);
        const avx3f tNear3 = min(t_min,t_max);
        const avx3f tFar3  = max(t_min,t_max);
        const avxf tNear = max(tNear3.x,tNear3.y,tNear3.z,avxf(ray.tnear));
        const avxf tFar  = min(tFar3.x,tFar3.y,tFar3.z,avxf(ray.tfar));
        const avxb vmask = tNear <= tFar;
        return any(vmask);
      }
      
      /*static __forceinline avxf shift_left1(avxf r, float a) 
        {
        #if defined(__AVX2__)
        avxf r0 = _mm256_castsi256_ps(_mm256_slli_si256(_mm256_castps_si256(r),4));
        r0[15] = a;
        return r0;
        #else
        avxf v[2] = { r, avxf(a) };
        return *(avxf*)&v[0][1];
        #endif
        }*/
      
#if 0
      static __forceinline void intersect2(const Precalculations& pre, Ray& ray, const Bezier1i& curve_in, const void* geom)
      {
        /* load bezier curve control points */
        STAT3(normal.trav_prims,1,1,1);
        const Vec3fa& v0 = curve_in.p[0];
        const Vec3fa& v1 = curve_in.p[1];
        const Vec3fa& v2 = curve_in.p[2];
        const Vec3fa& v3 = curve_in.p[3];
        
        /* transform control points into ray space */
        Vec3fa w0 = xfmVector(pre.ray_space,v0-ray.org); w0.w = v0.w;
        Vec3fa w1 = xfmVector(pre.ray_space,v1-ray.org); w1.w = v1.w;
        Vec3fa w2 = xfmVector(pre.ray_space,v2-ray.org); w2.w = v2.w;
        Vec3fa w3 = xfmVector(pre.ray_space,v3-ray.org); w3.w = v3.w;
        BezierCurve3D curve2D(w0,w1,w2,w3,0.0f,1.0f,4);
        
        /* subdivide 3 levels at once */ 
        const avx4f p0 = curve2D.eval(coeff_P0[0],coeff_P0[1],coeff_P0[2],coeff_P0[3]);
        const avx4f p1 = curve2D.eval(coeff_P1[0],coeff_P1[1],coeff_P1[2],coeff_P1[3]);
        //const avx4f p1(shift_left1(p0.x,w3.x),shift_left1(p0.y,w3.y),shift_left1(p0.z,w3.z),shift_left1(p0.w,w3.w));
        
        /* approximative intersection with cone */
        const avx4f v = p1-p0;
        const avx4f w = -p0;
        const avxf d0 = w.x*v.x + w.y*v.y;
        const avxf d1 = v.x*v.x + v.y*v.y;
        const avxf u = clamp(d0*rcp(d1),avxf(zero),avxf(one));
        const avx4f p = p0 + u*v;
        const avxf t = p.z;
        const avxf d2 = p.x*p.x + p.y*p.y; 
        const avxf r = p.w; //max(p.w,ray.org.w+ray.dir.w*t);
        const avxf r2 = r*r;
        avxb valid = d2 <= r2 & avxf(ray.tnear) < t & t < avxf(ray.tfar);
      retry:
        if (unlikely(none(valid))) return;
        const float one_over_8 = 1.0f/8.0f;
        size_t i = select_min(valid,t);
        STAT3(normal.trav_prim_hits,1,1,1);
        
        /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
        int geomID = curve_in.geomID;
        Geometry* geometry = ((Scene*)geom)->get(geomID);
        if (!likely(geometry->hasIntersectionFilter1())) 
        {
#endif
          /* update hit information */
          const float uu = (float(i)+u[i])*one_over_8; // FIXME: correct u range for subdivided segments
          const BezierCurve3D curve3D(v0,v1,v2,v3,0.0f,1.0f,0);
          Vec3fa P,T; curve3D.eval(uu,P,T);
          if (T == Vec3fa(zero)) { valid[i] = 0; goto retry; } // ignore denormalized curves
          ray.u = uu;
          ray.v = 0.0f;
          ray.tfar = t[i];
          ray.Ng = T;
          ray.geomID = curve_in.geomID;
          ray.primID = curve_in.primID;
#if defined(RTCORE_INTERSECTION_FILTER)
          return;
        }
        
        while (true) 
        {
          const float uu = (float(i)+u[i])*one_over_8;
          const BezierCurve3D curve3D(v0,v1,v2,v3,0.0f,1.0f,0);
          Vec3fa P,T; curve3D.eval(uu,P,T);
          if (T != Vec3fa(zero))
            if (runIntersectionFilter1(geometry,ray,uu,0.0f,t[i],T,geomID,curve_in.primID)) return;
          valid[i] = 0;
          if (none(valid)) return;
          i = select_min(valid,t);
          STAT3(normal.trav_prim_hits,1,1,1);
        }
#endif
      }
#endif
      
      static void intersect_recursive(const float max_radius, Ray& ray, const BezierCurve3D &curve2D, const void* geom, int geomID, int primID)
      {
        const BBox3fa bounds = curve2D.bounds();
#if 1
        if (bounds.lower.x >  max_radius ||
            bounds.upper.x < -max_radius ||
            bounds.lower.y >  max_radius ||
            bounds.upper.y < -max_radius) return;
#endif
        if (curve2D.depth == 0)
        {
          const Vec3fa &p0 = curve2D.v0;
          const Vec3fa &p1 = curve2D.v3;
          
          /* approximative intersection with cone */
          const Vec3fa v = p1-p0;
          const Vec3fa w = -p0;
          const float d0 = w.x*v.x + w.y*v.y;
          const float d1 = v.x*v.x + v.y*v.y;
          const float u = clamp(d0*rcp(d1),0.0f,1.0f);
          const Vec3fa p = p0 + u*v;
          const float t = p.z;
          const float d2 = p.x*p.x + p.y*p.y; 
          const float r = p.w; //max(p.w,ray.org.w+ray.dir.w*t);
          const float r2 = r*r;
          
          STAT3(normal.trav_prim_hits,1,1,1);
          
          if (unlikely(!( d2 <= r2 & ray.tnear < t & t < ray.tfar ))) return;
          
          const float uu = curve2D.t0; // FIXME: correct u range for subdivided segments
          
          Vec3fa P,T; curve2D.eval(0.5f,P,T); 
          if (T == Vec3fa(zero)) return; // { valid[i] = 0; goto retry; } // ignore denormalized curves */
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER) && 0
          Geometry* geometry = ((Scene*)geom)->get(geomID);
          if (!likely(geometry->hasIntersectionFilter1())) 
          {
            runIntersectionFilter1(geometry,ray,uu,0.0f,t,T,geomID,primID);
            return;
          }
#endif
          ray.u = uu;
          ray.v = 0.0f;
          ray.tfar = t;
          ray.Ng = T;
          ray.geomID = geomID;
          ray.primID = primID;
        }
        else
        {
          BezierCurve3D left,right;
          curve2D.subdivide(left,right);
          intersect_recursive(max_radius,ray,left,geom,geomID,primID);
          intersect_recursive(max_radius,ray,right,geom,geomID,primID);
        }
      }
      
      static __forceinline void intersect_recursive(const Precalculations& pre, Ray& ray, const Bezier1i& curve_in, const void* geom)
      {
        /* load bezier curve control points */
        STAT3(normal.trav_prims,1,1,1);
        const Vec3fa& v0 = curve_in.p[0];
        const Vec3fa& v1 = curve_in.p[1];
        const Vec3fa& v2 = curve_in.p[2];
        const Vec3fa& v3 = curve_in.p[3];
        
        /* transform control points into ray space */
        Vec3fa w0 = xfmVector(pre.ray_space,v0-ray.org); w0.w = v0.w;
        Vec3fa w1 = xfmVector(pre.ray_space,v1-ray.org); w1.w = v1.w;
        Vec3fa w2 = xfmVector(pre.ray_space,v2-ray.org); w2.w = v2.w;
        Vec3fa w3 = xfmVector(pre.ray_space,v3-ray.org); w3.w = v3.w;
        BezierCurve3D curve2D(w0,w1,w2,w3,0.0f,1.0f,3);
        
        const float max_radius = max(w0.w,w1.w,w2.w,w3.w);
        intersect_recursive(max_radius,ray,curve2D,geom,curve_in.geomID,curve_in.primID);
      }
      
      
      static __forceinline void intersect(const Precalculations& pre, Ray& ray, const Bezier1i& curve_in, const void* geom)
      {
        /* load bezier curve control points */
        STAT3(normal.trav_prims,1,1,1);
        const Vec3fa& v0 = curve_in.p[0];
        const Vec3fa& v1 = curve_in.p[1];
        const Vec3fa& v2 = curve_in.p[2];
        const Vec3fa& v3 = curve_in.p[3];
        
        //if (!intersectCylinder(ray,v0,v1,v2,v3)) return;
        //if (!intersectBoxes(ray,v0,v1,v2,v3)) return;
        
        /* transform control points into ray space */
        Vec3fa w0 = xfmVector(pre.ray_space,v0-ray.org); w0.w = v0.w;
        Vec3fa w1 = xfmVector(pre.ray_space,v1-ray.org); w1.w = v1.w;
        Vec3fa w2 = xfmVector(pre.ray_space,v2-ray.org); w2.w = v2.w;
        Vec3fa w3 = xfmVector(pre.ray_space,v3-ray.org); w3.w = v3.w;
        BezierCurve3D curve2D(w0,w1,w2,w3,0.0f,1.0f,4);
        
        /* subdivide 3 levels at once */ 
        const avx4f p0 = curve2D.eval(coeff0[0],coeff0[1],coeff0[2],coeff0[3]);
        const avx4f p1 = curve2D.eval(coeff1[0],coeff1[1],coeff1[2],coeff1[3]); // FIXME: can be calculated from p0 by shifting
        //const avx4f p1(shift_left1(p0.x,w3.x),shift_left1(p0.y,w3.y),shift_left1(p0.z,w3.z),shift_left1(p0.w,w3.w));
        
        /* approximative intersection with cone */
        const avx4f v = p1-p0;
        const avx4f w = -p0;
        const avxf d0 = w.x*v.x + w.y*v.y;
        const avxf d1 = v.x*v.x + v.y*v.y;
        const avxf u = clamp(d0*rcp(d1),avxf(zero),avxf(one));
        const avx4f p = p0 + u*v;
        const avxf t = p.z;
        const avxf d2 = p.x*p.x + p.y*p.y; 
        const avxf r = p.w; //max(p.w,ray.org.w+ray.dir.w*t);
        const avxf r2 = r*r;
        avxb valid = d2 <= r2 & avxf(ray.tnear) < t & t < avxf(ray.tfar);
      retry:
        if (unlikely(none(valid))) return;
        const float one_over_8 = 1.0f/8.0f;
        size_t i = select_min(valid,t);
        STAT3(normal.trav_prim_hits,1,1,1);
        
        /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
        int geomID = curve_in.geomID;
        Geometry* geometry = ((Scene*)geom)->get(geomID);
        if (!likely(geometry->hasIntersectionFilter1())) 
        {
#endif
          /* update hit information */
          const float uu = (float(i)+u[i])*one_over_8; // FIXME: correct u range for subdivided segments
          const BezierCurve3D curve3D(v0,v1,v2,v3,0.0f,1.0f,0);
          Vec3fa P,T; curve3D.eval(uu,P,T);
          if (T == Vec3fa(zero)) { valid[i] = 0; goto retry; } // ignore denormalized curves
          ray.u = uu;
          ray.v = 0.0f;
          ray.tfar = t[i];
          ray.Ng = T;
          ray.geomID = curve_in.geomID;
          ray.primID = curve_in.primID;
#if defined(RTCORE_INTERSECTION_FILTER)
          return;
        }
        
        while (true) 
        {
          const float uu = (float(i)+u[i])*one_over_8;
          const BezierCurve3D curve3D(v0,v1,v2,v3,0.0f,1.0f,0);
          Vec3fa P,T; curve3D.eval(uu,P,T);
          if (T != Vec3fa(zero))
            if (runIntersectionFilter1(geometry,ray,uu,0.0f,t[i],T,geomID,curve_in.primID)) return;
          valid[i] = 0;
          if (none(valid)) return;
          i = select_min(valid,t);
          STAT3(normal.trav_prim_hits,1,1,1);
        }
#endif
      }
      
      static __forceinline void intersect(const Precalculations& pre, Ray& ray, const Bezier1i* curves, size_t num, void* geom)
      {
        for (size_t i=0; i<num; i++)
        {
//#if defined(PRE_SUBDIVISION_HACK)
          if (unlikely(any(pre.mbox.ids == curves[i].primID))) continue; // FIXME: works only for single hair set
//#endif
          
#if USE_RECURSIVE_INTERSECTION == 0
          intersect(pre,ray,curves[i],geom);
#else
          intersect_recursive(pre,ray,curves[i],geom);
#endif
          
//#if defined(PRE_SUBDIVISION_HACK)
          *(unsigned int*)&pre.mbox.ids[pre.mbox.index] = curves[i].primID; // FIXME: works only for single hair set
          *(unsigned int*)&pre.mbox.index = (pre.mbox.index + 1 ) % 8;
//#endif
        }
      }
      
      static __forceinline bool occluded(const Precalculations& pre, Ray& ray, const Bezier1i& curve_in, const void* geom) 
      {
        /* load bezier curve control points */
        STAT3(shadow.trav_prims,1,1,1);
        const Vec3fa v0 = curve_in.p[0];
        const Vec3fa v1 = curve_in.p[1];
        const Vec3fa v2 = curve_in.p[2];
        const Vec3fa v3 = curve_in.p[3];
        
        /* transform control points into ray space */
        Vec3fa w0 = xfmVector(pre.ray_space,v0-ray.org); w0.w = v0.w;
        Vec3fa w1 = xfmVector(pre.ray_space,v1-ray.org); w1.w = v1.w;
        Vec3fa w2 = xfmVector(pre.ray_space,v2-ray.org); w2.w = v2.w;
        Vec3fa w3 = xfmVector(pre.ray_space,v3-ray.org); w3.w = v3.w;
        BezierCurve3D curve2D(w0,w1,w2,w3,0.0f,1.0f,4);
        
        /* subdivide 3 levels at once */ 
        const avx4f p0 = curve2D.eval(coeff0[0],coeff0[1],coeff0[2],coeff0[3]);
        const avx4f p1 = curve2D.eval(coeff1[0],coeff1[1],coeff1[2],coeff1[3]);
        
        /* approximative intersection with cone */
        const avx4f v = p1-p0;
        const avx4f w = -p0;
        const avxf d0 = w.x*v.x + w.y*v.y;
        const avxf d1 = v.x*v.x + v.y*v.y;
        const avxf u = clamp(d0*rcp(d1),avxf(zero),avxf(one));
        const avx4f p = p0 + u*v;
        const avxf t = p.z;
        const avxf d2 = p.x*p.x + p.y*p.y; 
        const avxf r = p.w; //max(p.w,ray.org.w+ray.dir.w*t);
        const avxf r2 = r*r;
        avxb valid = d2 <= r2 & avxf(ray.tnear) < t & t < avxf(ray.tfar);
        if (none(valid)) return false;
        STAT3(shadow.trav_prim_hits,1,1,1);
        
        /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER)
        
        size_t i = select_min(valid,t);
        int geomID = curve_in.geomID;
        Geometry* geometry = ((Scene*)geom)->get(geomID);
        if (likely(!geometry->hasOcclusionFilter1())) return true;
        const float one_over_8 = 1.0f/8.0f;
        
        while (true) 
        {
          /* calculate hit information */
          const float uu = (float(i)+u[i])*one_over_8;
          const BezierCurve3D curve3D(v0,v1,v2,v3,0.0f,1.0f,0);
          Vec3fa P,T; curve3D.eval(uu,P,T);
          if (T != Vec3fa(zero))
            if (runOcclusionFilter1(geometry,ray,uu,0.0f,t[i],T,geomID,curve_in.primID)) break;
          valid[i] = 0;
          if (none(valid)) return false;
          i = select_min(valid,t);
          STAT3(shadow.trav_prim_hits,1,1,1);
        }
#endif
        return true;
      }
      
      
      static bool occluded_recursive(const float max_radius, Ray& ray, const BezierCurve3D &curve2D, const void* geom, int geomID, int primID)
      {
        const BBox3fa bounds = curve2D.bounds();
        
#if 1
        if (bounds.lower.x >  max_radius ||
            bounds.upper.x < -max_radius ||
            bounds.lower.y >  max_radius ||
            bounds.upper.y < -max_radius) return false;
#endif
        
        if (curve2D.depth == 0)
        {
          const Vec3fa &p0 = curve2D.v0;
          const Vec3fa &p1 = curve2D.v3;
          
          /* approximative intersection with cone */
          const Vec3fa v = p1-p0;
          const Vec3fa w = -p0;
          const float d0 = w.x*v.x + w.y*v.y;
          const float d1 = v.x*v.x + v.y*v.y;
          const float u = clamp(d0*rcp(d1),0.0f,1.0f);
          const Vec3fa p = p0 + u*v;
          const float t = p.z;
          const float d2 = p.x*p.x + p.y*p.y; 
          const float r = p.w; //max(p.w,ray.org.w+ray.dir.w*t);
          const float r2 = r*r;
          
          STAT3(normal.trav_prim_hits,1,1,1);
          
          if (unlikely(!( d2 <= r2 & ray.tnear < t & t < ray.tfar ))) return false;
          
          const float uu = curve2D.t0; // FIXME: correct u range for subdivided segments
          
          Vec3fa P,T; curve2D.eval(0.5f,P,T); 
          if (T == Vec3fa(zero)) return false; // { valid[i] = 0; goto retry; } // ignore denormalized curves */
          
          /* intersection filter test */
#if defined(RTCORE_INTERSECTION_FILTER) && 0
          Geometry* geometry = ((Scene*)geom)->get(geomID);
          if (!likely(geometry->hasOcclusionFilter1())) 
          {
            return runOcclusionFilter1(geometry,ray,uu,0.0f,t,T,geomID,primID);
          }
#endif
          return true;
        }
        else
        {
          BezierCurve3D left,right;
          curve2D.subdivide(left,right);
          return                                                        \
            occluded_recursive(max_radius,ray,left,geom,geomID,primID) | 
            occluded_recursive(max_radius,ray,right,geom,geomID,primID);
        }
      }
      
      
      static __forceinline bool occluded_recursive(const Precalculations& pre, Ray& ray, const Bezier1i& curve_in, const void* geom)
      {
        /* load bezier curve control points */
        STAT3(normal.trav_prims,1,1,1);
        const Vec3fa& v0 = curve_in.p[0];
        const Vec3fa& v1 = curve_in.p[1];
        const Vec3fa& v2 = curve_in.p[2];
        const Vec3fa& v3 = curve_in.p[3];
        
        /* transform control points into ray space */
        Vec3fa w0 = xfmVector(pre.ray_space,v0-ray.org); w0.w = v0.w;
        Vec3fa w1 = xfmVector(pre.ray_space,v1-ray.org); w1.w = v1.w;
        Vec3fa w2 = xfmVector(pre.ray_space,v2-ray.org); w2.w = v2.w;
        Vec3fa w3 = xfmVector(pre.ray_space,v3-ray.org); w3.w = v3.w;
        BezierCurve3D curve2D(w0,w1,w2,w3,0.0f,1.0f,3);
        
        const float max_radius = max(w0.w,w1.w,w2.w,w3.w);
        return occluded_recursive(max_radius,ray,curve2D,geom,curve_in.geomID,curve_in.primID);
      }
      
      static __forceinline bool occluded(const Precalculations& pre, Ray& ray, const Bezier1i* curves, size_t num, void* geom) 
      {
        for (size_t i=0; i<num; i++) 
#if USE_RECURSIVE_INTERSECTION == 0
          if (occluded(pre,ray,curves[i],geom))
#else
            if (occluded_recursive(pre,ray,curves[i],geom))
#endif
              return true;
        
        return false;
      }
    };
  }
}
