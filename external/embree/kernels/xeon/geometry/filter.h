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

#include "common/geometry.h"

#include "common/ray.h"

#if defined(__SSE__)
#include "common/ray4.h"
#endif

#if defined(__AVX__)
#include "common/ray8.h"
#endif

#if defined(__MIC__)
#include "common/ray16.h"
#endif

namespace embree
{
  namespace isa
  {
    __forceinline bool runIntersectionFilter1(const Geometry* const geometry, Ray& ray, 
                                              const float& u, const float& v, const float& t, const Vec3fa& Ng, const int geomID, const int primID)
    {
      /* temporarily update hit information */
      const float  ray_tfar = ray.tfar;
      const Vec3fa ray_Ng   = ray.Ng;
      const ssef   ray_uv_ids = *(ssef*)&ray.u;
      ray.u = u;
      ray.v = v;
      ray.tfar = t;
      ray.geomID = geomID;
      ray.primID = primID;
      ray.Ng = Ng;
      
      /* invoke filter function */
      AVX_ZERO_UPPER();
      geometry->intersectionFilter1(geometry->userPtr,(RTCRay&)ray);
      
      /* restore hit if filter not passed */
      if (unlikely(ray.geomID == -1)) 
      {
        ray.tfar = ray_tfar;
        ray.Ng = ray_Ng;
        *(ssef*)&ray.u = ray_uv_ids;
        return false;
      }
      return true;
    }
    
    __forceinline bool runOcclusionFilter1(const Geometry* const geometry, Ray& ray, 
                                           const float& u, const float& v, const float& t, const Vec3fa& Ng, const int geomID, const int primID)
    {
      /* temporarily update hit information */
      const float ray_tfar = ray.tfar;
      const int   ray_geomID = ray.geomID;
      ray.u = u;
      ray.v = v;
      ray.tfar = t;
      ray.geomID = geomID;
      ray.primID = primID;
      ray.Ng = Ng;
      
      /* invoke filter function */
      AVX_ZERO_UPPER();
      geometry->occlusionFilter1(geometry->userPtr,(RTCRay&)ray);
      
      /* restore hit if filter not passed */
      if (unlikely(ray.geomID == -1)) 
      {
        ray.tfar = ray_tfar;
        ray.geomID = ray_geomID;
        return false;
      }
      return true;
    }
    
    __forceinline sseb runIntersectionFilter4(const sseb& valid, const Geometry* const geometry, Ray4& ray, 
                                              const ssef& u, const ssef& v, const ssef& t, const sse3f& Ng, const int geomID, const int primID)
    {
      /* temporarily update hit information */
      const ssef ray_u = ray.u;           store4f(valid,&ray.u,u);
      const ssef ray_v = ray.v;           store4f(valid,&ray.v,v);
      const ssef ray_tfar = ray.tfar;     store4f(valid,&ray.tfar,t);
      const ssei ray_geomID = ray.geomID; store4i(valid,&ray.geomID,geomID);
      const ssei ray_primID = ray.primID; store4i(valid,&ray.primID,primID);
      const ssef ray_Ng_x = ray.Ng.x;     store4f(valid,&ray.Ng.x,Ng.x);
      const ssef ray_Ng_y = ray.Ng.y;     store4f(valid,&ray.Ng.y,Ng.y);
      const ssef ray_Ng_z = ray.Ng.z;     store4f(valid,&ray.Ng.z,Ng.z);
      
      /* invoke filter function */
      RTCFilterFunc4  filter4     = (RTCFilterFunc4)  geometry->intersectionFilter4;
      ISPCFilterFunc4 ispcFilter4 = (ISPCFilterFunc4) geometry->ispcIntersectionFilter4;
      AVX_ZERO_UPPER();
      if (ispcFilter4) ispcFilter4(geometry->userPtr,(RTCRay4&)ray,valid);
      else { const sseb valid_temp = valid; filter4(&valid_temp,geometry->userPtr,(RTCRay4&)ray); }
      const sseb valid_failed = valid & (ray.geomID == ssei(-1));
      const sseb valid_passed = valid & (ray.geomID != ssei(-1));
      
      /* restore hit if filter not passed */
      if (unlikely(any(valid_failed))) 
      {
        store4f(valid_failed,&ray.u,ray_u);
        store4f(valid_failed,&ray.v,ray_v);
        store4f(valid_failed,&ray.tfar,ray_tfar);
        store4i(valid_failed,&ray.geomID,ray_geomID);
        store4i(valid_failed,&ray.primID,ray_primID);
        store4f(valid_failed,&ray.Ng.x,ray_Ng_x);
        store4f(valid_failed,&ray.Ng.y,ray_Ng_y);
        store4f(valid_failed,&ray.Ng.z,ray_Ng_z);
      }
      return valid_passed;
    }
    
    __forceinline sseb runOcclusionFilter4(const sseb& valid, const Geometry* const geometry, Ray4& ray, 
                                           const ssef& u, const ssef& v, const ssef& t, const sse3f& Ng, const int geomID, const int primID)
    {
      /* temporarily update hit information */
      const ssef ray_tfar = ray.tfar; 
      const ssei ray_geomID = ray.geomID;
      store4f(valid,&ray.u,u);
      store4f(valid,&ray.v,v);
      store4f(valid,&ray.tfar,t);
      store4i(valid,&ray.geomID,geomID);
      store4i(valid,&ray.primID,primID);
      store4f(valid,&ray.Ng.x,Ng.x);
      store4f(valid,&ray.Ng.y,Ng.y);
      store4f(valid,&ray.Ng.z,Ng.z);
      
      /* invoke filter function */
      RTCFilterFunc4  filter4     = (RTCFilterFunc4)  geometry->occlusionFilter4;
      ISPCFilterFunc4 ispcFilter4 = (ISPCFilterFunc4) geometry->ispcOcclusionFilter4;
      AVX_ZERO_UPPER();
      if (ispcFilter4) ispcFilter4(geometry->userPtr,(RTCRay4&)ray,valid);
      else { const sseb valid_temp = valid; filter4(&valid_temp,geometry->userPtr,(RTCRay4&)ray); }
      const sseb valid_failed = valid & (ray.geomID == ssei(-1));
      const sseb valid_passed = valid & (ray.geomID != ssei(-1));
      
      /* restore hit if filter not passed */
      store4f(valid_failed,&ray.tfar,ray_tfar);
      store4i(valid_failed,&ray.geomID,ray_geomID);
      return valid_passed;
    }
    
    __forceinline bool runIntersectionFilter4(const Geometry* const geometry, Ray4& ray, const size_t k,
                                              const float& u, const float& v, const float& t, const Vec3fa& Ng, const int geomID, const int primID)
    {
      /* temporarily update hit information */
      const ssef ray_u = ray.u;           ray.u[k] = u;
      const ssef ray_v = ray.v;           ray.v[k] = v;
      const ssef ray_tfar = ray.tfar;     ray.tfar[k] = t;
      const ssei ray_geomID = ray.geomID; ray.geomID[k] = geomID;
      const ssei ray_primID = ray.primID; ray.primID[k] = primID;
      const ssef ray_Ng_x = ray.Ng.x;     ray.Ng.x[k] = Ng.x;
      const ssef ray_Ng_y = ray.Ng.y;     ray.Ng.y[k] = Ng.y;
      const ssef ray_Ng_z = ray.Ng.z;     ray.Ng.z[k] = Ng.z;
      
      /* invoke filter function */
      const sseb valid(1 << k);
      RTCFilterFunc4  filter4     = (RTCFilterFunc4)  geometry->intersectionFilter4;
      ISPCFilterFunc4 ispcFilter4 = (ISPCFilterFunc4) geometry->ispcIntersectionFilter4;
      AVX_ZERO_UPPER();
      if (ispcFilter4) ispcFilter4(geometry->userPtr,(RTCRay4&)ray,valid);
      else { const sseb valid_temp = valid; filter4(&valid_temp,geometry->userPtr,(RTCRay4&)ray); }
      const bool passed = ray.geomID[k] != -1;
      
      /* restore hit if filter not passed */
      if (unlikely(!passed)) {
        store4f(&ray.u,ray_u);
        store4f(&ray.v,ray_v);
        store4f(&ray.tfar,ray_tfar);
        store4i(&ray.geomID,ray_geomID);
        store4i(&ray.primID,ray_primID);
        store4f(&ray.Ng.x,ray_Ng_x);
        store4f(&ray.Ng.y,ray_Ng_y);
        store4f(&ray.Ng.z,ray_Ng_z);
      }
      return passed;
    }
    
    __forceinline bool runOcclusionFilter4(const Geometry* const geometry, Ray4& ray, const size_t k,
                                           const float& u, const float& v, const float& t, const Vec3fa& Ng, const int geomID, const int primID)
    {
      /* temporarily update hit information */
      const ssef ray_tfar = ray.tfar; 
      const ssei ray_geomID = ray.geomID;
      ray.u[k] = u;
      ray.v[k] = v;
      ray.tfar[k] = t;
      ray.geomID[k] = geomID;
      ray.primID[k] = primID;
      ray.Ng.x[k] = Ng.x;
      ray.Ng.y[k] = Ng.y;
      ray.Ng.z[k] = Ng.z;
      
      /* invoke filter function */
      const sseb valid(1 << k);
      RTCFilterFunc4  filter4     = (RTCFilterFunc4)  geometry->occlusionFilter4;
      ISPCFilterFunc4 ispcFilter4 = (ISPCFilterFunc4) geometry->ispcOcclusionFilter4;
      AVX_ZERO_UPPER();
      if (ispcFilter4) ispcFilter4(geometry->userPtr,(RTCRay4&)ray,valid);
      else { const sseb valid_temp = valid; filter4(&valid_temp,geometry->userPtr,(RTCRay4&)ray); }
      const bool passed = ray.geomID[k] != -1;
      
      /* restore hit if filter not passed */
      if (unlikely(!passed)) {
        store4f(&ray.tfar,ray_tfar);
        store4i(&ray.geomID,ray_geomID);
      }
      return passed;
    }
    
#if defined(__AVX__)
    __forceinline avxb runIntersectionFilter8(const avxb& valid, const Geometry* const geometry, Ray8& ray, 
                                              const avxf& u, const avxf& v, const avxf& t, const avx3f& Ng, const int geomID, const int primID)
    {
      /* temporarily update hit information */
      const avxf ray_u = ray.u;           store8f(valid,&ray.u,u);
      const avxf ray_v = ray.v;           store8f(valid,&ray.v,v);
      const avxf ray_tfar = ray.tfar;     store8f(valid,&ray.tfar,t);
      const avxi ray_geomID = ray.geomID; store8i(valid,&ray.geomID,geomID);
      const avxi ray_primID = ray.primID; store8i(valid,&ray.primID,primID);
      const avxf ray_Ng_x = ray.Ng.x;     store8f(valid,&ray.Ng.x,Ng.x);
      const avxf ray_Ng_y = ray.Ng.y;     store8f(valid,&ray.Ng.y,Ng.y);
      const avxf ray_Ng_z = ray.Ng.z;     store8f(valid,&ray.Ng.z,Ng.z);
      
      /* invoke filter function */
      RTCFilterFunc8  filter8     = (RTCFilterFunc8)  geometry->intersectionFilter8;
      ISPCFilterFunc8 ispcFilter8 = (ISPCFilterFunc8) geometry->ispcIntersectionFilter8;
      if (ispcFilter8) ispcFilter8(geometry->userPtr,(RTCRay8&)ray,valid);
      else { const avxb valid_temp = valid; filter8(&valid_temp,geometry->userPtr,(RTCRay8&)ray); }
      const avxb valid_failed = valid & (ray.geomID == avxi(-1));
      const avxb valid_passed = valid & (ray.geomID != avxi(-1));
      
      /* restore hit if filter not passed */
      if (unlikely(any(valid_failed))) 
      {
        store8f(valid_failed,&ray.u,ray_u);
        store8f(valid_failed,&ray.v,ray_v);
        store8f(valid_failed,&ray.tfar,ray_tfar);
        store8i(valid_failed,&ray.geomID,ray_geomID);
        store8i(valid_failed,&ray.primID,ray_primID);
        store8f(valid_failed,&ray.Ng.x,ray_Ng_x);
        store8f(valid_failed,&ray.Ng.y,ray_Ng_y);
        store8f(valid_failed,&ray.Ng.z,ray_Ng_z);
      }
      return valid_passed;
    }
    
    __forceinline avxb runOcclusionFilter8(const avxb& valid, const Geometry* const geometry, Ray8& ray, 
                                           const avxf& u, const avxf& v, const avxf& t, const avx3f& Ng, const int geomID, const int primID)
    {
      /* temporarily update hit information */
      const avxf ray_tfar = ray.tfar; 
      const avxi ray_geomID = ray.geomID;
      store8f(valid,&ray.u,u);
      store8f(valid,&ray.v,v);
      store8f(valid,&ray.tfar,t);
      store8i(valid,&ray.geomID,geomID);
      store8i(valid,&ray.primID,primID);
      store8f(valid,&ray.Ng.x,Ng.x);
      store8f(valid,&ray.Ng.y,Ng.y);
      store8f(valid,&ray.Ng.z,Ng.z);
      
      /* invoke filter function */
      RTCFilterFunc8  filter8     = (RTCFilterFunc8)  geometry->occlusionFilter8;
      ISPCFilterFunc8 ispcFilter8 = (ISPCFilterFunc8) geometry->ispcOcclusionFilter8;
      if (ispcFilter8) ispcFilter8(geometry->userPtr,(RTCRay8&)ray,valid);
      else { const avxb valid_temp = valid; filter8(&valid_temp,geometry->userPtr,(RTCRay8&)ray); }
      const avxb valid_failed = valid & (ray.geomID == avxi(-1));
      const avxb valid_passed = valid & (ray.geomID != avxi(-1));
      
      /* restore hit if filter not passed */
      store8f(valid_failed,&ray.tfar,ray_tfar);
      store8i(valid_failed,&ray.geomID,ray_geomID);
      return valid_passed;
    }
    
    __forceinline bool runIntersectionFilter8(const Geometry* const geometry, Ray8& ray, const size_t k,
                                              const float& u, const float& v, const float& t, const Vec3fa& Ng, const int geomID, const int primID)
    {
      /* temporarily update hit information */
      const avxf ray_u = ray.u;           ray.u[k] = u;
      const avxf ray_v = ray.v;           ray.v[k] = v;
      const avxf ray_tfar = ray.tfar;     ray.tfar[k] = t;
      const avxi ray_geomID = ray.geomID; ray.geomID[k] = geomID;
      const avxi ray_primID = ray.primID; ray.primID[k] = primID;
      const avxf ray_Ng_x = ray.Ng.x;     ray.Ng.x[k] = Ng.x;
      const avxf ray_Ng_y = ray.Ng.y;     ray.Ng.y[k] = Ng.y;
      const avxf ray_Ng_z = ray.Ng.z;     ray.Ng.z[k] = Ng.z;
      
      /* invoke filter function */
      const avxb valid(1 << k);
      RTCFilterFunc8  filter8     = (RTCFilterFunc8)  geometry->intersectionFilter8;
      ISPCFilterFunc8 ispcFilter8 = (ISPCFilterFunc8) geometry->ispcIntersectionFilter8;
      if (ispcFilter8) ispcFilter8(geometry->userPtr,(RTCRay8&)ray,valid);
      else filter8(&valid,geometry->userPtr,(RTCRay8&)ray);
      const bool passed = ray.geomID[k] != -1;
      
      /* restore hit if filter not passed */
      if (unlikely(!passed)) {
        store8f(&ray.u,ray_u);
        store8f(&ray.v,ray_v);
        store8f(&ray.tfar,ray_tfar);
        store8i(&ray.geomID,ray_geomID);
        store8i(&ray.primID,ray_primID);
        store8f(&ray.Ng.x,ray_Ng_x);
        store8f(&ray.Ng.y,ray_Ng_y);
        store8f(&ray.Ng.z,ray_Ng_z);
      }
      return passed;
    }
    
    __forceinline bool runOcclusionFilter8(const Geometry* const geometry, Ray8& ray, const size_t k,
                                           const float& u, const float& v, const float& t, const Vec3fa& Ng, const int geomID, const int primID)
    {
      /* temporarily update hit information */
      const avxf ray_tfar = ray.tfar; 
      const avxi ray_geomID = ray.geomID;
      ray.u[k] = u;
      ray.v[k] = v;
      ray.tfar[k] = t;
      ray.geomID[k] = geomID;
      ray.primID[k] = primID;
      ray.Ng.x[k] = Ng.x;
      ray.Ng.y[k] = Ng.y;
      ray.Ng.z[k] = Ng.z;
      
      /* invoke filter function */
      const avxb valid(1 << k);
      RTCFilterFunc8  filter8     = (RTCFilterFunc8)  geometry->occlusionFilter8;
      ISPCFilterFunc8 ispcFilter8 = (ISPCFilterFunc8) geometry->ispcOcclusionFilter8;
      if (ispcFilter8) ispcFilter8(geometry->userPtr,(RTCRay8&)ray,valid);
      else filter8(&valid,geometry->userPtr,(RTCRay8&)ray);
      const bool passed = ray.geomID[k] != -1;
      
      /* restore hit if filter not passed */
      if (unlikely(!passed)) {
        store8f(&ray.tfar,ray_tfar);
        store8i(&ray.geomID,ray_geomID);
      }
      return passed;
    }
    
#endif
  }
}
