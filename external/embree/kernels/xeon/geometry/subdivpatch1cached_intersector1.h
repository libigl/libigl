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

#include "common/ray.h"
#include "common/scene_subdiv_mesh.h"
#include "geometry/filter.h"
#include "bvh4/bvh4.h"
#include "common/subdiv/tessellation.h"
#include "common/subdiv/tessellation_cache.h"
#include "geometry/subdivpatch1cached.h"

/* returns u,v based on individual triangles instead relative to original patch */
#define FORCE_TRIANGLE_UV 0

namespace embree
{
  namespace isa
  {
    class SubdivPatch1CachedIntersector1
    {
    public:
      typedef SubdivPatch1Cached Primitive;
      
      /*! Per thread tessellation cache */
      static __thread TessellationCache *thread_cache;
      
      /*! Creates per thread tessellation cache */
      static void createTessellationCache();
      
      /*! Precalculations for subdiv patch intersection */
      class Precalculations {
	  public:
        Vec3fa ray_rdir;
        Vec3fa ray_org_rdir;
        SubdivPatch1Cached *current_patch;
        SubdivPatch1Cached *hit_patch;
        TessellationCache *local_cache;      
        Ray &r;
        
        
        __forceinline Precalculations (Ray& ray) : r(ray) 
        {
          ray_rdir      = rcp_safe(ray.dir);
          ray_org_rdir  = ray.org*ray_rdir;
          current_patch = NULL;
          hit_patch     = NULL;

          /*! Initialize per thread tessellation cache */
          if (unlikely(!thread_cache))
            createTessellationCache();
          
          local_cache = thread_cache;
          assert(local_cache != NULL);
        }

          /*! Final per ray computations like smooth normal, patch u,v, etc. */        
        __forceinline ~Precalculations() 
        {
          if (unlikely(hit_patch != NULL))
          {

#if defined(RTCORE_RETURN_SUBDIV_NORMAL)
	    if (likely(!hit_patch->hasDisplacement()))
	      {		 
		Vec3fa normal = hit_patch->normal(r.u,r.v);
		r.Ng = normal;
	      }
#endif
	    const Vec2f uv0 = hit_patch->getUV(0);
	    const Vec2f uv1 = hit_patch->getUV(1);
	    const Vec2f uv2 = hit_patch->getUV(2);
	    const Vec2f uv3 = hit_patch->getUV(3);
	    
	    const float patch_u = bilinear_interpolate(uv0.x,uv1.x,uv2.x,uv3.x,r.u,r.v);
	    const float patch_v = bilinear_interpolate(uv0.y,uv1.y,uv2.y,uv3.y,r.u,r.v);

	    r.u      = patch_u;
	    r.v      = patch_v;
            r.geomID = hit_patch->geom;
            r.primID = hit_patch->prim;
          }
        }
        
      };
      
      
      /* intersect ray with Quad2x2 structure => 1 ray vs. 8 triangles */
      template<class M, class T>
        static __forceinline void intersect1_precise(Ray& ray,
                                                     const Quad2x2 &qquad,
                                                     const void* geom,
                                                     Precalculations &pre,
                                                     const size_t delta = 0)
      {
        const Vec3<T> v0_org = qquad.getVtx( 0, delta);
        const Vec3<T> v1_org = qquad.getVtx( 1, delta);
        const Vec3<T> v2_org = qquad.getVtx( 2, delta);
        
        const Vec3<T> O = ray.org;
        const Vec3<T> D = ray.dir;
        
        const Vec3<T> v0 = v0_org - O;
        const Vec3<T> v1 = v1_org - O;
        const Vec3<T> v2 = v2_org - O;
        
        const Vec3<T> e0 = v2 - v0;
        const Vec3<T> e1 = v0 - v1;	     
        const Vec3<T> e2 = v1 - v2;	     
        
        /* calculate geometry normal and denominator */
        const Vec3<T> Ng1 = cross(e1,e0);
        const Vec3<T> Ng = Ng1+Ng1;
        const T den = dot(Ng,D);
        const T absDen = abs(den);
        const T sgnDen = signmsk(den);
        
        M valid ( true );
        /* perform edge tests */
        const T U = dot(Vec3<T>(cross(v2+v0,e0)),D) ^ sgnDen;
        valid &= U >= 0.0f;
        if (likely(none(valid))) return;
        const T V = dot(Vec3<T>(cross(v0+v1,e1)),D) ^ sgnDen;
        valid &= V >= 0.0f;
        if (likely(none(valid))) return;
        const T W = dot(Vec3<T>(cross(v1+v2,e2)),D) ^ sgnDen;
        valid &= W >= 0.0f;
        if (likely(none(valid))) return;
        
        /* perform depth test */
        const T _t = dot(v0,Ng) ^ sgnDen;
        valid &= (_t >= absDen*ray.tnear) & (absDen*ray.tfar >= _t);
        if (unlikely(none(valid))) return;
        
        /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
        valid &= den > T(zero);
        if (unlikely(none(valid))) return;
#else
        valid &= den != T(zero);
        if (unlikely(none(valid))) return;
#endif
        
        /* calculate hit information */
        const T rcpAbsDen = rcp(absDen);
        const T u =  U*rcpAbsDen;
        const T v =  V*rcpAbsDen;
        const T t = _t*rcpAbsDen;
        
#if FORCE_TRIANGLE_UV == 0
        const Vec2<T> uv0 = qquad.getUV( 0, delta );
        const Vec2<T> uv1 = qquad.getUV( 1, delta );
        const Vec2<T> uv2 = qquad.getUV( 2, delta );
        
        const Vec2<T> uv = u * uv1 + v * uv2 + (1.0f-u-v) * uv0;
        
        const T u_final = uv[0];
        const T v_final = uv[1];
        
#else
        const T u_final = u;
        const T v_final = v;
#endif
        
        size_t i = select_min(valid,t);

        
        /* update hit information */
        pre.hit_patch = pre.current_patch;

        ray.u         = u_final[i];
        ray.v         = v_final[i];
        ray.tfar      = t[i];
	if (i % 2)
	  {
	    ray.Ng.x      = Ng.x[i];
	    ray.Ng.y      = Ng.y[i];
	    ray.Ng.z      = Ng.z[i];
	  }
	else
	  {
	    ray.Ng.x      = -Ng.x[i];
	    ray.Ng.y      = -Ng.y[i];
	    ray.Ng.z      = -Ng.z[i];	    
	  }
      };
      
      
      /*! intersect ray with Quad2x2 structure => 1 ray vs. 8 triangles */
      template<class M, class T>
        static __forceinline bool occluded1_precise(Ray& ray,
                                                    const Quad2x2 &qquad,
                                                    const void* geom,
						    const size_t delta = 0)
      {
        const Vec3<T> v0_org = qquad.getVtx( 0, delta );
        const Vec3<T> v1_org = qquad.getVtx( 1, delta );
        const Vec3<T> v2_org = qquad.getVtx( 2, delta );
        
        const Vec3<T> O = ray.org;
        const Vec3<T> D = ray.dir;
        
        const Vec3<T> v0 = v0_org - O;
        const Vec3<T> v1 = v1_org - O;
        const Vec3<T> v2 = v2_org - O;
        
        const Vec3<T> e0 = v2 - v0;
        const Vec3<T> e1 = v0 - v1;	     
        const Vec3<T> e2 = v1 - v2;	     
        
        /* calculate geometry normal and denominator */
        const Vec3<T> Ng1 = cross(e1,e0);
        const Vec3<T> Ng = Ng1+Ng1;
        const T den = dot(Ng,D);
        const T absDen = abs(den);
        const T sgnDen = signmsk(den);
        
        M valid ( true );
        /* perform edge tests */
        const T U = dot(Vec3<T>(cross(v2+v0,e0)),D) ^ sgnDen;
        valid &= U >= 0.0f;
        if (likely(none(valid))) return false;
        const T V = dot(Vec3<T>(cross(v0+v1,e1)),D) ^ sgnDen;
        valid &= V >= 0.0f;
        if (likely(none(valid))) return false;
        const T W = dot(Vec3<T>(cross(v1+v2,e2)),D) ^ sgnDen;
        valid &= W >= 0.0f;
        if (likely(none(valid))) return false;
        
        /* perform depth test */
        const T _t = dot(v0,Ng) ^ sgnDen;
        valid &= (_t >= absDen*ray.tnear) & (absDen*ray.tfar >= _t);
        if (unlikely(none(valid))) return false;
        
        /* perform backface culling */
#if defined(RTCORE_BACKFACE_CULLING)
        valid &= den > T(zero);
        if (unlikely(none(valid))) return false;
#else
        valid &= den != T(zero);
        if (unlikely(none(valid))) return false;
#endif
        return true;
      };
      

      /* eval grid over patch and stich edges when required */      
      static __forceinline void evalGrid(const SubdivPatch1Cached &patch,
                                         float *__restrict__ const grid_x,
                                         float *__restrict__ const grid_y,
                                         float *__restrict__ const grid_z,
                                         float *__restrict__ const grid_u,
                                         float *__restrict__ const grid_v,
                                         const SubdivMesh* const geom)
      {
        gridUVTessellator(patch.level,
                          patch.grid_u_res,
                          patch.grid_v_res,
                          grid_u,
                          grid_v);
        
        if (unlikely(patch.needsStiching()))
          stichUVGrid(patch.level,patch.grid_u_res,patch.grid_v_res,grid_u,grid_v);
        
        const Vec2f uv0 = patch.getUV(0);
        const Vec2f uv1 = patch.getUV(1);
        const Vec2f uv2 = patch.getUV(2);
        const Vec2f uv3 = patch.getUV(3);
        
#if defined(__AVX__)
        for (size_t i=0;i<patch.grid_size_simd_blocks;i++)
        {
          avxf uu = load8f(&grid_u[8*i]);
          avxf vv = load8f(&grid_v[8*i]);
          avx3f vtx = patch.eval8(uu,vv);
          
          
          if (unlikely(((SubdivMesh*)geom)->displFunc != NULL))
          {
            avx3f normal = patch.normal8(uu,vv);
            normal = normalize_safe(normal) ;
            
	    const avxf patch_uu = bilinear_interpolate(uv0.x,uv1.x,uv2.x,uv3.x,uu,vv);
	    const avxf patch_vv = bilinear_interpolate(uv0.y,uv1.y,uv2.y,uv3.y,uu,vv);
            
            ((SubdivMesh*)geom)->displFunc(((SubdivMesh*)geom)->userPtr,
                                           patch.geom,
                                           patch.prim,
                                           (const float*)&patch_uu,
                                           (const float*)&patch_vv,
                                           (const float*)&normal.x,
                                           (const float*)&normal.y,
                                           (const float*)&normal.z,
                                           (float*)&vtx.x,
                                           (float*)&vtx.y,
                                           (float*)&vtx.z,
                                           8);
          }
          
          *(avxf*)&grid_x[8*i] = vtx.x;
          *(avxf*)&grid_y[8*i] = vtx.y;
          *(avxf*)&grid_z[8*i] = vtx.z;        
          *(avxf*)&grid_u[8*i] = uu;
          *(avxf*)&grid_v[8*i] = vv;
        }
#else
        for (size_t i=0;i<patch.grid_size_simd_blocks*2;i++) // 4-wide blocks for SSE
        {
          ssef uu = load4f(&grid_u[4*i]);
          ssef vv = load4f(&grid_v[4*i]);
          sse3f vtx = patch.eval4(uu,vv);
          
          
          if (unlikely(((SubdivMesh*)geom)->displFunc != NULL))
          {
            sse3f normal = patch.normal4(uu,vv);
            normal = normalize_safe(normal);

	    const ssef patch_uu = bilinear_interpolate(uv0.x,uv1.x,uv2.x,uv3.x,uu,vv);
	    const ssef patch_vv = bilinear_interpolate(uv0.y,uv1.y,uv2.y,uv3.y,uu,vv);
            
            ((SubdivMesh*)geom)->displFunc(((SubdivMesh*)geom)->userPtr,
                                           patch.geom,
                                           patch.prim,
                                           (const float*)&patch_uu,
                                           (const float*)&patch_vv,
                                           (const float*)&normal.x,
                                           (const float*)&normal.y,
                                           (const float*)&normal.z,
                                           (float*)&vtx.x,
                                           (float*)&vtx.y,
                                           (float*)&vtx.z,
                                           4);
          }
          
          *(ssef*)&grid_x[4*i] = vtx.x;
          *(ssef*)&grid_y[4*i] = vtx.y;
          *(ssef*)&grid_z[4*i] = vtx.z;        
          *(ssef*)&grid_u[4*i] = uu;
          *(ssef*)&grid_v[4*i] = vv;
        }
        
#endif
      }
      
      
      /*! Returns BVH4 node reference for subtree over patch grid */
      static __forceinline size_t getSubtreeRootNode(TessellationCache *local_cache, const SubdivPatch1Cached* const subdiv_patch, const void* geom)
      {
        const unsigned int commitCounter = ((Scene*)geom)->commitCounter;
        TessellationCache::InputTagType tag = (TessellationCache::InputTagType)subdiv_patch;
        
        BVH4::NodeRef root = local_cache->lookup(tag,commitCounter);
        root.prefetch(0);
        if (unlikely(root == (size_t)-1))
        {
          subdiv_patch->prefetchData();
          const unsigned int blocks = subdiv_patch->grid_subtree_size_64b_blocks;
          
          TessellationCache::CacheTag &t = local_cache->request(tag,commitCounter,blocks);
          BVH4::Node* node = (BVH4::Node*)local_cache->getCacheMemoryPtr(t);
          prefetchL1(((float*)node + 0*16));
          prefetchL1(((float*)node + 1*16));
          prefetchL1(((float*)node + 2*16));
          prefetchL1(((float*)node + 3*16));
          
          size_t new_root = (size_t)buildSubdivPatchTree(*subdiv_patch,node,((Scene*)geom)->getSubdivMesh(subdiv_patch->geom));
          assert( new_root != BVH4::invalidNode);
          
          local_cache->updateRootRef(t,new_root);
          
          assert( (size_t)local_cache->getPtr() + (size_t)t.getRootRef() == new_root );
          
          return new_root;
        }
        
        
        return root;   
      }
      
      
      /*! Evaluates grid over patch and builds BVH4 tree over the grid. */
      static BVH4::NodeRef buildSubdivPatchTree(const SubdivPatch1Cached &patch,
                                                void *const lazymem,
                                                const SubdivMesh* const geom);
      
      /*! Create BVH4 tree over grid. */
      static BBox3fa createSubTree(BVH4::NodeRef &curNode,
                                   float *const lazymem,
                                   const SubdivPatch1Cached &patch,
                                   const float *const grid_x_array,
                                   const float *const grid_y_array,
                                   const float *const grid_z_array,
                                   const float *const grid_u_array,
                                   const float *const grid_v_array,
                                   const GridRange &range,
                                   unsigned int &localCounter,
                                   const SubdivMesh* const geom);
      
      
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      
      /*! Intersect a ray with the primitive. */
      static __forceinline void intersect(Precalculations& pre, Ray& ray, const Primitive* prim, size_t ty, const void* geom, size_t& lazy_node) 
      {
        STAT3(normal.trav_prims,1,1,1);
        
        if (likely(ty == 2))
        {
#if defined(__AVX__)
          intersect1_precise<avxb,avxf>( ray, *(Quad2x2*)prim, (SubdivMesh*)geom,pre);
#else
          intersect1_precise<sseb,ssef>( ray, *(Quad2x2*)prim, (SubdivMesh*)geom,pre,0);
          intersect1_precise<sseb,ssef>( ray, *(Quad2x2*)prim, (SubdivMesh*)geom,pre,6);
#endif
        }
        else 
        {
          lazy_node = getSubtreeRootNode(pre.local_cache, prim, geom);
          pre.current_patch = (SubdivPatch1Cached*)prim;
        }             
        
      }
      
      /*! Test if the ray is occluded by the primitive */
      static __forceinline bool occluded(Precalculations& pre, Ray& ray, const Primitive* prim, size_t ty, const void* geom, size_t& lazy_node) 
      {
        STAT3(shadow.trav_prims,1,1,1);
        
        if (likely(ty == 2))
        {
#if defined(__AVX__)
          return occluded1_precise<avxb,avxf>( ray, *(Quad2x2*)prim, (SubdivMesh*)geom);
#else
          if (occluded1_precise<sseb,ssef>( ray, *(Quad2x2*)prim, (SubdivMesh*)geom,0)) return true;
          if (occluded1_precise<sseb,ssef>( ray, *(Quad2x2*)prim, (SubdivMesh*)geom,6)) return true;
#endif
        }
        else 
        {
          lazy_node = getSubtreeRootNode(pre.local_cache, prim, geom);        
          pre.current_patch = (SubdivPatch1Cached*)prim;
        }             
        
        
        return false;
      }
      
      
    };
  }
}
