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

#include "common/default.h"
#include "common/scene.h"
#include "builders/primrefblock.h"
#include <typeinfo>
#include "common/ray.h"

#if defined(__SSE__)
#include "common/ray4.h"
#endif

#if defined(__AVX__)
#include "common/ray8.h"
#endif

namespace embree
{
  struct PrimitiveType
  {
    /*! constructs the primitive type */
    PrimitiveType (const std::string& name, size_t bytes, size_t blockSize, bool needVertices, int intCost) 
      : name(name), bytes(bytes), blockSize(blockSize), needVertices(needVertices), intCost(intCost) {}

    /*! Computes the number of blocks required to store a number of triangles. */
    virtual size_t blocks(size_t x) const = 0; // FIXME: are these still required

    /*! Returns the number of stored primitives in a block. */
    virtual size_t size(const char* This) const = 0;

    /*! Returns a hash number for the leaf. */
    virtual size_t hash(const char* This, size_t num) const { return 0; }

    /*! Updates all primitives stored in a leaf */
    virtual BBox3fa update(char* prim, size_t num, void* geom) const { return BBox3fa(empty); } // FIXME: remove

    /*! Updates all primitives stored in a leaf */
    virtual std::pair<BBox3fa,BBox3fa> update2(char* prim, size_t num, void* geom) const { return std::pair<BBox3fa,BBox3fa>(empty,empty); } // FIXME: remove

  public:
    std::string name;       //!< name of this primitive type
    size_t bytes;           //!< number of bytes of the triangle data
    size_t blockSize;       //!< block size
    bool   needVertices;    //!< determines if we need the vertex array
    int    intCost;         //!< cost of one ray/primitive intersection
  };

  //template<typename Primitive1, typename Primitive2>
    struct PrimitiveType2 : public PrimitiveType
  {
    PrimitiveType2 () 
      // : PrimitiveType(std::string(typeid(Primitive1).name()) + ", " + std::string(typeid(Primitive2).name()), 0, 0, 0, 0) {}
      : PrimitiveType("unknown", 0, 0, 0, 0) {}

    size_t blocks(size_t x) const { return x; }
    size_t size(const char* This) const { return 0; } // FIXME

    static PrimitiveType2 type;
  };

    namespace isa // FIMXE: move to separate file
  {
    template<typename Intersector>
      struct ListIntersector1
      {
        typedef typename Intersector::Primitive Primitive;
        typedef typename Intersector::Precalculations Precalculations;
        
        static __forceinline void intersect(Precalculations& pre, Ray& ray, const Primitive* prim, size_t num, Scene* scene, size_t& lazy_node)
        {
          while (true) {
            Intersector::intersect(pre,ray,*prim,scene);
            if (prim->last()) break;
            prim++;
          }
        }
        
        static __forceinline bool occluded(Precalculations& pre, Ray& ray, const Primitive* prim, size_t num, Scene* scene, size_t& lazy_node) 
        {
          while (true) {
            if (Intersector::occluded(pre,ray,*prim,scene))
              return true;
            if (prim->last()) break;
            prim++;
          }
          return false;
        }
      };
    
    template<typename Intersector>
      struct ArrayIntersector1
      {
        typedef typename Intersector::Primitive Primitive;
        typedef typename Intersector::Precalculations Precalculations;
        
        static __forceinline void intersect(Precalculations& pre, Ray& ray, const Primitive* prim, size_t num, Scene* scene, size_t& lazy_node)
        {
          for (size_t i=0; i<num; i++)
            Intersector::intersect(pre,ray,prim[i],scene);
        }
        
        static __forceinline bool occluded(Precalculations& pre, Ray& ray, const Primitive* prim, size_t num, Scene* scene, size_t& lazy_node) 
        {
          for (size_t i=0; i<num; i++) {
            if (Intersector::occluded(pre,ray,prim[i],scene))
              return true;
          }
          return false;
        }
      };
    
    template<typename Intersector1, typename Intersector2>
      struct Switch2Intersector1
      {
        typedef typename Intersector1::Primitive Primitive1;
        typedef typename Intersector2::Primitive Primitive2;
        typedef void Primitive;
        
        struct Precalculations 
        {
          __forceinline Precalculations (const Ray& ray) 
            : pre1(ray), pre2(ray) {}
          
          typename Intersector1::Precalculations pre1;
          typename Intersector2::Precalculations pre2;
        };
        
        static __forceinline void intersect(Precalculations& pre, Ray& ray, Primitive* prim, size_t ty, Scene* scene, size_t& lazy_node)
        {
          if (ty == 0) Intersector1::intersect(pre.pre1,ray,*(Primitive1*)prim,scene,lazy_node);
          else         Intersector2::intersect(pre.pre2,ray,*(Primitive2*)prim,scene,lazy_node);
        }
        
        static __forceinline bool occluded(Precalculations& pre, Ray& ray, Primitive* prim, size_t ty, Scene* scene, size_t& lazy_node) 
        {
          if (ty == 0) return Intersector1::occluded(pre.pre1,ray,*(Primitive1*)prim,scene,lazy_node);
          else         return Intersector2::occluded(pre.pre2,ray,*(Primitive2*)prim,scene,lazy_node);
        }
      };
    
#if defined __SSE__
    
    template<typename Intersector>
      struct ListIntersector4
      {
        typedef typename Intersector::Primitive Primitive;
        typedef typename Intersector::Precalculations Precalculations;
        
        static __forceinline void intersect(const sseb& valid, Precalculations& pre, Ray4& ray, const Primitive* prim, size_t num, Scene* scene)
        {
          while (true) {
            Intersector::intersect(valid,pre,ray,*prim,scene);
            if (prim->last()) break;
            prim++;
          }
        }
        
        static __forceinline sseb occluded(const sseb& valid, Precalculations& pre, Ray4& ray, const Primitive* prim, size_t num, Scene* scene) 
        {
          sseb valid0 = valid;
          while (true) {
            valid0 &= !Intersector::occluded(valid0,pre,ray,*prim,scene);
            if (none(valid0)) break;
            if (prim->last()) break;
            prim++;
          }
          return !valid0;
        }
      };
    
    template<typename Intersector>
      struct ArrayIntersector4
      {
        typedef typename Intersector::Primitive Primitive;
        typedef typename Intersector::Precalculations Precalculations;
        
        static __forceinline void intersect(const sseb& valid, Precalculations& pre, Ray4& ray, const Primitive* prim, size_t num, Scene* scene)
        {
          for (size_t i=0; i<num; i++) {
            Intersector::intersect(valid,pre,ray,prim[i],scene);
          }
        }
        
        static __forceinline sseb occluded(const sseb& valid, Precalculations& pre, Ray4& ray, const Primitive* prim, size_t num, Scene* scene) 
        {
          sseb valid0 = valid;
          for (size_t i=0; i<num; i++) {
            valid0 &= !Intersector::occluded(valid0,pre,ray,prim[i],scene);
            if (none(valid0)) break;
          }
          return !valid0;
        }
      };
    
    
    template<typename Intersector>
      struct ListIntersector4_1
      {
        typedef typename Intersector::Primitive Primitive;
        typedef typename Intersector::Precalculations Precalculations;
        
        static __forceinline void intersect(const sseb& valid, Precalculations& pre, Ray4& ray, const Primitive* prim, size_t num, Scene* scene)
        {
          while (true) {
            Intersector::intersect(valid,pre,ray,*prim,scene);
            if (prim->last()) break;
            prim++;
          }
        }
        
        static __forceinline sseb occluded(const sseb& valid, Precalculations& pre, Ray4& ray, const Primitive* prim, size_t num, Scene* scene) 
        {
          sseb valid0 = valid;
          while (true) {
            valid0 &= !Intersector::occluded(valid0,pre,ray,*prim,scene);
            if (none(valid0)) break;
            if (prim->last()) break;
            prim++;
          }
          return !valid0;
        }
        
        static __forceinline void intersect(Precalculations& pre, Ray4& ray, size_t k, const Primitive* prim, size_t num, Scene* scene)
        {
          while (true) {
            Intersector::intersect(pre,ray,k,*prim,scene);
            if (prim->last()) break;
            prim++;
          }
        }
        
        static __forceinline bool occluded(Precalculations& pre, Ray4& ray, size_t k, const Primitive* prim, size_t num, Scene* scene) 
        {
          while (true) {
            if (Intersector::occluded(pre,ray,k,*prim,scene))
              return true;
            if (prim->last()) break;
            prim++;
          }
          return false;
        }
      };
    
    template<typename Intersector>
      struct ArrayIntersector4_1
      {
        typedef typename Intersector::Primitive Primitive;
        typedef typename Intersector::Precalculations Precalculations;
        
        static __forceinline void intersect(const sseb& valid, Precalculations& pre, Ray4& ray, const Primitive* prim, size_t num, Scene* scene)
        {
          for (size_t i=0; i<num; i++) {
            Intersector::intersect(valid,pre,ray,prim[i],scene);
          }
        }
        
        static __forceinline sseb occluded(const sseb& valid, Precalculations& pre, Ray4& ray, const Primitive* prim, size_t num, Scene* scene) 
        {
          sseb valid0 = valid;
          for (size_t i=0; i<num; i++) {
            valid0 &= !Intersector::occluded(valid0,pre,ray,prim[i],scene);
            if (none(valid0)) break;
          }
          return !valid0;
        }
        
        static __forceinline void intersect(Precalculations& pre, Ray4& ray, size_t k, const Primitive* prim, size_t num, Scene* scene)
        {
          for (size_t i=0; i<num; i++) {
            Intersector::intersect(pre,ray,k,prim[i],scene);
          }
        }
        
        static __forceinline bool occluded(Precalculations& pre, Ray4& ray, size_t k, const Primitive* prim, size_t num, Scene* scene) 
        {
          for (size_t i=0; i<num; i++) {
            if (Intersector::occluded(pre,ray,k,prim[i],scene))
              return true;
          }
          return false;
        }
      };
    
#endif
    
#if defined __AVX__
    
    template<typename Intersector>
      struct ListIntersector8
      {
        typedef typename Intersector::Primitive Primitive;
        typedef typename Intersector::Precalculations Precalculations;
        
        static __forceinline void intersect(const avxb& valid, Precalculations& pre, Ray8& ray, const Primitive* prim, size_t num, Scene* scene)
        {
          while (true) {
            Intersector::intersect(valid,pre,ray,*prim,scene);
            if (prim->last()) break;
            prim++;
          }
        }
        
        static __forceinline avxb occluded(const avxb& valid, Precalculations& pre, Ray8& ray, const Primitive* prim, size_t num, Scene* scene) 
        {
          avxb valid0 = valid;
          while (true) {
            valid0 &= !Intersector::occluded(valid0,pre,ray,*prim,scene);
            if (none(valid0)) break;
            if (prim->last()) break;
            prim++;
          }
          return !valid0;
        }
      };
    
    template<typename Intersector>
      struct ArrayIntersector8
      {
        typedef typename Intersector::Primitive Primitive;
        typedef typename Intersector::Precalculations Precalculations;
        
        static __forceinline void intersect(const avxb& valid, Precalculations& pre, Ray8& ray, const Primitive* prim, size_t num, Scene* scene)
        {
          for (size_t i=0; i<num; i++) {
            Intersector::intersect(valid,pre,ray,prim[i],scene);
          }
        }
        
        static __forceinline avxb occluded(const avxb& valid, Precalculations& pre, Ray8& ray, const Primitive* prim, size_t num, Scene* scene) 
        {
          avxb valid0 = valid;
          for (size_t i=0; i<num; i++) {
            valid0 &= !Intersector::occluded(valid0,pre,ray,prim[i],scene);
            if (none(valid0)) break;
          }
          return !valid0;
        }
      };
    
    
    template<typename Intersector>
      struct ListIntersector8_1
      {
        typedef typename Intersector::Primitive Primitive;
        typedef typename Intersector::Precalculations Precalculations;
        
        static __forceinline void intersect(const avxb& valid, Precalculations& pre, Ray8& ray, const Primitive* prim, size_t num, Scene* scene)
        {
          while (true) {
            Intersector::intersect(valid,pre,ray,*prim,scene);
            if (prim->last()) break;
            prim++;
          }
        }
        
        static __forceinline avxb occluded(const avxb& valid, Precalculations& pre, Ray8& ray, const Primitive* prim, size_t num, Scene* scene) 
        {
          avxb valid0 = valid;
          while (true) {
            valid0 &= !Intersector::occluded(valid0,pre,ray,*prim,scene);
            if (none(valid0)) break;
            if (prim->last()) break;
            prim++;
          }
          return !valid0;
        }
        
        static __forceinline void intersect(Precalculations& pre, Ray8& ray, size_t k, const Primitive* prim, size_t num, Scene* scene)
        {
          while (true) {
            Intersector::intersect(pre,ray,k,*prim,scene);
            if (prim->last()) break;
            prim++;
          }
        }
        
        static __forceinline bool occluded(Precalculations& pre, Ray8& ray, size_t k, const Primitive* prim, size_t num, Scene* scene) 
        {
          while (true) {
            if (Intersector::occluded(pre,ray,k,*prim,scene))
              return true;
            if (prim->last()) break;
            prim++;
          }
          return false;
        }
      };
    
    template<typename Intersector>
      struct ArrayIntersector8_1
      {
        typedef typename Intersector::Primitive Primitive;
        typedef typename Intersector::Precalculations Precalculations;
        
        static __forceinline void intersect(const avxb& valid, Precalculations& pre, Ray8& ray, const Primitive* prim, size_t num, Scene* scene)
        {
          for (size_t i=0; i<num; i++) {
            Intersector::intersect(valid,pre,ray,prim[i],scene);
          }
        }
        
        static __forceinline avxb occluded(const avxb& valid, Precalculations& pre, Ray8& ray, const Primitive* prim, size_t num, Scene* scene) 
        {
          avxb valid0 = valid;
          for (size_t i=0; i<num; i++) {
            valid0 &= !Intersector::occluded(valid0,pre,ray,prim[i],scene);
            if (none(valid0)) break;
          }
          return !valid0;
        }
        
        static __forceinline void intersect(Precalculations& pre, Ray8& ray, size_t k, const Primitive* prim, size_t num, Scene* scene)
        {
          for (size_t i=0; i<num; i++) {
            Intersector::intersect(pre,ray,k,prim[i],scene);
          }
        }
        
        static __forceinline bool occluded(Precalculations& pre, Ray8& ray, size_t k, const Primitive* prim, size_t num, Scene* scene) 
        {
          for (size_t i=0; i<num; i++) {
            if (Intersector::occluded(pre,ray,k,prim[i],scene))
              return true;
          }
          return false;
        }
      };
    
#endif
  }
}
