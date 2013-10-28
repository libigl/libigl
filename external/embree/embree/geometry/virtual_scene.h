// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
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

#ifndef __EMBREE_VIRTUAL_SCENE_H__
#define __EMBREE_VIRTUAL_SCENE_H__

#include "geometry.h"
#include "../common/registry_accel.h"
#include "../common/registry_builder.h"
#include "../common/registry_intersector.h"

namespace embree
{
  struct VirtualScene : public RTCGeometry
  {
    /*! acceleration structure registry */
    static AccelRegistry accels; 

    /*! builder registry */
    static BuilderRegistry builders; 

    /*! intersector registrys */
    static IntersectorRegistry<Intersector1> intersectors1;
    
#if defined(__SSE__)
    static IntersectorRegistry<Intersector4> intersectors4;
#endif
    
#if defined(__AVX__)
    static IntersectorRegistry<Intersector8> intersectors8;
#endif
    
#if defined(__MIC__)
    static IntersectorRegistry<Intersector16> intersectors16;
#endif
    
    /*! Top level objects. */
    struct Object
    {
      Object () 
      : id0(0x7FFFFFFF), id1(0), mask(-1), hasTransform(false), 
        localBounds(empty), worldBounds(empty), 
        local2world(one), world2local(one)
      {
        intersector1 =  NULL;
#if defined(__SSE__)
        intersector4 = NULL;
#endif
#if defined(__AVX__)
        intersector8 = NULL;
#endif
#if defined(__MIC__)
        intersector16 = NULL;
#endif
      }

      void calculateWorldData ()
      {
        world2local = rcp(local2world);
        if (isEmpty(localBounds)) worldBounds = empty;
        else worldBounds = xfmBounds(local2world,localBounds);
      }
      
    public:
      int id0;                    //!< 1st user ID
      int id1;                    //!< 2nd user ID
      int mask;                   //!< for masking out this object during traversal
      bool hasTransform;          //!< true if tranformation is not the identity transform
      BBox3f localBounds;         //!< local object bounding box
      BBox3f worldBounds;         //!< tranformed world space bounding box
      AffineSpace3f local2world;  //!< transforms from local space to world space
      AffineSpace3f world2local;  //!< transforms from world space to local space

    public:
      const Intersector1* intersector1;
#if defined(__SSE__)
      const Intersector4* intersector4;
#endif
#if defined(__AVX__)
      const Intersector8* intersector8;
#endif
#if defined(__MIC__)
      const Intersector16* intersector16;
#endif
    };

    /*! Construction */
    VirtualScene (size_t numObjects, const char* accelTy);
    
    /*! Destruction */
    ~VirtualScene ();

    /*! clear registries */
    static void clearRegistry ();

    /*! Returns number of objects. */
    size_t size() const;
    
    /*! Returns bounds of ith object. */
    BBox3f bounds(size_t i) const;

    /*! builds acceleration structure */
    void build(TaskScheduler::Event* event, std::string builderName);

    /*! Deletes temporary data. */
    void freeze();

    /*! returns intersector */
    Intersector1* intersector1(std::string travName) const;

#if defined(__SSE__)
    Intersector4* intersector4(std::string travName) const;
#endif

#if defined(__AVX__)
    Intersector8* intersector8(std::string travName) const;
#endif

#if defined(__MIC__)
    Intersector16* intersector16(std::string travName) const;
#endif
    
    /*! returns the ith object */
    __forceinline Object& get(size_t i) {
      assert(i < numObjects);
      return objects[i];
    }

    /*! returns the ith object */
    __forceinline const Object& get(size_t i) const {
      assert(i < numObjects);
      return objects[i];
    }

    /*! returns pointer to vertex array */
    void* getVertices() { 
      return objects; 
    }

  public:
    Object* objects;
    size_t numObjects;
  };
}

#endif
