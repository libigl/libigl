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

#ifndef __EMBREE_BACKEND_SCENE_FLAT_H__
#define __EMBREE_BACKEND_SCENE_FLAT_H__

#include "scene.h"
#include "embree/include/embree.h"

namespace embree
{
  /*! Flat scene, no support for instancing, best render performance. */
  class BackendSceneFlat : public BackendScene
  {
  public:

    struct Primitive : public RefCount
    {
    public:
      Primitive (const Ref<Shape>& shape,
                 const Ref<Light>& light,
                 const Ref<Material>& material,
                 const light_mask_t illumMask,
                 const light_mask_t shadowMask)
        : shape(shape), light(light), material(material),
          illumMask(illumMask), shadowMask(shadowMask) {}

      __forceinline void postIntersect(const Ray& ray, DifferentialGeometry& dg) const 
      {
        dg.material   = material.ptr;
        dg.light      = (AreaLight*) light.ptr;
        dg.illumMask  = illumMask;
        dg.shadowMask = shadowMask;
        shape->postIntersect(ray,dg);
      }
      
    public:
      Ref<Shape> shape;
      Ref<Light> light;
      Ref<Material> material;
      light_mask_t illumMask;  /*! which light masks we receive illum from */
      light_mask_t shadowMask; /*! which light masks we cast shadows to */
    };

    /*! API handle that manages user actions. */
    class Handle : public BackendScene::Handle {
      ALIGNED_CLASS;
    public:
      
      void setPrimitive(size_t slot, Ref<PrimitiveHandle> prim) 
      {
        if (slot >= prims.size()) prims.resize(slot+1);
        Ref<Shape> shape = prim->getShapeInstance();
        Ref<Light> light = prim->getLightInstance();
        if (light) shape = light->shape();
        if (shape) shape = shape->transform(prim->transform);
        if (light) light = light->transform(prim->transform,prim->illumMask,prim->shadowMask);
        prims[slot] = new Primitive(shape,light,prim->getMaterialInstance(),prim->illumMask,prim->shadowMask);
      }
      
      void create() 
      {
        /* count number of vertices and triangles */
        size_t numAllocatedTriangles = 0, numAllocatedVertices = 0;
        for (size_t i=0; i<prims.size(); i++) {
          if (prims[i] && prims[i]->shape) {
            numAllocatedVertices  += prims[i]->shape->numVertices();
            numAllocatedTriangles += prims[i]->shape->numTriangles();
          }
        }
        
        /*! extact static geometry and build acceleration structure */
        size_t numTriangles = 0, numVertices = 0;
        RTCGeometry* mesh = rtcNewTriangleMesh(numAllocatedTriangles, numAllocatedVertices, accelTy.c_str());
        RTCTriangle* triangles = (RTCTriangle*) rtcMapTriangleBuffer(mesh);
        Vec3fa*        positions = (Vec3fa*       ) rtcMapPositionBuffer(mesh);
        for (size_t i=0; i<prims.size(); i++) {
          if (prims[i] && prims[i]->shape)
            prims[i]->shape->extract(i,triangles,numTriangles,positions,numVertices);
        }
        if (numTriangles > numAllocatedTriangles) throw std::runtime_error("internal error");
        if (numVertices  > numAllocatedVertices ) throw std::runtime_error("internal error");
        rtcUnmapTriangleBuffer(mesh);
        rtcUnmapPositionBuffer(mesh);

        rtcBuildAccel(mesh, builderTy.c_str());
        RTCIntersector1* intersector1 = rtcQueryIntersector1(mesh,traverserTy.c_str());
        
        /* create new scene */
        instance = new BackendSceneFlat(prims,mesh,intersector1);
      }
      
    public:
      std::vector<Ref<Primitive> > prims;
    };
        
    /*! Construction of scene. */
    BackendSceneFlat (const std::vector<Ref<Primitive> >& geometry, RTCGeometry* accel, RTCIntersector1* intersector1)
      : BackendScene(accel,intersector1), geometry(geometry)
    {
      for (size_t i=0; i<geometry.size(); i++) {
        const Ref<Primitive>& prim = geometry[i];
        if (prim && prim->light) add(prim->light);
      }
    }

    /*! Helper to call the post intersector of the shape instance,
     *  which will call the post intersector of the shape. */
    void postIntersect(const Ray& ray, DifferentialGeometry& dg) const {
      if (ray) geometry[ray.id0]->postIntersect(ray,dg);
    }

  private:
    std::vector<Ref<Primitive> > geometry;  //!< Geometry of the scene
  };
}

#endif
