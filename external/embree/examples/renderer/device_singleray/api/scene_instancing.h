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

#ifndef __EMBREE_BACKEND_SCENE_INSTANCING_H__
#define __EMBREE_BACKEND_SCENE_INSTANCING_H__

#include "scene.h"
#include "embree/include/embree.h"

namespace embree
{
  /*! Scene that supports instancing of geometry. */
  class BackendSceneInstancing : public BackendScene
  {
  public:

    struct Primitive : public RefCount
    {
      ALIGNED_CLASS;
    public:
      Primitive (const bool bstatic,
                 const Ref<Shape>& shape,
                 const Ref<Light>& light,
                 const Ref<Material>& material,
                 const AffineSpace3f& local2world,
                 const light_mask_t illumMask,
                 const light_mask_t shadowMask)
        : bstatic(bstatic), shape(shape), light(light), material(material), 
          local2world(local2world), normal2world(rcp(local2world.l).transposed()),
          illumMask(illumMask), shadowMask(shadowMask) 
      {
        hasTransform = local2world != AffineSpace3f(one);
      }

      __forceinline void postIntersect(const Ray& ray, DifferentialGeometry& dg) const 
      {
        dg.material   = material.ptr;
        dg.light      = (AreaLight*)light.ptr;
        dg.illumMask  = illumMask;
        dg.shadowMask = shadowMask;
        shape->postIntersect(ray,dg);
        if (hasTransform) {
          dg.Tx   = xfmVector(local2world,dg.Tx); 
          dg.Ty   = xfmVector(local2world,dg.Ty);
          dg.Ng = xfmVector(normal2world,dg.Ng);
          dg.Ns = normalize(xfmVector(normal2world,dg.Ns));
        }
      }
      
    public:
      bool bstatic;
      bool hasTransform;
      Ref<Shape> shape;
      Ref<Light> light;
      Ref<Material> material;
      AffineSpace3f local2world;
      LinearSpace3f normal2world;
      light_mask_t illumMask; 
      light_mask_t shadowMask;
    };

    /*! API handle that manages user actions. */
    class Handle : public BackendScene::Handle {
      ALIGNED_CLASS;
    public:

      Handle () : geom_static (NULL), geom_static_intersector1(NULL) {}

      ~Handle() {
        if (geom_static             ) rtcDeleteGeometry    (geom_static             ); geom_static             = NULL;
        if (geom_static_intersector1) rtcDeleteIntersector1(geom_static_intersector1); geom_static_intersector1 = NULL;
      }
             
      void setPrimitive(size_t slot, Ref<PrimitiveHandle> prim)
      {
        if (slot >= prims.size()) prims.resize(slot+1);
        if (prims[slot] && prims[slot]->bstatic) {
          if (geom_static             ) rtcDeleteGeometry    (geom_static             ); geom_static             = NULL;
          if (geom_static_intersector1) rtcDeleteIntersector1(geom_static_intersector1); geom_static_intersector1 = NULL;
        }

        Ref<Shape> shape = prim->getShapeInstance();
        Ref<Light> light = prim->getLightInstance();
        Ref<Material> material = prim->getMaterialInstance();
        AffineSpace3f transform = prim->transform;

        if (light) shape = light->shape();
        if (prim->bstatic) {
          if (geom_static             ) rtcDeleteGeometry    (geom_static             ); geom_static             = NULL;
          if (geom_static_intersector1) rtcDeleteIntersector1(geom_static_intersector1); geom_static_intersector1 = NULL;
          transform = AffineSpace3f(one);
          if (shape) shape = shape->transform(prim->transform);
        }
        if (light) light = light->transform(prim->transform,prim->illumMask,prim->shadowMask);
        prims[slot] = new Primitive(prim->bstatic,shape,light,material,transform,prim->illumMask,prim->shadowMask);
      }
      
      void create() 
      {
        /*! create static accel */
        if (!geom_static)
        {
          /* count number of static vertices and triangles */
          size_t numAllocatedTriangles = 0, numAllocatedVertices = 0;
          for (size_t i=0; i<prims.size(); i++) {
            if (prims[i] && prims[i]->shape && prims[i]->bstatic) {
              numAllocatedVertices  += prims[i]->shape->numVertices();
              numAllocatedTriangles += prims[i]->shape->numTriangles();
            }
          }

          /*! skip if we have no static geometry */
          if (numAllocatedTriangles != 0) 
          {
            /*! extact static geometry and build acceleration structure */
            size_t numTriangles = 0, numVertices = 0;
            geom_static = rtcNewTriangleMesh(numAllocatedTriangles, numAllocatedVertices, accelTy.c_str());
            RTCTriangle* triangles = (RTCTriangle*) rtcMapTriangleBuffer(geom_static);
            Vec3fa*      positions = (Vec3fa*     ) rtcMapPositionBuffer(geom_static);
            for (size_t i=0; i<prims.size(); i++) {
              if (prims[i] && prims[i]->shape && prims[i]->bstatic)
                prims[i]->shape->extract(i,triangles,numTriangles,positions,numVertices);
            }
            if (numTriangles > numAllocatedTriangles) throw std::runtime_error("internal error");
            if (numVertices  > numAllocatedVertices ) throw std::runtime_error("internal error");
            rtcUnmapTriangleBuffer(geom_static);
            rtcUnmapPositionBuffer(geom_static);
            rtcBuildAccel(geom_static, builderTy.c_str());
            geom_static_intersector1 = rtcQueryIntersector1(geom_static,traverserTy.c_str());
          }
        }

        /*! create toplevel accel */
        RTCGeometry* objs = rtcNewVirtualGeometry(prims.size()+1,"default.virtual");
        for (size_t i=0; i<prims.size(); i++) {
          Ref<Primitive>& prim = prims[i];
          if (prim && prim->shape && !prim->bstatic) {
            RTCIntersector1* intersector1 = NULL; 
            RTCGeometry* geom = prim->shape->getAccel(intersector1);
            BBox3f bounds; rtcGetBounds(geom,&bounds.lower.x,&bounds.upper.x);
            Array12f array = copyToArray(prim->local2world);
            rtcSetVirtualGeometryUserData (objs, i, i, 0);
            rtcSetVirtualGeometryBounds(objs, i, &bounds.lower.x, &bounds.upper.x, (RTCTransformation*) (float*) array);
            rtcSetVirtualGeometryIntersector1(objs,i,intersector1);
          }
        }
        if (geom_static) {
          size_t i = prims.size();
          BBox3f bounds; rtcGetBounds(geom_static,&bounds.lower.x,&bounds.upper.x);
          rtcSetVirtualGeometryUserData (objs, i, 0x7FFFFFFF, 0);
          rtcSetVirtualGeometryBounds(objs, i, &bounds.lower.x, &bounds.upper.x);
          rtcSetVirtualGeometryIntersector1(objs,i,geom_static_intersector1);
        }
        rtcBuildAccel(objs, "objectsplit");
        RTCIntersector1* intersector1 = rtcQueryIntersector1(objs,"default.default");
        
        /* create new scene */
        instance = new BackendSceneInstancing(prims,objs,intersector1);
      }
      
    public:
      RTCGeometry* geom_static;
      RTCIntersector1* geom_static_intersector1;
      std::vector<Ref<Primitive> > prims; //!< total geometry and lights
    };

    /*! Construction of scene. */
    BackendSceneInstancing (const std::vector<Ref<Primitive> >& geometry, RTCGeometry* accel, RTCIntersector1* intersector1)
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
    std::vector<Ref<Primitive> > geometry; //!< Geometry of the scene
  };
}

#endif
