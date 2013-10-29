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

#ifndef __EMBREE_SHAPE_H__
#define __EMBREE_SHAPE_H__

#include "../default.h"
#include "../api/parms.h"
#include "../api/data.h"
#include "embree/include/embree.h"
#include "embree/common/ray.h"
#include "differentialgeometry.h"

namespace embree
{
  struct AccelType
  {
    AccelType (const Parms& parms)
    : accel(parms.getString("accel","default")), 
      builder(parms.getString("builder","default")),
      traverser(parms.getString("traverser","default")) {}
    
    AccelType (const std::string& accel, const std::string& builder, const std::string& traverser)
    : accel(accel), builder(builder), traverser(traverser) {}
    
  public:
    std::string accel;
    std::string builder;
    std::string traverser;
  };

  /*! Interface to different shapes. A shape is the smallest geometric
   *  primitive materials and area lights can be assigned to. */
  class Shape : public RefCount {
    ALIGNED_CLASS
  public:

    Shape (const AccelType& ty) : ty(ty), accel(NULL), intersector1(NULL) {}
    Shape (const Parms& parms ) : ty(parms), accel(NULL), intersector1(NULL) {}

    Shape (const std::string& accelTy, const std::string& builderTy, const std::string& traverserTy) 
      : ty(accelTy,builderTy,traverserTy), accel(NULL), intersector1(NULL) {}

    /*! Shape virtual destructor. */
    virtual ~Shape() {
      if (accel) rtcDeleteGeometry(accel);
      if (intersector1) rtcDeleteIntersector1(intersector1);
    }
    
    /*! Instantiates a new shape by transforming this shape to a different location. */
    virtual Ref<Shape> transform(const AffineSpace3f& xfm) const = 0;

    /*! Builds acceleration structure for this shape */
    virtual RTCGeometry* getAccel(RTCIntersector1*& intersector1_o) 
    { 
      if (!accel) 
      {
        /* allocate mesh */
        size_t numAllocatedTriangles = numTriangles();
        size_t numAllocatedVertices = numVertices();
        accel = rtcNewTriangleMesh (numAllocatedTriangles, numAllocatedVertices, ty.accel.c_str());
        RTCTriangle* triangles = (RTCTriangle*) rtcMapTriangleBuffer(accel);
        Vec3fa*      positions = (Vec3fa*     ) rtcMapPositionBuffer(accel);
        
        /* extract all primitives */
        size_t numTriangles = 0, numVertices = 0;
        BBox3f bounds = extract(0x7FFFFFFF,triangles,numTriangles,positions,numVertices);
        rtcSetApproxBounds(accel, (float*)&bounds.lower, (float*)&bounds.upper);
        if (numTriangles > numAllocatedTriangles) throw std::runtime_error("internal error");
        if (numVertices  > numAllocatedVertices ) throw std::runtime_error("internal error");
        rtcUnmapTriangleBuffer(accel);
        rtcUnmapPositionBuffer(accel);
        
        /* build acceleration structure */
        rtcBuildAccel(accel, ty.builder.c_str());
        intersector1 = rtcQueryIntersector1(accel, ty.traverser.c_str());
        rtcCleanupGeometry(accel);
      }
      intersector1_o = intersector1;
      return accel; 
    }

    /*! Counts the number of triangles required for extraction. */
    virtual size_t numTriangles() const = 0;

    /*! Counts the number of vertices required for extraction. */
    virtual size_t numVertices() const = 0;

    /*! Extracts triangles for spatial index structure. */
    virtual BBox3f extract(size_t id, RTCTriangle* triangles, size_t& numTriangles, Vec3fa* positions, size_t& numVertices) const = 0;

    /*! Performs interpolation of shading vertex parameters. */
    virtual void postIntersect(const Ray& ray, DifferentialGeometry& dg) const = 0;

  public:
    AccelType ty;
    RTCGeometry* accel;            //!< acceleration structure for the shape
    RTCIntersector1* intersector1; //!< intersector for acceleration structure
  };
}

#endif
