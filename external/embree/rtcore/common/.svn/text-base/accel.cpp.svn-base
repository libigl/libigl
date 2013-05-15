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

/* include interfaces */
#include "accel.h"
#include "intersector.h"

/* include BVH2 */
#include "bvh2/bvh2.h"
#include "bvh2/bvh2_builder.h"

/* include BVH4 */
#include "bvh4/bvh4.h"
#include "bvh4/bvh4_builder.h"

/* include BVH4MB */
#include "bvh4mb/bvh4mb.h"
#include "bvh4mb/bvh4mb_builder.h"

/* include all build heuristics */
#include "heuristics.h"

/* include all triangle representations */
#include "triangle/triangles.h"

namespace embree
{
  /*! interface names */
  const char* const Intersector ::name = "Intersector";
  
  /*! triangle types */
  const Triangle1i::Type Triangle1i::type;
  const Triangle4i::Type Triangle4i::type;
  const Triangle1v::Type Triangle1v::type;
  const Triangle4v::Type Triangle4v::type;
  const Triangle1 ::Type Triangle1 ::type;
  const Triangle4 ::Type Triangle4 ::type;
  const Triangle8 ::Type Triangle8 ::type;
 
  void* rtcMalloc(size_t bytes) {
    return alignedMalloc(bytes);
  }

  void rtcFreeMemory() {
    Alloc::global.clear();
  }

  template<typename Builder>
  Ref<Accel> build(const TriangleType& trity, const std::string& intTy,
                   const BuildTriangle* triangles, size_t numTriangles, 
                   const Vec3fa* vertices, size_t numVertices, 
                   const BBox3f& bounds, bool freeArrays)
  {
#if 0 // FIXME: enable to measure build performance
    double dt = 0;
    Ref<typename Builder::Type> bvh = null;
    for (int i=0; i<7; i++) {
      bvh = null;
      double t0 = getSeconds();
      bvh = Builder(trity,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeArrays && (i == 6)).bvh;
      double t1 = getSeconds();
      std::cout << "build iteration " << i << ": " << (t1-t0)*1000.0f << " ms, " << double(numTriangles)/(t1-t0)*1E-6 << " Mtris/s" << std::endl;
      if (i>=2) dt += t1-t0;
    }
    dt/=5.0f;
#else
    double t0 = getSeconds();
    Ref<typename Builder::Type> bvh = Builder(trity,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeArrays).bvh;
    double dt = getSeconds()-t0;
#endif

    if (freeArrays) 
    {
      /*! free vertices if no longer required */
      if (!trity.needVertices) alignedFree(vertices);
      
      /*! free triangles */
      alignedFree(triangles);
    }

    /* output statistics */
    std::ostringstream stream;
    std::ios::fmtflags flags = stream.flags();
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(0);
    stream << "triangles = " << numTriangles << std::endl;
    stream << "build time = " << dt*1000.0f << " ms" << std::endl;
    stream.precision(3);
    stream << "build performance = " << double(numTriangles)/dt*1E-6 << " Mtris/s" << std::endl;
    stream << "build_stat " 
           << scheduler->getNumThreads() << " " 
           << numTriangles << " " << dt*1000.0f << " " 
           << double(numTriangles)/dt*1E-6 << std::endl;
    stream.setf(flags);
    bvh->print(stream);
    std::cout << stream.str();

    return bvh.ptr;
  }

  Ref<Accel> rtcCreateAccel(const char* accelTy_i, const char* triTy_i, 
                            const BuildTriangle* triangles , size_t numTriangles, 
                            const BuildVertex  * vertices_i, size_t numVertices,
                            const BBox3f& bounds, bool freeData)
  {
    const Vec3fa* vertices = (const Vec3fa*) vertices_i;

    /* parse accel.builder string */
    std::string accelTy = accelTy_i;
    
    /* parse triangle.type string */
    std::string triTy = triTy_i, intTy = "default";
    {
      size_t pos = triTy.find_first_of('.');
      if (pos != std::string::npos) {
        intTy = triTy.substr(pos+1);
        triTy.resize(pos);
      }
    }

    /* BVH2 with object split builder */
    if (accelTy == "bvh2.objectsplit" || accelTy == "bvh2") 
    {
      if (triTy == "default")
#ifdef __AVX__
        return build<BVH2Builder<HeuristicBinning<Triangle8::logBlockSize> > >(Triangle8::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
#else
        return build<BVH2Builder<HeuristicBinning<Triangle4::logBlockSize> > >(Triangle4::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
#endif
      else if (triTy == "triangle1i")
        return build<BVH2Builder<HeuristicBinning<Triangle1i::logBlockSize> > >(Triangle1i::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4i")
        return build<BVH2Builder<HeuristicBinning<Triangle4i::logBlockSize> > >(Triangle4i::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle1v")
        return build<BVH2Builder<HeuristicBinning<Triangle1v::logBlockSize> > >(Triangle1v::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4v")
        return build<BVH2Builder<HeuristicBinning<Triangle4v::logBlockSize> > >(Triangle4v::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle1")
        return build<BVH2Builder<HeuristicBinning<Triangle1::logBlockSize> > >(Triangle1::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4")
        return build<BVH2Builder<HeuristicBinning<Triangle4::logBlockSize> > >(Triangle4::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle8")
        return build<BVH2Builder<HeuristicBinning<Triangle8::logBlockSize> > >(Triangle8::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else {
        throw std::runtime_error("invalid triangle type for bvh2: "+std::string(triTy));
        return null;
      }
    }

    /* BVH2 with spatial split builder */
    if (accelTy == "bvh2.spatialsplit") 
    {
      if (triTy == "default")
#ifdef __AVX__
        return build<BVH2Builder<HeuristicSpatial<Triangle8::logBlockSize> > >(Triangle8::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
#else
        return build<BVH2Builder<HeuristicSpatial<Triangle4::logBlockSize> > >(Triangle4::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
#endif
      else if (triTy == "triangle1i")
        return build<BVH2Builder<HeuristicSpatial<Triangle1i::logBlockSize> > >(Triangle1i::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4i")
        return build<BVH2Builder<HeuristicSpatial<Triangle4i::logBlockSize> > >(Triangle4i::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle1v")
        return build<BVH2Builder<HeuristicSpatial<Triangle1v::logBlockSize> > >(Triangle1v::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4v")
        return build<BVH2Builder<HeuristicSpatial<Triangle4v::logBlockSize> > >(Triangle4v::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle1")
        return build<BVH2Builder<HeuristicSpatial<Triangle1::logBlockSize> > >(Triangle1::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4")
        return build<BVH2Builder<HeuristicSpatial<Triangle4::logBlockSize> > >(Triangle4::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle8")
        return build<BVH2Builder<HeuristicSpatial<Triangle8::logBlockSize> > >(Triangle8::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else {
        throw std::runtime_error("invalid triangle type for bvh2: "+std::string(triTy));
        return null;
      }
    }

    /* BVH4 with object split builder */
    if (accelTy == "bvh4.objectsplit" || accelTy == "bvh4" || accelTy == "default") 
    {
      if (triTy == "default")
#ifdef __AVX__
        return build<BVH4Builder<HeuristicBinning<Triangle8::logBlockSize> > >(Triangle8::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
#else
        return build<BVH4Builder<HeuristicBinning<Triangle4::logBlockSize> > >(Triangle4::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
#endif
      else if (triTy == "triangle1i")
        return build<BVH4Builder<HeuristicBinning<Triangle1i::logBlockSize> > >(Triangle1i::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4i")
        return build<BVH4Builder<HeuristicBinning<Triangle4i::logBlockSize> > >(Triangle4i::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle1v")
        return build<BVH4Builder<HeuristicBinning<Triangle1v::logBlockSize> > >(Triangle1v::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4v")
        return build<BVH4Builder<HeuristicBinning<Triangle4v::logBlockSize> > >(Triangle4v::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle1")
        return build<BVH4Builder<HeuristicBinning<Triangle1::logBlockSize> > >(Triangle1::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4")
        return build<BVH4Builder<HeuristicBinning<Triangle4::logBlockSize> > >(Triangle4::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle8")
         return build<BVH4Builder<HeuristicBinning<Triangle8::logBlockSize> > >(Triangle8::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else {
        throw std::runtime_error("invalid triangle type for bvh4: "+std::string(triTy));
        return null;
      }
    }

    /* BVH4 with spatial split builder */
    if (accelTy == "bvh4.spatialsplit") 
    {
      if (triTy == "default")
#ifdef __AVX__
        return build<BVH4Builder<HeuristicSpatial<Triangle8::logBlockSize> > >(Triangle8::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
#else
        return build<BVH4Builder<HeuristicSpatial<Triangle4::logBlockSize> > >(Triangle4::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
#endif
      else if (triTy == "triangle1i")
        return build<BVH4Builder<HeuristicSpatial<Triangle1i::logBlockSize> > >(Triangle1i::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4i")
        return build<BVH4Builder<HeuristicSpatial<Triangle4i::logBlockSize> > >(Triangle4i::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle1v")
        return build<BVH4Builder<HeuristicSpatial<Triangle1v::logBlockSize> > >(Triangle1v::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4v")
        return build<BVH4Builder<HeuristicSpatial<Triangle4v::logBlockSize> > >(Triangle4v::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle1")
        return build<BVH4Builder<HeuristicSpatial<Triangle1::logBlockSize> > >(Triangle1::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4")
        return build<BVH4Builder<HeuristicSpatial<Triangle4::logBlockSize> > >(Triangle4::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle8")
        return build<BVH4Builder<HeuristicSpatial<Triangle8::logBlockSize> > >(Triangle8::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else {
        throw std::runtime_error("invalid triangle type for bvh4: "+std::string(triTy));
        return null;
      }
    }

    /* BVH4MB with object split builder */
    if (accelTy == "bvh4mb.objectsplit" || accelTy == "bvh4mb") 
    {
      if (triTy == "default")
        return build<BVH4MBBuilder<HeuristicBinning<Triangle4::logBlockSize> > >(Triangle4i::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else if (triTy == "triangle4i")
        return build<BVH4MBBuilder<HeuristicBinning<Triangle4i::logBlockSize> > >(Triangle4i::type,intTy,triangles,numTriangles,vertices,numVertices,bounds,freeData);
      else {
        throw std::runtime_error("invalid triangle type for bvh4mb: "+std::string(triTy));
        return null;
      }
    }
    
    else {
      throw std::runtime_error("invalid acceleration structure \""+std::string(accelTy)+"\"");
      return null;
    }
  }
}
