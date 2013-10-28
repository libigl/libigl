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

#ifndef __EMBREE_BUILD_GEOMETRY_H__
#define __EMBREE_BUILD_GEOMETRY_H__

#include "../common/default.h"
#include "../include/embree.h"
#include "../common/accel.h"

namespace embree
{
  /*! Virtual interface to geometry. */
  struct RTCGeometry : public RefCount
  {
    ALIGNED_CLASS;
  public:

    /*! default constructor*/
    RTCGeometry (const BBox3f& bounds) : approx(bounds), accel(NULL) {}

    /*! this class own the acceleration structure, thus we have to delete it */
    ~RTCGeometry () { if (accel) delete accel; accel = NULL; }

    /*! returns the number of primitives */
    virtual size_t size() const = 0;
    
    /*! calculates the bounds of the ith primitive */
    virtual BBox3f bounds(size_t i) const = 0;

    /*! invokes an acceleration structure builder for the geometry */
    virtual void build(TaskScheduler::Event* event, std::string builderName) = 0;

    /*! frees data not required for static objects, e.g. for a triangle mesh 
        the list of triangle indices and in some cases also the vertex array 
        can be freed. */
    virtual void freeze() = 0;

    /*! return an intersector for single rays */
    virtual Intersector1* intersector1(std::string travName) const = 0;

    /*! return an intersector for ray packets of size 4 */
#if defined(__SSE__)
    virtual Intersector4* intersector4(std::string travName) const  = 0;
#endif

    /*! return an intersector for ray packets of size 8 */
#if defined(__AVX__)
    virtual Intersector8* intersector8(std::string travName) const  = 0;
#endif

    /*! return an intersector for ray packets of size 16 */
#if defined(__MIC__)
    virtual Intersector16* intersector16(std::string travName) const  = 0;
#endif

    /*! returns vertex array pointer of a triangle mesh */
    virtual void* getVertices() { return NULL; } 

    /*! returns number of vertices of a triangle mesh */
    virtual size_t getNumVertices() const { return 0; }

    /*! returns the bounding box of the geometry */
    __forceinline BBox3f bounds() const { return accel->bounds; }

  public:
    BBox3f approx; //!< optional approximate bounds
    Accel* accel;  //!< acceleration structure
  };
}

#endif
