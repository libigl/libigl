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

#include "bvh4i/bvh4i.h"
#include "geometry/triangle1.h"

namespace embree
{
  /*! Multi BVH with 4 children. Each node stores the bounding box of
   * it's 4 children as well as a 4 child indices. */
  class BVH4mb : public BVH4i
  {
  public:

    /*! BVH4mb Node */

    struct __aligned(64) Node
    {
    public:
      struct NodeStruct {
        float x,y,z;           // x,y, and z coordinates of bounds
        NodeRef child;         // encodes 1 is-leaf bit, 25 offset bits, and 6 num-items bits
      } lower[4], upper[4];    // lower and upper bounds of all 4 children

      struct Delta {
        float x,y,z,w;         
      } lower_t1[4], upper_t1[4];    // lower and upper bounds of all 4 children

      /*! Returns bounds of specified child. */
      __forceinline BBox3fa bounds(size_t i) const {
        Vec3fa l = *(Vec3fa*)&lower[i];
        Vec3fa u = *(Vec3fa*)&upper[i];
        return BBox3fa(l,u);
      }

      __forceinline BBox3fa bounds_t1(size_t i) const {
        Vec3fa l = *(Vec3fa*)&lower_t1[i];
        Vec3fa u = *(Vec3fa*)&upper_t1[i];
        return BBox3fa(l,u);
      }

      __forceinline void setInvalid(size_t i)
      {
	lower[i].x = pos_inf;
	lower[i].y = pos_inf;
	lower[i].z = pos_inf;
	lower[i].child = NodeRef(leaf_mask);

	upper[i].x = neg_inf;
	upper[i].y = neg_inf;
	upper[i].z = neg_inf;
	upper[i].child = NodeRef(0);
      }
      /*! Returns reference to specified child */
      __forceinline       NodeRef& child(size_t i)       { return lower[i].child; }
      __forceinline const NodeRef& child(size_t i) const { return lower[i].child; }

    };

    /*! BVH4mb Triangle01 */

    struct __aligned(64) Triangle01
    {
    public:
      Triangle1 t0;
      Triangle1 t1;
    };



  public:

    /*! BVH4 default constructor. */
    BVH4mb (const PrimitiveType& primTy, void* geometry = NULL) : BVH4i(primTy,geometry)
    {
    }


    static Accel* BVH4mbTriangle1ObjectSplitBinnedSAH(Scene* scene);

  };

  __forceinline std::ostream &operator<<(std::ostream &o, const BVH4mb::Node &v)
    {
      o << std::endl;
      o << "lower: ";
      for (size_t i=0;i<4;i++) o << "[" << v.lower[i].x << "," << v.lower[i].y << "," << v.lower[i].z << "," << v.lower[i].child <<"] ";
      o << std::endl;
      o << "upper: ";
      for (size_t i=0;i<4;i++) o << "[" << v.upper[i].x << "," << v.upper[i].y << "," << v.upper[i].z << "," << v.upper[i].child <<"] ";

      o << "lower_t1: ";
      for (size_t i=0;i<4;i++) o << "[" << v.lower_t1[i].x << "," << v.lower_t1[i].y << "," << v.lower_t1[i].z << "," << v.lower_t1[i].w <<"] ";
      o << std::endl;
      o << "upper_t1: ";
      for (size_t i=0;i<4;i++) o << "[" << v.upper_t1[i].x << "," << v.upper_t1[i].y << "," << v.upper_t1[i].z << "," << v.upper_t1[i].w <<"] ";
      o << std::endl;
      return o;
    } 

  __forceinline mic_i getTriMasks(const BVH4mb::Triangle01 * __restrict__ const tptr)
  {
    return swDDDD(gather16i_4i_align(&tptr[0].t0.v2,&tptr[1].t0.v2,&tptr[2].t0.v2,&tptr[3].t0.v2));
  }

}
