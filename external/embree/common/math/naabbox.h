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

#include "bbox.h"
#include "linearspace3.h"

namespace embree
{
  /*! Non-axis aligned bounds */
  template<typename T>
    struct NAABBox
  {
  public:
    
    __forceinline NAABBox () {}
    
    __forceinline NAABBox (EmptyTy) 
      : space(one), bounds(empty) {}
    
    __forceinline NAABBox (const BBox<T>& bounds) 
      : space(one), bounds(bounds) {}
      
    __forceinline NAABBox (const LinearSpace3<T>& space, const BBox<T>& bounds) 
      : space(space), bounds(bounds) {}
    
    friend std::ostream& operator<<(std::ostream& cout, const NAABBox& p) {
      return std::cout << "{ space = " << p.space << ", bounds = " << p.bounds << "}";
    }
    
  public:
    LinearSpace3<T> space; //!< orthonormal transformation
    BBox<T> bounds;       //!< bounds in transformed space
  };

  typedef NAABBox<Vec3f> NAABBox3f;
  typedef NAABBox<Vec3fa> NAABBox3fa;
}
