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

namespace embree
{
  struct PrimitiveType
  {
    /*! constructs the triangle type */
    PrimitiveType (const char* name, size_t bytes, size_t blockSize, bool needVertices, int intCost) 
      : name(name), bytes(bytes), blockSize(blockSize), needVertices(needVertices), intCost(intCost) {}

    /*! Computes the number of blocks required to store a number of triangles. */
    virtual size_t blocks(size_t x) const = 0;

    /*! Returns the number of stored primitives in a block. */
    virtual size_t size(const char* This) const = 0;

  public:
    std::string name;       //!< name of this primitive type
    size_t bytes;           //!< number of bytes of the triangle data
    size_t blockSize;       //!< block size
    bool   needVertices;    //!< determines if we need the vertex array
    int    intCost;         //!< cost of one ray/primitive intersection
  };
}
