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

#ifndef __EMBREE_BUILDER_REGISTRY_H__
#define __EMBREE_BUILDER_REGISTRY_H__

#include "default.h"
#include "../geometry/triangle.h"
#include "../include/embree.h"
#include "accel.h"

namespace embree
{
  /*! Registry for acceleration structure builders */
  class BuilderRegistry 
  {
  public:
    typedef void (*Builder)(TaskScheduler::Event* group, Accel* accel, const size_t maxLeafSize, const size_t forceLeafSize);
   
    struct BuildClosure 
    {
      BuildClosure () {}
      BuildClosure (Builder builder, const size_t minLeafSize, const size_t maxLeafSize)
      : builder(builder), minLeafSize(minLeafSize), maxLeafSize(maxLeafSize) {}

    public:
      Builder builder;
      size_t minLeafSize;
      size_t maxLeafSize;
    };
        
  public:

    /*! sets the default builder for an acceleration structure */
    void setDefaultBuilder(const std::string& accelName, const std::string& builderName);

    /*! adds a new builder to the registry */
    void add(const std::string& accelName, const std::string& triName, const std::string& builderName, 
             Builder builder, const size_t minLeafSize = 1, const size_t maxLeafSize = inf);
      
    /*! builds an acceleration structure using the specified builder */
    void build(std::string builderName, TaskScheduler::Event* group, Accel* accel);
    
    /*! clears the registry */
    void clear();

  private:
    std::map<std::string, std::string> defaultBuilder;
    std::map<std::string, BuildClosure> table; 
  }; 
}

#endif

