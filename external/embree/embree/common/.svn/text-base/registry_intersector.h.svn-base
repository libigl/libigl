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

#ifndef __EMBREE_INTERSECTOR_REGISTRY_H__
#define __EMBREE_INTERSECTOR_REGISTRY_H__

#include "default.h"
#include "../geometry/triangle.h"
#include "../include/embree.h"

namespace embree
{
  /*! Registry for intersectors */
  template<typename T>
    class IntersectorRegistry 
  {
  public:

    typedef T* (*Constructor) (const Accel* accel);

    /*! sets the default intersector for an acceleration structure */
    void setAccelDefaultTraverser(const std::string& accelName, const std::string& travName);
    
    /*! adds a new intersector */
    void add(const std::string& accelName, const std::string& triName, const std::string& travName, const std::string& intName, bool def, const Constructor constructor);
    
    /*! get registered intersector by name */
    T* get(std::string accelName, std::string travName, Accel* accel);
    
    /*! clears the registry */
    void clear();

  private:
    std::map<std::string, std::string> defaultTraverser; 
    std::map<std::string, std::string> defaultIntersector; 
    std::map<std::string, Constructor> table; 
  }; 
}

#endif

