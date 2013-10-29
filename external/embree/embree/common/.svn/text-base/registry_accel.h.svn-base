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

#ifndef __EMBREE_ACCEL_REGISTRY_H__
#define __EMBREE_ACCEL_REGISTRY_H__

#include "default.h"
#include "../geometry/triangle.h"
#include "../include/embree.h"

namespace embree
{
  /*! Registry for acceleration structures and triangles */
  class AccelRegistry 
  {
  public:
    typedef Accel* (*Constructor)(RTCGeometry* geom, const TriangleType& trity);
   
  public:

    /*! sets the default acceleration structure */
    void setDefaultAccel(const std::string& accelName);

    /*! sets the default triangle for an acceleration structure */
    void setDefaultTriangle(const std::string& accelName, const std::string& triangleName);

    /*! adds a new acceleration structure to the registry */
    void add(const std::string& accelName, Constructor constructor);

    /*! adds a new acceleration structure to the registry */
    void add(const std::string& triName, const TriangleType& trity);
      
    /*! constructs an empty acceleration structure */
    Accel* create(std::string accelName, RTCGeometry* geom);

    /*! clears the registry */
    void clear();
    
  private:
    std::string defaultAccel;
    std::map<std::string, std::string> defaultTriangle;
    std::map<std::string, Constructor> accelTable; 
    std::map<std::string, const TriangleType*> triangleTable; 
  }; 
}

#endif

