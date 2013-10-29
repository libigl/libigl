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

#include "registry_accel.h"
#include "accel.h"

namespace embree
{
  void AccelRegistry::setDefaultAccel(const std::string& accelName) {
    defaultAccel = accelName; 
  }
  
  void AccelRegistry::setDefaultTriangle(const std::string& accelName, const std::string& triangleName) {
    defaultTriangle[accelName] = triangleName; 
  }

  void AccelRegistry::add(const std::string& accelName, Constructor constructor) {
    accelTable[accelName] = constructor;
  }
  
  void AccelRegistry::add(const std::string& triName, const TriangleType& trity) {
    triangleTable[triName] = &trity;
  }
  
  Accel* AccelRegistry::create(std::string accelName, RTCGeometry* geom)
  {
    /*! split accelName into accel and triangle */
    std::string triangleName = "default";
    size_t pos = accelName.find_first_of('.');
    if (pos != std::string::npos) {
      triangleName = accelName.substr(pos+1);
      accelName.resize(pos);
    }
    
    /*! set default acceleration structure */
    if (accelName == "default") {
      accelName = defaultAccel;
      if (accelName == "") throw std::runtime_error("No default acceleration structure defined.");
    }

    /*! check if acceleration structure is registered */
    std::map<std::string, Constructor>::iterator accel = accelTable.find(accelName);
    if (accel == accelTable.end())
      throw std::runtime_error("Acceleration structure \""+accelName+"\" is not registered.");
    
    /*! set default triangle */
    if (triangleName == "default") {
      triangleName = defaultTriangle[accelName];
      if (triangleName == "") throw std::runtime_error("Acceleration structure \""+accelName+"\" has no default triangle.");
    }

    /*! check if triangle is registered */
    std::map<std::string, const TriangleType*>::iterator triangle = triangleTable.find(triangleName);
    if (triangle == triangleTable.end())
      throw std::runtime_error("Triangle \""+triangleName+"\" is not registered.");

    /*! construct empty acceleration structure */
    return (*accel->second)(geom,(const TriangleType&)(*triangle->second));
  }

  void AccelRegistry::clear()
  {
    defaultAccel = "";
    defaultTriangle.clear();
    accelTable.clear();
    triangleTable.clear();
  }
}
