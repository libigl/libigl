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

#include "registry_builder.h"

namespace embree
{
  void BuilderRegistry::setDefaultBuilder(const std::string& accelName, const std::string& builderName) {
    defaultBuilder[accelName] = builderName; 
  }

  void BuilderRegistry::add(const std::string& accelName, const std::string& triName, const std::string& builderName, 
                            Builder builder, const size_t minLeafSize, const size_t maxLeafSize) 
  { 
    table[accelName+"."+triName+"."+builderName] = BuildClosure(builder,minLeafSize,maxLeafSize);
  }

  void BuilderRegistry::build(std::string builderName, TaskScheduler::Event* group, Accel* accel) 
  {
    std::string accelType = accel->name();

    /*! split accelName into accel and triangle */
    std::string accelName = accelType;
    std::string triangleName = "default";
    size_t pos = accelName.find_first_of('.');
    if (pos != std::string::npos) {
      triangleName = accelName.substr(pos+1);
      accelName.resize(pos);
    }

    /*! select default builder */
    if (builderName == "default") {
      builderName = defaultBuilder[accelName];
      if (builderName == "") throw std::runtime_error("Acceleration structure \""+accelName+"\" has no default builder.");
    }
     
    /*! find specified builder */
    std::string key = accelType + "." + builderName;
    if (table.find(key) == table.end()) 
      throw std::runtime_error("Acceleration structure \""+accel->name()+"\" has no builder \""+builderName+"\".");
    
    /*! build acceleration structure */
    const BuildClosure& closure = table[key];
    closure.builder(group,accel,closure.minLeafSize,closure.maxLeafSize);
  }

  void BuilderRegistry::clear()
  {
    defaultBuilder.clear();
    table.clear();
  }
}

