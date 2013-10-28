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

#include "registry_intersector.h"

namespace embree
{
  template<typename T>
  void IntersectorRegistry<T>::setAccelDefaultTraverser(const std::string& accelName, const std::string& travName) {
    defaultTraverser[accelName] = travName;
  }
    
  template<typename T>
  void IntersectorRegistry<T>::add(const std::string& accelName, const std::string& triName, 
                                   const std::string& travName, const std::string& intName, bool def, 
                                   const Constructor constructor) 
  { 
    std::string accelKey = accelName+"."+triName;
    std::string travKey = travName+"."+intName;
    std::string key = accelKey+"."+travKey;
    if (def) defaultIntersector[accelKey+"."+travName] = intName;
    table[key] = constructor;
  }
    
  template<typename T>
  T* IntersectorRegistry<T>::get(std::string accelName, std::string travName, Accel* accel) 
  {
    /*! split accelName into accel and triangle */
    std::string triangleName = "default";
    size_t pos = accelName.find_first_of('.');
    if (pos != std::string::npos) {
      triangleName = accelName.substr(pos+1);
      accelName.resize(pos);
    }
    
    /*! split travName into traverser and triangle intersector */
    std::string intName = "default";
    pos = travName.find_first_of('.');
    if (pos != std::string::npos) {
      intName = travName.substr(pos+1);
      travName.resize(pos);
    }

    /*! set defaults */
    if (travName == "default") {
      travName = defaultTraverser[accelName+"."+triangleName];
      if (travName == "") {
        travName = defaultTraverser[accelName];
        if (travName == "") throw std::runtime_error("Acceleration structure \""+accelName+"."+triangleName+"\" has no default "+T::name);
      }
    }
    if (intName  == "default") {
      intName  = defaultIntersector[accelName+"."+triangleName+"."+travName];
      if (travName == "") throw std::runtime_error("Acceleration structure \""+accelName+"."+triangleName+"\" has no default triangle for \""+travName+"\"");
    }
    
    /*! find specified intersector */
    const std::string key = accelName+"."+triangleName+"."+travName+"."+intName;
    if (table.find(key) == table.end()) 
      throw std::runtime_error("Acceleration structure \""+accelName+"."+triangleName+"\" has no "+T::name+" \""+travName+"."+intName+"\"");
    
    /*! return intersector */
    if (g_verbose > 0)
      std::cout << "using intersector " << travName << "." << intName << " for " << accelName << "." << triangleName << std::endl;
    return table[key](accel);
  }

  template<typename T>
  void IntersectorRegistry<T>::clear()
  {
    defaultTraverser.clear();
    defaultIntersector.clear();
    table.clear();
  }

  /*! explicit template instantiations */
  template class IntersectorRegistry<Intersector1>;

#if defined(__SSE__)
  template class IntersectorRegistry<Intersector4>;
#endif

#if defined(__AVX__)
  template class IntersectorRegistry<Intersector8>;
#endif

#if defined(__MIC__)
  template class IntersectorRegistry<Intersector16>;
#endif
}

