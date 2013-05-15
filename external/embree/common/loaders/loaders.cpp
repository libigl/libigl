// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "loaders.h"
#include "sys/stl/string.h"
#include <map>

namespace embree
{
  Device::RTImage loadImage(const FileName &fileName, Device *device) 
  {
    static std::map<std::string, Device::RTImage> image_map;
    
    if (image_map.find(fileName.str()) != image_map.end()) 
      return(image_map[fileName.str()]);
    
    return(image_map[fileName.str()] = device->rtNewImageFromFile(fileName.c_str()));
  }

  Device::RTTexture loadTexture(const FileName &fileName, Device *device) 
  {
    static std::map<std::string, Device::RTTexture> texture_map;
    
    if (texture_map.find(fileName.str()) != texture_map.end()) 
      return(texture_map[fileName.str()]);
    
    Device::RTTexture texture = device->rtNewTexture("nearest");
    device->rtSetImage(texture, "image", loadImage(fileName, device));
    device->rtCommit(texture);
    
    return(texture_map[fileName.str()] = texture);
  }

  std::vector<Device::RTPrimitive> loadScene(const FileName &fileName, Device *device) 
  {
    std::string ext = strlwr( fileName.ext() );
    if (ext == "obj") return loadOBJ(fileName, device);
    if (ext == "xml") return loadXML(fileName, device);
    throw std::runtime_error("file format " + ext + " not supported");
  }
}
