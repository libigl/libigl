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

#ifndef __EMBREE_LOADERS_H__
#define __EMBREE_LOADERS_H__

#include "obj_loader.h"
#include "xml_loader.h"

namespace embree
{
  Device::RTImage                  loadImage  (const FileName& fileName, Device* device);
  Device::RTTexture                loadTexture(const FileName& fileName, Device* device);
  std::vector<Device::RTPrimitive> loadScene  (const FileName& fileName, Device* device);
}

#endif
