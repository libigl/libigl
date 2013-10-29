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

#ifndef __EMBREE_OBJ_LOADER_H__
#define __EMBREE_OBJ_LOADER_H__

#include <vector>
#include "sys/filename.h"
#include "device/device.h"
#include "device/handle.h"

namespace embree
{
  std::vector<Handle<Device::RTPrimitive> > loadOBJ(const FileName &fileName);
}

#endif // __EMBREE_OBJ_LOADER_H__

