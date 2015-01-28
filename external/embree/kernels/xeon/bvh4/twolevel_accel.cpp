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

#include "twolevel_accel.h"
#include "virtual_accel.h"
#include "common/scene.h"

#include "bvh4/bvh4.h"
#include "bvh4i/bvh4i.h"

namespace embree
{
  TwoLevelAccel::TwoLevelAccel (const std::string topAccel, Scene* scene)
    : scene(scene), accel(new VirtualAccel(topAccel,enabled_accels)) {}

  TwoLevelAccel::~TwoLevelAccel () {
    delete accel;
  }

  void TwoLevelAccel::buildUserGeometryAccels(size_t threadIndex, size_t threadCount)
  {
    if (scene->numUserGeometries == 0) 
      return;

    size_t N = scene->size();
    for (size_t i=0; i<N; i++) 
    {
      Geometry* geom = scene->getUserGeometrySafe(i);
      if (geom == NULL) continue;
      AccelSet* accel = scene->getUserGeometrySafe(i);

      if (geom->isEnabled()) {
        enabled_accels.push_back(accel);
      }

      if (geom->isModified()) {
        accel->build(threadIndex,threadCount);
        geom->state = Geometry::ENABLED;
      }
    }
  }

  void TwoLevelAccel::build(size_t threadIndex, size_t threadCount)
  {
    /* build all object accels */
    enabled_accels.clear();
    buildUserGeometryAccels(threadIndex,threadCount);

    /* build toplevel accel */
    accel->build(threadIndex,threadCount);
    bounds = accel->bounds;
    intersectors = accel->intersectors;
  }
}
