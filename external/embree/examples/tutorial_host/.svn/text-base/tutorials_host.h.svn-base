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

#include "math/vec3.h"
#include "../tutorials/obj_loader.h"

namespace embree
{
  /* initialize renderer */
  void init(int verbose);

  /* set scene */
  void setScene(OBJMesh& mesh, 
                const Vec3f pointLightPosition,
                const Vec3f pointLightIntensity,
                const Vec3f ambientLightIntensity);

  /* resize framebuffer */
  void resize(int32 width, int32 height);

  /* render frame and map framebuffer */
  int* render(const float time, const Vec3f& vx, const Vec3f& vy, const Vec3f& vz, const Vec3f& p);
  
  /* unmap framebuffer */
  void unmap ();

  /* cleanup renderer */
  void cleanup();
}
