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

#include "../tutorial_host/tutorials_host.h"
#include "../tutorial_device/tutorials_device.h"
#include "../tutorials/obj_loader.h"

namespace embree
{
  /* scene */
  extern OBJMesh g_mesh;
  extern Vec3f pointLightPosition;
  extern Vec3f pointLightIntensity;
  extern Vec3f ambientLightIntensity;
  ispc::Scene* g_scene = NULL;

  /* framebuffer */
  int* g_pixels = NULL;
  int g_width = -1;
  int g_height = -1;

  void init(int verbose) {
    ispc::init(verbose);
  }

  void setScene(OBJMesh& g_mesh, 
                const Vec3f pointLightPosition,
                const Vec3f pointLightIntensity,
                const Vec3f ambientLightIntensity)
  {
    g_scene = new ispc::Scene;
    g_scene->materials = g_mesh.material.size() ? &g_mesh.material[0] : NULL;
    g_scene->positions = g_mesh.v.size() ? &g_mesh.v[0] : NULL;
    g_scene->normals = g_mesh.vn.size() ? &g_mesh.vn[0] : NULL;
    g_scene->texcoords = g_mesh.vt.size() ? &g_mesh.vt[0] : NULL;
    g_scene->triangles = g_mesh.triangles.size() ? &g_mesh.triangles[0] : NULL;
    g_scene->numMaterials = g_mesh.material.size();
    g_scene->numVertices = g_mesh.v.size();
    g_scene->numTriangles = g_mesh.triangles.size();
    g_scene->pointLightPosition = pointLightPosition;
    g_scene->pointLightIntensity = pointLightIntensity;
    g_scene->ambientLightIntensity = ambientLightIntensity;
    ispc::set_scene(g_scene);
  }

  void resize(int32 width, int32 height)
  {
    if (width == g_width && height == g_height)
      return;

    if (g_pixels) alignedFree(g_pixels);
    g_width = width;
    g_height = height;
    g_pixels = (int*) alignedMalloc(g_width*g_height*sizeof(int),64);
  }

  int* render(const float time, const Vec3f& vx, const Vec3f& vy, const Vec3f& vz, const Vec3f& p)
  {
    ispc::render(g_pixels,g_width,g_height,time, 
                 (ispc::vec3f&)vx, 
                 (ispc::vec3f&)vy, 
                 (ispc::vec3f&)vz,
                 (ispc::vec3f&)p);

    return g_pixels;
  }
  
  void unmap ()
  {
  }

  void cleanup()
  {
    ispc::cleanup();
    delete g_scene; g_scene = NULL;
    alignedFree(g_pixels); g_pixels = NULL;
  }
}
