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

#include "transport_host.h"
#include "transport_device.h"
#include "tutorial/obj_loader.h"

extern "C" int64 get_tsc() {
  return __rdtsc();
}

namespace embree
{
  /* framebuffer */
  int* g_pixels = NULL;
  int g_width = -1;
  int g_height = -1;

  /* ISPC compatible mesh */
  struct ISPCMesh
  {
    Vec3fa* positions;    //!< vertex position array
    Vec3fa* normals;       //!< vertex normal array
    Vec2f* texcoords;     //!< vertex texcoord array
    OBJScene::Triangle* triangles;  //!< list of triangles
    int numVertices;
    int numTriangles;
    Vec3f dir;
    float offset;
  };

  /* ISPC compatible scene */
  struct ISPCScene
  {
    ISPCMesh** meshes;
    OBJScene::Material* materials;  //!< material list
    int numMeshes;
    int numMaterials;
    bool animate;
  };

  /* scene */
  extern "C" ISPCScene* g_ispc_scene = NULL;
  
  ISPCMesh* convertMesh (OBJScene::Mesh* in)
  {
    ISPCMesh* out = new ISPCMesh;
    out->positions = in->v.size() ? &in->v[0] : NULL;
    out->normals = in->vn.size() ? &in->vn[0] : NULL;
    out->texcoords = in->vt.size() ? &in->vt[0] : NULL;
    out->triangles = in->triangles.size() ? &in->triangles[0] : NULL;
    out->numVertices = in->v.size();
    out->numTriangles = in->triangles.size();
    out->dir = normalize(Vec3f(drand48(),drand48(),drand48())-Vec3f(0.5f));
    out->offset = 5.0f*drand48();
    return out;
  }

  void init(const char* cfg) {
    device_init(cfg);
  }

  void key_pressed (int32 key) {
    device_key_pressed(key);
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

  void set_scene (OBJScene* in) 
  {
    ISPCScene* out = new ISPCScene;
    out->materials = in->materials.size() ? &in->materials[0] : NULL;
    out->meshes = new ISPCMesh*[in->meshes.size()];
    for (size_t i=0; i<in->meshes.size(); i++) out->meshes[i] = convertMesh(in->meshes[i]);
    out->numMaterials = in->materials.size();
    out->numMeshes = in->meshes.size();
    out->animate = false; //g_animate;
    g_ispc_scene = out;
  }

  bool pick(const float x, const float y, const Vec3f& vx, const Vec3f& vy, const Vec3f& vz, const Vec3f& p, Vec3f& hitPos) {
    return device_pick(x,y,vx,vy,vz,p,hitPos);
  }

  void render(const float time, const Vec3f& vx, const Vec3f& vy, const Vec3f& vz, const Vec3f& p) {
    device_render(g_pixels,g_width,g_height,time,vx,vy,vz,p);
  }

  int* map () {
    return g_pixels;
  }
  
  void unmap () {
  }

  void cleanup()
  {
    device_cleanup();
    alignedFree(g_pixels); g_pixels = NULL;
  }
}
