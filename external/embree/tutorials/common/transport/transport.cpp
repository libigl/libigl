// ======================================================================== //
// Copyright 2009-2014 Intel Corporation                                    //
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
#include "tutorial/tutorial_device.h"
#include "tutorial/scene_device.h"

extern "C" int64 get_tsc() {
  return read_tsc();
}

namespace embree
{
  /* framebuffer */
  int* g_pixels = NULL;
  int g_width = -1;
  int g_height = -1;

  /* scene */
  extern "C" ISPCScene* g_ispc_scene = NULL;

  ISPCHairSet* convertHair (OBJScene::HairSet* in)
  {
    ISPCHairSet* out = new ISPCHairSet;
    out->v = in->v.size() ? &in->v[0] : NULL;
    out->v2 = in->v2.size() ? &in->v2[0] : NULL;
    out->hairs = (ISPCHair*) (in->hairs.size() ? &in->hairs[0] : NULL);
    out->numVertices = in->v.size();
    out->numHairs = in->hairs.size();
    return out;
  }
  
  ISPCMesh* convertMesh (OBJScene::Mesh* in)
  {
    ISPCMesh* out = new ISPCMesh;
    out->positions = in->v.size() ? &in->v[0] : NULL;
    out->positions2 = in->v2.size() ? &in->v2[0] : NULL;
    out->normals = in->vn.size() ? &in->vn[0] : NULL;
    out->texcoords = in->vt.size() ? &in->vt[0] : NULL;
    out->triangles = (ISPCTriangle*) (in->triangles.size() ? &in->triangles[0] : NULL);
    out->quads = (ISPCQuad*) (in->quads.size() ? &in->quads[0] : NULL);
    out->numVertices = in->v.size();
    out->numTriangles = in->triangles.size();
    out->numQuads = in->quads.size();   
    out->geomID = -1;
    return out;
  }

  ISPCSubdivMesh* convertSubdivMesh (OBJScene::SubdivMesh* in)
  {
    ISPCSubdivMesh* out = new ISPCSubdivMesh;
    out->positions = in->positions.size() ? &in->positions[0] : NULL;
    out->normals = in->normals.size() ? &in->normals[0] : NULL;
    out->texcoords = in->texcoords.size() ? &in->texcoords[0] : NULL;
    out->position_indices = in->position_indices.size()   ? &in->position_indices[0] : NULL;
    out->normal_indices = in->normal_indices.size()   ? &in->normal_indices[0] : NULL;
    out->texcoord_indices = in->texcoord_indices.size()   ? &in->texcoord_indices[0] : NULL;
    out->verticesPerFace = in->verticesPerFace.size() ? &in->verticesPerFace[0] : NULL;
    out->holes = in->holes.size() ? &in->holes[0] : NULL;
    out->edge_creases = in->edge_creases.size() ? &in->edge_creases[0] : NULL;
    out->edge_crease_weights = in->edge_crease_weights.size() ? &in->edge_crease_weights[0] : NULL;
    out->vertex_creases = in->vertex_creases.size() ? &in->vertex_creases[0] : NULL;
    out->vertex_crease_weights = in->vertex_crease_weights.size() ? &in->vertex_crease_weights[0] : NULL;
    out->numVertices = in->positions.size();
    out->numFaces = in->verticesPerFace.size();
    out->numEdges = in->position_indices.size();   
    out->numEdgeCreases = in->edge_creases.size();
    out->numVertexCreases = in->vertex_creases.size();
    out->numHoles = in->holes.size();
    out->materialID = in->materialID;
    out->geomID = -1;

    size_t numEdges = in->position_indices.size();
    out->subdivlevel = new float[numEdges];
    for (size_t i=0; i<numEdges; i++) out->subdivlevel[i] = 1.0f;
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

    out->meshes = new ISPCMesh*[in->meshes.size()];
    for (size_t i=0; i<in->meshes.size(); i++) out->meshes[i] = convertMesh(in->meshes[i]);
    out->numMeshes = in->meshes.size();

    out->materials = (ISPCMaterial*) (in->materials.size() ? &in->materials[0] : NULL);
    out->numMaterials = in->materials.size();

    out->hairs = new ISPCHairSet*[in->hairsets.size()];
    for (size_t i=0; i<in->hairsets.size(); i++) out->hairs[i] = convertHair(in->hairsets[i]);
    out->numHairSets = in->hairsets.size();

    out->ambientLights = (ISPCAmbientLight*) (in->ambientLights.size() ? &in->ambientLights[0] : NULL);
    out->numAmbientLights = in->ambientLights.size();

    out->pointLights = (ISPCPointLight*) (in->pointLights.size() ? &in->pointLights[0] : NULL);
    out->numPointLights = in->pointLights.size();

    out->dirLights = (ISPCDirectionalLight*) (in->directionalLights.size() ? &in->directionalLights[0] : NULL);
    out->numDirectionalLights = in->directionalLights.size();

    out->distantLights = (ISPCDistantLight*) (in->distantLights.size() ? &in->distantLights[0] : NULL);
    out->numDistantLights = in->distantLights.size();

    out->subdiv = new ISPCSubdivMesh*[in->subdiv.size()];
    for (size_t i=0; i<in->subdiv.size(); i++) out->subdiv[i] = convertSubdivMesh(in->subdiv[i]);
    out->numSubdivMeshes = in->subdiv.size();

    g_ispc_scene = out;
  }

  bool pick(const float x, const float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p, Vec3fa& hitPos) {
    return device_pick(x,y,vx,vy,vz,p,hitPos);
  }

  void render(const float time, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p) {
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
    alignedFree(g_pixels); 
    g_pixels = NULL;
  }
}
