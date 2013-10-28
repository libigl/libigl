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

#include "trianglemesh.h"
#include "trianglemesh_ispc.h"

namespace embree
{
  void* ISPCTriangleMesh::create (const Parms& parms)
  {
    int numVertices = 0;
    Vec3fa* position = NULL;      //!< Position array.
    Vec3fa* motion = NULL;        //!< Motion array.
    Vec3fa* normal = NULL;        //!< Normal array (can be empty).
    Vec2f* texcoord = NULL;      //!< Texture coordinates array (can be empty).
    int numTriangles = 0;
    Vec4i* triangles = NULL;  //!< Triangle indices array.
    
    if (Variant v = parms.getData("positions")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong position format");
      position = new Vec3fa[v.data->size()];
      for (size_t i=0; i<v.data->size(); i++) position[i] = v.data->getVector3f(i);
      numVertices = v.data->size();
    }
    if (Variant v = parms.getData("motions")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong motion vector format");
      motion = new Vec3fa[v.data->size()];
      for (size_t i=0; i<v.data->size(); i++) motion[i] = v.data->getVector3f(i);
    }
    if (Variant v = parms.getData("normals")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong normal format");
      normal = new Vec3fa[v.data->size()];
      for (size_t i=0; i<v.data->size(); i++) normal[i] = v.data->getVector3f(i);
    }
    if (Variant v = parms.getData("texcoords")) {
      if (!v.data || v.type != Variant::FLOAT2) throw std::runtime_error("wrong texcoords0 format");
      texcoord = new Vec2f[v.data->size()];
      for (size_t i=0; i<v.data->size(); i++) texcoord[i] = v.data->getVec2f(i);
    }
    if (Variant v = parms.getData("texcoords0")) {
      if (!v.data || v.type != Variant::FLOAT2) throw std::runtime_error("wrong texcoords0 format");
      if (texcoord) delete[] texcoord;
      texcoord = new Vec2f[v.data->size()];
      for (size_t i=0; i<v.data->size(); i++) texcoord[i] = v.data->getVec2f(i);
    }
    if (Variant v = parms.getData("indices")) {
      if (!v.data || v.type != Variant::INT3) throw std::runtime_error("wrong triangle format");
      triangles = new Vec4i[v.data->size()];
      for (size_t i=0; i<v.data->size(); i++) {
        Vector3i t = v.data->getVector3i(i);
        triangles[i] = Vec4i(t.x,t.y,t.z,0);
      }
      numTriangles = v.data->size();
    }

    void* mesh = ispc::TriangleMesh__new((ispc::vec3fa*)position,
                                    (ispc::vec3fa*)motion,
                                    (ispc::vec3fa*)normal,
                                    (ispc::vec2f*)texcoord,
                                    numVertices,
                                    (ispc::vec4i*)triangles,
                                    numTriangles);

    delete position;
    delete motion;
    delete normal;
    delete texcoord;
    delete triangles;
	return mesh;
  }
}
