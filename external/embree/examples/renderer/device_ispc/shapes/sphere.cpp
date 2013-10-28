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

#include "sphere.h"
#include "trianglemesh_ispc.h"

namespace embree
{
  Vector3f eval(float theta, float phi) {
    return Vector3f(sinf(theta)*cosf(phi),cosf(theta),sinf(theta)*sinf(phi));
  }
  
  void* Sphere::create (const Parms& parms)
  {
    Vector3f P    = parms.getVector3f("P");
    float r       = parms.getFloat("r");
    size_t numTheta = parms.getInt("numTheta");
    size_t numPhi   = parms.getInt("numPhi");

    size_t allocatedVertices = (numTheta+1)*numPhi;
    size_t allocatedTriangles = 2*numTheta*numPhi;

    size_t numVertices = 0;
    Vec3fa* position = new Vec3fa[allocatedVertices];
    Vec3fa* motion = NULL;        //!< Motion array.
    Vec3fa* normal = new Vec3fa[allocatedVertices];
    Vec2f* texcoord = new Vec2f[allocatedVertices];
    size_t numTriangles = 0;
    Vec4i* triangles = new Vec4i[allocatedTriangles];

    /* triangulate sphere */
    for (size_t theta=0; theta<=numTheta; theta++)
    {
      float rcpNumTheta = rcp(float(numTheta));
      for (size_t phi=0; phi<numPhi; phi++)
      {
        float rcpNumPhi = rcp(float(numPhi));
        
        Vector3f d = eval(theta*float(pi)*rcpNumTheta,phi*2.0f*float(pi)*rcpNumPhi);
        //Vector3f dpdu = eval((theta+0.001f)*float(pi)*rcpNumTheta,phi*2.0f*float(pi)*rcpNumPhi)-p;
        //Vector3f dpdv = eval(theta*float(pi)*rcpNumTheta,(phi+0.001f)*2.0f*float(pi)*rcpNumPhi)-p;
        Vector3f p = r*d+P;
        
        position[numVertices] = p;
        //if (dPdt != Vector3f(zero)) motion.push_back(dPdt);
        normal[numVertices] = normalize(d);//cross(dpdv,dpdu));
        texcoord[numVertices] = Vec2f(phi*rcpNumPhi,theta*rcpNumTheta);
        numVertices++;
      }
      if (theta == 0) continue;
      for (size_t phi=1; phi<=numPhi; phi++) {
        size_t p00 = (theta-1)*numPhi+phi-1;
        size_t p01 = (theta-1)*numPhi+phi%numPhi;
        size_t p10 = theta*numPhi+phi-1;
        size_t p11 = theta*numPhi+phi%numPhi;
        if (theta > 1) triangles[numTriangles++] = Vec4i((int)p10,(int)p00,(int)p01,0);
        if (theta < numTheta) triangles[numTriangles++] = Vec4i((int)p11,(int)p10,(int)p01,0);
      }
    }

    if (numVertices  > allocatedVertices ) throw std::runtime_error("internal error");
    if (numTriangles > allocatedTriangles) throw std::runtime_error("internal error");
    
    return  ispc::TriangleMesh__new((ispc::vec3fa*)position,
                                    (ispc::vec3fa*)motion,
                                    (ispc::vec3fa*)normal,
                                    (ispc::vec2f*)texcoord,
                                    numVertices,
                                    (ispc::vec4i*)triangles,
                                    numTriangles);
  }
}
