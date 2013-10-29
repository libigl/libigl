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

#ifndef __EMBREE_SPHERE_H__
#define __EMBREE_SPHERE_H__

#include "../shapes/trianglemesh_full.h"

namespace embree
{
  /*! Implements a sphere shape. */
  class Sphere : public TriangleMeshFull
  {
  public:

    /*! Construction from acceleration structure type. */
    Sphere (const AccelType& ty)
      : TriangleMeshFull(ty) {}

    /*! Construction from parameter container. */
    Sphere (const Parms& parms)
      : TriangleMeshFull(AccelType(parms))
    {
      P    = parms.getVector3f("P");
      dPdt = parms.getVector3f("dPdt");
      r    = parms.getFloat("r");
      numTheta = parms.getInt("numTheta");
      numPhi   = parms.getInt("numPhi");
      triangulate();
    }

  private:

    Vector3f eval(float theta, float phi) const {
      return Vector3f(sinf(theta)*cosf(phi),cosf(theta),sinf(theta)*sinf(phi));
    }

    void triangulate()
    {
      /* triangulate sphere */
      for (size_t theta=0; theta<=numTheta; theta++)
      {
        float rcpNumTheta = rcp(float(numTheta));
        for (size_t phi=0; phi<numPhi; phi++)
        {
          float rcpNumPhi = rcp(float(numPhi));

          Vector3f p = eval(theta*float(pi)*rcpNumTheta,phi*2.0f*float(pi)*rcpNumPhi);
          Vector3f dpdu = eval((theta+0.001f)*float(pi)*rcpNumTheta,phi*2.0f*float(pi)*rcpNumPhi)-p;
          Vector3f dpdv = eval(theta*float(pi)*rcpNumTheta,(phi+0.001f)*2.0f*float(pi)*rcpNumPhi)-p;
          p = r*p+P;

          position.push_back(p);
          if (dPdt != Vector3f(zero)) motion.push_back(dPdt);
          normal.push_back(normalize(cross(dpdv,dpdu)));
          texcoord.push_back(Vec2f(phi*rcpNumPhi,theta*rcpNumTheta));
        }
        if (theta == 0) continue;
        for (size_t phi=1; phi<=numPhi; phi++) {
          size_t p00 = (theta-1)*numPhi+phi-1;
          size_t p01 = (theta-1)*numPhi+phi%numPhi;
          size_t p10 = theta*numPhi+phi-1;
          size_t p11 = theta*numPhi+phi%numPhi;
          if (theta > 1) triangles.push_back(Vector3i((int)p10,(int)p00,(int)p01));
          if (theta < numTheta) triangles.push_back(Vector3i((int)p11,(int)p10,(int)p01));
        }
      }
    }

  public:
    Vector3f P;           //!< center of sphere
    Vector3f dPdt;        //!< motion of sphere with time
    float r;           //!< radius of sphere
    size_t numTheta;   //!< triangulation amount from north to south pole
    size_t numPhi;     //!< triangulation amount around the sphere
  };
}

#endif
