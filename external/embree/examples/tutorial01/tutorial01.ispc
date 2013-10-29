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

#include "../tutorials/tutorials.isph"

const uniform int numPhi = 120; 
//const uniform int numPhi = 400; 
const uniform int numTheta = 2*numPhi;

/* scene data */
uniform RTCGeometry* uniform mesh = NULL;
uniform RTCIntersector* uniform intersector = NULL;

/* called by the C++ code for initialization */
export void init (uniform int verbose)
{
  /* initialize ray tracing core */
  rtcInit();
  rtcStartThreads();
  rtcSetVerbose(verbose);

  /* create a triangulated sphere */
  mesh = rtcNewTriangleMesh (2*numTheta*(numPhi-1), numTheta*(numPhi+1));
  intersector = rtcQueryIntersector (mesh);

  /* set triangles and vertices */
  uniform RTCVertex* uniform vertices = rtcMapPositionBuffer(mesh); 
  uniform RTCTriangle* uniform triangles = rtcMapTriangleBuffer(mesh);
  const uniform float rcpNumTheta = rcp(numTheta);
  const uniform float rcpNumPhi   = rcp(numPhi);
  
  /* create sphere geometry */
  uniform int tri = 0;
  for (uniform int phi=0; phi<=numPhi; phi++)
  {
    for (uniform int theta=0; theta<numTheta; theta++)
    {
      const uniform float phif   = phi*pi*rcpNumPhi;
      const uniform float thetaf = theta*2.0f*pi*rcpNumTheta;
      uniform RTCVertex& v = vertices[phi*numTheta+theta];
      v.x = sin(phif)*sin(thetaf);
      v.y = cos(phif);
      v.z = sin(phif)*cos(thetaf);
    }
    if (phi == 0) continue;

    for (uniform int theta=1; theta<=numTheta; theta++) 
    {
      uniform int p00 = (phi-1)*numTheta+theta-1;
      uniform int p01 = (phi-1)*numTheta+theta%numTheta;
      uniform int p10 = phi*numTheta+theta-1;
      uniform int p11 = phi*numTheta+theta%numTheta;

      if (phi > 1) {
        triangles[tri].v0 = p10; 
        triangles[tri].v1 = p00; 
        triangles[tri].v2 = p01; 
        triangles[tri].id0 = 0; 
        triangles[tri].id1 = tri++;
      }

      if (phi < numPhi) {
        triangles[tri].v0 = p11; 
        triangles[tri].v1 = p10;
        triangles[tri].v2 = p01; 
        triangles[tri].id0 = 0; 
        triangles[tri].id1 = tri++;
      }
    }
  }

  rtcUnmapPositionBuffer(mesh); 
  rtcUnmapTriangleBuffer(mesh); 

  /* build acceleration structure */
  launch rtcBuildAccel (mesh); sync;
}

/* called by the C++ code to set scene */
export void set_scene (uniform Scene* uniform scene) {
}

/* animates the sphere */
task void animateSphere (uniform RTCVertex* uniform vertices, 
                         const uniform float rcpNumTheta,
                         const uniform float rcpNumPhi,
                         const uniform float f)
{
  uniform int phi = taskIndex;
  foreach (theta = 0 ... numTheta)
  {
    uniform RTCVertex* v = &vertices[phi*numTheta+theta];
    const float phif   = phi*pi*rcpNumPhi;
    const float thetaf = theta*2.0f*pi*rcpNumTheta;
    v->x = sin(f*phif)*sin(thetaf);
    v->y = cos(phif);
    v->z = sin(f*phif)*cos(thetaf);
  }
}

/* task that renders a single screen tile */
task void renderTile(uniform int* uniform pixels,
                     const uniform int width,
                     const uniform int height, 
                     const uniform float time,
                     const uniform vec3f& vx, 
                     const uniform vec3f& vy, 
                     const uniform vec3f& vz, 
                     const uniform vec3f& p,
                     const uniform int numTilesX, 
                     const uniform int numTilesY)
{
  const uniform int tileY = taskIndex / numTilesX;
  const uniform int tileX = taskIndex - tileY * numTilesX;
  const uniform int x0 = tileX * TILE_SIZE_X;
  const uniform int x1 = min(x0+TILE_SIZE_X,width);
  const uniform int y0 = tileY * TILE_SIZE_Y;
  const uniform int y1 = min(y0+TILE_SIZE_Y,height);

  foreach (y = y0 ... y1, x = x0 ... x1)
  {
    /* initialize ray */
    Ray ray;
    ray.org = p;
    ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
    ray.tnear = 0.0f;
    ray.tfar = inf;
    ray.id0 = -1;
    ray.id1 = -1;
    ray.mask = -1;
    ray.time = 0;

    /* intersect ray with scene */
    intersector->intersect(intersector,ray);

    /* shade pixels */
    if (ray.id0 != -1) {
      vec3f c = make_vec3f(abs(dot(normalize(ray.Ng),ray.dir)));
      unsigned int r = (unsigned int) (255.0f * c.x);
      unsigned int g = (unsigned int) (255.0f * c.y);
      unsigned int b = (unsigned int) (255.0f * c.z);
      pixels[y*width+x] = (b << 16) + (g << 8) + r;
    }
    else pixels[y*width+x] = 0;
  }
}

/* called by the C++ code to render */
export void render (uniform int* uniform pixels,
                    const uniform int width,
                    const uniform int height, 
                    const uniform float time,
                    const uniform vec3f& vx, 
                    const uniform vec3f& vy, 
                    const uniform vec3f& vz, 
                    const uniform vec3f& p)
{
  /* animate vertices */
  uniform RTCVertex* uniform vertices = rtcMapPositionBuffer(mesh); 
  const uniform float rcpNumTheta = rcp(numTheta);
  const uniform float rcpNumPhi   = rcp(numPhi);
  const uniform float f = 2.0f*(1.0f+0.5f*sin(time));

  /* loop over all vertices */
#if 1
  launch[numPhi+1] animateSphere(vertices,rcpNumTheta,rcpNumPhi,f); sync;
#else
  foreach (phi = 0 ... numPhi+1, theta = 0 ... numTheta)
  {
    uniform RTCVertex* v = &vertices[phi*numTheta+theta];
    const float phif   = phi*pi*rcpNumPhi;
    const float thetaf = theta*2.0f*pi*rcpNumTheta;
    v->x = sin(f*phif)*sin(thetaf);
    v->y = cos(phif);
    v->z = sin(f*phif)*cos(thetaf);
  }
#endif
  rtcUnmapPositionBuffer(mesh); 

  /* build spatial index structure and get intersector */
  launch rtcBuildAccel (mesh); sync;
 
  /* render all pixels */
  const uniform int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const uniform int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  launch[numTilesX*numTilesY] renderTile(pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY); sync;
}

/* called by the C++ code for cleanup */
export void cleanup ()
{
  rtcDeleteIntersector (intersector);
  rtcDeleteGeometry (mesh);
  rtcStopThreads();
  rtcExit();
}
