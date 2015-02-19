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

#include "../common/tutorial/tutorial_device.h"



#if 0
const int numSpheres = 1000;
const int numPhi = 5; 
#else
const int numSpheres = 20;
const int numPhi = 120; 
//const int numPhi = 400; 
#endif
const int numTheta = 2*numPhi;

/* scene data */
RTCScene g_scene = NULL;
Vec3fa position[numSpheres];
Vec3fa colors[numSpheres+1];
float radius[numSpheres];
int disabledID = -1;

/* render function to use */
renderPixelFunc renderPixel;

/* error reporting function */
void error_handler(const RTCError code, const int8* str)
{
  printf("Embree: ");
  switch (code) {
  case RTC_UNKNOWN_ERROR    : printf("RTC_UNKNOWN_ERROR"); break;
  case RTC_INVALID_ARGUMENT : printf("RTC_INVALID_ARGUMENT"); break;
  case RTC_INVALID_OPERATION: printf("RTC_INVALID_OPERATION"); break;
  case RTC_OUT_OF_MEMORY    : printf("RTC_OUT_OF_MEMORY"); break;
  case RTC_UNSUPPORTED_CPU  : printf("RTC_UNSUPPORTED_CPU"); break;
  default                   : printf("invalid error code"); break;
  }
  if (str) { 
    printf(" ("); 
    while (*str) putchar(*str++); 
    printf(")\n"); 
  }
  abort();
}

/* adds a sphere to the scene */
unsigned int createSphere (RTCGeometryFlags flags, const Vec3fa& pos, const float r)
{
  /* create a triangulated sphere */
  unsigned int mesh = rtcNewTriangleMesh (g_scene, flags, 2*numTheta*(numPhi-1), numTheta*(numPhi+1));

  /* map triangle and vertex buffer */
  Vertex*   vertices  = (Vertex*  ) rtcMapBuffer(g_scene,mesh,RTC_VERTEX_BUFFER); 
  Triangle* triangles = (Triangle*) rtcMapBuffer(g_scene,mesh,RTC_INDEX_BUFFER);
  
  /* create sphere geometry */
  int tri = 0;
  const float rcpNumTheta = rcp((float)numTheta);
  const float rcpNumPhi   = rcp((float)numPhi);
  for (int phi=0; phi<=numPhi; phi++)
  {
    for (int theta=0; theta<numTheta; theta++)
    {
      const float phif   = phi*float(pi)*rcpNumPhi;
      const float thetaf = theta*2.0f*float(pi)*rcpNumTheta;
      Vertex& v = vertices[phi*numTheta+theta];
      v.x = pos.x + r*sin(phif)*sin(thetaf);
      v.y = pos.y + r*cos(phif);
      v.z = pos.z + r*sin(phif)*cos(thetaf);
    }
    if (phi == 0) continue;

    for (int theta=1; theta<=numTheta; theta++) 
    {
      int p00 = (phi-1)*numTheta+theta-1;
      int p01 = (phi-1)*numTheta+theta%numTheta;
      int p10 = phi*numTheta+theta-1;
      int p11 = phi*numTheta+theta%numTheta;

      if (phi > 1) {
        triangles[tri].v0 = p10; 
        triangles[tri].v1 = p00; 
        triangles[tri].v2 = p01; 
        tri++;
      }

      if (phi < numPhi) {
        triangles[tri].v0 = p11; 
        triangles[tri].v1 = p10;
        triangles[tri].v2 = p01; 
        tri++;
      }
    }
  }
  rtcUnmapBuffer(g_scene,mesh,RTC_VERTEX_BUFFER); 
  rtcUnmapBuffer(g_scene,mesh,RTC_INDEX_BUFFER);

  return mesh;
}

/* rtcCommitThread called by all ISPC worker threads to enable parallel build */
#if defined(PARALLEL_COMMIT)
task void parallelCommit(RTCScene scene) {
  rtcCommitThread (scene,threadIndex,threadCount); 
}
#endif

/* adds a ground plane to the scene */
unsigned int addGroundPlane (RTCScene scene_i)
{
  /* create a triangulated plane with 2 triangles and 4 vertices */
  unsigned int mesh = rtcNewTriangleMesh (scene_i, RTC_GEOMETRY_STATIC, 2, 4);

  /* set vertices */
  Vertex* vertices = (Vertex*) rtcMapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 
  vertices[0].x = -10; vertices[0].y = -2; vertices[0].z = -10; 
  vertices[1].x = -10; vertices[1].y = -2; vertices[1].z = +10; 
  vertices[2].x = +10; vertices[2].y = -2; vertices[2].z = -10; 
  vertices[3].x = +10; vertices[3].y = -2; vertices[3].z = +10;
  rtcUnmapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 

  /* set triangles */
  Triangle* triangles = (Triangle*) rtcMapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);
  triangles[0].v0 = 0; triangles[0].v1 = 2; triangles[0].v2 = 1;
  triangles[1].v0 = 1; triangles[1].v1 = 2; triangles[1].v2 = 3;
  rtcUnmapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);

  return mesh;
}

/* called by the C++ code for initialization */
extern "C" void device_init (int8* cfg)
{
  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set error handler */
  rtcSetErrorFunction(error_handler);

  /* create scene */
  g_scene = rtcNewScene(RTC_SCENE_DYNAMIC,RTC_INTERSECT1);

  /* create some triangulated spheres */
  for (int i=0; i<numSpheres; i++)
  {
    const float phi = i*2.0f*float(pi)/numSpheres;
    const float r = 2.0f*float(pi)/numSpheres;
    const Vec3fa p = 2.0f*Vec3fa(sin(phi),0.0f,-cos(phi));
    //RTCGeometryFlags flags = i%3 == 0 ? RTC_GEOMETRY_STATIC : i%3 == 1 ? RTC_GEOMETRY_DEFORMABLE : RTC_GEOMETRY_DYNAMIC;
    RTCGeometryFlags flags = i%2 ? RTC_GEOMETRY_DEFORMABLE : RTC_GEOMETRY_DYNAMIC;
    //RTCGeometryFlags flags = RTC_GEOMETRY_DEFORMABLE;
    int id = createSphere(flags,p,r);
    position[id] = p;
    radius[id] = r;
    colors[id].x = (i%16+1)/17.0f;
    colors[id].y = (i%8+1)/9.0f;
    colors[id].z = (i%4+1)/5.0f;
  }

  /* add ground plane to scene */
  int id = addGroundPlane(g_scene);
  colors[id] = Vec3fa(1.0f,1.0f,1.0f);

  /* commit changes to scene */
#if !defined(PARALLEL_COMMIT)
  rtcCommit (g_scene);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene); 
#endif

  /* set start render mode */
  renderPixel = renderPixelStandard;
}

/* animates the sphere */
void animateSphere (int taskIndex, Vertex* vertices, 
                         const float rcpNumTheta,
                         const float rcpNumPhi,
                         const Vec3fa& pos, 
                         const float r,
                         const float f)
{
  int phi = taskIndex;
  for (int theta = 0; theta<numTheta; theta++)
  {
    Vertex* v = &vertices[phi*numTheta+theta];
    const float phif   = phi*float(pi)*rcpNumPhi;
    const float thetaf = theta*2.0f*float(pi)*rcpNumTheta;
    v->x = pos.x + r*sin(f*phif)*sin(thetaf);
    v->y = pos.y + r*cos(phif);
    v->z = pos.z + r*sin(f*phif)*cos(thetaf);
  }
}

/* task that renders a single screen tile */
Vec3fa renderPixelStandard(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  /* initialize ray */
  RTCRay ray;
  ray.org = p;
  ray.dir = normalize(x*vx + y*vy + vz);
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = 0;
  
  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);
  
  /* shade pixels */
  Vec3fa color = Vec3fa(0.0f);
  if (ray.geomID != RTC_INVALID_GEOMETRY_ID) 
  {
    Vec3fa diffuse = colors[ray.geomID];
    color = color + diffuse*0.1f; // FIXME: +=
    Vec3fa lightDir = normalize(Vec3fa(-1,-1,-1));
    
    /* initialize shadow ray */
    RTCRay shadow;
    shadow.org = ray.org + ray.tfar*ray.dir;
    shadow.dir = neg(lightDir);
    shadow.tnear = 0.001f;
    shadow.tfar = inf;
    shadow.geomID = 1;
    shadow.primID = 0;
    shadow.mask = -1;
    shadow.time = 0;
    
    /* trace shadow ray */
    rtcOccluded(g_scene,shadow);
    
    /* add light contribution */
    if (shadow.geomID)
      color = color + diffuse*clamp(-dot(lightDir,normalize(ray.Ng)),0.0f,1.0f); // FIXME: +=
  }
  return color;
}

/* task that renders a single screen tile */
void renderTile(int taskIndex, int* pixels,
                     const int width,
                     const int height, 
                     const float time,
                     const Vec3fa& vx, 
                     const Vec3fa& vy, 
                     const Vec3fa& vz, 
                     const Vec3fa& p,
                     const int numTilesX, 
                     const int numTilesY)
{
  const int tileY = taskIndex / numTilesX;
  const int tileX = taskIndex - tileY * numTilesX;
  const int x0 = tileX * TILE_SIZE_X;
  const int x1 = min(x0+TILE_SIZE_X,width);
  const int y0 = tileY * TILE_SIZE_Y;
  const int y1 = min(y0+TILE_SIZE_Y,height);

  for (int y = y0; y<y1; y++) for (int x = x0; x<x1; x++)
  {
    /* calculate pixel color */
    Vec3fa color = renderPixel(x,y,vx,vy,vz,p);

    /* write color to framebuffer */
    unsigned int r = (unsigned int) (255.0f * clamp(color.x,0.0f,1.0f));
    unsigned int g = (unsigned int) (255.0f * clamp(color.y,0.0f,1.0f));
    unsigned int b = (unsigned int) (255.0f * clamp(color.z,0.0f,1.0f));
    pixels[y*width+x] = (b << 16) + (g << 8) + r;
  }
}

/* animates a sphere */
void animateSphere (int id, float time)
{
  /* animate vertices */
  Vertex* vertices = (Vertex*) rtcMapBuffer(g_scene,id,RTC_VERTEX_BUFFER); 
  const float rcpNumTheta = rcp((float)numTheta);
  const float rcpNumPhi   = rcp((float)numPhi);
  const Vec3fa pos = position[id];
  const float r = radius[id];
  const float f = 2.0f*(1.0f+0.5f*sin(time));

  /* loop over all vertices */
#if 1 // enables parallel execution
  launch_animateSphere(animateSphere,numPhi+1,vertices,rcpNumTheta,rcpNumPhi,pos,r,f); 
#else
  for (int phi = 0; phi <numPhi+1; phi++) for (int theta = 0; theta<numTheta; theta++)
  {
    Vertex* v = &vertices[phi*numTheta+theta];
    const float phif   = phi*float(pi)*rcpNumPhi;
    const float thetaf = theta*2.0f*float(pi)*rcpNumTheta;
    v->x = pos.x+r*sin(f*phif)*sin(thetaf);
    v->y = pos.y+r*cos(phif);
    v->z = pos.z+r*sin(f*phif)*cos(thetaf);
  }
#endif

  rtcUnmapBuffer(g_scene,id,RTC_VERTEX_BUFFER); 

  /* update mesh */
  rtcUpdate (g_scene,id);
}

/* called by the C++ code to render */
extern "C" void device_render (int* pixels,
                           const int width,
                           const int height, 
                           const float time,
                           const Vec3fa& vx, 
                           const Vec3fa& vy, 
                           const Vec3fa& vz, 
                           const Vec3fa& p)
{
  /* animate sphere */
  for (int i=0; i<numSpheres; i++)
    animateSphere(i,time+i);

  /* commit changes to scene */
#if !defined(PARALLEL_COMMIT)
  rtcCommit (g_scene);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene); 
#endif

 
  /* render all pixels */
  const int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  launch_renderTile(numTilesX*numTilesY,pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY); 
  rtcDebug();
}

/* called by the C++ code for cleanup */
extern "C" void device_cleanup ()
{
  rtcDeleteScene (g_scene);
  rtcExit();
}
