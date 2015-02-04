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

const int numPhi = 5;
const int numTheta = 2*numPhi;

//

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

// ======================================================================== //
//                         User defined instancing                          //
// ======================================================================== //

struct Instance 
{
  ALIGNED_STRUCT
  unsigned int geometry;
  RTCScene object;
  int userID;
  AffineSpace3f local2world;
  AffineSpace3f world2local;
  Vec3fa lower;
  Vec3fa upper;
};

void instanceBoundsFunc(const Instance* instance, size_t item, RTCBounds* bounds_o)
{
  Vec3fa l = instance->lower;
  Vec3fa u = instance->upper;
  Vec3fa p000 = xfmPoint(instance->local2world,Vec3fa(l.x,l.y,l.z));
  Vec3fa p001 = xfmPoint(instance->local2world,Vec3fa(l.x,l.y,u.z));
  Vec3fa p010 = xfmPoint(instance->local2world,Vec3fa(l.x,u.y,l.z));
  Vec3fa p011 = xfmPoint(instance->local2world,Vec3fa(l.x,u.y,u.z));
  Vec3fa p100 = xfmPoint(instance->local2world,Vec3fa(u.x,l.y,l.z));
  Vec3fa p101 = xfmPoint(instance->local2world,Vec3fa(u.x,l.y,u.z));
  Vec3fa p110 = xfmPoint(instance->local2world,Vec3fa(u.x,u.y,l.z));
  Vec3fa p111 = xfmPoint(instance->local2world,Vec3fa(u.x,u.y,u.z));
  Vec3fa lower = min(min(min(p000,p001),min(p010,p011)),min(min(p100,p101),min(p110,p111)));
  Vec3fa upper = max(max(max(p000,p001),max(p010,p011)),max(max(p100,p101),max(p110,p111)));
  bounds_o->lower_x = lower.x;
  bounds_o->lower_y = lower.y;
  bounds_o->lower_z = lower.z;
  bounds_o->upper_x = upper.x;
  bounds_o->upper_y = upper.y;
  bounds_o->upper_z = upper.z;
}

void instanceIntersectFunc(const Instance* instance, RTCRay& ray, size_t item)
{
  const Vec3fa ray_org = ray.org;
  const Vec3fa ray_dir = ray.dir;
  const int geomID = ray.geomID;
  ray.org = xfmPoint (instance->world2local,ray_org);
  ray.dir = xfmVector(instance->world2local,ray_dir);
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  rtcIntersect(instance->object,ray);
  ray.org = ray_org;
  ray.dir = ray_dir;
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) ray.geomID = geomID;
  else ray.instID = instance->userID;
}

void instanceOccludedFunc(const Instance* instance, RTCRay& ray, size_t item)
{
  const Vec3fa ray_org = ray.org;
  const Vec3fa ray_dir = ray.dir;
  ray.org = xfmPoint (instance->world2local,ray_org);
  ray.dir = xfmVector(instance->world2local,ray_dir);
  rtcOccluded(instance->object,ray);
  ray.org = ray_org;
  ray.dir = ray_dir;
}

Instance* createInstance (RTCScene scene, RTCScene object, int userID, const Vec3fa& lower, const Vec3fa& upper)
{
  Instance* instance = new Instance;
  instance->object = object;
  instance->userID = userID;
  instance->lower = lower;
  instance->upper = upper;
  instance->local2world.l.vx = Vec3fa(1,0,0);
  instance->local2world.l.vy = Vec3fa(0,1,0);
  instance->local2world.l.vz = Vec3fa(0,0,1);
  instance->local2world.p    = Vec3fa(0,0,0);
  instance->geometry = rtcNewUserGeometry(scene,1);
  rtcSetUserData(scene,instance->geometry,instance);
  rtcSetBoundsFunction(scene,instance->geometry,(RTCBoundsFunc)&instanceBoundsFunc);
  rtcSetIntersectFunction(scene,instance->geometry,(RTCIntersectFunc)&instanceIntersectFunc);
  rtcSetOccludedFunction (scene,instance->geometry,(RTCOccludedFunc )&instanceOccludedFunc);
  return instance;
}

void updateInstance (RTCScene scene, Instance* instance)
{
  unsigned int geometry = instance->geometry;
  instance->world2local = rcp(instance->local2world);
  rtcUpdate(scene,instance->geometry);
}

// ======================================================================== //
//                     User defined sphere geometry                         //
// ======================================================================== //

struct Sphere
{
  ALIGNED_STRUCT
  Vec3fa p;                      //!< position of the sphere
  float r;                      //!< radius of the sphere
  unsigned int geomID;
};

void sphereBoundsFunc(const Sphere* spheres, size_t item, RTCBounds* bounds_o)
{
  const Sphere& sphere = spheres[item];
  bounds_o->lower_x = sphere.p.x-sphere.r;
  bounds_o->lower_y = sphere.p.y-sphere.r;
  bounds_o->lower_z = sphere.p.z-sphere.r;
  bounds_o->upper_x = sphere.p.x+sphere.r;
  bounds_o->upper_y = sphere.p.y+sphere.r;
  bounds_o->upper_z = sphere.p.z+sphere.r;
}

void sphereIntersectFunc(const Sphere* spheres, RTCRay& ray, size_t item)
{
  const Sphere& sphere = spheres[item];
  const Vec3fa v = ray.org-sphere.p;
  const float A = dot(ray.dir,ray.dir);
  const float B = 2.0f*dot(v,ray.dir);
  const float C = dot(v,v) - sqr(sphere.r);
  const float D = B*B - 4.0f*A*C;
  if (D < 0.0f) return;
  const float Q = sqrt(D);
  const float rcpA = rcp(A);
  const float t0 = 0.5f*rcpA*(-B-Q);
  const float t1 = 0.5f*rcpA*(-B+Q);
  if ((ray.tnear < t0) & (t0 < ray.tfar)) {
    ray.u = 0.0f;
    ray.v = 0.0f;
    ray.tfar = t0;
    ray.geomID = sphere.geomID;
    ray.primID = item;
    ray.Ng = ray.org+t0*ray.dir-sphere.p;
  }
  if ((ray.tnear < t1) & (t1 < ray.tfar)) {
    ray.u = 0.0f;
    ray.v = 0.0f;
    ray.tfar = t1;
    ray.geomID = sphere.geomID;
    ray.primID = item;
    ray.Ng = ray.org+t1*ray.dir-sphere.p;
  }
}

void sphereOccludedFunc(const Sphere* spheres, RTCRay& ray, size_t item)
{
  const Sphere& sphere = spheres[item];
  const Vec3fa v = ray.org-sphere.p;
  const float A = dot(ray.dir,ray.dir);
  const float B = 2.0f*dot(v,ray.dir);
  const float C = dot(v,v) - sqr(sphere.r);
  const float D = B*B - 4.0f*A*C;
  if (D < 0.0f) return;
  const float Q = sqrt(D);
  const float rcpA = rcp(A);
  const float t0 = 0.5f*rcpA*(-B-Q);
  const float t1 = 0.5f*rcpA*(-B+Q);
  if ((ray.tnear < t0) & (t0 < ray.tfar)) {
    ray.geomID = 0;
  }
  if ((ray.tnear < t1) & (t1 < ray.tfar)) {
    ray.geomID = 0;
  }
}

Sphere* createAnalyticalSphere (RTCScene scene, const Vec3fa& p, float r)
{
  unsigned int geomID = rtcNewUserGeometry(scene,1);
  Sphere* sphere = new Sphere;
  sphere->p = p;
  sphere->r = r;
  sphere->geomID = geomID;
  rtcSetUserData(scene,geomID,sphere);
  rtcSetBoundsFunction(scene,geomID,(RTCBoundsFunc)&sphereBoundsFunc);
  rtcSetIntersectFunction(scene,geomID,(RTCIntersectFunc)&sphereIntersectFunc);
  rtcSetOccludedFunction (scene,geomID,(RTCOccludedFunc )&sphereOccludedFunc);
  return sphere;
}

Sphere* createAnalyticalSpheres (RTCScene scene, size_t N)
{
  unsigned int geomID = rtcNewUserGeometry(scene,N);
  Sphere* spheres = new Sphere[N];
  for (int i=0; i<N; i++) spheres[i].geomID = geomID;
  rtcSetUserData(scene,geomID,spheres);
  rtcSetBoundsFunction(scene,geomID,(RTCBoundsFunc)&sphereBoundsFunc);
  rtcSetIntersectFunction(scene,geomID,(RTCIntersectFunc)&sphereIntersectFunc);
  rtcSetOccludedFunction (scene,geomID,(RTCOccludedFunc )&sphereOccludedFunc);
  return spheres;
}

// ======================================================================== //
//                      Triangular sphere geometry                          //
// ======================================================================== //

unsigned int createTriangulatedSphere (RTCScene scene, const Vec3fa& p, float r)
{
  /* create triangle mesh */
  unsigned int mesh = rtcNewTriangleMesh (scene, RTC_GEOMETRY_STATIC, 2*numTheta*(numPhi-1), numTheta*(numPhi+1));
  
  /* map triangle and vertex buffers */
  Vertex* vertices = (Vertex*) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
  Triangle* triangles = (Triangle*) rtcMapBuffer(scene,mesh,RTC_INDEX_BUFFER);

  /* create sphere */
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
      v.x = p.x + r*sin(phif)*sin(thetaf);
      v.y = p.y + r*cos(phif);
      v.z = p.z + r*sin(phif)*cos(thetaf);
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
  rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
  rtcUnmapBuffer(scene,mesh,RTC_INDEX_BUFFER);
  return mesh;
}

/* creates a ground plane */
unsigned int createGroundPlane (RTCScene scene)
{
  /* create a triangulated plane with 2 triangles and 4 vertices */
  unsigned int mesh = rtcNewTriangleMesh (scene, RTC_GEOMETRY_STATIC, 2, 4);

  /* set vertices */
  Vertex* vertices = (Vertex*) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
  vertices[0].x = -10; vertices[0].y = -2; vertices[0].z = -10; 
  vertices[1].x = -10; vertices[1].y = -2; vertices[1].z = +10; 
  vertices[2].x = +10; vertices[2].y = -2; vertices[2].z = -10; 
  vertices[3].x = +10; vertices[3].y = -2; vertices[3].z = +10;
  rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 

  /* set triangles */
  Triangle* triangles = (Triangle*) rtcMapBuffer(scene,mesh,RTC_INDEX_BUFFER);
  triangles[0].v0 = 0; triangles[0].v1 = 2; triangles[0].v2 = 1;
  triangles[1].v0 = 1; triangles[1].v1 = 2; triangles[1].v2 = 3;
  rtcUnmapBuffer(scene,mesh,RTC_INDEX_BUFFER);

  return mesh;
}

/* scene data */
RTCScene g_scene  = NULL;
RTCScene g_scene0 = NULL;
RTCScene g_scene1 = NULL;
RTCScene g_scene2 = NULL;

Instance* g_instance0 = NULL;
Instance* g_instance1 = NULL;
Instance* g_instance2 = NULL;
Instance* g_instance3 = NULL;

Vec3fa colors[5][4];

/* rtcCommitThread called by all ISPC worker threads to enable parallel build */
#if defined(PARALLEL_COMMIT)
task void parallelCommit(RTCScene scene) {
  rtcCommitThread (scene,threadIndex,threadCount); 
}
#endif

/* called by the C++ code for initialization */
extern "C" void device_init (int8* cfg)
{
  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set error handler */
  rtcSetErrorFunction(error_handler);

  /* create scene */
  g_scene = rtcNewScene(RTC_SCENE_DYNAMIC,RTC_INTERSECT1);

  /* create scene with 4 analytical spheres */
  g_scene0 = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT1);
  Sphere* spheres = createAnalyticalSpheres(g_scene0,4);
  spheres[0].p = Vec3fa( 0, 0,+1); spheres[0].r = 0.5f;
  spheres[1].p = Vec3fa(+1, 0, 0); spheres[1].r = 0.5f;
  spheres[2].p = Vec3fa( 0, 0,-1); spheres[2].r = 0.5f;
  spheres[3].p = Vec3fa(-1, 0, 0); spheres[3].r = 0.5f;
#if !defined(PARALLEL_COMMIT)
  rtcCommit(g_scene0);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene0); 
#endif

  /* create scene with 4 triangulated spheres */
  g_scene1 = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT1);
  createTriangulatedSphere(g_scene1,Vec3fa( 0, 0,+1),0.5);
  createTriangulatedSphere(g_scene1,Vec3fa(+1, 0, 0),0.5);
  createTriangulatedSphere(g_scene1,Vec3fa( 0, 0,-1),0.5);
  createTriangulatedSphere(g_scene1,Vec3fa(-1, 0, 0),0.5);
#if !defined(PARALLEL_COMMIT)
  rtcCommit(g_scene1);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene1); 
#endif

  /* create scene with 2 triangulated and 2 analytical spheres */
  g_scene2 = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT1);
  createTriangulatedSphere(g_scene2,Vec3fa( 0, 0,+1),0.5);
  createAnalyticalSphere  (g_scene2,Vec3fa(+1, 0, 0),0.5);
  createTriangulatedSphere(g_scene2,Vec3fa( 0, 0,-1),0.5);
  createAnalyticalSphere  (g_scene2,Vec3fa(-1, 0, 0),0.5);
#if !defined(PARALLEL_COMMIT)
  rtcCommit(g_scene2);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene2); 
#endif

  /* instantiate geometry */
  createGroundPlane(g_scene);
  g_instance0 = createInstance(g_scene,g_scene0,0,Vec3fa(-2,-2,-2),Vec3fa(+2,+2,+2));
  g_instance1 = createInstance(g_scene,g_scene1,1,Vec3fa(-2,-2,-2),Vec3fa(+2,+2,+2));
  g_instance2 = createInstance(g_scene,g_scene2,2,Vec3fa(-2,-2,-2),Vec3fa(+2,+2,+2));
  g_instance3 = createInstance(g_scene,g_scene2,3,Vec3fa(-2,-2,-2),Vec3fa(+2,+2,+2));

  /* set all colors */
  colors[0][0] = Vec3fa(0.25,0,0);
  colors[0][1] = Vec3fa(0.50,0,0);
  colors[0][2] = Vec3fa(0.75,0,0);
  colors[0][3] = Vec3fa(1.00,0,0);

  colors[1][0] = Vec3fa(0,0.25,0);
  colors[1][1] = Vec3fa(0,0.50,0);
  colors[1][2] = Vec3fa(0,0.75,0);
  colors[1][3] = Vec3fa(0,1.00,0);

  colors[2][0] = Vec3fa(0,0,0.25);
  colors[2][1] = Vec3fa(0,0,0.50);
  colors[2][2] = Vec3fa(0,0,0.75);
  colors[2][3] = Vec3fa(0,0,1.00);

  colors[3][0] = Vec3fa(0.25,0.25,0);
  colors[3][1] = Vec3fa(0.50,0.50,0);
  colors[3][2] = Vec3fa(0.75,0.75,0);
  colors[3][3] = Vec3fa(1.00,1.00,0);

  colors[4][0] = Vec3fa(1.0,1.0,1.0);
  colors[4][1] = Vec3fa(1.0,1.0,1.0);
  colors[4][2] = Vec3fa(1.0,1.0,1.0);
  colors[4][3] = Vec3fa(1.0,1.0,1.0);

  /* set start render mode */
  renderPixel = renderPixelStandard;
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
  ray.instID = 4; // set default instance ID
  ray.mask = -1;
  ray.time = 0;
  
  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);
  
  /* shade pixels */
  Vec3fa color = Vec3fa(0.0f);
  if (ray.geomID != RTC_INVALID_GEOMETRY_ID) 
  {
    Vec3fa diffuse = Vec3fa(0.0f);
    if (ray.instID == 0) diffuse = colors[ray.instID][ray.primID];
    else                 diffuse = colors[ray.instID][ray.geomID];
    color = color + diffuse*0.5; // FIXME: +=
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
  /* move instances */
  float t = 0.7f*time;
  g_instance0->local2world.p = 2.0f*Vec3fa(+cos(t),0.0f,+sin(t));
  g_instance1->local2world.p = 2.0f*Vec3fa(-cos(t),0.0f,-sin(t));
  g_instance2->local2world.p = 2.0f*Vec3fa(-sin(t),0.0f,+cos(t));
  g_instance3->local2world.p = 2.0f*Vec3fa(+sin(t),0.0f,-cos(t));
  updateInstance(g_scene,g_instance0);
  updateInstance(g_scene,g_instance1);
  updateInstance(g_scene,g_instance2);
  updateInstance(g_scene,g_instance3);
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
  rtcDeleteScene (g_scene0);
  rtcDeleteScene (g_scene1);
  rtcDeleteScene (g_scene2);
  rtcExit();
}
