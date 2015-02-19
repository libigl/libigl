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

#include "../common/tutorial/tutorial_device.isph"

const uniform int numPhi = 5;
const uniform int numTheta = 2*numPhi;

//#define PARALLEL_COMMIT

/* render function to use */
renderPixelFunc renderPixel;

/* error reporting function */
void error_handler(const uniform RTCError code, const uniform int8* uniform str)
{
  print("Embree: ");
  switch (code) {
  case RTC_UNKNOWN_ERROR    : print("RTC_UNKNOWN_ERROR"); break;
  case RTC_INVALID_ARGUMENT : print("RTC_INVALID_ARGUMENT"); break;
  case RTC_INVALID_OPERATION: print("RTC_INVALID_OPERATION"); break;
  case RTC_OUT_OF_MEMORY    : print("RTC_OUT_OF_MEMORY"); break;
  case RTC_UNSUPPORTED_CPU  : print("RTC_UNSUPPORTED_CPU"); break;
  default                   : print("invalid error code"); break;
  }
  if (str) { 
    print(" ("); 
    while (*str) putchar(*str++); 
    print(")\n"); 
  }
  abort();
}

// ======================================================================== //
//                         User defined instancing                          //
// ======================================================================== //

struct Instance 
{
  ALIGNED_STRUCT
  uniform unsigned int geometry;
  RTCScene object;
  uniform int userID;
  AffineSpace3f local2world;
  AffineSpace3f world2local;
  Vec3f lower;
  Vec3f upper;
};

void instanceBoundsFunc(const uniform Instance* uniform instance, uniform size_t item, uniform RTCBounds* uniform bounds_o)
{
  uniform Vec3f l = instance->lower;
  uniform Vec3f u = instance->upper;
  uniform Vec3f p000 = xfmPoint(instance->local2world,make_Vec3f(l.x,l.y,l.z));
  uniform Vec3f p001 = xfmPoint(instance->local2world,make_Vec3f(l.x,l.y,u.z));
  uniform Vec3f p010 = xfmPoint(instance->local2world,make_Vec3f(l.x,u.y,l.z));
  uniform Vec3f p011 = xfmPoint(instance->local2world,make_Vec3f(l.x,u.y,u.z));
  uniform Vec3f p100 = xfmPoint(instance->local2world,make_Vec3f(u.x,l.y,l.z));
  uniform Vec3f p101 = xfmPoint(instance->local2world,make_Vec3f(u.x,l.y,u.z));
  uniform Vec3f p110 = xfmPoint(instance->local2world,make_Vec3f(u.x,u.y,l.z));
  uniform Vec3f p111 = xfmPoint(instance->local2world,make_Vec3f(u.x,u.y,u.z));
  uniform Vec3f lower = min(min(min(p000,p001),min(p010,p011)),min(min(p100,p101),min(p110,p111)));
  uniform Vec3f upper = max(max(max(p000,p001),max(p010,p011)),max(max(p100,p101),max(p110,p111)));
  bounds_o->lower_x = lower.x;
  bounds_o->lower_y = lower.y;
  bounds_o->lower_z = lower.z;
  bounds_o->upper_x = upper.x;
  bounds_o->upper_y = upper.y;
  bounds_o->upper_z = upper.z;
}

void instanceIntersectFunc(const uniform Instance* uniform instance, varying RTCRay& ray, uniform size_t item)
{
  const Vec3f ray_org = ray.org;
  const Vec3f ray_dir = ray.dir;
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

void instanceOccludedFunc(const uniform Instance* uniform instance, varying RTCRay& ray, uniform size_t item)
{
  const Vec3f ray_org = ray.org;
  const Vec3f ray_dir = ray.dir;
  ray.org = xfmPoint (instance->world2local,ray_org);
  ray.dir = xfmVector(instance->world2local,ray_dir);
  rtcOccluded(instance->object,ray);
  ray.org = ray_org;
  ray.dir = ray_dir;
}

uniform Instance* uniform createInstance (RTCScene scene, RTCScene object, uniform int userID, const uniform Vec3f& lower, const uniform Vec3f& upper)
{
  uniform Instance* uniform instance = uniform new uniform Instance;
  instance->object = object;
  instance->userID = userID;
  instance->lower = lower;
  instance->upper = upper;
  instance->local2world.l.vx = make_Vec3f(1,0,0);
  instance->local2world.l.vy = make_Vec3f(0,1,0);
  instance->local2world.l.vz = make_Vec3f(0,0,1);
  instance->local2world.p    = make_Vec3f(0,0,0);
  instance->geometry = rtcNewUserGeometry(scene,1);
  rtcSetUserData(scene,instance->geometry,instance);
  rtcSetBoundsFunction(scene,instance->geometry,(RTCBoundsFunc)&instanceBoundsFunc);
  rtcSetIntersectFunction(scene,instance->geometry,(RTCIntersectFuncVarying)&instanceIntersectFunc);
  rtcSetOccludedFunction (scene,instance->geometry,(RTCOccludedFuncVarying )&instanceOccludedFunc);
  return instance;
}

void updateInstance (RTCScene scene, uniform Instance* uniform instance)
{
  uniform unsigned int geometry = instance->geometry;
  instance->world2local = rcp(instance->local2world);
  rtcUpdate(scene,instance->geometry);
}

// ======================================================================== //
//                     User defined sphere geometry                         //
// ======================================================================== //

struct Sphere
{
  ALIGNED_STRUCT
  Vec3f p;                      //!< position of the sphere
  float r;                      //!< radius of the sphere
  unsigned int geomID;
};

void sphereBoundsFunc(const uniform Sphere* uniform spheres, uniform size_t item, uniform RTCBounds* uniform bounds_o)
{
  const uniform Sphere& sphere = spheres[item];
  bounds_o->lower_x = sphere.p.x-sphere.r;
  bounds_o->lower_y = sphere.p.y-sphere.r;
  bounds_o->lower_z = sphere.p.z-sphere.r;
  bounds_o->upper_x = sphere.p.x+sphere.r;
  bounds_o->upper_y = sphere.p.y+sphere.r;
  bounds_o->upper_z = sphere.p.z+sphere.r;
}

void sphereIntersectFunc(const uniform Sphere* uniform spheres, varying RTCRay& ray, uniform size_t item)
{
  const uniform Sphere& sphere = spheres[item];
  const Vec3f v = ray.org-sphere.p;
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

void sphereOccludedFunc(const uniform Sphere* uniform spheres, varying RTCRay& ray, uniform size_t item)
{
  const uniform Sphere& sphere = spheres[item];
  const Vec3f v = ray.org-sphere.p;
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

uniform Sphere* uniform createAnalyticalSphere (RTCScene scene, const uniform Vec3f& p, uniform float r)
{
  uniform unsigned int geomID = rtcNewUserGeometry(scene,1);
  uniform Sphere* uniform sphere = uniform new uniform Sphere;
  sphere->p = p;
  sphere->r = r;
  sphere->geomID = geomID;
  rtcSetUserData(scene,geomID,sphere);
  rtcSetBoundsFunction(scene,geomID,(RTCBoundsFunc)&sphereBoundsFunc);
  rtcSetIntersectFunction(scene,geomID,(RTCIntersectFuncVarying)&sphereIntersectFunc);
  rtcSetOccludedFunction (scene,geomID,(RTCOccludedFuncVarying )&sphereOccludedFunc);
  return sphere;
}

uniform Sphere* uniform createAnalyticalSpheres (RTCScene scene, uniform size_t N)
{
  uniform unsigned int geomID = rtcNewUserGeometry(scene,N);
  uniform Sphere* uniform spheres = uniform new uniform Sphere[N];
  for (uniform int i=0; i<N; i++) spheres[i].geomID = geomID;
  rtcSetUserData(scene,geomID,spheres);
  rtcSetBoundsFunction(scene,geomID,(RTCBoundsFunc)&sphereBoundsFunc);
  rtcSetIntersectFunction(scene,geomID,(RTCIntersectFuncVarying)&sphereIntersectFunc);
  rtcSetOccludedFunction (scene,geomID,(RTCOccludedFuncVarying )&sphereOccludedFunc);
  return spheres;
}

// ======================================================================== //
//                      Triangular sphere geometry                          //
// ======================================================================== //

uniform unsigned int createTriangulatedSphere (RTCScene scene, const uniform Vec3f& p, uniform float r)
{
  /* create triangle mesh */
  uniform unsigned int mesh = rtcNewTriangleMesh (scene, RTC_GEOMETRY_STATIC, 2*numTheta*(numPhi-1), numTheta*(numPhi+1));
  
  /* map triangle and vertex buffers */
  uniform Vertex* uniform vertices = (uniform Vertex* uniform) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
  uniform Triangle* uniform triangles = (uniform Triangle* uniform) rtcMapBuffer(scene,mesh,RTC_INDEX_BUFFER);

  /* create sphere */
  uniform int tri = 0;
  const uniform float rcpNumTheta = rcp((uniform float)numTheta);
  const uniform float rcpNumPhi   = rcp((uniform float)numPhi);
  for (uniform int phi=0; phi<=numPhi; phi++)
  {
    for (uniform int theta=0; theta<numTheta; theta++)
    {
      const uniform float phif   = phi*pi*rcpNumPhi;
      const uniform float thetaf = theta*2.0f*pi*rcpNumTheta;

      uniform Vertex& v = vertices[phi*numTheta+theta];
      v.x = p.x + r*sin(phif)*sin(thetaf);
      v.y = p.y + r*cos(phif);
      v.z = p.z + r*sin(phif)*cos(thetaf);
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
uniform unsigned int createGroundPlane (RTCScene scene)
{
  /* create a triangulated plane with 2 triangles and 4 vertices */
  uniform unsigned int mesh = rtcNewTriangleMesh (scene, RTC_GEOMETRY_STATIC, 2, 4);

  /* set vertices */
  uniform Vertex* uniform vertices = (uniform Vertex* uniform) rtcMapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 
  vertices[0].x = -10; vertices[0].y = -2; vertices[0].z = -10; 
  vertices[1].x = -10; vertices[1].y = -2; vertices[1].z = +10; 
  vertices[2].x = +10; vertices[2].y = -2; vertices[2].z = -10; 
  vertices[3].x = +10; vertices[3].y = -2; vertices[3].z = +10;
  rtcUnmapBuffer(scene,mesh,RTC_VERTEX_BUFFER); 

  /* set triangles */
  uniform Triangle* uniform triangles = (uniform Triangle* uniform) rtcMapBuffer(scene,mesh,RTC_INDEX_BUFFER);
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

uniform Instance* uniform g_instance0 = NULL;
uniform Instance* uniform g_instance1 = NULL;
uniform Instance* uniform g_instance2 = NULL;
uniform Instance* uniform g_instance3 = NULL;

uniform Vec3f colors[5][4];

/* rtcCommitThread called by all ISPC worker threads to enable parallel build */
#if defined(PARALLEL_COMMIT)
task void parallelCommit(RTCScene scene) {
  rtcCommitThread (scene,threadIndex,threadCount); 
}
#endif

/* called by the C++ code for initialization */
export void device_init (uniform int8* uniform cfg)
{
  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set error handler */
  rtcSetErrorFunction(error_handler);

  /* create scene */
  g_scene = rtcNewScene(RTC_SCENE_DYNAMIC,RTC_INTERSECT_UNIFORM | RTC_INTERSECT_VARYING);

  /* create scene with 4 analytical spheres */
  g_scene0 = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT_UNIFORM | RTC_INTERSECT_VARYING);
  uniform Sphere* uniform spheres = createAnalyticalSpheres(g_scene0,4);
  spheres[0].p = make_Vec3f( 0, 0,+1); spheres[0].r = 0.5f;
  spheres[1].p = make_Vec3f(+1, 0, 0); spheres[1].r = 0.5f;
  spheres[2].p = make_Vec3f( 0, 0,-1); spheres[2].r = 0.5f;
  spheres[3].p = make_Vec3f(-1, 0, 0); spheres[3].r = 0.5f;
#if !defined(PARALLEL_COMMIT)
  rtcCommit(g_scene0);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene0); sync;
#endif

  /* create scene with 4 triangulated spheres */
  g_scene1 = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT_UNIFORM | RTC_INTERSECT_VARYING);
  createTriangulatedSphere(g_scene1,make_Vec3f( 0, 0,+1),0.5);
  createTriangulatedSphere(g_scene1,make_Vec3f(+1, 0, 0),0.5);
  createTriangulatedSphere(g_scene1,make_Vec3f( 0, 0,-1),0.5);
  createTriangulatedSphere(g_scene1,make_Vec3f(-1, 0, 0),0.5);
#if !defined(PARALLEL_COMMIT)
  rtcCommit(g_scene1);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene1); sync;
#endif

  /* create scene with 2 triangulated and 2 analytical spheres */
  g_scene2 = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT_UNIFORM | RTC_INTERSECT_VARYING);
  createTriangulatedSphere(g_scene2,make_Vec3f( 0, 0,+1),0.5);
  createAnalyticalSphere  (g_scene2,make_Vec3f(+1, 0, 0),0.5);
  createTriangulatedSphere(g_scene2,make_Vec3f( 0, 0,-1),0.5);
  createAnalyticalSphere  (g_scene2,make_Vec3f(-1, 0, 0),0.5);
#if !defined(PARALLEL_COMMIT)
  rtcCommit(g_scene2);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene2); sync;
#endif

  /* instantiate geometry */
  createGroundPlane(g_scene);
  g_instance0 = createInstance(g_scene,g_scene0,0,make_Vec3f(-2,-2,-2),make_Vec3f(+2,+2,+2));
  g_instance1 = createInstance(g_scene,g_scene1,1,make_Vec3f(-2,-2,-2),make_Vec3f(+2,+2,+2));
  g_instance2 = createInstance(g_scene,g_scene2,2,make_Vec3f(-2,-2,-2),make_Vec3f(+2,+2,+2));
  g_instance3 = createInstance(g_scene,g_scene2,3,make_Vec3f(-2,-2,-2),make_Vec3f(+2,+2,+2));

  /* set all colors */
  colors[0][0] = make_Vec3f(0.25,0,0);
  colors[0][1] = make_Vec3f(0.50,0,0);
  colors[0][2] = make_Vec3f(0.75,0,0);
  colors[0][3] = make_Vec3f(1.00,0,0);

  colors[1][0] = make_Vec3f(0,0.25,0);
  colors[1][1] = make_Vec3f(0,0.50,0);
  colors[1][2] = make_Vec3f(0,0.75,0);
  colors[1][3] = make_Vec3f(0,1.00,0);

  colors[2][0] = make_Vec3f(0,0,0.25);
  colors[2][1] = make_Vec3f(0,0,0.50);
  colors[2][2] = make_Vec3f(0,0,0.75);
  colors[2][3] = make_Vec3f(0,0,1.00);

  colors[3][0] = make_Vec3f(0.25,0.25,0);
  colors[3][1] = make_Vec3f(0.50,0.50,0);
  colors[3][2] = make_Vec3f(0.75,0.75,0);
  colors[3][3] = make_Vec3f(1.00,1.00,0);

  colors[4][0] = make_Vec3f(1.0,1.0,1.0);
  colors[4][1] = make_Vec3f(1.0,1.0,1.0);
  colors[4][2] = make_Vec3f(1.0,1.0,1.0);
  colors[4][3] = make_Vec3f(1.0,1.0,1.0);

  /* set start render mode */
  renderPixel = renderPixelStandard;
}

/* task that renders a single screen tile */
Vec3f renderPixelStandard(float x, float y, const uniform Vec3f& vx, const uniform Vec3f& vy, const uniform Vec3f& vz, const uniform Vec3f& p)
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
  Vec3f color = make_Vec3f(0.0f);
  if (ray.geomID != RTC_INVALID_GEOMETRY_ID) 
  {
    Vec3f diffuse = make_Vec3f(0.0f);
    if (ray.instID == 0) diffuse = colors[ray.instID][ray.primID];
    else                 diffuse = colors[ray.instID][ray.geomID];
    color = color + diffuse*0.5; // FIXME: +=
    Vec3f lightDir = normalize(make_Vec3f(-1,-1,-1));
    
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
task void renderTile(uniform int* uniform pixels,
                     const uniform int width,
                     const uniform int height, 
                     const uniform float time,
                     const uniform Vec3f& vx, 
                     const uniform Vec3f& vy, 
                     const uniform Vec3f& vz, 
                     const uniform Vec3f& p,
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
    /* calculate pixel color */
    Vec3f color = renderPixel(x,y,vx,vy,vz,p);

    /* write color to framebuffer */
    unsigned int r = (unsigned int) (255.0f * clamp(color.x,0.0f,1.0f));
    unsigned int g = (unsigned int) (255.0f * clamp(color.y,0.0f,1.0f));
    unsigned int b = (unsigned int) (255.0f * clamp(color.z,0.0f,1.0f));
    pixels[y*width+x] = (b << 16) + (g << 8) + r;
  }
}

/* called by the C++ code to render */
export void device_render (uniform int* uniform pixels,
                           const uniform int width,
                           const uniform int height, 
                           const uniform float time,
                           const uniform Vec3f& vx, 
                           const uniform Vec3f& vy, 
                           const uniform Vec3f& vz, 
                           const uniform Vec3f& p)
{
  /* move instances */
  uniform float t = 0.7f*time;
  g_instance0->local2world.p = 2.0f*make_Vec3f(+cos(t),0.0f,+sin(t));
  g_instance1->local2world.p = 2.0f*make_Vec3f(-cos(t),0.0f,-sin(t));
  g_instance2->local2world.p = 2.0f*make_Vec3f(-sin(t),0.0f,+cos(t));
  g_instance3->local2world.p = 2.0f*make_Vec3f(+sin(t),0.0f,-cos(t));
  updateInstance(g_scene,g_instance0);
  updateInstance(g_scene,g_instance1);
  updateInstance(g_scene,g_instance2);
  updateInstance(g_scene,g_instance3);
#if !defined(PARALLEL_COMMIT)
  rtcCommit (g_scene);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene); sync;
#endif

  /* render all pixels */
  const uniform int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const uniform int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  launch[numTilesX*numTilesY] renderTile(pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY); sync;
  rtcDebug();
}

/* called by the C++ code for cleanup */
export void device_cleanup ()
{
  rtcDeleteScene (g_scene);
  rtcDeleteScene (g_scene0);
  rtcDeleteScene (g_scene1);
  rtcDeleteScene (g_scene2);
  rtcExit();
}
