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

const int numPhi = 5;
const int numTheta = 2*numPhi;

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
RTCScene g_scene1 = NULL;

unsigned int g_instance0 = -1;
unsigned int g_instance1 = -1;
unsigned int g_instance2 = -1;
unsigned int g_instance3 = -1;

Vec3fa colors[4][4];

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

  /* create scene with 4 triangulated spheres */
  g_scene1 = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT1);
  createTriangulatedSphere(g_scene1,Vec3fa( 0, 0,+1),0.5);
  createTriangulatedSphere(g_scene1,Vec3fa(+1, 0, 0),0.5);
  createTriangulatedSphere(g_scene1,Vec3fa( 0, 0,-1),0.5);
  createTriangulatedSphere(g_scene1,Vec3fa(-1, 0, 0),0.5);

#if !defined(PARALLEL_COMMIT)
  rtcCommit (g_scene1);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene1); 
#endif

  /* instantiate geometry */
  g_instance0 = rtcNewInstance(g_scene,g_scene1);
  g_instance1 = rtcNewInstance(g_scene,g_scene1);
  g_instance2 = rtcNewInstance(g_scene,g_scene1);
  g_instance3 = rtcNewInstance(g_scene,g_scene1);
  createGroundPlane(g_scene);

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
  ray.instID = -1;
  ray.mask = -1;
  ray.time = 0;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);
  
  /* shade pixels */
  Vec3fa color = Vec3fa(0.0f);
  if (ray.geomID != RTC_INVALID_GEOMETRY_ID) 
  {
    Vec3fa diffuse = Vec3fa(1,1,1);
    if (ray.instID != -1)
      diffuse = colors[ray.instID][ray.geomID];
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
  /* create identity matrix */
  AffineSpace3f xfm;
  xfm.l.vx = Vec3fa(1,0,0);
  xfm.l.vy = Vec3fa(0,1,0);
  xfm.l.vz = Vec3fa(0,0,1);
  xfm.p    = Vec3fa(0,0,0);
  float t = 0.7f*time;

  /* move instances */
  xfm.p = 2.0f*Vec3fa(+cos(t),0.0f,+sin(t));
  rtcSetTransform(g_scene,g_instance0,RTC_MATRIX_COLUMN_MAJOR,(float*)&xfm);
  xfm.p = 2.0f*Vec3fa(-cos(t),0.0f,-sin(t));
  rtcSetTransform(g_scene,g_instance1,RTC_MATRIX_COLUMN_MAJOR,(float*)&xfm);
  xfm.p = 2.0f*Vec3fa(-sin(t),0.0f,+cos(t));
  rtcSetTransform(g_scene,g_instance2,RTC_MATRIX_COLUMN_MAJOR,(float*)&xfm);
  xfm.p = 2.0f*Vec3fa(+sin(t),0.0f,-cos(t));
  rtcSetTransform(g_scene,g_instance3,RTC_MATRIX_COLUMN_MAJOR,(float*)&xfm);

  /* update scene */
  rtcUpdate(g_scene,g_instance0);
  rtcUpdate(g_scene,g_instance1);
  rtcUpdate(g_scene,g_instance2);
  rtcUpdate(g_scene,g_instance3);
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
  rtcDeleteScene (g_scene1);
  rtcExit();
}
