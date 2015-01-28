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

#include "tutorial_device.h"

/* the scene to render */
extern RTCScene g_scene;

/* intensity scaling for traversal cost visualization */
float scale = 0.001f;

/* stores pointer to currently used rendePixel function */
extern renderPixelFunc renderPixel;

/* standard rendering function for each tutorial */
Vec3fa renderPixelStandard(int x, int y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p);

/* renders a single pixel with eyelight shading */
Vec3fa renderPixelEyeLight(int x, int y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  /* initialize ray */
  RTCRay ray;
  ray.org = p;
  ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = 0;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);
  else return Vec3fa(embree::abs(dot(ray.dir,normalize(ray.Ng))));
}

/* renders a single pixel with ambient occlusion shading */
Vec3fa renderPixelAmbientOcclusion(int x, int y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  /* initialize ray */
  RTCRay ray;
  ray.org = p;
  ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = 0;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* return black if nothing hit */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);

  /* calculate hit point */
  float intensity = 0;
  Vec3fa hitPos = add(ray.org,mul(ray.tfar,ray.dir));

  /* trace some ambient occlusion rays */
  int seed = 34*x+12*y;
  for (int i=0; i<32; i++) 
  {
    Vec3fa dir; 
    const float oneOver10000f = 1.f/10000.f;
    seed = 1103515245 * seed + 12345;
    dir.x = (seed%10000)*oneOver10000f;
    seed = 1103515245 * seed + 12345;
    dir.y = (seed%10000)*oneOver10000f;
    seed = 1103515245 * seed + 12345;
    dir.z = (seed%10000)*oneOver10000f;
    
    /* initialize shadow ray */
    RTCRay shadow;
    shadow.org = hitPos;
    shadow.dir = dir;
    shadow.tnear = 0.001f;
    shadow.tfar = inf;
    shadow.geomID = 1;
    shadow.primID = 0;
    shadow.mask = -1;
    shadow.time = 0;
    
    /* trace shadow ray */
    rtcOccluded(g_scene,shadow);
    
    /* add light contribution */
    intensity += shadow.geomID;
  }
  intensity *= 1.0f/32.0f;

  /* shade pixel */
  return Vec3fa(intensity);
}

/* renders a single pixel with UV shading */
Vec3fa renderPixelUV(int x, int y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  /* initialize ray */
  RTCRay ray;
  ray.org = p;
  ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = 0;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);
  else return Vec3fa(ray.u,ray.v,1.0f-ray.u-ray.v);
}

/* renders a single pixel with geometry normal shading */
Vec3fa renderPixelNg(int x, int y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  /* initialize ray */
  RTCRay ray;
  ray.org = p;
  ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = 0;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);
  else return abs(Vec3fa(ray.Ng.x,ray.Ng.y,ray.Ng.z));
}

Vec3fa randomColor(const int ID) 
{
  int r = ((ID+13)*17*23) & 255;
  int g = ((ID+15)*11*13) & 255;
  int b = ((ID+17)* 7*19) & 255;
  const float oneOver255f = 1.f/255.f;
  return Vec3fa(r*oneOver255f,g*oneOver255f,b*oneOver255f);
}

/* geometry ID shading */
Vec3fa renderPixelGeomID(int x, int y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  /* initialize ray */
  RTCRay ray;
  ray.org = p;
  ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = 0;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);
  else return randomColor(ray.geomID);
}

/* geometry ID and primitive ID shading */
Vec3fa renderPixelGeomIDPrimID(int x, int y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  /* initialize ray */
  RTCRay ray;
  ray.org = p;
  ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = 0;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);
  else return randomColor(ray.geomID ^ ray.primID);
}

/* vizualizes the traversal cost of a pixel */
Vec3fa renderPixelCycles(int x, int y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  /* initialize ray */
  RTCRay ray;
  ray.org = p;
  ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = 0;

  /* intersect ray with scene */
  int64 c0 = get_tsc();
  rtcIntersect(g_scene,ray);
  int64 c1 = get_tsc();
  
  /* shade pixel */
  return Vec3fa((float)(c1-c0)*scale,0.0f,0.0f);
}

/* renders a single pixel with UV shading */
Vec3fa renderPixelUV16(int x, int y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  /* initialize ray */
  RTCRay ray;
  ray.org = p;
  ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = 0;

  /* intersect ray with scene */
  for (int i=0; i<16; i++) {
    ray.tfar = inf;
    rtcIntersect(g_scene,ray);
  }

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);
  else return Vec3fa(ray.u,ray.v,1.0f-ray.u-ray.v);
}

/* returns the point seen through specified pixel */
extern "C" bool device_pick(const float x,
                            const float y, 
                            const Vec3f& vx, 
                            const Vec3f& vy, 
                            const Vec3f& vz, 
                            const Vec3f& p,
                            Vec3f& hitPos)
{
  /* initialize ray */
  RTCRay ray;
  ray.org = p;
  ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
  ray.tnear = 0.0f;
  ray.tfar = inf;
  ray.geomID = RTC_INVALID_GEOMETRY_ID;
  ray.primID = RTC_INVALID_GEOMETRY_ID;
  ray.mask = -1;
  ray.time = 0;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) {
    hitPos = Vec3fa(0.0f,0.0f,0.0f);
    return false;
  }
  else {
    hitPos = add(ray.org,mul(ray.tfar,ray.dir));
    return true;
  }
}

/* called when a key is pressed */
extern "C" void device_key_pressed(int key)
{
  if      (key == GLUT_KEY_F1) renderPixel = renderPixelStandard;
  else if (key == GLUT_KEY_F2) renderPixel = renderPixelEyeLight;
  else if (key == GLUT_KEY_F3) renderPixel = renderPixelAmbientOcclusion;
  else if (key == GLUT_KEY_F4) renderPixel = renderPixelUV;
  else if (key == GLUT_KEY_F5) renderPixel = renderPixelNg;
  else if (key == GLUT_KEY_F6) renderPixel = renderPixelGeomID;
  else if (key == GLUT_KEY_F7) renderPixel = renderPixelGeomIDPrimID;
  else if (key == GLUT_KEY_F8) renderPixel = renderPixelUV16;
  else if (key == GLUT_KEY_F9) {
    if (renderPixel == renderPixelCycles) scale *= 1.1f;
    renderPixel = renderPixelCycles;
  }
  else if (key == GLUT_KEY_F10) {
    if (renderPixel == renderPixelCycles) scale *= 0.9f;
    renderPixel = renderPixelCycles;
  }
}

void renderTile(int taskIndex,
                int* pixels,
                const int width,
                const int height, 
                const float time,
                const Vec3f& vx, 
                const Vec3f& vy, 
                const Vec3f& vz, 
                const Vec3f& p,
                const int numTilesX, 
                const int numTilesY);

struct RenderTileTask
{
  RenderTileTask (int* pixels, const int width, const int height, const float time, 
                  const Vec3f& vx, const Vec3f& vy, const Vec3f& vz, const Vec3f& p, const int numTilesX, const int numTilesY)
    : pixels(pixels), width(width), height(height), time(time), vx(vx), vy(vy), vz(vz), p(p), numTilesX(numTilesX), numTilesY(numTilesY) {}

public:
  int* pixels;
  const int width;
  const int height;
  const float time;
  const Vec3f vx;
  const Vec3f vy;
  const Vec3f vz;
  const Vec3f p;
  const int numTilesX;
  const int numTilesY;
};

void renderTile_parallel(RenderTileTask* task, size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) {
  renderTile(taskIndex,task->pixels,task->width,task->height,task->time,task->vx,task->vy,task->vz,task->p,task->numTilesX,task->numTilesY);
}

void launch_renderTile (int numTiles, 
                        int* pixels, const int width, const int height, const float time, 
                        const Vec3f& vx, const Vec3f& vy, const Vec3f& vz, const Vec3f& p, const int numTilesX, const int numTilesY)
{
  TaskScheduler::EventSync event;
  RenderTileTask parms(pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY);
  TaskScheduler::Task task(&event,(TaskScheduler::runFunction)renderTile_parallel,&parms,numTiles,NULL,NULL,"render");
  TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_FRONT,&task);
  event.sync();
}
