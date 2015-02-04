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

#include "tutorial_device.h"

/* the scene to render */
extern RTCScene g_scene;

/* global subdivision level for subdivision geometry */
unsigned int g_subdivision_levels = 0;

/* intensity scaling for traversal cost visualization */
float scale = 0.001f;
extern "C" bool g_changed = false;

/* stores pointer to currently used rendePixel function */
extern renderPixelFunc renderPixel;

extern float g_debug;

/* standard rendering function for each tutorial */
Vec3fa renderPixelStandard(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p);

/* renders a single pixel with eyelight shading */
Vec3fa renderPixelEyeLight(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
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
  ray.time = g_debug;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);
  else return Vec3fa(embree::abs(dot(ray.dir,normalize(ray.Ng))));
}

__noinline void setray(RTCRay& ray)
{
  ray.u = ray.v = 0.001f;
  ray.Ng = Vec3fa(0,1,0);
  ray.geomID = 0;
  ray.primID = 0;
}

/* renders a single pixel with wireframe shading */
Vec3fa renderPixelWireframe(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
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
  ray.time = g_debug;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* return black if nothing hit */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);

  /* calculate wireframe around triangles */
  const float border = 0.05f;
  Vec3fa color = Vec3fa(1.0f);
  if (ray.u < border) color = Vec3fa(0.0f);
  if (ray.v < border) color = Vec3fa(0.0f);
  if (1.0f-ray.u-ray.v < border) color = Vec3fa(0.0f);

  /* perform eyelight shading */
  return color*Vec3fa(embree::abs(dot(ray.dir,normalize(ray.Ng))));
}

/* renders a single pixel with UV shading */
Vec3fa renderPixelUV(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
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
  ray.time = g_debug;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);
  else return Vec3fa(ray.u,ray.v,1.0f-ray.u-ray.v);
}

/* renders a single pixel with geometry normal shading */
Vec3fa renderPixelNg(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
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
  ray.time = g_debug;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);
  else return abs(normalize(Vec3fa(ray.Ng.x,ray.Ng.y,ray.Ng.z)));
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
Vec3fa renderPixelGeomID(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
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
  ray.time = g_debug;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);
  else return randomColor(ray.geomID);
}

/* geometry ID and primitive ID shading */
Vec3fa renderPixelGeomIDPrimID(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
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
  ray.time = g_debug;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f);
  else return randomColor(ray.geomID ^ ray.primID);
}

/* vizualizes the traversal cost of a pixel */
Vec3fa renderPixelCycles(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
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
  ray.time = g_debug;

  /* intersect ray with scene */
  int64 c0 = get_tsc();
  rtcIntersect(g_scene,ray);
  int64 c1 = get_tsc();
  
  /* shade pixel */
  return Vec3fa((float)(c1-c0)*scale,0.0f,0.0f);
}

/* renders a single pixel with UV shading */
Vec3fa renderPixelUV16(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
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
  ray.time = g_debug;

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
                            const Vec3fa& vx, 
                            const Vec3fa& vy, 
                            const Vec3fa& vz, 
                            const Vec3fa& p,
                            Vec3fa& hitPos)
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
  ray.time = g_debug;

  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);
  PRINT(ray.geomID);
  PRINT(ray.primID);

  /* shade pixel */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) {
    hitPos = Vec3fa(0.0f,0.0f,0.0f);
    return false;
  }
  else {
    hitPos = ray.org + ray.tfar*ray.dir;
    return true;
  }
}

/* called when a key is pressed */
extern "C" void device_key_pressed(int key)
{
  if (key == GLUT_KEY_F1) {
    renderPixel = renderPixelStandard;
    g_changed = true;
  }
  else if (key == GLUT_KEY_F2) {
    renderPixel = renderPixelEyeLight;
    g_changed = true;
  }    
  else if (key == GLUT_KEY_F3) {
    renderPixel = renderPixelWireframe;
    g_changed = true;
  }
  else if (key == GLUT_KEY_F4) {
    renderPixel = renderPixelUV;
    g_changed = true;
  }
  else if (key == GLUT_KEY_F5) {
    renderPixel = renderPixelNg;
    g_changed = true;
  }
  else if (key == GLUT_KEY_F6) {
    renderPixel = renderPixelGeomID;
    g_changed = true;
  }
  else if (key == GLUT_KEY_F7) {
    renderPixel = renderPixelGeomIDPrimID;
    g_changed = true;
  }
  else if (key == GLUT_KEY_F8) {
    renderPixel = renderPixelUV16;
    g_changed = true;
  }
  else if (key == GLUT_KEY_F9) {
    if (renderPixel == renderPixelCycles) scale *= 1.1f;
    renderPixel = renderPixelCycles;
    g_changed = true;
  }
  else if (key == GLUT_KEY_F10) {
    if (renderPixel == renderPixelCycles) scale *= 0.9f;
    renderPixel = renderPixelCycles;
    g_changed = true;
  }
  else if (key == GLUT_KEY_F11) {
    if (g_subdivision_levels > 0)	
      g_subdivision_levels--;
    PRINT(g_subdivision_levels);
    g_changed = true;
  }
  else if (key == GLUT_KEY_F12) {
    g_subdivision_levels++;
    PRINT(g_subdivision_levels);
    g_changed = true;
  }

}

void renderTile(int taskIndex,
                int* pixels,
                const int width,
                const int height, 
                const float time,
                const Vec3fa& vx, 
                const Vec3fa& vy, 
                const Vec3fa& vz, 
                const Vec3fa& p,
                const int numTilesX, 
                const int numTilesY);

struct RenderTileTask
{
  RenderTileTask (int* pixels, const int width, const int height, const float time, 
                  const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p, const int numTilesX, const int numTilesY)
    : pixels(pixels), width(width), height(height), time(time), vx(vx), vy(vy), vz(vz), p(p), numTilesX(numTilesX), numTilesY(numTilesY) {}

public:
  int* pixels;
  const int width;
  const int height;
  const float time;
  const Vec3fa vx;
  const Vec3fa vy;
  const Vec3fa vz;
  const Vec3fa p;
  const int numTilesX;
  const int numTilesY;
};

void renderTile_parallel(RenderTileTask* task, size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) {
  renderTile(taskIndex,task->pixels,task->width,task->height,task->time,task->vx,task->vy,task->vz,task->p,task->numTilesX,task->numTilesY);
}

void launch_renderTile (int numTiles, 
                        int* pixels, const int width, const int height, const float time, 
                        const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p, const int numTilesX, const int numTilesY)
{
  TaskScheduler::EventSync event;
  RenderTileTask parms(pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY);
  TaskScheduler::Task task(&event,(TaskScheduler::runFunction)renderTile_parallel,&parms,numTiles,NULL,NULL,"render");
  TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_FRONT,&task);
  event.sync();
}

typedef void (*animateSphereFunc) (int taskIndex, Vertex* vertices, 
				   const float rcpNumTheta,
				   const float rcpNumPhi,
				   const Vec3fa& pos, 
				   const float r,
				   const float f);

struct AnimateSphereTask
{
  AnimateSphereTask (animateSphereFunc func,
		     Vertex* vertices, 
		     const float rcpNumTheta,
		     const float rcpNumPhi,
		     const Vec3fa& pos, 
		     const float r,
		     const float f)
    : func(func), vertices(vertices), rcpNumTheta(rcpNumTheta), rcpNumPhi(rcpNumPhi), pos(pos), r(r), f(f) {}

public:
  animateSphereFunc func;
  Vertex* vertices;
  const float rcpNumTheta;
  const float rcpNumPhi;
  const Vec3fa pos;
  const float r;
  const float f;
};

void animateSphere_parallel(AnimateSphereTask* task, size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) {
  task->func(taskIndex,task->vertices,task->rcpNumTheta,task->rcpNumPhi,task->pos,task->r,task->f);
}

void launch_animateSphere(animateSphereFunc func,
			  int taskSize,
			  Vertex* vertices, 
			  const float rcpNumTheta,
			  const float rcpNumPhi,
			  const Vec3fa& pos, 
			  const float r,
			  const float f)
{
  TaskScheduler::EventSync event;
  AnimateSphereTask parms(func,vertices,rcpNumTheta,rcpNumPhi,pos,r,f);
  TaskScheduler::Task task(&event,(TaskScheduler::runFunction)animateSphere_parallel,&parms,taskSize,NULL,NULL,"render");
  TaskScheduler::addTask(-1,TaskScheduler::GLOBAL_FRONT,&task);
  event.sync();
}

static int p[513] = { 
  151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
  190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
  88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
  77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
  102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
  135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
  5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
  223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
  129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
  251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
  49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
  138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
  151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
  190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
  88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
  77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
  102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
  135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
  5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
  223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
  129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
  251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
  49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
  138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
  151
};

static float g1[128] = {
  -0.20707, 0.680971, -0.293328, -0.106833, -0.362614, 0.772857, -0.968834, 0.16818, -0.681263, -0.232568, 0.382009, -0.882282, 0.799709, -0.672908, -0.681857, 0.0661294, 0.208288, 0.165398, -0.460058, -0.219044, -0.413199, 0.484755, -0.402949, -0.848924, -0.190035, 0.714756, 0.883937, 0.325661, 0.692952, -0.99449, -0.0752415, 0.065192, 0.575753, -0.468776, 0.965505, -0.38643, 0.20171, 0.217431, -0.575122, 0.77179, -0.390686, -0.69628, -0.324676, -0.225046, 0.28722, 0.507107, 0.207232, 0.0632565, -0.0812794, 0.304977, -0.345638, 0.892741, -0.26392, 0.887781, -0.985143, 0.0331999, -0.454458, -0.951402, 0.183909, -0.590073, 0.755387, -0.881263, -0.478315, -0.394342, 0.78299, -0.00360388, 0.420051, -0.427172, 0.729847, 0.351081, -0.0830201, 0.919271, 0.549351, -0.246897, -0.542722, -0.290932, -0.399364, 0.339532, 0.437933, 0.131909, 0.648931, -0.218776, 0.637533, 0.688017, -0.639064, 0.886792, -0.150226, 0.0413316, -0.868712, 0.827016, 0.765169, 0.522728, -0.202155, 0.376514, 0.523097, -0.189982, -0.749498, -0.0307322, -0.555075, 0.746242, 0.0576438, -0.997172, 0.721028, -0.962605, 0.629784, -0.514232, -0.370856, 0.931465, 0.87112, 0.618863, -0.0157817, -0.559729, 0.152707, -0.421942, -0.357866, -0.477353, -0.652024, -0.996365, -0.910432, -0.517651, -0.169098, 0.403249, -0.556309, 0.00782069, -0.86594, -0.213873, -0.0410469, -0.563716
};

static float g2[128*2] = {
0.605609, 0.538399, 0.796519, -0.944204, 0.908294, 0.756016, 0.0977536, -0.863638, 0.842196, -0.744751, -0.932081, 0.932392, -0.588525, 0.516884, 0.841188, -0.978497, -0.608649, -0.868011, 0.992137, -0.772425, 0.963049, -0.0478757, 0.953878, 0.889467, 0.562174, 0.624644, -0.356598, -0.520726, -0.821833, 0.99985, 0.234183, -0.9791, -0.971815, -0.0979374, -0.108159, -0.34927, -0.592124, -0.775632, 0.97228, 0.753819, 0.941608, 0.578291, 0.852108, -0.760312, -0.784772, 0.0223242, -0.606013, -0.980319, 0.252581, -0.575064, 0.884701, 0.943763, 0.737344, 0.938496, 0.0466562, -0.994566, 0.989782, 0.988368, -0.546155, 0.279211, -0.69504, 0.931229, 0.99768, -0.325874, -0.630157, -0.999936, -0.968623, -0.226805, -0.750428, -0.450961, 0.257868, 0.968011, -0.988005, -0.713965, 0.991007, -0.61059, 0.950437, -0.483042, -0.98105, -0.915356, -0.892527, -0.772958, -0.9081, 0.55692, 0.906075, 0.937419, 0.454624, -0.991582, 0.400857, 0.855933, -0.672619, 0.0713424, 0.593249, -0.378286, -0.997369, -0.827112, 0.708222, -0.995343, 0.985069, 0.698711, -0.180105, 0.999961, -0.768451, 0.993107, -0.918024, 0.0446961, 0.91882, 0.97691, -0.393915, 0.364803, 0.0495592, 0.186545, -0.461553, -0.242776, 0.901952, -0.0710866, 0.888101, 0.999935, 0.277688, 0.0554235, 0.506599, -0.299293, 0.984394, -0.999698, 0.408822, -0.782639, 0.128596, 0.198834, 0.981707, 0.864566, 0.808197, 0.352335, 0.970484, -0.667503, 0.330243, 0.208392, 0.191539, -0.938943, 0.895002, 0.910575, -0.537691, -0.98548, -0.721635, -0.335382, -0.424701, -0.960452, 0.595047, 0.783579, -0.937749, 0.529096, -0.997906, -0.581313, -0.899828, -0.88461, 0.989469, 0.91872, -0.850793, 0.955954, 0.715768, -0.736686, 0.80392, -0.717276, -0.788579, 0.987003, -0.839648, 0.885176, -0.998929, -0.0376033, -0.578371, -0.718771, 0.906081, 0.239947, -0.803563, -0.00826282, 0.991011, -0.0057943, -0.349232, 0.65319, 0.992067, -0.953535, 0.893781, 0.661689, 0.957253, -0.425442, -0.866609, 0.712892, -0.807777, 0.89632, -0.595147, -0.0224999, -0.643786, 0.545815, -0.870124, -0.696306, -0.99902, 0.773648, -0.806008, -0.931319, -0.780114, -0.552154, -0.933812, -0.563108, -0.619909, 0.966532, 0.692454, 0.993284, 0.338885, -0.75104, 0.237272, -0.713619, -0.160187, -0.199242, -0.371265, -0.781439, -0.914125, -0.944104, 0.169525, -0.984403, 0.976056, -0.265228, 0.94232, 0.993906, -0.877517, -0.89618, 0.611817, -0.106758, 0.680403, 0.163329, -0.325386, -0.0687362, -0.901164, 0.460314, 0.999981, -0.0408026, 0.850356, -0.763343, -0.170806, -0.102919, 0.581564, 0.688634, 0.284368, -0.276419, 0.616641, -0.929771, 0.927865, 0.440373, 0.153446, 0.840456, 0.996966, 0.867209, -0.135077, -0.493238, -0.577193, 0.0588088, 0.715215, 0.0143633
};

static float g3[128*4] = {
  -0.582745, 0.443494, -0.680971, 0, -0.601153, 0.791961, 0.106833, 0, -0.265466, 0.576385, -0.772857, 0, 0.981035, 0.0963612, -0.16818, 0, 0.524388, 0.819103, 0.232568, 0, -0.170518, -0.43875, 0.882282, 0, 0.598053, -0.435348, 0.672908, 0, 0.53956, 0.839346, -0.0661294, 0, -0.782511, -0.600267, -0.165398, 0, -0.122114, 0.968043, 0.219044, 0, -0.235567, 0.842331, -0.484755, 0, -0.158657, 0.504139, 0.848924, 0, -0.578396, 0.39317, -0.714756, 0, 0.883328, -0.337159, -0.325661, 0, 0.0597264, -0.0861552, 0.99449, 0, -0.970124, 0.233685, -0.0651921, 0, 0.208238, -0.858421, 0.468776, 0, 0.916908, -0.0997567, 0.38643, 0, -0.786568, -0.577957, -0.217431, 0, 0.14868, 0.618251, -0.77179, 0, -0.24168, 0.675858, 0.69628, 0, -0.50994, 0.83025, 0.225046, 0, -0.534183, -0.676382, -0.507107, 0, -0.793861, -0.6048, -0.0632565, 0, -0.92148, 0.240548, -0.304977, 0, -0.210037, 0.39862, -0.892741, 0, -0.310918, 0.339375, -0.887781, 0, 0.99836, 0.0466305, -0.0331999, 0, -0.0439099, 0.304806, 0.951402, 0, -0.676304, -0.440938, 0.590073, 0, 0.339805, -0.328495, 0.881263, 0, -0.0625568, 0.916832, 0.394342, 0, 0.776463, -0.630153, 0.00360388, 0, -0.224717, -0.8758, 0.427172, 0, 0.618879, -0.70266, -0.351081, 0, -0.380313, 0.101503, -0.919271, 0, 0.149639, -0.957418, 0.246897, 0, 0.128024, 0.948139, 0.290932, 0, -0.292448, 0.893976, -0.339532, 0, -0.192062, -0.972477, -0.131909, 0, 0.44007, -0.870905, 0.218776, 0, 0.303887, -0.659003, -0.688017, 0, 0.195552, 0.41876, -0.886792, 0, -0.889922, 0.454236, -0.0413315, 0, 0.515034, 0.225353, -0.827016, 0, 0.63084, -0.573408, -0.522728, 0, -0.745779, 0.549592, -0.376514, 0, 0.0711763, -0.979204, 0.189982, 0, 0.705657, 0.707887, 0.0307322, 0, 0.114603, 0.655735, -0.746242, 0, -0.0739232, -0.0135353, 0.997172, 0, 0.173356, -0.20818, 0.962605, 0, 0.34008, -0.787344, 0.514232, 0, -0.143596, 0.334295, -0.931465, 0, 0.721989, -0.30942, -0.618863, 0, -0.827657, 0.0410685, 0.559729, 0, -0.804277, -0.418454, 0.421942, 0, -0.379459, 0.792556, 0.477353, 0, 0.0391537, 0.0756503, 0.996365, 0, 0.821943, 0.237588, 0.517651, 0, -0.788974, 0.463584, -0.403249, 0, 0.175972, 0.984364, -0.00782073, 0, 0.891497, 0.399363, 0.213873, 0, -0.819111, 0.106216, 0.563716, 0, 0.105511, 0.544028, -0.832406, 0, -0.464551, 0.63753, 0.614612, 0, 0.232387, 0.935154, -0.267363, 0, 0.777619, 0.272068, -0.566823, 0, 0.975331, 0.190338, 0.111807, 0, 0.224313, 0.450072, -0.86436, 0, 0.841897, -0.536898, 0.0543103, 0, 0.637123, -0.664145, -0.391135, 0, 0.901675, -0.422984, 0.0898189, 0, -0.496241, 0.367413, -0.786608, 0, -0.255468, -0.689763, -0.677469, 0, -0.0616459, -0.951141, -0.302539, 0, -0.431011, -0.889035, -0.154425, 0, -0.0711688, 0.486502, -0.870776, 0, -0.223359, -0.36162, 0.905175, 0, -0.678546, 0.695482, -0.23639, 0, 0.576553, 0.77934, 0.245389, 0, -0.194568, -0.24951, 0.948624, 0, 0.28962, -0.447736, 0.845962, 0, -0.0403821, -0.871893, 0.488028, 0, 0.790972, -0.560788, 0.244705, 0, -0.34553, 0.739953, 0.57713, 0, -0.516376, -0.697122, 0.49737, 0, 0.115998, 0.859293, 0.498156, 0, 0.643831, -0.239955, 0.72657, 0, -0.125114, 0.987348, -0.0974144, 0, -0.306452, 0.610699, -0.73016, 0, -0.269845, 0.893027, -0.360119, 0, 0.328563, -0.570628, -0.752615, 0, -0.306918, -0.42057, 0.853769, 0, 0.699245, -0.51785, 0.492837, 0, -0.558362, -0.469763, -0.68378, 0, 0.476563, -0.841398, 0.254826, 0, 0.0276172, -0.623206, 0.78157, 0, 0.587723, -0.800313, -0.118659, 0, 0.594035, -0.740708, 0.313806, 0, -0.340185, -0.887929, 0.309605, 0, 0.312245, -0.246681, -0.917416, 0, 0.194206, 0.186398, -0.963089, 0, 0.915704, 0.329835, -0.229553, 0, 0.94133, 0.229917, 0.247055, 0, -0.888253, -0.144148, 0.436152, 0, -0.906917, -0.362625, -0.214486, 0, 0.403108, -0.908884, 0.10693, 0, 0.983963, 0.169256, 0.056292, 0, -0.197949, 0.888236, 0.414553, 0, 0.0879741, 0.247673, 0.964841, 0, 0.474384, -0.868071, -0.146331, 0, 0.699884, 0.541342, -0.465953, 0, 0.610965, 0.567249, 0.552223, 0, 0.830508, -0.285788, -0.478103, 0, 0.328573, -0.683076, -0.652263, 0, -0.00537775, 0.873381, 0.487009, 0, -0.51289, 0.828835, 0.223557, 0, -0.871168, -0.15102, 0.467182, 0, -0.545561, 0.390016, -0.741789, 0, 0.874063, 0.259258, 0.410852, 0, -0.781555, 0.612184, -0.120005, 0, -0.284928, 0.708938, -0.645154, 0, -0.568809, 0.0883274, 0.817713, 0, -0.0429388, 0.549957, -0.834088, 0, 0.933296, -0.127233, 0.335813, 0, 0.698149, -0.493464, 0.51873, 0, -0.603413, 0.617495, -0.504572, 0 
};

__forceinline float fade(float t) {
  return (t * t * t) * (t * (t * 6 - 15) + 10); 
}

/*__forceinline float fade(float t) {
  return t * t * (3 - t * 2); 
  }*/

__forceinline float lerp(float t, float a, float b) { 
  return a + t * (b - a); 
}

__forceinline float grad(int hash, float x) {
  return x*g1[hash&127];
}

__forceinline float grad(int hash, float x, float y) {
  int h = hash&127;
  return x*g2[2*h+0]+y*g2[2*h+1];
}

__forceinline float grad(int hash, float x, float y, float z) {
  int h = hash&127;
  return x*g3[4*h+0]+y*g3[4*h+1]+z*g3[4*h+2];
}

float noise(float x) 
{
  float fx = floorf(x);
  int X = (int)fx & 255;
  x -= fx;
  float u = fade(x);
  float g00 = grad(p[X  ],x  );
  float g10 = grad(p[X+1],x-1);
  return lerp(u,g00,g10);
}

float noise(float x, float y) 
{
  float fx = floorf(x);
  float fy = floorf(y);
  
  int X = (int)fx & 255;
  int Y = (int)fy & 255;
  
  x -= fx;
  y -= fy;
  
  float u = fade(x);
  float v = fade(y);
  
  int p00  = p[X  ]+Y;
  int p10  = p[X+1]+Y;
  int p01  = p[X  ]+Y+1;
  int p11  = p[X+1]+Y+1;
  
  float g00 = grad(p[p00],x  ,y  );
  float g10 = grad(p[p10],x-1,y  );
  float g01 = grad(p[p01],x  ,y-1);
  float g11 = grad(p[p11],x-1,y-1);
  
  return lerp(v,lerp(u,g00,g10),lerp(u,g01,g11));
}

float noise(const Vec3fa& pos)  
{
  float x = pos.x;
  float y = pos.y;
  float z = pos.z;

  float fx = floorf(x);
  float fy = floorf(y);
  float fz = floorf(z);
  
  int X = (int)fx & 255;
  int Y = (int)fy & 255;
  int Z = (int)fz & 255;
  
  x -= fx;
  y -= fy;
  z -= fz;
  
  float u = fade(x);
  float v = fade(y);
  float w = fade(z);
  
  int p00  = p[X]+Y;
  int p000 = p[p00]+Z;
  int p010 = p[p00+1]+Z;
  int p001 = p000+1; 
  int p011 = p010+1;
  int p10  = p[X+1]+Y;
  int p100 = p[p10]+Z;
  int p110 = p[p10+1]+Z;
  int p101 = p100+1;
  int p111 = p110+1;
  
  float g000 = grad(p[p000],x  ,y  ,z  );
  float g100 = grad(p[p100],x-1,y  ,z  );
  float g010 = grad(p[p010],x  ,y-1,z  );
  float g110 = grad(p[p110],x-1,y-1,z  );
  float g001 = grad(p[p001],x  ,y  ,z-1);
  float g101 = grad(p[p101],x-1,y  ,z-1);
  float g011 = grad(p[p011],x  ,y-1,z-1);
  float g111 = grad(p[p111],x-1,y-1,z-1);
  
  return lerp(w,
              lerp(v,lerp(u,g000,g100),lerp(u,g010,g110)),
              lerp(v,lerp(u,g001,g101),lerp(u,g011,g111)));
}

Vec3fa noise3D(const Vec3fa& p)
{
  float x = noise(p.x+128.0f);
  float y = noise(p.y+64.0f);
  float z = noise(p.z+192.0f);
  return Vec3fa(x,y,z);
}
