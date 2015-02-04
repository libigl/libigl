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

#define PARALLEL_COMMIT

/* scene data */
RTCScene g_scene = NULL;
uniform Vec3f* uniform colors = NULL;

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

/* rtcCommitThread called by all ISPC worker threads to enable parallel build */
#if defined(PARALLEL_COMMIT)
task void parallelCommit(RTCScene scene) {
  rtcCommitThread (scene,threadIndex,threadCount); 
}
#endif

/* extended ray structure that includes total transparency along the ray */
struct RTCRay2
{
  Vec3f org;     //!< Ray origin
  Vec3f dir;     //!< Ray direction
  float tnear;   //!< Start of ray segment
  float tfar;    //!< End of ray segment
  float time;    //!< Time of this ray for motion blur.
  int mask;      //!< used to mask out objects during traversal
  Vec3f Ng;      //!< Geometric normal.
  float u;       //!< Barycentric u coordinate of hit
  float v;       //!< Barycentric v coordinate of hit
  int geomID;    //!< geometry ID
  int primID;    //!< primitive ID
  int instID;    //!< instance ID

  // ray extensions
  float transparency; //!< accumulated transparency value
};

/* 3D procedural transparency */
inline float transparencyFunction(RTCRay2& ray)
{
  Vec3f h = ray.org + ray.dir*ray.tfar;
  float v = abs(sin(4.0f*h.x)*cos(4.0f*h.y)*sin(4.0f*h.z));
  float T = clamp((v-0.1f)*3.0f,0.0f,1.0f);
  return T;
}

/* intersection filter function */
void intersectionFilter(void* uniform ptr, RTCRay2& ray)
{
  float T = transparencyFunction(ray);
  if (T >= 1.0f) ray.geomID = RTC_INVALID_GEOMETRY_ID;
  else ray.transparency = T;
}

/* occlusion filter function */
void occlusionFilter(void* uniform ptr, RTCRay2& ray)
{
  float T = transparencyFunction(ray);
  T *= ray.transparency;
  ray.transparency = T;
  if (T != 0.0f) ray.geomID = RTC_INVALID_GEOMETRY_ID;
}

/* adds a cube to the scene */
uniform unsigned int addCube (RTCScene scene_i)
{
  /* create a triangulated cube with 12 triangles and 8 vertices */
  uniform unsigned int mesh = rtcNewTriangleMesh (scene_i, RTC_GEOMETRY_STATIC, 12, 8);

  /* set vertices */
  uniform Vertex* uniform vertices = (uniform Vertex* uniform) rtcMapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 
  vertices[0].x = -1; vertices[0].y = -1; vertices[0].z = -1; 
  vertices[1].x = -1; vertices[1].y = -1; vertices[1].z = +1; 
  vertices[2].x = -1; vertices[2].y = +1; vertices[2].z = -1; 
  vertices[3].x = -1; vertices[3].y = +1; vertices[3].z = +1; 
  vertices[4].x = +1; vertices[4].y = -1; vertices[4].z = -1; 
  vertices[5].x = +1; vertices[5].y = -1; vertices[5].z = +1; 
  vertices[6].x = +1; vertices[6].y = +1; vertices[6].z = -1; 
  vertices[7].x = +1; vertices[7].y = +1; vertices[7].z = +1; 
  rtcUnmapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 

  /* create triangle color array */
  colors = uniform new uniform Vec3f[12];

  /* set triangles and colors */
  uniform int tri = 0;
  uniform Triangle* uniform triangles = (uniform Triangle* uniform) rtcMapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);
  
  // left side
  colors[tri] = make_Vec3f(1,0,0); triangles[tri].v0 = 0; triangles[tri].v1 = 2; triangles[tri].v2 = 1; tri++;
  colors[tri] = make_Vec3f(1,0,0); triangles[tri].v0 = 1; triangles[tri].v1 = 2; triangles[tri].v2 = 3; tri++;

  // right side
  colors[tri] = make_Vec3f(0,1,0); triangles[tri].v0 = 4; triangles[tri].v1 = 5; triangles[tri].v2 = 6; tri++;
  colors[tri] = make_Vec3f(0,1,0); triangles[tri].v0 = 5; triangles[tri].v1 = 7; triangles[tri].v2 = 6; tri++;

  // bottom side
  colors[tri] = make_Vec3f(0.5f);  triangles[tri].v0 = 0; triangles[tri].v1 = 1; triangles[tri].v2 = 4; tri++;
  colors[tri] = make_Vec3f(0.5f);  triangles[tri].v0 = 1; triangles[tri].v1 = 5; triangles[tri].v2 = 4; tri++;

  // top side
  colors[tri] = make_Vec3f(1.0f);  triangles[tri].v0 = 2; triangles[tri].v1 = 6; triangles[tri].v2 = 3; tri++;
  colors[tri] = make_Vec3f(1.0f);  triangles[tri].v0 = 3; triangles[tri].v1 = 6; triangles[tri].v2 = 7; tri++;

  // front side
  colors[tri] = make_Vec3f(0,0,1); triangles[tri].v0 = 0; triangles[tri].v1 = 4; triangles[tri].v2 = 2; tri++;
  colors[tri] = make_Vec3f(0,0,1); triangles[tri].v0 = 2; triangles[tri].v1 = 4; triangles[tri].v2 = 6; tri++;

  // back side
  colors[tri] = make_Vec3f(1,1,0); triangles[tri].v0 = 1; triangles[tri].v1 = 3; triangles[tri].v2 = 5; tri++;
  colors[tri] = make_Vec3f(1,1,0); triangles[tri].v0 = 3; triangles[tri].v1 = 7; triangles[tri].v2 = 5; tri++;

  rtcUnmapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);

  /* set intersection filter for the cube */
  rtcSetIntersectionFilterFunction(scene_i,mesh,(RTCFilterFuncVarying)&intersectionFilter);
  rtcSetOcclusionFilterFunction   (scene_i,mesh,(RTCFilterFuncVarying)&occlusionFilter);

  return mesh;
}

/* adds a ground plane to the scene */
uniform unsigned int addGroundPlane (RTCScene scene_i)
{
  /* create a triangulated plane with 2 triangles and 4 vertices */
  uniform unsigned int mesh = rtcNewTriangleMesh (scene_i, RTC_GEOMETRY_STATIC, 2, 4);

  /* set vertices */
  uniform Vertex* uniform vertices = (uniform Vertex* uniform) rtcMapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 
  vertices[0].x = -10; vertices[0].y = -2; vertices[0].z = -10; 
  vertices[1].x = -10; vertices[1].y = -2; vertices[1].z = +10; 
  vertices[2].x = +10; vertices[2].y = -2; vertices[2].z = -10; 
  vertices[3].x = +10; vertices[3].y = -2; vertices[3].z = +10;
  rtcUnmapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 

  /* set triangles */
  uniform Triangle* uniform triangles = (uniform Triangle* uniform) rtcMapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);
  triangles[0].v0 = 0; triangles[0].v1 = 2; triangles[0].v2 = 1;
  triangles[1].v0 = 1; triangles[1].v1 = 2; triangles[1].v2 = 3;
  rtcUnmapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);

  return mesh;
}

/* called by the C++ code for initialization */
export void device_init (uniform int8* uniform cfg)
{
  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set error handler */
  rtcSetErrorFunction(error_handler);

  /* create scene */
  g_scene = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT_UNIFORM | RTC_INTERSECT_VARYING);

  /* add cube */
  addCube(g_scene);

  /* add ground plane */
  addGroundPlane(g_scene);

  /* commit changes to scene */
#if !defined(PARALLEL_COMMIT)
  rtcCommit (g_scene);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene); sync;
#endif

  /* set start render mode */
  renderPixel = renderPixelStandard;
}

/* task that renders a single screen tile */
Vec3f renderPixelStandard(float x, float y, const uniform Vec3f& vx, const uniform Vec3f& vy, const uniform Vec3f& vz, const uniform Vec3f& p)
{
  float weight = 1.0f;
  Vec3f color = make_Vec3f(0.0f);

  /* initialize ray */
  RTCRay2 primary;
  primary.org = p;
  primary.dir = normalize(x*vx + y*vy + vz);
  primary.tnear = 0.0f;
  primary.tfar = inf;
  primary.geomID = RTC_INVALID_GEOMETRY_ID;
  primary.primID = RTC_INVALID_GEOMETRY_ID;
  primary.mask = -1;
  primary.time = 0;
  primary.transparency = 0.0f;

  while (true)
  {
    /* intersect ray with scene */
    rtcIntersect(g_scene,*((varying RTCRay*)&primary)); // FIXME: use (RTCRay&) cast
    
    /* shade pixels */
    if (primary.geomID == RTC_INVALID_GEOMETRY_ID) 
      break;

    float opacity = 1.0f-primary.transparency;
    Vec3f diffuse = colors[primary.primID];
    Vec3f La = diffuse*0.5f;
    color = color + weight*opacity*La; // FIXME: +=
    Vec3f lightDir = normalize(make_Vec3f(-1,-1,-1));
      
    /* initialize shadow ray */
    RTCRay2 shadow;
    shadow.org = primary.org + primary.tfar*primary.dir;
    shadow.dir = neg(lightDir);
    shadow.tnear = 0.001f;
    shadow.tfar = inf;
    shadow.geomID = RTC_INVALID_GEOMETRY_ID;
    shadow.primID = RTC_INVALID_GEOMETRY_ID;
    shadow.mask = -1;
    shadow.time = 0;
    shadow.transparency = 1.0f;
    
    /* trace shadow ray */
    rtcOccluded(g_scene,*((varying RTCRay*)&shadow)); // FIXME: use (RTCRay&) cast
    
    /* add light contribution */
    if (shadow.geomID) {
      Vec3f Ll = diffuse*shadow.transparency*clamp(-dot(lightDir,normalize(primary.Ng)),0.0f,1.0f);
      color = color + weight*opacity*Ll; // FIXME: +=
    }

    /* shoot transmission ray */
    weight *= primary.transparency;
    primary.tnear = 1.001f*primary.tfar;
    primary.tfar = inf;
    primary.geomID = RTC_INVALID_GEOMETRY_ID;
    primary.primID = RTC_INVALID_GEOMETRY_ID;
    primary.transparency = 0.0f;
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
  const uniform int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const uniform int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  launch[numTilesX*numTilesY] renderTile(pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY); sync;
  rtcDebug();
}

/* called by the C++ code for cleanup */
export void device_cleanup ()
{
  rtcDeleteScene (g_scene);
  delete[] colors;
  rtcExit();
}

