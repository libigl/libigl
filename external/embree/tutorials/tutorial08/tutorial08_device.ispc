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

/* configuration */
#define PARALLEL_COMMIT
#define MIN_EDGE_LEVEL 2.0f
#define MAX_EDGE_LEVEL 64.0f
#define LEVEL_FACTOR 64.0f

/* scene data */
RTCScene g_scene = NULL;

/* render function to use */
renderPixelFunc renderPixel;

/* previous camera position */
uniform Vec3f old_p; 

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

uniform float cube_vertices[8][4] = 
{
  { -1.0f, -1.0f, -1.0f, 0.0f },
  {  1.0f, -1.0f, -1.0f, 0.0f },
  {  1.0f, -1.0f,  1.0f, 0.0f },
  { -1.0f, -1.0f,  1.0f, 0.0f },
  { -1.0f,  1.0f, -1.0f, 0.0f },
  {  1.0f,  1.0f, -1.0f, 0.0f },
  {  1.0f,  1.0f,  1.0f, 0.0f },
  { -1.0f,  1.0f,  1.0f, 0.0f }
};

#if 1

#define NUM_INDICES 24
#define NUM_FACES 6
#define FACE_SIZE 4

uniform unsigned int cube_indices[24] = { 
  0, 1, 5, 4, 
  1, 2, 6, 5, 
  2, 3, 7, 6, 
  0, 4, 7, 3, 
  4, 5, 6, 7, 
  0, 3, 2, 1, 
};

uniform unsigned int cube_faces[6] = { 
  4, 4, 4, 4, 4, 4 
};

#else

#define NUM_INDICES 36
#define NUM_FACES 12
#define FACE_SIZE 3

uniform unsigned int cube_indices[36] = { 
  1, 5, 4,  0, 1, 4,   
  2, 6, 5,  1, 2, 5,
  3, 7, 6,  2, 3, 6,  
  4, 7, 3,  0, 4, 3,
  5, 6, 7,  4, 5, 7,    
  3, 2, 1,  0, 3, 1 
};

uniform unsigned int cube_faces[12] = { 
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
};

#endif

/* adds a cube to the scene */
uniform unsigned int addCube (RTCScene scene_i)
{
  /* create a triangulated cube with 6 quads and 8 vertices */
  uniform unsigned int geomID = rtcNewSubdivisionMesh(scene_i, RTC_GEOMETRY_STATIC, NUM_FACES, NUM_INDICES, 8, 0, 0, 0);

  rtcSetBuffer(scene_i, geomID, RTC_VERTEX_BUFFER, cube_vertices, 0, sizeof(uniform Vec3fa  ));
  rtcSetBuffer(scene_i, geomID, RTC_INDEX_BUFFER,  cube_indices , 0, sizeof(uniform unsigned int));
  rtcSetBuffer(scene_i, geomID, RTC_FACE_BUFFER,   cube_faces,    0, sizeof(uniform unsigned int));

  return geomID;
}

/*! updates the tessellation level for each edge */
void updateEdgeLevelBuffer( RTCScene scene_i, uniform unsigned geomID, const uniform Vec3f& cam_pos )
{
  uniform float*  uniform level    = (uniform float*  uniform) rtcMapBuffer(scene_i, geomID, RTC_LEVEL_BUFFER);
  uniform int*    uniform faces    = (uniform int*    uniform) rtcMapBuffer(scene_i, geomID, RTC_INDEX_BUFFER);
  uniform Vec3fa* uniform vertices = (uniform Vec3fa* uniform) rtcMapBuffer(scene_i, geomID, RTC_VERTEX_BUFFER);
  
  for (uniform size_t f=0; f<NUM_FACES; f++) 
  {
    for (uniform size_t i=0; i<FACE_SIZE; i++) {
      const uniform Vec3f v0 = make_Vec3f(vertices[faces[FACE_SIZE*f+(i+0)%FACE_SIZE]]);
      const uniform Vec3f v1 = make_Vec3f(vertices[faces[FACE_SIZE*f+(i+1)%FACE_SIZE]]);
      const uniform float l  = LEVEL_FACTOR*length(v1-v0)/length(cam_pos-0.5f*(v1+v0));
      level[FACE_SIZE*f+i] = max(min(l,MAX_EDGE_LEVEL),MIN_EDGE_LEVEL);
    }
  }
  rtcUnmapBuffer(scene_i, geomID, RTC_VERTEX_BUFFER);
  rtcUnmapBuffer(scene_i, geomID, RTC_INDEX_BUFFER);
  rtcUnmapBuffer(scene_i, geomID, RTC_LEVEL_BUFFER);

  rtcUpdateBuffer(scene_i,geomID,RTC_LEVEL_BUFFER);
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
  g_scene = rtcNewScene(RTC_SCENE_DYNAMIC | RTC_SCENE_ROBUST,RTC_INTERSECT_UNIFORM | RTC_INTERSECT_VARYING);

  /* add cube */
  addCube(g_scene);

  /* add ground plane */
  addGroundPlane(g_scene);

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
  ray.mask = -1;
  ray.time = 0;
  
  /* intersect ray with scene */
  rtcIntersect(g_scene,ray);
  
  /* shade pixels */
  Vec3f color = make_Vec3f(0.0f);
  if (ray.geomID != RTC_INVALID_GEOMETRY_ID) 
  {
    Vec3f diffuse = ray.geomID == 0 ? make_Vec3f(0.9f,0.6f,0.5f) : make_Vec3f(0.8f,0.0f,0.0f);
    color = color + diffuse*0.5f; // FIXME: +=
    Vec3f lightDir = normalize(make_Vec3f(-1,-1,-1));
    
    /* initialize shadow ray */
    RTCRay shadow;
    shadow.org = ray.org + ray.tfar*ray.dir;
    shadow.dir = neg(lightDir);
    shadow.tnear = 0.001f;
    shadow.tfar = inf;
    shadow.geomID = RTC_INVALID_GEOMETRY_ID;
    shadow.primID = RTC_INVALID_GEOMETRY_ID;
    shadow.mask = -1;
    shadow.time = 0;
    
    /* trace shadow ray */
     rtcOccluded(g_scene,shadow);
             
    /* add light contribution */
    if (shadow.geomID == RTC_INVALID_GEOMETRY_ID)
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
  /* recompute levels */
  updateEdgeLevelBuffer(g_scene,0,p);
    
  /* rebuild scene */
#if !defined(PARALLEL_COMMIT)
  rtcCommit (g_scene);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene); sync;
#endif
  
  /* render image */
  const uniform int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const uniform int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  launch[numTilesX*numTilesY] renderTile(pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY); sync;

  rtcDebug();
}

/* called by the C++ code for cleanup */
export void device_cleanup ()
{
  rtcDeleteScene (g_scene);
  rtcExit();
}
