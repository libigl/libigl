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

#include "../common/tutorial/tutorial_device.h"

struct ISPCTriangle 
{
  int v0;                /*< first triangle vertex */
  int v1;                /*< second triangle vertex */
  int v2;                /*< third triangle vertex */
  int materialID;        /*< material of triangle */
};

struct ISPCMaterial
{
  int illum;             /*< illumination model */
  
  float d;               /*< dissolve factor, 1=opaque, 0=transparent */
  float Ns;              /*< specular exponent */
  float Ni;              /*< optical density for the surface (index of refraction) */
  
  Vec3f Ka;              /*< ambient reflectivity */
  Vec3f Kd;              /*< diffuse reflectivity */
  Vec3f Ks;              /*< specular reflectivity */
  Vec3f Tf;              /*< transmission filter */
};

struct ISPCMesh
{
  Vec3fa* positions;    //!< vertex position array
  Vec3fa* normals;       //!< vertex normal array
  Vec2f* texcoords;     //!< vertex texcoord array
  ISPCTriangle* triangles;  //!< list of triangles
  int numVertices;
  int numTriangles;
};

struct ISPCScene
{
  ISPCMesh** meshes;         //!< list of meshes
  ISPCMaterial* materials;  //!< material list
  int numMeshes;
  int numMaterials;
};

/* scene data */
extern "C" ISPCScene* g_ispc_scene;
RTCScene g_scene = NULL;

/* render function to use */
renderPixelFunc renderPixel;

/* light */
Vec3f AmbientLight__L;

/* called by the C++ code for initialization */
extern "C" void device_init (int8* cfg)
{
  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set start render mode */
  renderPixel = renderPixelStandard;

  /* set light */
  AmbientLight__L = Vec3f(1,1,1);
}

RTCScene convertScene(ISPCScene* scene_in)
{
  /* create scene */
  RTCScene scene_out = rtcNewScene(RTC_SCENE_STATIC | RTC_SCENE_INCOHERENT,RTC_INTERSECT1);

  /* add all meshes to the scene */
  for (int i=0; i<scene_in->numMeshes; i++)
  {
    /* get ith mesh */
    ISPCMesh* mesh = scene_in->meshes[i];

    /* create a triangle mesh */
    unsigned int geometry = rtcNewTriangleMesh (scene_out, RTC_GEOMETRY_STATIC, mesh->numTriangles, mesh->numVertices);
    
    /* set vertices */
    Vertex* vertices = (Vertex*) rtcMapBuffer(scene_out,geometry,RTC_VERTEX_BUFFER); 
    for (int j=0; j<mesh->numVertices; j++) {
      vertices[j].x = mesh->positions[j].x;
      vertices[j].y = mesh->positions[j].y;
      vertices[j].z = mesh->positions[j].z;
    }

    /* set triangles */
    Triangle* triangles = (Triangle*) rtcMapBuffer(scene_out,geometry,RTC_INDEX_BUFFER);
    for (int j=0; j<mesh->numTriangles; j++) {
      triangles[j].v0 = mesh->triangles[j].v0;
      triangles[j].v1 = mesh->triangles[j].v1;
      triangles[j].v2 = mesh->triangles[j].v2;
    }
    rtcUnmapBuffer(scene_out,geometry,RTC_VERTEX_BUFFER); 
    rtcUnmapBuffer(scene_out,geometry,RTC_INDEX_BUFFER);
  }

  /* commit changes to scene */
  rtcCommit (scene_out);
  return scene_out;
}

/*! Cosine weighted hemisphere sampling. Up direction is the z direction. */
inline Vec3f cosineSampleHemisphere(float& pdf, const float u, const float v) 
{
  const float phi = 2.0f * (float)pi * u;
  const float cosTheta = sqrt(v), sinTheta = sqrt(1.0f - v);
  pdf = cosTheta*(1.0f/(float)pi);
  return Vec3f(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

/*! Cosine weighted hemisphere sampling. Up direction is provided as argument. */
inline Vec3f cosineSampleHemisphere(float& pdf, const float& u, const float& v, const Vec3f& N) {
  return mul(frame(N),cosineSampleHemisphere(pdf,u,v));
}

inline Vec3f AmbientLight__eval(const Vec3f& Ns, const Vec3f& wi) {
  return AmbientLight__L;
}

inline Vec3f AmbientLight__sample(const Vec3f& Ns, 
                                  Vec3f& wi,
                                  float& tMax,
                                  const Vec2f& s) 
{
  float pdf; wi = cosineSampleHemisphere(pdf,s.x,s.y,Ns);
  tMax = 1e20f;
  return div(AmbientLight__L,pdf);
}

inline Vec3f Matte__eval(const int& materialID, const Vec3f& wo, const Vec3f& Ns, const Vec3f& wi) 
{
  ISPCMaterial* material = &g_ispc_scene->materials[materialID];
  Vec3f diffuse = material->Kd;
  return mul(diffuse, (1.0f/(float)pi) * clamp(dot(wi,Ns),0.0f,1.0f));
}

inline Vec3f Matte__sample(const int& materialID, const Vec3f& wo, const Vec3f& Ns, Vec3f& wi, const Vec2f& s) 
{
  float pdf; wi = cosineSampleHemisphere(pdf,s.x,s.y,Ns);
  return div(Matte__eval(materialID, wo, Ns, wi),pdf);
}

inline float frand(int& seed) {
  seed = 1103515245 * seed + 12345;
  seed = 235543534 * seed + 2341233;
  seed = 43565 * seed + 2332443;
  return (seed & 0xFFFF)/(float)0xFFFF;
}

inline Vec3f face_forward(const Vec3fa& dir, const Vec3fa& Ng) {
  return dot(dir,Ng) < 0.0f ? Ng : neg(Ng);
}


Vec3f renderPixelSeed(int x, int y, int& seed, const Vec3f& vx, const Vec3f& vy, const Vec3f& vz, const Vec3f& p)
{
  /* radiance accumulator and weight */
  Vec3f L = Vec3f(0.0f);
  Vec3f Lw = Vec3f(1.0f);

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

  /* iterative path tracer loop */
  for (int i=0; i<10; i++)
  {
    /* terminate if contribution too low */
    if (max(Lw.x,max(Lw.y,Lw.z)) < 0.01f)
      break;

    /* intersect ray with scene */ 
    rtcIntersect(g_scene,ray);
    Vec3f Ns = face_forward(ray.dir,normalize(ray.Ng));
    Vec3f Ph = add(ray.org,mul(ray.tfar,ray.dir));

    /* shade background with ambient light */
    if (ray.geomID == RTC_INVALID_GEOMETRY_ID) {
      Vec3f La = AmbientLight__eval(Ns,neg(ray.dir));
      L = add(L,mul(Lw,La));
      break;
    }
        
    /* shade all rays that hit something */
#if 1 // FIXME: pointer gather not implemented on ISPC for Xeon Phi
    int materialID = g_ispc_scene->meshes[ray.geomID]->triangles[ray.primID].materialID; 
#else
    int materialID = 0;
    foreach_unique (geomID in ray.geomID) {
      if (geomID >= 0 && geomID < g_ispc_scene->numMeshes) { // FIXME: workaround for ISPC bug
        ISPCMesh* mesh = g_ispc_scene->meshes[geomID];
        materialID = mesh->triangles[ray.primID].materialID;
      }
    }
#endif

    /* sample ambient light */
    Vec3f wi; float tMax;
    Vec2f s = Vec2f(frand(seed),frand(seed));
    Vec3f Ll = AmbientLight__sample(Ns,wi,tMax,s);
        
    /* initialize shadow ray */
    RTCRay shadow;
    shadow.org = Ph;
    shadow.dir = wi;
    shadow.tnear = 0.001f;
    shadow.tfar = inf;
    shadow.geomID = RTC_INVALID_GEOMETRY_ID;
    shadow.primID = RTC_INVALID_GEOMETRY_ID;
    shadow.mask = -1;
    shadow.time = 0;
    
    /* trace shadow ray */
    rtcOccluded(g_scene,shadow);
    
    /* add light contribution */
    if (shadow.geomID == RTC_INVALID_GEOMETRY_ID) {
      Vec3f Lm = Matte__eval(materialID,neg(ray.dir),Ns,wi);
      L = add(L,mul(Lw,mul(Ll,Lm)));
    }

    /* calculate diffuce bounce */
    s = Vec2f(frand(seed),frand(seed));
    Vec3f c = Matte__sample(materialID,neg(ray.dir), Ns, wi, s);
    Lw = mul(Lw,c);

    /* setup secondary ray */
    ray.org = Ph;
    ray.dir = normalize(wi);
    ray.tnear = 0.001f;
    ray.tfar = inf;
    ray.geomID = RTC_INVALID_GEOMETRY_ID;
    ray.primID = RTC_INVALID_GEOMETRY_ID;
    ray.mask = -1;
    ray.time = 0;
  }
  return L;
}

/* task that renders a single screen tile */
Vec3fa renderPixelStandard(int x, int y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  int seed = x*233+y*234234+237;
  Vec3f L = Vec3f(0.0f,0.0f,0.0f);
  //for (int i=0; i<16; i++) {
    L = add(L,renderPixelSeed(x,y,seed,vx,vy,vz,p));
  //}
  //L = mul(L,1.0f/16.0f);
  return L;
}
  
/* task that renders a single screen tile */
void renderTile(int taskIndex, int* pixels,
                     const int width,
                     const int height, 
                     const float time,
                     const Vec3f& vx, 
                     const Vec3f& vy, 
                     const Vec3f& vz, 
                     const Vec3f& p,
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
    //Vec3f color = Vec3f(0.0f,0.0f,0.0f);

    Vec3f color = renderPixel(x,y,vx,vy,vz,p);

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
                           const Vec3f& vx, 
                           const Vec3f& vy, 
                           const Vec3f& vz, 
                           const Vec3f& p)
{
  /* create scene */
  if (g_scene == NULL)
    g_scene = convertScene(g_ispc_scene);

  /* render image */
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
