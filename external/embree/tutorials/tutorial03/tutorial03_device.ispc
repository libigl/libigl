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

#include "../common/tutorial/tutorial_device.isph"

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
  uniform Vec3fa* positions;    //!< vertex position array
  uniform Vec3fa* normals;       //!< vertex normal array
  uniform Vec2f* texcoords;     //!< vertex texcoord array
  uniform ISPCTriangle* triangles;  //!< list of triangles
  int numVertices;
  int numTriangles;
};

struct ISPCScene
{
  uniform ISPCMesh* uniform* meshes;         //!< list of meshes
  uniform ISPCMaterial* materials;  //!< material list
  int numMeshes;
  int numMaterials;
};

/* scene data */
extern uniform ISPCScene* uniform g_ispc_scene;
RTCScene g_scene = NULL;

/* render function to use */
renderPixelFunc renderPixel;

/* called by the C++ code for initialization */
export void device_init (uniform int8* uniform cfg)
{
  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set start render mode */
  renderPixel = renderPixelStandard;
}

RTCScene convertScene(uniform ISPCScene* uniform scene_in)
{
  /* create scene */
  RTCScene scene_out = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT_UNIFORM | RTC_INTERSECT_VARYING);

  /* add all meshes to the scene */
  for (uniform int i=0; i<scene_in->numMeshes; i++)
  {
    /* get ith mesh */
    uniform ISPCMesh* uniform mesh = scene_in->meshes[i];

    /* create a triangle mesh */
    uniform unsigned int geometry = rtcNewTriangleMesh (scene_out, RTC_GEOMETRY_STATIC, mesh->numTriangles, mesh->numVertices);

#if 1
    /* share vertex buffer */
    rtcSetBuffer(scene_out, geometry, RTC_VERTEX_BUFFER, mesh->positions, 0, sizeof(uniform Vec3fa      ));
    //rtcSetBuffer(scene_out, geometry, RTC_INDEX_BUFFER,  mesh->triangles, 0, sizeof(uniform ISPCTriangle));

    /* set triangles */
    uniform Triangle* uniform triangles = (uniform Triangle* uniform) rtcMapBuffer(scene_out,geometry,RTC_INDEX_BUFFER);
    for (uniform int j=0; j<mesh->numTriangles; j++) {
      triangles[j].v0 = mesh->triangles[j].v0;
      triangles[j].v1 = mesh->triangles[j].v1;
      triangles[j].v2 = mesh->triangles[j].v2;
    }
    rtcUnmapBuffer(scene_out,geometry,RTC_INDEX_BUFFER);

#elif 0
    uniform Vec3f* uniform positions = uniform new uniform Vec3f[mesh->numVertices+1];
    for (uniform int j=0; j<mesh->numVertices; j++) {
      positions[j].x = mesh->positions[j].x;
      positions[j].y = mesh->positions[j].y;
      positions[j].z = mesh->positions[j].z;
    }
    rtcSetBuffer(scene_out, geometry, RTC_VERTEX_BUFFER, positions, 0, sizeof(uniform Vec3f));
    rtcSetBuffer(scene_out, geometry, RTC_INDEX_BUFFER,  mesh->triangles, 0, sizeof(uniform ISPCTriangle));
#else

    /* set vertices */
    uniform Vertex* uniform vertices = (uniform Vertex* uniform) rtcMapBuffer(scene_out,geometry,RTC_VERTEX_BUFFER); 
    for (uniform int j=0; j<mesh->numVertices; j++) {
      vertices[j].x = mesh->positions[j].x;
      vertices[j].y = mesh->positions[j].y;
      vertices[j].z = mesh->positions[j].z;
    }
    rtcUnmapBuffer(scene_out,geometry,RTC_VERTEX_BUFFER); 

    /* set triangles */
    uniform Triangle* uniform triangles = (uniform Triangle* uniform) rtcMapBuffer(scene_out,geometry,RTC_INDEX_BUFFER);
    for (uniform int j=0; j<mesh->numTriangles; j++) {
      triangles[j].v0 = mesh->triangles[j].v0;
      triangles[j].v1 = mesh->triangles[j].v1;
      triangles[j].v2 = mesh->triangles[j].v2;
    }
    rtcUnmapBuffer(scene_out,geometry,RTC_INDEX_BUFFER);
#endif
  }

  /* commit changes to scene */
  rtcCommit (scene_out);
  return scene_out;
}

/* task that renders a single screen tile */
Vec3f renderPixelStandard(int x, int y, const uniform Vec3f& vx, const uniform Vec3f& vy, const uniform Vec3f& vz, const uniform Vec3f& p)
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
  
  /* shade background black */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return make_Vec3f(0.0f);
  
  /* shade all rays that hit something */
  Vec3f color = make_Vec3f(0.0f);
#if 0 // FIXME: pointer gather not implemented on ISPC for Xeon Phi
  uniform ISPCMesh* varying mesh = g_ispc_scene->meshes[ray.geomID];
  uniform ISPCTriangle* varying tri = &mesh->triangles[ray.primID];

  /* load material ID */
  int materialID = tri->materialID;

  /* interpolate shading normal */
  Vec3f n0 = make_Vec3f(mesh->normals[tri->v0]);
  Vec3f n1 = make_Vec3f(mesh->normals[tri->v1]);
  Vec3f n2 = make_Vec3f(mesh->normals[tri->v2]);
  float u = ray.u, v = ray.v, w = 1.0f-ray.u-ray.v;
  Vec3f Ns = add(add(mul(w,n0),mul(u,n1)),mul(v,n2));
  Ns = normalize(Ns);

#else

  int materialID = 0;
  Vec3f Ns = make_Vec3f(0.0f);
  foreach_unique (geomID in ray.geomID) 
  {
    uniform ISPCMesh* uniform mesh = g_ispc_scene->meshes[geomID];
    
    foreach_unique (primID in ray.primID) 
    {
      uniform ISPCTriangle* uniform tri = &mesh->triangles[primID];
      
      /* load material ID */
      materialID = tri->materialID;

      /* interpolate shading normal */
      Vec3f n0 = make_Vec3f(mesh->normals[tri->v0]);
      Vec3f n1 = make_Vec3f(mesh->normals[tri->v1]);
      Vec3f n2 = make_Vec3f(mesh->normals[tri->v2]);
      float u = ray.u, v = ray.v, w = 1.0f-ray.u-ray.v;
      Ns = add(add(mul(w,n0),mul(u,n1)),mul(v,n2));
    }
  }
  Ns = normalize(Ns);
#endif
  uniform ISPCMaterial* material = &g_ispc_scene->materials[materialID];
  color = material->Kd;

  /* apply ambient light */
  Vec3f Nf = faceforward(Ns,neg(ray.dir),ray.Ng);
  //Vec3f Ng = normalize(ray.Ng);
  //Vec3f Nf = dot(ray.dir,Ng) < 0.0f ? Ng : neg(Ng);
  color = mul(color,dot(ray.dir,Nf));    
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
  /* create scene */
  if (g_scene == NULL)
    g_scene = convertScene(g_ispc_scene);

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
