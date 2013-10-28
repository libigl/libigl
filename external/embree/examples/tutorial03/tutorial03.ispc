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

#include "../tutorials/tutorials.isph"

/* scene data */
uniform Scene* uniform scene = NULL;
uniform RTCGeometry* uniform geometry = NULL;
uniform RTCIntersector* uniform intersector = NULL;

/* called by the C++ code for initialization */
export void init (uniform int verbose)
{
  /* initialize ray tracing core */
  rtcInit();
  rtcStartThreads();
  rtcSetVerbose(verbose);
}

/* called by the C++ code to set scene */
export void set_scene (uniform Scene* uniform scene_i) 
{
  /* set scene */
  scene = scene_i;

  /* create a triangle mesh */
  geometry = rtcNewTriangleMesh (scene->numTriangles, scene->numVertices);

  /* set vertices */
  uniform RTCVertex* uniform vertices = rtcMapPositionBuffer(geometry); 
  memcpy(vertices,scene->positions,scene->numVertices*sizeof(uniform RTCVertex));
  rtcUnmapPositionBuffer(geometry); 

  /* set triangles and colors */
  uniform RTCTriangle* uniform triangles = rtcMapTriangleBuffer(geometry);
  memcpy(triangles,scene->triangles,scene->numTriangles*sizeof(uniform RTCTriangle));
  rtcUnmapTriangleBuffer(geometry); 

  /* build spatial index structure and mark as static */
  launch rtcBuildAccel (geometry); sync;
  rtcCleanupGeometry (geometry);

  /* get intersector for the geometry */
  intersector = rtcQueryIntersector (geometry);
}

/* renders a single pixel */
inline vec3f renderPixel (const uniform vec3f& vx, 
                          const uniform vec3f& vy, 
                          const uniform vec3f& vz, 
                          const uniform vec3f& p,
                          const varying int x,
                          const varying int y)
{
    /* initialize ray */
    Ray ray;
    ray.org = p;
    ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
    ray.tnear = 0.0f;
    ray.tfar = inf;
    ray.id0 = -1;
    ray.id1 = -1;
    ray.mask = -1;
    ray.time = 0;

    /* intersect ray with scene */
    intersector->intersect(intersector,ray);

    /* shade background black */
    if (ray.id0 == -1) return make_vec3f(0.0f);

    /* shade all rays that hit something */
    vec3f color = make_vec3f(0.0f);
    foreach_unique(materialID in ray.id1) {
      uniform Material* uniform material = &scene->materials[materialID];
      color = material->Kd;
    }

    /* apply ambient light */
    vec3f Ng = normalize(ray.Ng);
    vec3f Nf = dot(ray.dir,Ng) < 0.0f ? Ng : neg(Ng);
    color = mul(color,mul(scene->ambientLightIntensity,abs(dot(ray.dir,Ng))));

    /* trace shadow ray */
    uniform vec3f pointLightIntensity = scene->pointLightIntensity;
    if (ne(pointLightIntensity,make_vec3f(0.0f)))
    {
      vec3f hitPosition = add(ray.org,mul(ray.tfar,ray.dir));
      vec3f lightVector = sub(scene->pointLightPosition,hitPosition);
      Ray shadow;
      shadow.org = hitPosition;
      shadow.dir = lightVector;
      shadow.tnear = 0.001f;
      shadow.tfar  = 0.999f;
      shadow.id0 = -1;
      shadow.id1 = -1;
      shadow.mask = -1;
      shadow.time = 0;

      /* intersect ray with scene */
      if (!intersector->occluded(intersector,shadow))
      {
        lightVector = normalize(lightVector);
        color = add(color,mul(clamp(dot(Nf,lightVector)),scene->pointLightIntensity));
      }
    }
    return color;
}

/* task that renders a single screen tile */
task void renderTile(uniform int* uniform pixels,
                     const uniform int width,
                     const uniform int height, 
                     const uniform float time,
                     const uniform vec3f& vx, 
                     const uniform vec3f& vy, 
                     const uniform vec3f& vz, 
                     const uniform vec3f& p,
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
    /* render single pixel */
    vec3f color = renderPixel(vx,vy,vz,p,x,y);

    /* write color into framebuffer */
    unsigned int r = (unsigned int) (255.0f * clamp(color.x));
    unsigned int g = (unsigned int) (255.0f * clamp(color.y));
    unsigned int b = (unsigned int) (255.0f * clamp(color.z));
    pixels[y*width+x] = (b << 16) + (g << 8) + r;
  }
}

/* called by the C++ code to render */
export void render (uniform int* uniform pixels,
                    const uniform int width,
                    const uniform int height, 
                    const uniform float time,
                    const uniform vec3f& vx, 
                    const uniform vec3f& vy, 
                    const uniform vec3f& vz, 
                    const uniform vec3f& p)
{
  const uniform int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const uniform int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  launch[numTilesX*numTilesY] renderTile(pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY); sync;
}

/* called by the C++ code for cleanup */
export void cleanup ()
{
  rtcDeleteIntersector (intersector);
  rtcDeleteGeometry (geometry);
  rtcStopThreads();
  rtcExit();
}
