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
uniform RTCGeometry* uniform mesh = NULL;
uniform vec3f* uniform colors = NULL;
uniform RTCIntersector* uniform intersector = NULL;

/* called by the C++ code for initialization */
export void init (uniform int verbose)
{
  /* initialize ray tracing core */
  rtcInit();
  rtcStartThreads();
  rtcSetVerbose(verbose);

  /* create a triangulated cube with 12 triangles and 8 vertices */
  mesh = rtcNewTriangleMesh (12, 8);

  /* set vertices */
  uniform RTCVertex* uniform vertices = rtcMapPositionBuffer(mesh); 
  vertices[0].x = -1; vertices[0].y = -1; vertices[0].z = -1; 
  vertices[1].x = -1; vertices[1].y = -1; vertices[1].z = +1; 
  vertices[2].x = -1; vertices[2].y = +1; vertices[2].z = -1; 
  vertices[3].x = -1; vertices[3].y = +1; vertices[3].z = +1; 
  vertices[4].x = +1; vertices[4].y = -1; vertices[4].z = -1; 
  vertices[5].x = +1; vertices[5].y = -1; vertices[5].z = +1; 
  vertices[6].x = +1; vertices[6].y = +1; vertices[6].z = -1; 
  vertices[7].x = +1; vertices[7].y = +1; vertices[7].z = +1; 
  rtcUnmapPositionBuffer(mesh); 

  /* create triangle color array */
  colors = uniform new uniform vec3f[12];

  /* set triangles and colors */
  uniform int tri = 0;
  uniform RTCTriangle* uniform triangles = rtcMapTriangleBuffer(mesh);
  
  // left side
  colors[tri] = make_vec3f(1,0,0); triangles[tri].v0 = 0; triangles[tri].v1 = 2; triangles[tri].v2 = 1; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;
  colors[tri] = make_vec3f(1,0,0); triangles[tri].v0 = 1; triangles[tri].v1 = 2; triangles[tri].v2 = 3; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;

  // right side
  colors[tri] = make_vec3f(0,1,0); triangles[tri].v0 = 4; triangles[tri].v1 = 5; triangles[tri].v2 = 6; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;
  colors[tri] = make_vec3f(0,1,0); triangles[tri].v0 = 5; triangles[tri].v1 = 7; triangles[tri].v2 = 6; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;

  // bottom side
  colors[tri] = make_vec3f(0.5f);  triangles[tri].v0 = 0; triangles[tri].v1 = 1; triangles[tri].v2 = 4; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;
  colors[tri] = make_vec3f(0.5f);  triangles[tri].v0 = 1; triangles[tri].v1 = 5; triangles[tri].v2 = 4; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;

  // top side
  colors[tri] = make_vec3f(1.0f);  triangles[tri].v0 = 2; triangles[tri].v1 = 6; triangles[tri].v2 = 3; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;
  colors[tri] = make_vec3f(1.0f);  triangles[tri].v0 = 3; triangles[tri].v1 = 6; triangles[tri].v2 = 7; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;

  // front side
  colors[tri] = make_vec3f(0,0,1); triangles[tri].v0 = 0; triangles[tri].v1 = 4; triangles[tri].v2 = 2; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;
  colors[tri] = make_vec3f(0,0,1); triangles[tri].v0 = 2; triangles[tri].v1 = 4; triangles[tri].v2 = 6; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;

  // back side
  colors[tri] = make_vec3f(1,1,0); triangles[tri].v0 = 1; triangles[tri].v1 = 3; triangles[tri].v2 = 5; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;
  colors[tri] = make_vec3f(1,1,0); triangles[tri].v0 = 3; triangles[tri].v1 = 7; triangles[tri].v2 = 5; triangles[tri].id0 = 0; triangles[tri].id1 = tri++;

  rtcUnmapTriangleBuffer(mesh); 

  /* build spatial index structure and mark as static */
  launch rtcBuildAccel (mesh); sync;
  rtcCleanupGeometry (mesh);
  
  /* get intersector for the mesh */
  intersector = rtcQueryIntersector (mesh);
}

/* called by the C++ code to set scene */
export void set_scene (uniform Scene* uniform scene) {
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

    /* shade pixels */
    if (ray.id0 != -1) {
      vec3f c = colors[ray.id1];
      unsigned int r = (unsigned int) (255.0f * c.x);
      unsigned int g = (unsigned int) (255.0f * c.y);
      unsigned int b = (unsigned int) (255.0f * c.z);
      pixels[y*width+x] = (b << 16) + (g << 8) + r;
    }
    else pixels[y*width+x] = 0;
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
  rtcDeleteGeometry (mesh);
  delete[] colors;
  rtcStopThreads();
  rtcExit();
}
