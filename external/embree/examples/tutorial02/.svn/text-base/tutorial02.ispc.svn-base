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

const uniform int numPhi = 20;
const uniform int numTheta = 2*numPhi;

/* scene data */
struct Instance 
{
  uniform RTCGeometry* geometry;
  uniform RTCIntersector* intersector;
  AffineSpace3f local2world;
  vec3f lower;
  vec3f upper;
};
typedef uniform Instance* uniform uniformInstancePtr;
uniformInstancePtr* uniform g_instances = NULL;

/* creates a triangulated sphere */
uniform Instance* uniform createTriangulatedSphere ()
{
  /* create instance */
  uniform Instance* uniform instance = uniform new uniform Instance;
 
  /* create triangle mesh */
  uniform RTCGeometry* uniform mesh = rtcNewTriangleMesh (2*numTheta*(numPhi-1), numTheta*(numPhi+1));

  /* set triangles and vertices */
  uniform RTCVertex* uniform vertices = rtcMapPositionBuffer(mesh); 
  uniform RTCTriangle* uniform triangles = rtcMapTriangleBuffer(mesh);
  const uniform float rcpNumTheta = rcp(numTheta);
  const uniform float rcpNumPhi   = rcp(numPhi);
  
  uniform int tri = 0;
  for (uniform int phi=0; phi<=numPhi; phi++)
  {
    for (uniform int theta=0; theta<numTheta; theta++)
    {
      const uniform float phif   = phi*pi*rcpNumPhi;
      const uniform float thetaf = theta*2.0f*pi*rcpNumTheta;

      uniform RTCVertex& v = vertices[phi*numTheta+theta];
      v.x = sin(phif)*sin(thetaf);
      v.y = cos(phif);
      v.z = sin(phif)*cos(thetaf);
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
        triangles[tri].id0 = 0; 
        triangles[tri].id1 = tri++;
      }

      if (phi < numPhi) {
        triangles[tri].v0 = p11; 
        triangles[tri].v1 = p10;
        triangles[tri].v2 = p01; 
        triangles[tri].id0 = 0; 
        triangles[tri].id1 = tri++;
      }
    }
  }
  rtcUnmapPositionBuffer(mesh); 
  rtcUnmapTriangleBuffer(mesh); 
  launch rtcBuildAccel(mesh); sync;
  rtcCleanupGeometry(mesh);

  instance->geometry = mesh;
  instance->intersector = rtcQueryIntersector(mesh);
  rtcGetBounds (mesh, &instance->lower.x, &instance->upper.x);
  instance->local2world.l.vx = make_vec3f(1,0,0);
  instance->local2world.l.vy = make_vec3f(0,1,0);
  instance->local2world.l.vz = make_vec3f(0,0,1);
  instance->local2world.p    = make_vec3f(-1.5f,0,0);
  
  return instance;
}

/* Sphere structure */
struct Sphere
{
  RTCIntersector intersector;
  vec3f p;   //!< position of the sphere
  float r;   //!< radius of the sphere
};

/* Intersects the ray with the sphere. */
void intersectSphereFunc(const uniform Sphere* uniform sphere, varying Ray& ray)
{
  const vec3f v = sub(ray.org,sphere->p);
  const float A = dot(ray.dir,ray.dir);
  const float B = 2.0f*dot(v,ray.dir);
  const float C = dot(v,v) - sqr(sphere->r);
  const float D = B*B - 4.0f*A*C;
  if (D < 0.0f) return;
  const float Q = sqrt(D);
  const float rcpA = rcp(A);
  const float t0 = 0.5f*rcpA*(-B-Q);
  const float t1 = 0.5f*rcpA*(-B+Q);
  if (ray.tnear < t0 & t0 < ray.tfar) {
    ray.u = 0.0f;
    ray.v = 0.0f;
    ray.tfar = t0;
    ray.id0 = 0;
    ray.id1 = 0;
    ray.Ng = sub(add(ray.org,mul(t0,ray.dir)),sphere->p);
  }
  if (ray.tnear < t1 & t1 < ray.tfar) {
    ray.u = 0.0f;
    ray.v = 0.0f;
    ray.tfar = t1;
    ray.id0 = 0;
    ray.id1 = 0;
    ray.Ng = sub(add(ray.org,mul(t1,ray.dir)),sphere->p);
  }
}

/* Tests the ray for occlusion with the sphere. */
varying bool occludedSphereFunc(const uniform Sphere* uniform sphere, varying Ray& ray)
{
  const vec3f v = sub(ray.org,sphere->p);
  const float A = dot(ray.dir,ray.dir);
  const float B = 2.0f*dot(v,ray.dir);
  const float C = dot(v,v) - sqr(sphere->r);
  const float discr = B*B - 4.0f*A*C;
  return discr >= 0.0f;
}

/* creates an analytical sphere */
uniform Instance* uniform createAnalyticalSphere ()
{
  uniform Sphere* uniform sphere = uniform new uniform Sphere;
  sphere->intersector.intersect = (intersectFunc) &intersectSphereFunc;
  sphere->intersector.occluded  = (occludedFunc ) &occludedSphereFunc;
  sphere->p = make_vec3f(0,0,0);
  sphere->r = 1.0f;
  
  uniform Instance* uniform instance = uniform new uniform Instance;
  instance->geometry = NULL;
  instance->intersector = (uniform RTCIntersector* uniform) sphere;
  instance->lower = sub(sphere->p,make_vec3f(sphere->r));
  instance->upper = add(sphere->p,make_vec3f(sphere->r));
  instance->local2world.l.vx = make_vec3f(1,0,0);
  instance->local2world.l.vy = make_vec3f(0,1,0);
  instance->local2world.l.vz = make_vec3f(0,0,1);
  instance->local2world.p    = make_vec3f(1.5f,0,0);
  return instance;
}

/* create top level scene */
uniform RTCGeometry* uniform createScene (uniform Instance** uniform instances, uniform int numInstances)
{
  uniform RTCGeometry* uniform scene = rtcNewVirtualGeometry (numInstances);
  for (uniform int i=0; i<numInstances; i++) 
  {
    uniform int mask = i+1;
    uniform Instance* uniform instance = instances[i];
    rtcSetVirtualGeometryUserData    (scene, i, 0, 0, mask);
    rtcSetVirtualGeometryBounds      (scene, i, &instance->lower.x, &instance->upper.x, (uniform RTCTransformation* uniform) &instance->local2world);
    rtcSetVirtualGeometryIntersector (scene, i, instance->intersector);
  }
  launch rtcBuildAccel(scene); sync;
  return scene;
}

/* called by the C++ code for initialization */
export void init (uniform int verbose)
{
  /* initialize ray tracing core */
  rtcInit();
  rtcStartThreads();
  rtcSetVerbose(verbose);

  /* create scene */
  g_instances = uniform new uniformInstancePtr [2];
  g_instances[0] = createTriangulatedSphere();
  g_instances[1] = createAnalyticalSphere();
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
                     uniform RTCIntersector* uniform intersector,
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
    //ray.mask = 2; // masks out first sphere
    //ray.mask = 1; // masks out second sphere
    ray.time = 0;

    /* intersect ray with scene */
    intersector->intersect(intersector,ray);

    /* shade pixels */
    if (ray.id0 != -1) {
      vec3f c = make_vec3f(abs(dot(normalize(ray.Ng),ray.dir)));
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
  /* move instances */
  uniform float t = 0.7f*time;
  g_instances[0]->local2world.p = mul(1.5f,make_vec3f(+cos(t),0.0f,+sin(t)));
  g_instances[1]->local2world.p = mul(1.5f,make_vec3f(-cos(t),0.0f,-sin(t)));

  /* create scene */
  uniform RTCGeometry* uniform scene = createScene(g_instances,2);
  uniform RTCIntersector* uniform intersector = rtcQueryIntersector(scene);

  /* render all pixels */
  const uniform int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const uniform int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  launch[numTilesX*numTilesY] renderTile(pixels,width,height,time,vx,vy,vz,p,intersector,numTilesX,numTilesY); sync;

  /* cleanup */
  rtcDeleteIntersector(intersector);
  rtcDeleteGeometry(scene);
}

/* called by the C++ code for cleanup */
export void cleanup ()
{
  delete[] g_instances;
  rtcStopThreads();
  rtcExit();
}
