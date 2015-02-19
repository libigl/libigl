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
#include "../common/tutorial/scene_device.isph"

#define PARALLEL_COMMIT

//#define FORCE_FIXED_EDGE_TESSELLATION
#define FIXED_EDGE_TESSELLATION_VALUE 4

#define MAX_EDGE_LEVEL 64.0f
#define MIN_EDGE_LEVEL 2.0f
#define ENABLE_DISPLACEMENTS 0
#if ENABLE_DISPLACEMENTS
#  define LEVEL_FACTOR 256.0f
#else
#  define LEVEL_FACTOR 64.0f
#endif

/* scene data */
extern uniform ISPCScene* uniform g_ispc_scene;
RTCScene g_scene = NULL;
RTCScene g_embree_scene = NULL;
RTCScene g_osd_scene = NULL;

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

Vec3fa renderPixelEyeLight(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p);


/* called by the C++ code for initialization */
export void device_init (uniform int8* uniform cfg)
{
  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set error handler */
  rtcSetErrorFunction(error_handler);

  /* set start render mode */
  renderPixel = renderPixelStandard;
  //  renderPixel = renderPixelUV;	
}

void updateEdgeLevelBuffer( uniform ISPCMesh* uniform mesh, const uniform Vec3fa& cam_pos, uniform size_t startID, uniform size_t endID )
{
  for (uniform size_t f=startID; f<endID; f++) 
  {
    const uniform int N = 4;
    for (uniform size_t i=0; i<N; i++) {
      const uniform int * uniform index = (uniform int * uniform)&mesh->quads[f];	
      const uniform Vec3fa v0 = mesh->positions[index[i]];
      const uniform Vec3fa v1 = mesh->positions[index[(i+1)%N]];
      const uniform Vec3fa edge = v1-v0;
      const uniform Vec3fa P = 0.5f*(v1+v0);
      const uniform Vec3fa dist = cam_pos - P;
      mesh->edge_level[f*4+i] = max(min(LEVEL_FACTOR*(0.5f*length(edge)/length(dist)),MAX_EDGE_LEVEL),MIN_EDGE_LEVEL); // FIXME: map mesh subdivlevel buffer
      assert( mesh->edge_level[f*4+i] >= 1.0f );
    }
  }
}

#if defined(ISPC)
task void updateEdgeLevelBufferTask( uniform ISPCMesh* uniform mesh, const uniform Vec3fa& cam_pos )
{
  const uniform size_t size = mesh->numQuads;
  const uniform size_t startID = ((taskIndex+0)*size)/taskCount;
  const uniform size_t endID   = ((taskIndex+1)*size)/taskCount;
  updateEdgeLevelBuffer(mesh,cam_pos,startID,endID);
}
#endif

void updateEdgeLevels(uniform ISPCScene* uniform scene_in, const uniform Vec3fa& cam_pos)
{
  for (uniform int i=0; i<g_ispc_scene->numMeshes; i++)
  {
    uniform ISPCMesh* uniform mesh = g_ispc_scene->meshes[i];
    if (mesh->numQuads == 0) continue;

    if (mesh->edge_level) 
    {
      uniform unsigned int geomID = mesh->geomID;
#if defined(ISPC)
      launch[  getNumHWThreads() ] updateEdgeLevelBufferTask(mesh,cam_pos); sync;	           
#else
      updateEdgeLevelBuffer(mesh,cam_pos,0,mesh->numQuads);
#endif
      rtcUpdateBuffer(g_scene,geomID,RTC_LEVEL_BUFFER);
    }
  }

  for (uniform size_t g=0; g<scene_in->numSubdivMeshes; g++)
  {
    uniform ISPCSubdivMesh* uniform mesh = g_ispc_scene->subdiv[g];
    uniform unsigned int geomID = mesh->geomID;
    for (uniform size_t f=0, e=0; f<mesh->numFaces; e+=mesh->verticesPerFace[f++]) {
      uniform int N = mesh->verticesPerFace[f];
      for (uniform size_t i=0; i<N; i++) {
        const uniform Vec3fa v0 = mesh->positions[mesh->position_indices[e+(i+0)]];
        const uniform Vec3fa v1 = mesh->positions[mesh->position_indices[e+(i+1)%N]];
        const uniform Vec3fa edge = v1-v0;
        const uniform Vec3fa P = 0.5f*(v1+v0);
	const uniform Vec3fa dist = cam_pos - P;
        mesh->subdivlevel[e+i] = max(min(LEVEL_FACTOR*(0.5f*length(edge)/length(dist)),MAX_EDGE_LEVEL),MIN_EDGE_LEVEL);
      }
    }
  }
}

void displacementFunction(void* uniform ptr, uniform unsigned int geomID, uniform int unsigned primID, 
                      const uniform float* uniform u,      /*!< u coordinates (source) */
                      const uniform float* uniform v,      /*!< v coordinates (source) */
                      const uniform float* uniform nx,     /*!< x coordinates of normal at point to displace (source) */
                      const uniform float* uniform ny,     /*!< y coordinates of normal at point to displace (source) */
                      const uniform float* uniform nz,     /*!< z coordinates of normal at point to displace (source) */
                      uniform float* uniform px,           /*!< x coordinates of points to displace (source and target) */
                      uniform float* uniform py,           /*!< y coordinates of points to displace (source and target) */
                      uniform float* uniform pz,           /*!< z coordinates of points to displace (source and target) */
                      uniform size_t N)
{
#if 0
  foreach (i=0 ... N) {
    const Vec3fa dP = 0.02f*make_Vec3fa(sin(100.0f*px[i]+0.5f),sin(100.0f*pz[i]+1.5f),cos(100.0f*py[i]));
    px[i] += dP.x; py[i] += dP.y; pz[i] += dP.z;
  }
#else
  foreach (i=0 ... N) {
    const Vec3fa P = make_Vec3fa(px[i],py[i],pz[i]);
    const Vec3fa nor = make_Vec3fa(nx[i],ny[i],nz[i]);
    float dN = 0.0f;
    for (float freq = 1.0f; freq<40.0f; freq*= 2) {
      float n = abs(noise(freq*P));
      dN += 1.4f*n*n/freq;
    }
    const Vec3fa dP = dN*nor;
    px[i] += dP.x; py[i] += dP.y; pz[i] += dP.z;
  }
#endif
}

void convertScene(uniform ISPCScene* uniform scene_in, const uniform Vec3fa& p)
{
  /* add all meshes to the scene */
  for (uniform int i=0; i<scene_in->numMeshes; i++)
  {
    /* get ith mesh */
    uniform ISPCMesh* uniform mesh = scene_in->meshes[i];

    if (mesh->numQuads)
    {   
      mesh->edge_level  = uniform new uniform float[mesh->numQuads*4];
      for (uniform size_t i=0; i<4*mesh->numQuads; i++) 
	mesh->edge_level[i] = FIXED_EDGE_TESSELLATION_VALUE;
      
      /* create a triangle mesh */
      uniform unsigned int geomID = rtcNewSubdivisionMesh (g_scene, RTC_GEOMETRY_STATIC, mesh->numQuads, mesh->numQuads*4, mesh->numVertices, 0, 0, 0);
      mesh->geomID = geomID;

      uniform unsigned int* faces = (uniform unsigned int*) rtcMapBuffer(g_scene, geomID, RTC_FACE_BUFFER);
      for (uniform size_t i=0; i<mesh->numQuads; i++) faces[i] = 4;
      rtcUnmapBuffer(g_scene,geomID,RTC_FACE_BUFFER);
      
      rtcSetBuffer(g_scene, geomID, RTC_VERTEX_BUFFER, mesh->positions , 0, sizeof(uniform Vec3fa  ));
      rtcSetBuffer(g_scene, geomID, RTC_INDEX_BUFFER,  mesh->quads     , 0, sizeof(uniform unsigned int));
      rtcSetBuffer(g_scene, geomID, RTC_LEVEL_BUFFER,  mesh->edge_level, 0, sizeof(uniform float));
      
#if ENABLE_DISPLACEMENTS == 1
      rtcSetDisplacementFunction(g_scene,geomID,(RTCDisplacementFunc)&displacementFunction,NULL);
#endif
     }
     else
      mesh->edge_level = NULL;
  }
  
  /* add all subdiv meshes to the scene */
  for (uniform size_t i=0; i<scene_in->numSubdivMeshes; i++)
  {
    uniform ISPCSubdivMesh* uniform mesh = scene_in->subdiv[i];
    uniform unsigned int geomID = rtcNewSubdivisionMesh(g_scene, RTC_GEOMETRY_STATIC, mesh->numFaces, mesh->numEdges, mesh->numVertices, 
							mesh->numEdgeCreases, mesh->numVertexCreases, mesh->numHoles);
    mesh->geomID = geomID;

    rtcSetBuffer(g_scene, geomID, RTC_VERTEX_BUFFER, mesh->positions, 0, sizeof(uniform Vec3fa  ));
    rtcSetBuffer(g_scene, geomID, RTC_LEVEL_BUFFER,  mesh->subdivlevel, 0, sizeof(uniform float));
    rtcSetBuffer(g_scene, geomID, RTC_INDEX_BUFFER,  mesh->position_indices  , 0, sizeof(uniform unsigned int));
    rtcSetBuffer(g_scene, geomID, RTC_FACE_BUFFER,   mesh->verticesPerFace, 0, sizeof(uniform unsigned int));
    rtcSetBuffer(g_scene, geomID, RTC_HOLE_BUFFER,   mesh->holes, 0, sizeof(uniform unsigned int));
    rtcSetBuffer(g_scene, geomID, RTC_EDGE_CREASE_INDEX_BUFFER,    mesh->edge_creases,          0, 2*sizeof(uniform unsigned int));
    rtcSetBuffer(g_scene, geomID, RTC_EDGE_CREASE_WEIGHT_BUFFER,   mesh->edge_crease_weights,   0, sizeof(uniform float));
    rtcSetBuffer(g_scene, geomID, RTC_VERTEX_CREASE_INDEX_BUFFER,  mesh->vertex_creases,        0, sizeof(uniform unsigned int));
    rtcSetBuffer(g_scene, geomID, RTC_VERTEX_CREASE_WEIGHT_BUFFER, mesh->vertex_crease_weights, 0, sizeof(uniform float));

#if ENABLE_DISPLACEMENTS == 1
    //BBox3fa bounds(Vec3fa(-0.1f,-0.1f,-0.1f),Vec3fa(0.1f,0.1f,0.1f));
    rtcSetDisplacementFunction(g_scene, geomID, (RTCDisplacementFunc)&displacementFunction,NULL);
#endif
  }
}

#if defined(__USE_OPENSUBDIV__)

#include <opensubdiv/far/topologyRefinerFactory.h>
using namespace OpenSubdiv;

struct OSDVertex {

    // Minimal required interface ----------------------
    OSDVertex() { }

    OSDVertex(OSDVertex const & src) {
        _position[0] = src._position[0];
        _position[1] = src._position[1];
        _position[1] = src._position[1];
    }

    void Clear( void * =0 ) {
        _position[0]=_position[1]=_position[2]=0.0f;
    }

    void AddWithWeight(OSDVertex const & src, float weight) {
        _position[0]+=weight*src._position[0];
        _position[1]+=weight*src._position[1];
        _position[2]+=weight*src._position[2];
    }

    void AddVaryingWithWeight(OSDVertex const &, float) { }

    // Public interface ------------------------------------
    void SetPosition(float x, float y, float z) {
        _position[0]=x;
        _position[1]=y;
        _position[2]=z;
    }

    const float * GetPosition() const {
        return _position;
    }

private:
    float _position[3];
};

RTCScene constructSceneOpenSubdiv() 
{
  if (!g_ispc_scene) return NULL;

  typedef Far::TopologyRefinerFactoryBase::TopologyDescriptor Descriptor;

  Sdc::Options options;
  options.SetVVarBoundaryInterpolation(Sdc::Options::VVAR_BOUNDARY_EDGE_ONLY);
  options.SetCreasingMethod(Sdc::Options::CREASE_CHAIKIN);

  RTCScene scene = rtcNewScene(RTC_SCENE_STATIC, RTC_INTERSECT1);
  
  for (size_t i=0; i<g_ispc_scene->numSubdivMeshes; i++)
  {
    ISPCSubdivMesh* mesh = g_ispc_scene->subdiv[i];
    
    Descriptor desc;
    desc.numVertices  = mesh->numVertices;
    desc.numFaces     = mesh->numFaces;
    desc.vertsPerFace = mesh->verticesPerFace;
    desc.vertIndices  = mesh->position_indices;
    desc.numCreases   = mesh->numEdgeCreases;
    desc.creaseVertexIndexPairs = (int*) mesh->edge_creases;
    desc.creaseWeights = mesh->edge_crease_weights;
    desc.numCorners    = mesh->numVertexCreases;
    desc.cornerVertexIndices = mesh->vertex_creases;
    desc.cornerWeights = mesh->vertex_crease_weights;
    
    size_t maxlevel = 5;
    Far::TopologyRefiner* refiner = Far::TopologyRefinerFactory<Descriptor>::Create(OpenSubdiv::Sdc::TYPE_CATMARK, options, desc);
    refiner->RefineUniform(maxlevel);
    
    std::vector<OSDVertex> vbuffer(refiner->GetNumVerticesTotal());
    OSDVertex* verts = &vbuffer[0];
    
    for (int i=0; i<mesh->numVertices; ++i)
      verts[i].SetPosition(mesh->positions[i].x,mesh->positions[i].y,mesh->positions[i].z);
    
    refiner->Interpolate(verts, verts + mesh->numVertices);
    
    for (int level=0; level<maxlevel; ++level)
        verts += refiner->GetNumVertices(level);
    
    const size_t numVertices = refiner->GetNumVertices(maxlevel);
    const size_t numFaces    = refiner->GetNumFaces(maxlevel);
    
    unsigned int meshID = rtcNewTriangleMesh(scene, RTC_GEOMETRY_STATIC, 2*numFaces, numVertices);
    rtcSetBuffer(scene, meshID, RTC_VERTEX_BUFFER, verts, 0, sizeof(Vec3f));
    
    Vec3i* tris = (Vec3i*) rtcMapBuffer(scene, meshID, RTC_INDEX_BUFFER);
    for (size_t i=0; i<numFaces; i++) {
      Far::IndexArray fverts = refiner->GetFaceVertices(maxlevel, i);
      assert(fverts.size() == 4);
      tris[2*i+0] = Vec3i(fverts[0],fverts[1],fverts[2]);
      tris[2*i+1] = Vec3i(fverts[2],fverts[3],fverts[0]);
    }
    rtcUnmapBuffer(scene,meshID,RTC_INDEX_BUFFER);
  }
  rtcCommit(scene);
  return scene;
}

#endif

Vec3f rndColor(const int ID) 
{
  int r = ((ID+13)*17*23) & 255;
  int g = ((ID+15)*11*13) & 255;
  int b = ((ID+17)* 7*19) & 255;
  const float oneOver255f = 1.f/255.f;
  return make_Vec3f(r*oneOver255f,g*oneOver255f,b*oneOver255f);
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
  
  /* shade background black */
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return make_Vec3f(0.0f,0.0f,1.0f);
  
  /* shade all rays that hit something */
  Vec3f color = make_Vec3f(1.0f);
#if 0
    color = rndColor(ray.geomID);
#else
  /* apply ambient light */
  Vec3f Ng = normalize(ray.Ng);
  color = color*abs(dot(ray.dir,Ng));   
#endif
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

#if defined(PARALLEL_COMMIT)
task void parallelCommit(RTCScene scene) {
  rtcCommitThread (scene,threadIndex,threadCount); 
}
#endif

export void toggleOpenSubdiv(uniform int key, uniform int x, uniform int y)
{
#if defined(__USE_OPENSUBDIV__)
  if (g_osd_scene == NULL) {
    g_osd_scene = constructSceneOpenSubdiv();
    g_embree_scene = g_scene;
  }
  if (g_scene == g_embree_scene) g_scene = g_osd_scene;
  else                           g_scene = g_embree_scene;
#endif
}

uniform Vec3f old_p; 
 
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
   uniform Vec3fa cam_org = make_Vec3fa(p.x,p.y,p.z);

   /* create scene */
   if (g_scene == NULL)
  { 
    g_scene = rtcNewScene(RTC_SCENE_DYNAMIC,RTC_INTERSECT_UNIFORM | RTC_INTERSECT_VARYING);
    convertScene(g_ispc_scene,cam_org);

#if !defined(FORCE_FIXED_EDGE_TESSELLATION)
    updateEdgeLevels(g_ispc_scene, cam_org);
#endif

    old_p = p;

#if !defined(PARALLEL_COMMIT)
  rtcCommit (g_scene);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene); sync;
#endif
  }

#if !defined(FORCE_FIXED_EDGE_TESSELLATION)
  {
    if ((p.x != old_p.x | p.y != old_p.y | p.z != old_p.z))
    {
     old_p = p;

   updateEdgeLevels(g_ispc_scene, cam_org);

#if !defined(PARALLEL_COMMIT)
     rtcCommit (g_scene);
#else
     launch[ getNumHWThreads() ] parallelCommit(g_scene); sync;
#endif
    }
   }
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
