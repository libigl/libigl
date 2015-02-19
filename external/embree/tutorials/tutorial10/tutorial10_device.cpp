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

#include "../common/tutorial/tutorial_device.h"
#include "../common/tutorial/scene_device.h"



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
extern "C" ISPCScene* g_ispc_scene;
RTCScene g_scene = NULL;
RTCScene g_embree_scene = NULL;
RTCScene g_osd_scene = NULL;

/* render function to use */
renderPixelFunc renderPixel;

/* error reporting function */
void error_handler(const RTCError code, const int8* str)
{
  printf("Embree: ");
  switch (code) {
  case RTC_UNKNOWN_ERROR    : printf("RTC_UNKNOWN_ERROR"); break;
  case RTC_INVALID_ARGUMENT : printf("RTC_INVALID_ARGUMENT"); break;
  case RTC_INVALID_OPERATION: printf("RTC_INVALID_OPERATION"); break;
  case RTC_OUT_OF_MEMORY    : printf("RTC_OUT_OF_MEMORY"); break;
  case RTC_UNSUPPORTED_CPU  : printf("RTC_UNSUPPORTED_CPU"); break;
  default                   : printf("invalid error code"); break;
  }
  if (str) { 
    printf(" ("); 
    while (*str) putchar(*str++); 
    printf(")\n"); 
  }
  abort();
}

Vec3fa renderPixelEyeLight(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p);


/* called by the C++ code for initialization */
extern "C" void device_init (int8* cfg)
{
  /* initialize ray tracing core */
  rtcInit(cfg);

  /* set error handler */
  rtcSetErrorFunction(error_handler);

  /* set start render mode */
  renderPixel = renderPixelStandard;
  //  renderPixel = renderPixelUV;	
}

void updateEdgeLevelBuffer( ISPCMesh* mesh, const Vec3fa& cam_pos, size_t startID, size_t endID )
{
  for (size_t f=startID; f<endID; f++) 
  {
    const int N = 4;
    for (size_t i=0; i<N; i++) {
      const int * index = (int *)&mesh->quads[f];	
      const Vec3fa v0 = mesh->positions[index[i]];
      const Vec3fa v1 = mesh->positions[index[(i+1)%N]];
      const Vec3fa edge = v1-v0;
      const Vec3fa P = 0.5f*(v1+v0);
      const Vec3fa dist = cam_pos - P;
      mesh->edge_level[f*4+i] = max(min(LEVEL_FACTOR*(0.5f*length(edge)/length(dist)),MAX_EDGE_LEVEL),MIN_EDGE_LEVEL); // FIXME: map mesh subdivlevel buffer
      assert( mesh->edge_level[f*4+i] >= 1.0f );
    }
  }
}

#if defined(ISPC)
task void updateEdgeLevelBufferTask( ISPCMesh* mesh, const Vec3fa& cam_pos )
{
  const size_t size = mesh->numQuads;
  const size_t startID = ((taskIndex+0)*size)/taskCount;
  const size_t endID   = ((taskIndex+1)*size)/taskCount;
  updateEdgeLevelBuffer(mesh,cam_pos,startID,endID);
}
#endif

void updateEdgeLevels(ISPCScene* scene_in, const Vec3fa& cam_pos)
{
  for (int i=0; i<g_ispc_scene->numMeshes; i++)
  {
    ISPCMesh* mesh = g_ispc_scene->meshes[i];
    if (mesh->numQuads == 0) continue;

    if (mesh->edge_level) 
    {
      unsigned int geomID = mesh->geomID;
#if defined(ISPC)
      launch[  getNumHWThreads() ] updateEdgeLevelBufferTask(mesh,cam_pos); 	           
#else
      updateEdgeLevelBuffer(mesh,cam_pos,0,mesh->numQuads);
#endif
      rtcUpdateBuffer(g_scene,geomID,RTC_LEVEL_BUFFER);
    }
  }

  for (size_t g=0; g<scene_in->numSubdivMeshes; g++)
  {
    ISPCSubdivMesh* mesh = g_ispc_scene->subdiv[g];
    unsigned int geomID = mesh->geomID;
    for (size_t f=0, e=0; f<mesh->numFaces; e+=mesh->verticesPerFace[f++]) {
      int N = mesh->verticesPerFace[f];
      for (size_t i=0; i<N; i++) {
        const Vec3fa v0 = mesh->positions[mesh->position_indices[e+(i+0)]];
        const Vec3fa v1 = mesh->positions[mesh->position_indices[e+(i+1)%N]];
        const Vec3fa edge = v1-v0;
        const Vec3fa P = 0.5f*(v1+v0);
	const Vec3fa dist = cam_pos - P;
        mesh->subdivlevel[e+i] = max(min(LEVEL_FACTOR*(0.5f*length(edge)/length(dist)),MAX_EDGE_LEVEL),MIN_EDGE_LEVEL);
      }
    }
  }
}

void displacementFunction(void* ptr, unsigned int geomID, int unsigned primID, 
                      const float* u,      /*!< u coordinates (source) */
                      const float* v,      /*!< v coordinates (source) */
                      const float* nx,     /*!< x coordinates of normal at point to displace (source) */
                      const float* ny,     /*!< y coordinates of normal at point to displace (source) */
                      const float* nz,     /*!< z coordinates of normal at point to displace (source) */
                      float* px,           /*!< x coordinates of points to displace (source and target) */
                      float* py,           /*!< y coordinates of points to displace (source and target) */
                      float* pz,           /*!< z coordinates of points to displace (source and target) */
                      size_t N)
{
#if 0
  for (size_t i = 0; i<N; i++) {
    const Vec3fa dP = 0.02f*Vec3fa(sin(100.0f*px[i]+0.5f),sin(100.0f*pz[i]+1.5f),cos(100.0f*py[i]));
    px[i] += dP.x; py[i] += dP.y; pz[i] += dP.z;
  }
#else
  for (size_t i = 0; i<N; i++) {
    const Vec3fa P = Vec3fa(px[i],py[i],pz[i]);
    const Vec3fa nor = Vec3fa(nx[i],ny[i],nz[i]);
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

void convertScene(ISPCScene* scene_in, const Vec3fa& p)
{
  /* add all meshes to the scene */
  for (int i=0; i<scene_in->numMeshes; i++)
  {
    /* get ith mesh */
    ISPCMesh* mesh = scene_in->meshes[i];

    if (mesh->numQuads)
    {   
      mesh->edge_level  = new float[mesh->numQuads*4];
      for (size_t i=0; i<4*mesh->numQuads; i++) 
	mesh->edge_level[i] = FIXED_EDGE_TESSELLATION_VALUE;
      
      /* create a triangle mesh */
      unsigned int geomID = rtcNewSubdivisionMesh (g_scene, RTC_GEOMETRY_STATIC, mesh->numQuads, mesh->numQuads*4, mesh->numVertices, 0, 0, 0);
      mesh->geomID = geomID;

      unsigned int* faces = (unsigned int*) rtcMapBuffer(g_scene, geomID, RTC_FACE_BUFFER);
      for (size_t i=0; i<mesh->numQuads; i++) faces[i] = 4;
      rtcUnmapBuffer(g_scene,geomID,RTC_FACE_BUFFER);
      
      rtcSetBuffer(g_scene, geomID, RTC_VERTEX_BUFFER, mesh->positions , 0, sizeof(Vec3fa  ));
      rtcSetBuffer(g_scene, geomID, RTC_INDEX_BUFFER,  mesh->quads     , 0, sizeof(unsigned int));
      rtcSetBuffer(g_scene, geomID, RTC_LEVEL_BUFFER,  mesh->edge_level, 0, sizeof(float));
      
#if ENABLE_DISPLACEMENTS == 1
      rtcSetDisplacementFunction(g_scene,geomID,(RTCDisplacementFunc)&displacementFunction,NULL);
#endif
     }
     else
      mesh->edge_level = NULL;
  }
  
  /* add all subdiv meshes to the scene */
  for (size_t i=0; i<scene_in->numSubdivMeshes; i++)
  {
    ISPCSubdivMesh* mesh = scene_in->subdiv[i];
    unsigned int geomID = rtcNewSubdivisionMesh(g_scene, RTC_GEOMETRY_STATIC, mesh->numFaces, mesh->numEdges, mesh->numVertices, 
							mesh->numEdgeCreases, mesh->numVertexCreases, mesh->numHoles);
    mesh->geomID = geomID;

    rtcSetBuffer(g_scene, geomID, RTC_VERTEX_BUFFER, mesh->positions, 0, sizeof(Vec3fa  ));
    rtcSetBuffer(g_scene, geomID, RTC_LEVEL_BUFFER,  mesh->subdivlevel, 0, sizeof(float));
    rtcSetBuffer(g_scene, geomID, RTC_INDEX_BUFFER,  mesh->position_indices  , 0, sizeof(unsigned int));
    rtcSetBuffer(g_scene, geomID, RTC_FACE_BUFFER,   mesh->verticesPerFace, 0, sizeof(unsigned int));
    rtcSetBuffer(g_scene, geomID, RTC_HOLE_BUFFER,   mesh->holes, 0, sizeof(unsigned int));
    rtcSetBuffer(g_scene, geomID, RTC_EDGE_CREASE_INDEX_BUFFER,    mesh->edge_creases,          0, 2*sizeof(unsigned int));
    rtcSetBuffer(g_scene, geomID, RTC_EDGE_CREASE_WEIGHT_BUFFER,   mesh->edge_crease_weights,   0, sizeof(float));
    rtcSetBuffer(g_scene, geomID, RTC_VERTEX_CREASE_INDEX_BUFFER,  mesh->vertex_creases,        0, sizeof(unsigned int));
    rtcSetBuffer(g_scene, geomID, RTC_VERTEX_CREASE_WEIGHT_BUFFER, mesh->vertex_crease_weights, 0, sizeof(float));

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
    rtcSetBuffer(scene, meshID, RTC_VERTEX_BUFFER, verts, 0, sizeof(Vec3fa));
    
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

Vec3fa rndColor(const int ID) 
{
  int r = ((ID+13)*17*23) & 255;
  int g = ((ID+15)*11*13) & 255;
  int b = ((ID+17)* 7*19) & 255;
  const float oneOver255f = 1.f/255.f;
  return Vec3fa(r*oneOver255f,g*oneOver255f,b*oneOver255f);
}


/* task that renders a single screen tile */
Vec3fa renderPixelStandard(float x, float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
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
  if (ray.geomID == RTC_INVALID_GEOMETRY_ID) return Vec3fa(0.0f,0.0f,1.0f);
  
  /* shade all rays that hit something */
  Vec3fa color = Vec3fa(1.0f);
#if 0
    color = rndColor(ray.geomID);
#else
  /* apply ambient light */
  Vec3fa Ng = normalize(ray.Ng);
  color = color*abs(dot(ray.dir,Ng));   
#endif
  return color;
}

/* task that renders a single screen tile */
void renderTile(int taskIndex, int* pixels,
                     const int width,
                     const int height, 
                     const float time,
                     const Vec3fa& vx, 
                     const Vec3fa& vy, 
                     const Vec3fa& vz, 
                     const Vec3fa& p,
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
    Vec3fa color = renderPixel(x,y,vx,vy,vz,p);

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

extern "C" void toggleOpenSubdiv(int key, int x, int y)
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

Vec3fa old_p; 
 
/* called by the C++ code to render */
extern "C" void device_render (int* pixels,
                           const int width,
                           const int height, 
                           const float time,
                           const Vec3fa& vx, 
                           const Vec3fa& vy, 
                           const Vec3fa& vz, 
                           const Vec3fa& p)
{
   Vec3fa cam_org = Vec3fa(p.x,p.y,p.z);

   /* create scene */
   if (g_scene == NULL)
  { 
    g_scene = rtcNewScene(RTC_SCENE_DYNAMIC,RTC_INTERSECT1);
    convertScene(g_ispc_scene,cam_org);

#if !defined(FORCE_FIXED_EDGE_TESSELLATION)
    updateEdgeLevels(g_ispc_scene, cam_org);
#endif

    old_p = p;

#if !defined(PARALLEL_COMMIT)
  rtcCommit (g_scene);
#else
  launch[ getNumHWThreads() ] parallelCommit(g_scene); 
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
     launch[ getNumHWThreads() ] parallelCommit(g_scene); 
#endif
    }
   }
#endif


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
