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

#include "tutorial/obj_loader.h"
#include "transport/transport_host.h"
#include "transport/transport_device.h"
#include "transport_coi/common.h"

#include <sink/COIPipeline_sink.h>
#include <sink/COIProcess_sink.h>
#include <common/COIMacros_common.h>
#include <common/COISysInfo_common.h>
#include <common/COIEvent_common.h>


extern "C" int64 get_tsc() {
  return read_tsc();
}

namespace embree
{
  /* ISPC compatible mesh */
  struct ISPCMesh
  {
    ALIGNED_CLASS;
  public:
    ISPCMesh (int numTriangles, int numQuads, int numVertices) 
      : numTriangles(numTriangles), numQuads(numQuads), numVertices(numVertices),
        positions(NULL), positions2(NULL), normals(NULL), texcoords(NULL), triangles(NULL), quads(NULL), edge_level(NULL)
    {
      sizePositions = 0;
      sizeNormals   = 0;
      sizeTexCoords = 0;
      sizeTriangles = 0;
      sizeQuads     = 0;
    }

    ~ISPCMesh () {
      if (positions)  os_free(positions ,sizePositions);
      if (positions2) os_free(positions2,sizePositions);
      if (normals)    os_free(normals   ,sizeNormals);
      if (texcoords)  os_free(texcoords ,sizeTexCoords);
      if (triangles)  os_free(triangles ,sizeTriangles);
      if (quads)      os_free(quads     ,sizeQuads);

      positions = NULL;
      positions2 = NULL;
      normals   = NULL;
      texcoords = NULL;
      triangles = NULL;
      quads = NULL;
    }

  public:
    Vec3fa* positions;     //!< vertex position array
    Vec3fa* positions2;    //!< vertex position array
    Vec3fa* normals;       //!< vertex normal array
    Vec2f*  texcoords;     //!< vertex texcoord array
    OBJScene::Triangle* triangles;  //!< list of triangles
    OBJScene::Quad* quads;          //!< list of quads
    float *edge_level;

    int numVertices;
    int numTriangles;
    int numQuads;
    int geomID;
    size_t sizePositions;
    size_t sizeNormals;
    size_t sizeTexCoords;
    size_t sizeTriangles;
    size_t sizeQuads;
  };

  /* ISPC compatible scene */
  struct ISPCHair
  {
  public:
    ISPCHair () {}
    ISPCHair (int vertex, int id)
      : vertex(vertex), id(id) {}

    int vertex,id;  //!< index of first control point and hair ID
  };

  /*! Hair Set. */
  struct ISPCHairSet
  {
    ALIGNED_CLASS;
  public:
    Vec3fa *positions;   //!< hair control points (x,y,z,r)
    Vec3fa *positions2;   //!< hair control points (x,y,z,r)
    ISPCHair *hairs;    //!< list of hairs
    int numVertices;
    int numHairs;
    ISPCHairSet(int numHairs, int numVertices) 
      : numHairs(numHairs),numVertices(numVertices),positions(NULL),positions2(NULL),hairs(NULL) {}
    ~ISPCHairSet() {
      if (positions) free(positions);
      if (positions2) free(positions2);
      if (hairs) free(hairs);
    }
  };

  /* ISPC compatible scene */
  struct ISPCScene
  {
    ALIGNED_CLASS;
  public:
    ISPCScene (int numMeshes, int numHairSets, 
               void* materials_in, int numMaterials,
               void* ambientLights_in, int numAmbientLights,
               void* pointLights_in, int numPointLights,
               void* directionalLights_in, int numDirectionalLights,
               void* distantLights_in, int numDistantLights)

      : meshes(NULL), numMeshes(numMeshes), numHairSets(numHairSets), 
        materials(NULL), numMaterials(numMaterials),
        ambientLights(NULL), numAmbientLights(numAmbientLights),
        pointLights(NULL), numPointLights(numPointLights),
        directionalLights(NULL), numDirectionalLights(numDirectionalLights),
        distantLights(NULL), numDistantLights(numDistantLights),
	subdiv(NULL), numSubdivMeshes(0)
      {
        meshes = new ISPCMesh*[numMeshes];
        for (size_t i=0; i<numMeshes; i++)
          meshes[i] = NULL;

        hairsets = new ISPCHairSet*[numHairSets];
        for (size_t i=0; i<numHairSets; i++)
          hairsets[i] = NULL;
        
        materials = new OBJScene::Material[numMaterials];
        memcpy(materials,materials_in,numMaterials*sizeof(OBJScene::Material));

        ambientLights = new OBJScene::AmbientLight[numAmbientLights];
        memcpy(ambientLights,ambientLights_in,numAmbientLights*sizeof(OBJScene::AmbientLight));

        pointLights = new OBJScene::PointLight[numPointLights];
        memcpy(pointLights,pointLights_in,numPointLights*sizeof(OBJScene::PointLight));

        directionalLights = new OBJScene::DirectionalLight[numDirectionalLights];
        memcpy(directionalLights,directionalLights_in,numDirectionalLights*sizeof(OBJScene::DirectionalLight));

        distantLights = new OBJScene::DistantLight[numDistantLights];
        memcpy(distantLights,distantLights_in,numDistantLights*sizeof(OBJScene::DistantLight));
      }

    ~ISPCScene () 
    {
      delete[] materials;
      delete[] ambientLights;
      delete[] pointLights;
      delete[] directionalLights;
      delete[] distantLights;

      if (meshes) {
        for (size_t i=0; i<numMeshes; i++)
          if (meshes[i]) delete meshes[i];
	delete[] meshes;
	meshes = NULL;
      }
    }

  public:
    ISPCMesh** meshes;
    OBJScene::Material* materials;  //!< material list
    int numMeshes;
    int numMaterials;

    ISPCHairSet** hairsets;
    int numHairSets;

    OBJScene::AmbientLight* ambientLights;
    int numAmbientLights;

    OBJScene::PointLight* pointLights;
    int numPointLights;

    OBJScene::DirectionalLight* directionalLights;
    int numDirectionalLights;

    OBJScene::DistantLight* distantLights;
    int numDistantLights;

    //ISPCSubdivMesh** subdiv;
    void* subdiv;
    int numSubdivMeshes; 
  };

  /* scene */
  static size_t g_meshID = 0;
  static size_t g_hairsetID = 0;

  extern "C" ISPCScene* g_ispc_scene = NULL;

  extern "C" void run_init(uint32_t         in_BufferCount,
                           void**           in_ppBufferPointers,
                           uint64_t*        in_pBufferLengths,
                           InitData*        in_pMiscData,
                           uint16_t         in_MiscDataLength,
                           void*            in_pReturnValue,
                           uint16_t         in_ReturnValueLength)
  {
    device_init(in_pMiscData->cfg);
  }

  extern "C" void run_key_pressed(uint32_t         in_BufferCount,
                                  void**           in_ppBufferPointers,
                                  uint64_t*        in_pBufferLengths,
                                  KeyPressedData* in_pMiscData,
                                  uint16_t         in_MiscDataLength,
                                  void*            in_pReturnValue,
                                  uint16_t         in_ReturnValueLength)
  {
    device_key_pressed(in_pMiscData->key);
  }

  extern "C" void run_create_mesh(uint32_t         in_BufferCount,
                                  void**           in_ppBufferPointers,
                                  uint64_t*        in_pBufferLengths,
                                  CreateMeshData*  in_pMiscData,
                                  uint16_t         in_MiscDataLength,
                                  void*            in_pReturnValue,
                                  uint16_t         in_ReturnValueLength)
  {
    size_t meshID = g_meshID++;

#if 0
    DBG_PRINT( in_pMiscData->numTriangles );
    DBG_PRINT( in_pMiscData->numQuads );
    DBG_PRINT( in_pMiscData->numVertices );
#endif

    ISPCMesh* mesh = new ISPCMesh(in_pMiscData->numTriangles,in_pMiscData->numQuads,in_pMiscData->numVertices);
    assert( mesh );
    assert( in_pMiscData->numTriangles*sizeof(OBJScene::Triangle) == in_pBufferLengths[3] );
    assert( in_pMiscData->numQuads*sizeof(OBJScene::Quad) == in_pBufferLengths[4] );

    //assert( in_pMiscData->numVertices*sizeof(Vec3fa) == in_pBufferLengths[1] );

    mesh->positions = (Vec3fa*)os_malloc(in_pBufferLengths[0]);
    mesh->normals   = (Vec3fa*)os_malloc(in_pBufferLengths[1]);
    mesh->texcoords = (Vec2f* )os_malloc(in_pBufferLengths[2]);
    mesh->triangles = (OBJScene::Triangle*)os_malloc(in_pBufferLengths[3]);
    mesh->quads     = (OBJScene::Quad*)os_malloc(in_pBufferLengths[4]);

    memcpy(mesh->positions,in_ppBufferPointers[0],in_pBufferLengths[0]);
    memcpy(mesh->normals  ,in_ppBufferPointers[1],in_pBufferLengths[1]);
    memcpy(mesh->texcoords,in_ppBufferPointers[2],in_pBufferLengths[2]);
    memcpy(mesh->triangles,in_ppBufferPointers[3],in_pBufferLengths[3]);
    memcpy(mesh->quads    ,in_ppBufferPointers[4],in_pBufferLengths[4]);

    mesh->sizePositions = in_pBufferLengths[0];
    mesh->sizeNormals   = in_pBufferLengths[1];
    mesh->sizeTexCoords = in_pBufferLengths[2];
    mesh->sizeTriangles = in_pBufferLengths[3];
    mesh->sizeQuads     = in_pBufferLengths[4];
    
#if 0
    DBG_PRINT( mesh->sizePositions );
    DBG_PRINT( mesh->sizeNormals );
    DBG_PRINT( mesh->sizeTexCoords );
    DBG_PRINT( mesh->sizeTriangles );
    DBG_PRINT( mesh->sizeQuads );
#endif

    g_ispc_scene->meshes[meshID] = mesh;
  }

  extern "C" void run_create_hairset(uint32_t         in_BufferCount,
				     void**           in_ppBufferPointers,
				     uint64_t*        in_pBufferLengths,
				     CreateHairSetData*  in_pMiscData,
				     uint16_t         in_MiscDataLength,
				     void*            in_pReturnValue,
				     uint16_t         in_ReturnValueLength)
  {
    size_t hairsetID = g_hairsetID++;
    ISPCHairSet* hairset = new ISPCHairSet(in_pMiscData->numHairs,in_pMiscData->numVertices);
    memcpy(hairset->positions = (Vec3fa*)malloc(in_pBufferLengths[0]),in_ppBufferPointers[0],in_pBufferLengths[0]);
    memcpy(hairset->hairs = (ISPCHair*)malloc(in_pBufferLengths[1]),in_ppBufferPointers[1],in_pBufferLengths[1]);
    g_ispc_scene->hairsets[hairsetID] = hairset;
  }

  extern "C" void run_create_scene(uint32_t         in_BufferCount,
                                   void**           in_ppBufferPointers,
                                   uint64_t*        in_pBufferLengths,
                                   CreateSceneData* in_pMiscData,
                                   uint16_t         in_MiscDataLength,
                                   void*            in_pReturnValue,
                                   uint16_t         in_ReturnValueLength)
  {
    g_meshID = 0;
    g_ispc_scene = new ISPCScene(in_pMiscData->numMeshes,
                                 in_pMiscData->numHairSets,
                                 in_ppBufferPointers[0],in_pMiscData->numMaterials,
                                 in_ppBufferPointers[1],in_pMiscData->numAmbientLights,
                                 in_ppBufferPointers[2],in_pMiscData->numPointLights,
                                 in_ppBufferPointers[3],in_pMiscData->numDirectionalLights,
                                 in_ppBufferPointers[4],in_pMiscData->numDistantLights);
  }

  extern "C" void run_pick(uint32_t         in_BufferCount,
                           void**           in_ppBufferPointers,
                           uint64_t*        in_pBufferLengths,
                           PickDataSend*    in_pMiscData,
                           uint16_t         in_MiscDataLength,
                           PickDataReceive* in_pReturnValue,
                           uint16_t         in_ReturnValueLength)
  {
    Vec3fa hitPos = zero;
    bool hit = device_pick(in_pMiscData->x,
                           in_pMiscData->y,
                           in_pMiscData->vx,
                           in_pMiscData->vy,
                           in_pMiscData->vz,
                           in_pMiscData->p,
                           hitPos);
    in_pReturnValue->pos = hitPos;
    in_pReturnValue->hit = hit;
  }

  extern "C" void run_render(uint32_t         in_BufferCount,
                             void**           in_ppBufferPointers,
                             uint64_t*        in_pBufferLengths,
                             RenderData*      in_pMiscData,
                             uint16_t         in_MiscDataLength,
                             void*            in_pReturnValue,
                             uint16_t         in_ReturnValueLength)
  {
    //double t0 = getSeconds();
    device_render((int*)in_ppBufferPointers[0],
                  in_pMiscData->width,
                  in_pMiscData->height,
                  in_pMiscData->time, 
                  in_pMiscData->vx, 
                  in_pMiscData->vy, 
                  in_pMiscData->vz,
                  in_pMiscData->p);
    //double dt = getSeconds() - t0;
    //printf("render %3.2f fps, %.2f ms\n",1.0f/dt,dt*1000.0f); flush(std::cout);
  }

  extern "C" void run_cleanup(uint32_t         in_BufferCount,
                              void**           in_ppBufferPointers,
                              uint64_t*        in_pBufferLengths,
                              void*            in_pMiscData,
                              uint16_t         in_MiscDataLength,
                              void*            in_pReturnValue,
                              uint16_t         in_ReturnValueLength)
  {
    device_cleanup();
    if (g_ispc_scene) delete g_ispc_scene; g_ispc_scene = NULL;
  }
}

int main(int argc, char** argv) 
{
  UNUSED_ATTR COIRESULT result;
  UNREFERENCED_PARAM (argc);
  UNREFERENCED_PARAM (argv);

  /* enable wait to attach with debugger */
#if 0
  std::cout << "waiting for debugger to attach ..." << std::flush;
#if 0
  volatile int loop = 1;
  do {
    volatile int a = 1;
  } while (loop);
  
#else
  for (int i=0; i<20; i++) {
    sleep(1);
    std::cout << "." << std::flush;
  }
#endif
  std::cout << " [DONE]" << std::endl;
#endif
  
  // Functions enqueued on the sink side will not start executing until
  // you call COIPipelineStartExecutingRunFunctions(). This call is to
  // synchronize any initialization required on the sink side
  result = COIPipelineStartExecutingRunFunctions();
  assert(result == COI_SUCCESS);

  // This call will wait until COIProcessDestroy() gets called on the source
  // side. If COIProcessDestroy is called without force flag set, this call
  // will make sure all the functions enqueued are executed and does all
  // clean up required to exit gracefully.
  COIProcessWaitForShutdown();
  return 0;
}
