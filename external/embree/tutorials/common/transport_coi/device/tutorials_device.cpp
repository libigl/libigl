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
  return __rdtsc();
}

namespace embree
{
  /* ISPC compatible mesh */
  struct ISPCMesh
  {
    ALIGNED_CLASS;
  public:
    ISPCMesh (int numTriangles, int numVertices) 
      : numTriangles(numTriangles), numVertices(numVertices),
        positions(NULL), normals(NULL), texcoords(NULL), triangles(NULL), dir(zero), offset(zero) {}

    ~ISPCMesh () {
      if (positions) free(positions);
      if (normals) free(normals);
      if (texcoords) free(texcoords);
      if (triangles) free(triangles);
    }

  public:
    Vec3fa* positions;    //!< vertex position array
    Vec3fa* normals;       //!< vertex normal array
    Vec2f* texcoords;     //!< vertex texcoord array
    OBJScene::Triangle* triangles;  //!< list of triangles
    int numVertices;
    int numTriangles;
    Vec3f dir;
    float offset;
  };

  /* ISPC compatible scene */
  struct ISPCScene
  {
    ALIGNED_CLASS;
  public:
    ISPCScene (int numMeshes, int numMaterials, void* materials_in)
      : numMeshes(numMeshes), numMaterials(numMaterials),
        meshes(NULL), materials(NULL) 
    {
      meshes = new ISPCMesh*[numMeshes];
      for (size_t i=0; i<numMeshes; i++)
        meshes[i] = NULL;

      materials = new OBJScene::Material[numMaterials];
      memcpy(materials,materials_in,numMaterials*sizeof(OBJScene::Material));
    }

    ~ISPCScene () {
      delete[] materials;
      if (meshes) {
        for (size_t i=0; i<numMeshes; i++)
          if (meshes[i]) delete meshes[i];
        free(meshes);
      }
    }

  public:
    ISPCMesh** meshes;
    OBJScene::Material* materials;  //!< material list
    int numMeshes;
    int numMaterials;
    bool animate;
  };

  /* scene */
  static size_t g_meshID = 0;
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
    ISPCMesh* mesh = new ISPCMesh(in_pMiscData->numTriangles,in_pMiscData->numVertices);
    memcpy(mesh->positions = (Vec3fa*)malloc(in_pBufferLengths[0]),in_ppBufferPointers[0],in_pBufferLengths[0]);
    memcpy(mesh->normals   = (Vec3fa*)malloc(in_pBufferLengths[1]),in_ppBufferPointers[1],in_pBufferLengths[1]);
    memcpy(mesh->texcoords = (Vec2f* )malloc(in_pBufferLengths[2]),in_ppBufferPointers[2],in_pBufferLengths[2]);
    memcpy(mesh->triangles = (OBJScene::Triangle*)malloc(in_pBufferLengths[3]),in_ppBufferPointers[3],in_pBufferLengths[3]);
    g_ispc_scene->meshes[meshID] = mesh;
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
    g_ispc_scene = new ISPCScene(in_pMiscData->numMeshes,in_pMiscData->numMaterials,in_ppBufferPointers[0]);
    //g_ispc_scene->pointLightPosition = in_pMiscData->pointLightPosition;
    //g_ispc_scene->pointLightIntensity = in_pMiscData->pointLightIntensity;
    //g_ispc_scene->ambientLightIntensity = in_pMiscData->ambientLightIntensity;
  }

  extern "C" void run_pick(uint32_t         in_BufferCount,
                           void**           in_ppBufferPointers,
                           uint64_t*        in_pBufferLengths,
                           PickDataSend*    in_pMiscData,
                           uint16_t         in_MiscDataLength,
                           PickDataReceive* in_pReturnValue,
                           uint16_t         in_ReturnValueLength)
  {
    Vec3f hitPos = zero;
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
