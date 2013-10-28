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

#include "tutorials_device.h"
#include "../tutorials/tutorials_common.h"

#include <sink/COIPipeline_sink.h>
#include <sink/COIProcess_sink.h>
#include <common/COIMacros_common.h>
#include <common/COISysInfo_common.h>
#include <common/COIEvent_common.h>

namespace embree
{
  ispc::Scene* scene = NULL;

  extern "C" void run_init(uint32_t         in_BufferCount,
                           void**           in_ppBufferPointers,
                           uint64_t*        in_pBufferLengths,
                           InitData*        in_pMiscData,
                           uint16_t         in_MiscDataLength,
                           void*            in_pReturnValue,
                           uint16_t         in_ReturnValueLength)
  {
    ispc::init(in_pMiscData->verbose);
  }

  extern "C" void run_set_scene(uint32_t         in_BufferCount,
                                void**           in_ppBufferPointers,
                                uint64_t*        in_pBufferLengths,
                                SetSceneData*    in_pMiscData,
                                uint16_t         in_MiscDataLength,
                                void*            in_pReturnValue,
                                uint16_t         in_ReturnValueLength)
  {

    scene = new ispc::Scene;
    memcpy(scene->materials = malloc(in_pBufferLengths[0]),in_ppBufferPointers[0],in_pBufferLengths[0]);
    memcpy(scene->positions = malloc(in_pBufferLengths[1]),in_ppBufferPointers[1],in_pBufferLengths[1]);
    memcpy(scene->normals   = malloc(in_pBufferLengths[2]),in_ppBufferPointers[2],in_pBufferLengths[2]);
    memcpy(scene->texcoords = malloc(in_pBufferLengths[3]),in_ppBufferPointers[3],in_pBufferLengths[3]);
    memcpy(scene->triangles = malloc(in_pBufferLengths[4]),in_ppBufferPointers[4],in_pBufferLengths[4]);
    scene->numMaterials = in_pMiscData->numMaterials;
    scene->numVertices = in_pMiscData->numVertices;
    scene->numTriangles = in_pMiscData->numTriangles;
    scene->pointLightPosition = in_pMiscData->pointLightPosition;
    scene->pointLightIntensity = in_pMiscData->pointLightIntensity;
    scene->ambientLightIntensity = in_pMiscData->ambientLightIntensity;
    ispc::set_scene(scene);
  }

  extern "C" void run_render(uint32_t         in_BufferCount,
                             void**           in_ppBufferPointers,
                             uint64_t*        in_pBufferLengths,
                             RenderData*      in_pMiscData,
                             uint16_t         in_MiscDataLength,
                             void*            in_pReturnValue,
                             uint16_t         in_ReturnValueLength)
  {
    ispc::render(in_ppBufferPointers[0],
                 in_pMiscData->width,
                 in_pMiscData->height,
                 in_pMiscData->time, 
                 (ispc::vec3f&)in_pMiscData->vx, 
                 (ispc::vec3f&)in_pMiscData->vy, 
                 (ispc::vec3f&)in_pMiscData->vz,
                 (ispc::vec3f&)in_pMiscData->p);


  }

  extern "C" void run_cleanup(uint32_t         in_BufferCount,
                              void**           in_ppBufferPointers,
                              uint64_t*        in_pBufferLengths,
                              void*            in_pMiscData,
                              uint16_t         in_MiscDataLength,
                              void*            in_pReturnValue,
                              uint16_t         in_ReturnValueLength)
  {
    ispc::cleanup();
    if (scene) {
      free (scene->materials);
      free (scene->positions);
      free (scene->normals);
      free (scene->texcoords);
      free (scene->triangles);
      delete scene;
      scene = NULL;
    }
  }
}

int main(int argc, char** argv) 
{
  UNUSED_ATTR COIRESULT result;
  UNREFERENCED_PARAM (argc);
  UNREFERENCED_PARAM (argv);

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
