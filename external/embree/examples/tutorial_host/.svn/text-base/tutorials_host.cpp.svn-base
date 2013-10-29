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

#include "tutorials_host.h"
#include "../tutorials/tutorials_common.h"
#include "../tutorials/obj_loader.h"
#include "sys/stl/string.h"

#include <source/COIProcess_source.h>
#include <source/COIEngine_source.h>
#include <source/COIBuffer_source.h>
#include <source/COIPipeline_source.h>

namespace embree
{
  COIENGINE   engine;
  COIPROCESS  process;
  COIPIPELINE pipeline;
  
  COIFUNCTION runInit;
  COIFUNCTION runSetScene;
  COIFUNCTION runRender;
  COIFUNCTION runCleanup;

  COIBUFFER frameBuffer;
  COIMAPINSTANCE mapInst;
  int g_width = -1, g_height = -1;

  COIBUFFER materialBuffer;      //!< materials array
  COIBUFFER positionBuffer;      //!< vertex position array
  COIBUFFER normalBuffer;        //!< vertex normal array
  COIBUFFER texcoordBuffer;      //!< vertex texcoord array
  COIBUFFER triangleBuffer;      //!< list of triangles

  extern const char* tutorialName;
  
  void init(int verbose)
  {
    /* get number of Xeon Phi devices */
    uint32_t engines = 0;
    COIEngineGetCount( COI_ISA_MIC, &engines );
    if (engines == 0) throw std::runtime_error("No Xeon Phi device found.");

    /* get engine handle */
    COIRESULT result;
    result = COIEngineGetHandle( COI_ISA_MIC, 0, &engine );
    if (result != COI_SUCCESS)
      throw std::runtime_error("Failed to load engine number " + std::stringOf(0) + ": " + COIResultGetName(result));
    
    /* print info of engine */
    COI_ENGINE_INFO info;
    result = COIEngineGetInfo(engine,sizeof(info),&info);
    std::cout << "Found Xeon Phi device with " << info.NumCores << " cores and " << (info.PhysicalMemory/1024/1024) << "MB memory" << std::endl;
    
    /* create process */
    const std::string executable = std::string(tutorialName)+"_device_knc";
    result = COIProcessCreateFromFile
      (engine,
       executable.c_str(), // The local path to the sink side binary to launch.
       0, NULL,            // argc and argv for the sink process.
       false, NULL,        // Environment variables to set for the sink process.
       true, NULL,         // Enable the proxy but don't specify a proxy root path.
       0,                  // The amount of memory to reserve for COIBuffers.
       NULL,               // Path to search for dependencies
       &process            // The resulting process handle.
       );
    
    if (result != COI_SUCCESS) 
      throw std::runtime_error("Failed to create process " + std::string(executable) +": " + COIResultGetName(result));

    /* create pipeline */
    COI_CPU_MASK cpuMask;
    result = COIPipelineClearCPUMask(&cpuMask);
    if (result != COI_SUCCESS) 
      throw std::runtime_error(std::string("COIPipelineClearCPUMask failed: ") + COIResultGetName(result));

    result = COIPipelineSetCPUMask(process,info.NumCores-1,0,&cpuMask);
    result = COIPipelineSetCPUMask(process,info.NumCores-1,1,&cpuMask);
    result = COIPipelineSetCPUMask(process,info.NumCores-1,2,&cpuMask);
    result = COIPipelineSetCPUMask(process,info.NumCores-1,3,&cpuMask);
    if (result != COI_SUCCESS) 
      throw std::runtime_error(std::string("COIPipelineSetCPUMask failed: ") + COIResultGetName(result));

    result = COIPipelineCreate(process,cpuMask,0,&pipeline);
    if (result != COI_SUCCESS) 
      throw std::runtime_error(std::string("COIPipelineCreate failed: ") + COIResultGetName(result));

    /* get run functions */
    const char *fctNameArray[4] = { "run_init", "run_set_scene", "run_render", "run_cleanup" };
    result = COIProcessGetFunctionHandles (process, 4, fctNameArray, &runInit);
    if (result != COI_SUCCESS) 
      throw std::runtime_error("COIProcessGetFunctionHandles failed: "+std::string(COIResultGetName(result)));

    /* run init runfunction */
    InitData parms;
    parms.verbose = verbose;

    result = COIPipelineRunFunction (pipeline, runInit, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) 
      throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));

  }

  /* set scene */
  void setScene(OBJMesh& mesh, 
                const Vec3f pointLightPosition,
                const Vec3f pointLightIntensity,
                const Vec3f ambientLightIntensity)
  {
    COIRESULT result;

    /* send scene */
    size_t materialBytes = max(size_t(16),mesh.material.size()*sizeof(OBJMesh::Material));
    void* materialPtr = mesh.material.size() ? &mesh.material.front() : NULL;
    result = COIBufferCreate(materialBytes,COI_BUFFER_STREAMING_TO_SINK,0,materialPtr,1,&process,&materialBuffer);
    if (result != COI_SUCCESS) throw std::runtime_error("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    size_t positionBytes = max(size_t(16),mesh.v.size()*sizeof(Vec3fa));
    void* positionPtr = mesh.v.size() ? &mesh.v.front() : NULL;
    result = COIBufferCreate(positionBytes,COI_BUFFER_STREAMING_TO_SINK,0,positionPtr,1,&process,&positionBuffer);
    if (result != COI_SUCCESS) throw std::runtime_error("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    size_t normalBytes = max(size_t(16),mesh.vn.size()*sizeof(Vec3f));
    void* normalPtr = mesh.vn.size() ? &mesh.vn.front() : NULL;
    result = COIBufferCreate(normalBytes,COI_BUFFER_STREAMING_TO_SINK,0,normalPtr,1,&process,&normalBuffer);
    if (result != COI_SUCCESS) throw std::runtime_error("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    size_t texcoordBytes = max(size_t(16),mesh.vt.size()*sizeof(Vec2f));
    void* texcoordPtr = mesh.vt.size() ? &mesh.vt.front() : NULL;
    result = COIBufferCreate(texcoordBytes,COI_BUFFER_STREAMING_TO_SINK,0,texcoordPtr,1,&process,&texcoordBuffer);
    if (result != COI_SUCCESS) throw std::runtime_error("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    size_t triangleBytes = max(size_t(16),mesh.triangles.size()*sizeof(OBJMesh::Triangle));
    void* trianglePtr = mesh.triangles.size() ? &mesh.triangles.front() : NULL;
    result = COIBufferCreate(triangleBytes,COI_BUFFER_STREAMING_TO_SINK,0,trianglePtr,1,&process,&triangleBuffer);
    if (result != COI_SUCCESS) throw std::runtime_error("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    SetSceneData parms;
    parms.numMaterials = mesh.material.size();
    parms.numVertices = mesh.v.size();
    parms.numTriangles = mesh.triangles.size();
    parms.pointLightPosition = pointLightPosition;
    parms.pointLightIntensity = pointLightIntensity;
    parms.ambientLightIntensity = ambientLightIntensity;

    COI_ACCESS_FLAGS flags[5] = { COI_SINK_READ, COI_SINK_READ, COI_SINK_READ, COI_SINK_READ, COI_SINK_READ };

    /* run set scene runfunction */
    result = COIPipelineRunFunction (pipeline, runSetScene, 5, &materialBuffer, flags, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);

    if (result != COI_SUCCESS) 
      throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void resize(int32_t width, int32_t height)
  {
    COIRESULT result;
    if (g_width == width && g_height == height)
      return;

    /* destroy old framebuffer */
    if (g_width != -1 && g_height != -1) {
      result = COIBufferDestroy(frameBuffer);
      if (result != COI_SUCCESS)
        throw std::runtime_error("COIBufferDestroy failed: "+std::string(COIResultGetName(result)));
    }

    /* create new framebuffer */
    g_width  = width;
    g_height = height;
    result = COIBufferCreate (width*height*4, COI_BUFFER_STREAMING_TO_SOURCE, 0, NULL, 1, &process, &frameBuffer);

    if (result != COI_SUCCESS)  
      throw std::runtime_error("COIBufferCreate failed: " + std::string(COIResultGetName(result)));
  }

  int* render(const float time, const Vec3f& vx, const Vec3f& vy, const Vec3f& vz, const Vec3f& p)
  {
    /* set parameters for rendering */
    RenderData parms;
    parms.time = time;
    parms.vx = vx;
    parms.vy = vy;
    parms.vz = vz;
    parms.p  = p;
    parms.width = g_width;
    parms.height = g_height;
    COI_ACCESS_FLAGS flags = COI_SINK_WRITE_ENTIRE;

    /* run init runfunction */
    COIRESULT result = COIPipelineRunFunction (pipeline, runRender, 1, &frameBuffer, &flags, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);

    if (result != COI_SUCCESS) 
      throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));

    /* map the framebuffer */
    void* ptr = NULL;
    result = COIBufferMap(frameBuffer,0,g_width*g_height*4,COI_MAP_READ_ONLY,0,NULL,NULL,&mapInst,&ptr);
    if (result != COI_SUCCESS)
      throw std::runtime_error("COIBufferMap failed: "+std::string(COIResultGetName(result)));

    return (int*) ptr;
  }
  
  void unmap ()
  {
    COIRESULT result = COIBufferUnmap(mapInst,0,NULL,NULL);
    if (result != COI_SUCCESS)
      throw std::runtime_error("COIBufferUnmap failed: "+std::string(COIResultGetName(result)));
  }

  void cleanup()
  {
    /* run cleanup runfunction */
    COIRESULT result = COIPipelineRunFunction (pipeline, runCleanup, 0, NULL, NULL, 0, NULL, NULL, 0, NULL, 0, NULL);

    if (result != COI_SUCCESS) 
      throw std::runtime_error("Error launching runfunction: "+std::string(COIResultGetName(result)));

    /* destroy Xeon Phi process */
    result = COIProcessDestroy(process,-1,0,NULL,NULL);
    if (result != COI_SUCCESS) 
      throw std::runtime_error(std::string("Destroying COI process failed: ") + std::string(COIResultGetName(result)));
  }
}
