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

#include "transport/transport_host.h"
#include "transport_coi/common.h"
#include "tutorial/obj_loader.h"
#include "sys/stl/string.h"

#include <source/COIProcess_source.h>
#include <source/COIEngine_source.h>
#include <source/COIBuffer_source.h>
#include <source/COIPipeline_source.h>
#include <source/COIEvent_source.h>

namespace embree
{
  COIENGINE   engine;
  COIPROCESS  process;
  COIPIPELINE pipeline;
  
  COIFUNCTION runInit;
  COIFUNCTION runKeyPressed;
  COIFUNCTION runCreateMesh;
  COIFUNCTION runCreateHairSet;
  COIFUNCTION runCreateScene;
  COIFUNCTION runPick;
  COIFUNCTION runRender;
  COIFUNCTION runCleanup;

  COIBUFFER frameBuffer;
  COIMAPINSTANCE mapInst;
  int g_width = -1, g_height = -1;

  extern const char* tutorialName;

  void init(const char* cfg)
  {
    /* get number of Xeon Phi devices */
    uint32_t engines = 0;
    COIEngineGetCount( COI_ISA_MIC, &engines );
    if (engines == 0) THROW_RUNTIME_ERROR("No Xeon Phi device found.");

    /* get engine handle */
    COIRESULT result;
    result = COIEngineGetHandle( COI_ISA_MIC, 0, &engine );
    if (result != COI_SUCCESS)
      THROW_RUNTIME_ERROR("Failed to load engine number " + std::stringOf(0) + ": " + COIResultGetName(result));
    
    /* print info of engine */
    COI_ENGINE_INFO info;
    result = COIEngineGetInfo(engine,sizeof(info),&info);
    std::cout << "Found Xeon Phi device with " << info.NumCores << " cores and " << (info.PhysicalMemory/1024/1024) << "MB memory" << std::endl;
    
    /* create process */
    const std::string executable = std::string(tutorialName)+"_xeonphi_device";
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
      THROW_RUNTIME_ERROR("Failed to create process " + std::string(executable) +": " + COIResultGetName(result));

    /* create pipeline */
    COI_CPU_MASK cpuMask;
    result = COIPipelineClearCPUMask(&cpuMask);
    if (result != COI_SUCCESS) 
      THROW_RUNTIME_ERROR(std::string("COIPipelineClearCPUMask failed: ") + COIResultGetName(result));

    result = COIPipelineSetCPUMask(process,info.NumCores-1,0,&cpuMask);
    result = COIPipelineSetCPUMask(process,info.NumCores-1,1,&cpuMask);
    result = COIPipelineSetCPUMask(process,info.NumCores-1,2,&cpuMask);
    result = COIPipelineSetCPUMask(process,info.NumCores-1,3,&cpuMask);
    if (result != COI_SUCCESS) 
      THROW_RUNTIME_ERROR(std::string("COIPipelineSetCPUMask failed: ") + COIResultGetName(result));

    result = COIPipelineCreate(process,cpuMask,0,&pipeline);
    if (result != COI_SUCCESS) 
      THROW_RUNTIME_ERROR(std::string("COIPipelineCreate failed: ") + COIResultGetName(result));

    /* get run functions */
    const char *fctNameArray[8] = { "run_init", "run_key_pressed", "run_create_mesh", "run_create_hairset", "run_create_scene", "run_pick", "run_render", "run_cleanup" };
    result = COIProcessGetFunctionHandles (process, 8, fctNameArray, &runInit);
    if (result != COI_SUCCESS) 
      THROW_RUNTIME_ERROR("COIProcessGetFunctionHandles failed: "+std::string(COIResultGetName(result)));

    /* run init runfunction */
    InitData parms;
    strncpy(parms.cfg,cfg,sizeof(parms.cfg));
    result = COIPipelineRunFunction (pipeline, runInit, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) 
      THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));

  }

  void key_pressed(int32 key)
  {
    KeyPressedData parms;
    parms.key = key;
    COIRESULT result = COIPipelineRunFunction (pipeline, runKeyPressed, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) 
      THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void send_hairset (OBJScene::HairSet* hairset)
  {
    COIRESULT result;
    struct {
      COIBUFFER position;    //!< vertex position array
      COIBUFFER hairs;      //!< hair array
    } buffers;

    size_t positionBytes = max(size_t(16),hairset->v.size()*sizeof(Vec3fa));
    void* positionPtr = hairset->v.size() ? &hairset->v.front() : NULL;
    result = COIBufferCreate(positionBytes,COI_BUFFER_STREAMING_TO_SINK,0,positionPtr,1,&process,&buffers.position);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    size_t hairsBytes = max(size_t(16),hairset->hairs.size()*sizeof(OBJScene::Hair));
    void* hairsPtr = hairset->hairs.size() ? &hairset->hairs.front() : NULL;
    result = COIBufferCreate(hairsBytes,COI_BUFFER_STREAMING_TO_SINK,0,hairsPtr,1,&process,&buffers.hairs);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    CreateHairSetData parms;

    parms.numVertices = hairset->v.size();
    parms.numHairs    = hairset->hairs.size();
    COI_ACCESS_FLAGS flags[2] = { COI_SINK_READ, COI_SINK_READ};

    COIEVENT event;
    memset(&event,0,sizeof(event));

    /* run set scene runfunction */
    result = COIPipelineRunFunction (pipeline, runCreateHairSet, 2, &buffers.position, flags, 0, NULL, &parms, sizeof(parms), NULL, 0, &event);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
 
    result = COIEventWait(1,&event,-1,1,NULL,NULL);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIEventWait failed: "+std::string(COIResultGetName(result)));

    /* destroy buffers again */
    result = COIBufferDestroy(buffers.position);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
    result = COIBufferDestroy(buffers.hairs);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void send_mesh (OBJScene::Mesh* mesh)
  {
    COIRESULT result;
    struct {
      COIBUFFER position;      //!< vertex position array
      COIBUFFER normal;        //!< vertex normal array
      COIBUFFER texcoord;      //!< vertex texcoord array
      COIBUFFER triangle;      //!< list of triangles
      COIBUFFER quad;          //!< list of quads
    } buffers;

    assert( mesh->v.size() );

    if (mesh->triangles.size() == 0)
      {
	OBJScene::Triangle dummy(0,0,0,0);
	mesh->triangles.push_back(dummy);
      }

    if (mesh->vn.size() == 0)
      for (size_t i=0;i<4;i++)
	mesh->vn.push_back(Vec3f(0.0f,0.0f,0.0f));

    if (mesh->vt.size() == 0)
      for (size_t i=0;i<2;i++)
	mesh->vt.push_back(Vec2f(0.0f,0.0f));

    if (mesh->quads.size() == 0)
      {
	OBJScene::Quad dummy(0,0,0,0);
	mesh->quads.push_back(dummy);
      }
    
    assert( mesh->vn.size() );
    assert( mesh->vt.size() );
    assert( mesh->quads.size() );
    
    
    size_t positionBytes = max(size_t(16),mesh->v.size()*sizeof(Vec3fa));

    void* positionPtr = mesh->v.size() ? &mesh->v.front() : NULL;
    //result = COIBufferCreate(positionBytes,COI_BUFFER_STREAMING_TO_SINK,0,positionPtr,1,&process,&buffers.position);
    result = COIBufferCreateFromMemory(positionBytes,COI_BUFFER_NORMAL,0,positionPtr,1,&process,&buffers.position);

    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    size_t normalBytes = max(size_t(16),mesh->vn.size()*sizeof(Vec3fa));
    void* normalPtr = mesh->vn.size() ? &mesh->vn.front() : NULL;
    //result = COIBufferCreate(normalBytes,COI_BUFFER_STREAMING_TO_SINK,0,normalPtr,1,&process,&buffers.normal);
    result = COIBufferCreateFromMemory(normalBytes,COI_BUFFER_NORMAL,0,normalPtr,1,&process,&buffers.normal);

    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    size_t texcoordBytes = max(size_t(16),mesh->vt.size()*sizeof(Vec2f));
    void* texcoordPtr = mesh->vt.size() ? &mesh->vt.front() : NULL;
    //result = COIBufferCreate(texcoordBytes,COI_BUFFER_STREAMING_TO_SINK,0,texcoordPtr,1,&process,&buffers.texcoord);
    result = COIBufferCreateFromMemory(texcoordBytes,COI_BUFFER_NORMAL,0,texcoordPtr,1,&process,&buffers.texcoord);

    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    size_t triangleBytes = max(size_t(16),mesh->triangles.size()*sizeof(OBJScene::Triangle));
    void* trianglePtr = mesh->triangles.size() ? &mesh->triangles.front() : NULL;
    //result = COIBufferCreate(triangleBytes,COI_BUFFER_STREAMING_TO_SINK,0,trianglePtr,1,&process,&buffers.triangle);
    result = COIBufferCreateFromMemory(triangleBytes,COI_BUFFER_NORMAL,0,trianglePtr,1,&process,&buffers.triangle);

    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    size_t quadBytes = max(size_t(16),mesh->quads.size()*sizeof(OBJScene::Quad));
    void* quadPtr = mesh->quads.size() ? &mesh->quads.front() : NULL;
    //result = COIBufferCreate(quadBytes,COI_BUFFER_STREAMING_TO_SINK,0,quadPtr,1,&process,&buffers.quad);
    result = COIBufferCreateFromMemory(quadBytes,COI_BUFFER_NORMAL,0,quadPtr,1,&process,&buffers.quad);

    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    CreateMeshData parms;
    parms.numVertices = mesh->v.size();
    parms.numTriangles = mesh->triangles.size();
    parms.numQuads = mesh->quads.size();
    COI_ACCESS_FLAGS flags[5] = { COI_SINK_READ, COI_SINK_READ, COI_SINK_READ, COI_SINK_READ, COI_SINK_READ };

    COIEVENT event;
    memset(&event,0,sizeof(event));

    /* run set scene runfunction */
    result = COIPipelineRunFunction (pipeline, runCreateMesh, 5, &buffers.position, flags, 0, NULL, &parms, sizeof(parms), NULL, 0, &event);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
 
    result = COIEventWait(1,&event,-1,1,NULL,NULL);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIEventWait failed: "+std::string(COIResultGetName(result)));

    /* destroy buffers again */
    result = COIBufferDestroy(buffers.position);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
    result = COIBufferDestroy(buffers.normal);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
    result = COIBufferDestroy(buffers.texcoord);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
    result = COIBufferDestroy(buffers.triangle);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
    result = COIBufferDestroy(buffers.quad);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  /* set scene to use */
  void set_scene (OBJScene* scene)
  {
    COIRESULT result;
    COIBUFFER buffers[5];
    COI_ACCESS_FLAGS flags[5] = { COI_SINK_READ, COI_SINK_READ, COI_SINK_READ, COI_SINK_READ, COI_SINK_READ };

    /* send materials */
    size_t materialBytes = max(size_t(16),scene->materials.size()*sizeof(OBJScene::Material));
    void* materialPtr = scene->materials.size() ? &scene->materials.front() : NULL;
    result = COIBufferCreate(materialBytes,COI_BUFFER_STREAMING_TO_SINK,0,materialPtr,1,&process,&buffers[0]);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    /* send ambient lights */
    COIBUFFER ambientLightsBuffer;
    size_t ambientLightsBytes = max(size_t(16),scene->ambientLights.size()*sizeof(OBJScene::AmbientLight));
    void* ambientLightsPtr = scene->ambientLights.size() ? &scene->ambientLights.front() : NULL;
    result = COIBufferCreate(ambientLightsBytes,COI_BUFFER_STREAMING_TO_SINK,0,ambientLightsPtr,1,&process,&buffers[1]);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    /* send point lights */
    COIBUFFER pointLightsBuffer;
    size_t pointLightsBytes = max(size_t(16),scene->pointLights.size()*sizeof(OBJScene::PointLight));
    void* pointLightsPtr = scene->pointLights.size() ? &scene->pointLights.front() : NULL;
    result = COIBufferCreate(pointLightsBytes,COI_BUFFER_STREAMING_TO_SINK,0,pointLightsPtr,1,&process,&buffers[2]);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    /* send directional lights */
    COIBUFFER directionalLightsBuffer;
    size_t directionalLightsBytes = max(size_t(16),scene->directionalLights.size()*sizeof(OBJScene::DirectionalLight));
    void* directionalLightsPtr = scene->directionalLights.size() ? &scene->directionalLights.front() : NULL;
    result = COIBufferCreate(directionalLightsBytes,COI_BUFFER_STREAMING_TO_SINK,0,directionalLightsPtr,1,&process,&buffers[3]);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));

    /* send distant lights */
    COIBUFFER distantLightsBuffer;
    size_t distantLightsBytes = max(size_t(16),scene->distantLights.size()*sizeof(OBJScene::DistantLight));
    void* distantLightsPtr = scene->distantLights.size() ? &scene->distantLights.front() : NULL;
    result = COIBufferCreate(distantLightsBytes,COI_BUFFER_STREAMING_TO_SINK,0,distantLightsPtr,1,&process,&buffers[4]);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));
    
    CreateSceneData parms;
    parms.numMaterials = scene->materials.size();
    parms.numMeshes    = scene->meshes.size();
    parms.numHairSets  = scene->hairsets.size();
    parms.numAmbientLights = scene->ambientLights.size();
    parms.numPointLights = scene->pointLights.size();
    parms.numDirectionalLights = scene->directionalLights.size();
    parms.numDistantLights = scene->distantLights.size();

    COIEVENT event;
    memset(&event,0,sizeof(event));

    /* run set scene runfunction */
    result = COIPipelineRunFunction (pipeline, runCreateScene, 5, buffers, flags, 0, NULL, &parms, sizeof(parms), NULL, 0, &event);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));

    result = COIEventWait(1,&event,-1,1,NULL,NULL);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIEventWait failed: "+std::string(COIResultGetName(result)));

    /* destroy buffers again */
    // result = COIBufferDestroy(materialBuffer);
    // if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));

    /* send all meshes */
    for (size_t i=0; i<scene->meshes.size(); i++) 
      send_mesh(scene->meshes[i]);

    /* send all hairsets */
    for (size_t i=0; i<scene->hairsets.size(); i++) 
      send_hairset(scene->hairsets[i]);
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
        THROW_RUNTIME_ERROR("COIBufferDestroy failed: "+std::string(COIResultGetName(result)));
    }

    /* create new framebuffer */
    g_width  = width;
    g_height = height;
    result = COIBufferCreate (width*height*4, COI_BUFFER_NORMAL, COI_OPTIMIZE_SOURCE_READ|COI_OPTIMIZE_SINK_WRITE|COI_OPTIMIZE_HUGE_PAGE_SIZE, NULL, 1, &process, &frameBuffer);

    if (result != COI_SUCCESS)  
      THROW_RUNTIME_ERROR("COIBufferCreate failed: " + std::string(COIResultGetName(result)));
  }

  bool pick(const float x, const float y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p, Vec3fa& hitPos)
  {
    PickDataSend send;
    send.x = x; send.y = y;
    send.vx = vx; send.vy = vy; send.vz = vz; send.p = p;

    COIEVENT event;
    memset(&event,0,sizeof(event));

    PickDataReceive receive;
    COIRESULT result = COIPipelineRunFunction (pipeline, runPick, 0, NULL, NULL, 0, NULL, &send, sizeof(send), &receive, sizeof(receive), &event);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));

    result = COIEventWait(1,&event,-1,1,NULL,NULL);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIEventWait failed: "+std::string(COIResultGetName(result)));

    hitPos = receive.pos;
    return receive.hit;
  }

  void render(const float time, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
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

    COIEVENT event;
    memset(&event,0,sizeof(event));

    /* run init runfunction */
    COIRESULT result = COIPipelineRunFunction (pipeline, runRender, 1, &frameBuffer, &flags, 0, NULL, &parms, sizeof(parms), NULL, 0, &event);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));

    result = COIEventWait(1,&event,-1,1,NULL,NULL);
    if (result != COI_SUCCESS) THROW_RUNTIME_ERROR("COIEventWait failed: "+std::string(COIResultGetName(result)));
  }
  
  int* map ()
  {
    /* map the framebuffer */
    void* ptr = NULL;
    COIRESULT result = COIBufferMap(frameBuffer,0,g_width*g_height*4,COI_MAP_READ_ONLY,0,NULL,NULL,&mapInst,&ptr);
    if (result != COI_SUCCESS)
      THROW_RUNTIME_ERROR("COIBufferMap failed: "+std::string(COIResultGetName(result)));
    
    return (int*) ptr;
  }
  
  void unmap ()
  {
    COIRESULT result = COIBufferUnmap(mapInst,0,NULL,NULL);
    if (result != COI_SUCCESS)
      THROW_RUNTIME_ERROR("COIBufferUnmap failed: "+std::string(COIResultGetName(result)));
  }

  void cleanup()
  {
    /* run cleanup runfunction */
    COIRESULT result = COIPipelineRunFunction (pipeline, runCleanup, 0, NULL, NULL, 0, NULL, NULL, 0, NULL, 0, NULL);

    if (result != COI_SUCCESS) 
      THROW_RUNTIME_ERROR("Error launching runfunction: "+std::string(COIResultGetName(result)));

    /* destroy Xeon Phi process */
    result = COIProcessDestroy(process,-1,0,NULL,NULL);
    if (result != COI_SUCCESS) 
      THROW_RUNTIME_ERROR(std::string("Destroying COI process failed: ") + std::string(COIResultGetName(result)));
  }
}
