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

#include "coi_device.h"
#include "sys/stl/string.h"
#include "sys/sysinfo.h"
#include "sys/filename.h"
#include "image/image.h"

namespace embree
{
  COIDevice::COIProcess::COIProcess (int cardID, const char* executable, size_t numThreads, size_t verbose)
  {
    /* assume MIC executable to be in same folder as host executable */
    FileName fname = FileName(getExecutableFileName()).path() + FileName(executable);

    /* get engine handle */
    COIRESULT result;
    result = COIEngineGetHandle( COI_ISA_MIC, cardID, &engine );
    if (result != COI_SUCCESS)
      throw std::runtime_error("Failed to load engine number " + std::stringOf(cardID) + ": " + COIResultGetName(result));
    
    /* print info of engine */
    COI_ENGINE_INFO info;
    result = COIEngineGetInfo(engine,sizeof(info),&info);
    if (result != COI_SUCCESS) 
      throw std::runtime_error("COIEngineGetInfo failed: "+std::string(COIResultGetName(result)));
    
    std::cout << "Found Xeon Phi device with " << info.NumCores << " cores and " << (info.PhysicalMemory/1024/1024) << "MB memory" << std::endl;
    std::string strNumThreads = std::stringOf(numThreads);
    std::string strVerbose = std::stringOf(verbose);
    const char* argv[2] = { strNumThreads.c_str(), strVerbose.c_str() }; 

    /* create process */
    result = COIProcessCreateFromFile
      (engine,
       fname.c_str(),    // The local path to the sink side binary to launch.
       2, argv,          // argc and argv for the sink process.
       false, NULL,      // Environment variables to set for the sink process.
       true, NULL,       // Enable the proxy but don't specify a proxy root path.
       0,                // The amount of memory to reserve for COIBuffers.
       NULL,             // Path to search for dependencies
       &process          // The resulting process handle.
       );
    // if fail check loading by name
    if (result != COI_SUCCESS) {
      fname = FileName(executable);
      result = COIProcessCreateFromFile
        (engine,
         fname.c_str(),    // The local path to the sink side binary to launch.
         2, argv,          // argc and argv for the sink process.
         false, NULL,      // Environment variables to set for the sink process.
         true, NULL,       // Enable the proxy but don't specify a proxy root path.
         0,                // The amount of memory to reserve for COIBuffers.
         NULL,             // Path to search for dependencies
         &process          // The resulting process handle.
         );
    }
    
    if (result != COI_SUCCESS) 
      throw std::runtime_error("Failed to create process " + std::string(executable) +": " + COIResultGetName(result));
 
    /* create pipeline */
    COI_CPU_MASK cpuMask;
    COIPipelineClearCPUMask(&cpuMask);
    COIPipelineSetCPUMask(process,info.NumCores-1,0,&cpuMask);
    COIPipelineSetCPUMask(process,info.NumCores-1,1,&cpuMask);
    COIPipelineSetCPUMask(process,info.NumCores-1,2,&cpuMask);
    COIPipelineSetCPUMask(process,info.NumCores-1,3,&cpuMask);
    result = COIPipelineCreate(process,cpuMask,0,&pipeline);
    if (result != COI_SUCCESS) 
      throw std::runtime_error(std::string("Failed to create pipeline : ") + COIResultGetName(result));

    /* get run functions */
    const char *fctNameArray[] = { 
      "rtNewCamera", 
      "rtNewData", 
      "rtNewImage", 
      "rtNewTexture", 
      "rtNewMaterial", 
      "rtNewShape", 
      "rtNewLight", 
      "rtNewShapePrimitive", 
      "rtNewLightPrimitive", 
      "rtTransformPrimitive", 
      "rtNewScene", 
      "rtSetPrimitive", 
      "rtNewToneMapper", 
      "rtNewRenderer", 
      "rtNewFrameBuffer", 
      "rtSwapBuffers", 
      "rtIncRef", 
      "rtDecRef", 
      "rtSetBool1", 
      "rtSetBool2", 
      "rtSetBool3", 
      "rtSetBool4", 
      "rtSetInt1", 
      "rtSetInt2", 
      "rtSetInt3", 
      "rtSetInt4", 
      "rtSetFloat1", 
      "rtSetFloat2", 
      "rtSetFloat3", 
      "rtSetFloat4", 
      "rtSetArray", 
      "rtSetString", 
      "rtSetImage", 
      "rtSetTexture", 
      "rtSetTransform", 
      "rtClear", 
      "rtCommit", 
      "rtRenderFrame", 
      "rtPick",
      "rtNewDataStart",
      "rtNewDataSet",
      "rtNewDataEnd"
    };

    result = COIProcessGetFunctionHandles (process, sizeof(fctNameArray)/sizeof(char*), fctNameArray, &runNewCamera);
    if (result != COI_SUCCESS) 
      throw std::runtime_error("COIProcessGetFunctionHandles failed: "+std::string(COIResultGetName(result)));

    result = COIBufferCreate(STREAM_BUFFER_SIZE,COI_BUFFER_NORMAL,0,NULL,1,&process,&stream);
    if (result != COI_SUCCESS) throw std::runtime_error("COIBufferCreate failed: " + std::string(COIResultGetName(result)));
  }
  
  void COIDevice::COIProcess::loadLibrary (const char* library)
  {
    std::cout << "Loading library from file \"" << library << "\"" << std::endl;

    COILIBRARY lib;
    COIRESULT result = COIProcessLoadLibraryFromFile(process,library,library,NULL,&lib);
    if (result != COI_SUCCESS) 
      throw std::runtime_error(std::string("Failed to load libary: ") + COIResultGetName(result));
    
    libs.push_back(lib);
  }

  COIDevice::COIProcess::~COIProcess ()
  {
    for (size_t i=0; i<libs.size(); i++)
    {
      COIRESULT result = COIProcessUnloadLibrary(process,libs[i]);
      if (result != COI_SUCCESS) 
        throw std::runtime_error(std::string("Unloading library failed: ") + std::string(COIResultGetName(result)));
    }

    COIRESULT result = COIProcessDestroy(process,-1,0,NULL,NULL);
    if (result != COI_SUCCESS) 
      throw std::runtime_error(std::string("Destroying COI process failed: ") + std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::free (int id) {
    swapchains.erase(id);
  }

  COIRESULT COIPipelineRunFunctionSync(COIPIPELINE               in_Pipeline,
                                       COIFUNCTION               in_Function,
                                       uint32_t                  in_NumBuffers,
                                       const   COIBUFFER*        in_pBuffers,
                                       const   COI_ACCESS_FLAGS* in_pBufferAccessFlags,
                                       uint32_t                  in_NumDependencies,
                                       const   COIEVENT*         in_pDependencies,
                                       const   void*             in_pMiscData,
                                       uint16_t                  in_MiscDataLen,
                                       void*                     out_pAsyncReturnValue,
                                       uint16_t                  in_AsyncReturnValueLen)
  {
    COIEVENT completion;
    COIRESULT result = COIPipelineRunFunction (in_Pipeline,
                                               in_Function,
                                               in_NumBuffers,
                                               in_pBuffers,
                                               in_pBufferAccessFlags,
                                               in_NumDependencies,
                                               in_pDependencies,
                                               in_pMiscData,
                                               in_MiscDataLen,
                                               out_pAsyncReturnValue,
                                               in_AsyncReturnValueLen,
                                               &completion);
    
    COIEventWait(1,&completion,-1,true,NULL,NULL);
    return result;
  }

  /*******************************************************************
                      creation of objects
  *******************************************************************/

  void COIDevice::COIProcess::rtNewCamera(Device::RTCamera id, const char* type) 
  {
    parmsNewCamera parms;
    parms.id = (int) (long) id;
    strncpy(parms.type,type,sizeof(parms.type)-1);

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewCamera, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunctionSync failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtNewData(Device::RTData id, const char* type, size_t bytes, const void* data) 
  {
#if 0

    parmsNewData parms;
    parms.id = (int) (long) id;
    parms.bytes = bytes;

    COIBUFFER buffer;   
    COIRESULT result = COIBufferCreate(max(bytes,size_t(1)),COI_BUFFER_STREAMING_TO_SINK,0,data,1,&process,&buffer);
    if (result != COI_SUCCESS) throw std::runtime_error("COIBufferCreate failed: " + std::string(COIResultGetName(result)));
    
    COI_ACCESS_FLAGS flag = COI_SINK_READ;
    result = COIPipelineRunFunctionSync (pipeline, runNewData, 1, &buffer, &flag, 0, NULL, &parms, sizeof(parms), NULL, 0);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunctionSync failed: "+std::string(COIResultGetName(result)));    

    result = COIBufferDestroy(buffer);
    if (result != COI_SUCCESS) throw std::runtime_error("COIBufferDestroy failed: "+std::string(COIResultGetName(result)));

#else

    parmsNewDataStart parms0;
    parms0.id = (int) (long) id;
    parms0.bytes = bytes;

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewDataStart, 0, NULL, NULL, 0, NULL, &parms0, sizeof(parms0), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));  

    for (size_t i=0; i<bytes; i+=STREAM_BUFFER_SIZE) 
    {
      void* dst = NULL;
      size_t offset = i;
      size_t dbytes = min(size_t(STREAM_BUFFER_SIZE),bytes-i);

      COIEVENT completion;
      COIRESULT result = COIBufferWrite(stream,0,(char*)data+offset,dbytes,COI_COPY_USE_DMA,0,NULL,&completion);
      if (result != COI_SUCCESS) throw std::runtime_error("COIBufferWrite failed: "+std::string(COIResultGetName(result)));
      COIEventWait(1,&completion,-1,true,NULL,NULL);

      parmsNewDataSet parms1;
      parms1.offset = offset;
      parms1.bytes = dbytes;

      COI_ACCESS_FLAGS flag = COI_SINK_READ;
      result = COIPipelineRunFunctionSync (pipeline, runNewDataSet, 1, &stream, &flag, 0, NULL, &parms1, sizeof(parms1), NULL, 0);
      if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));  
    }

    COI_ACCESS_FLAGS flag = COI_SINK_READ;
    result = COIPipelineRunFunction (pipeline, runNewDataEnd, 0, NULL, NULL, 0, NULL, &parms0, sizeof(parms0), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));  

#endif
  }

  void COIDevice::COIProcess::rtNewImage(Device::RTImage id, const char* type, size_t width, size_t height, const void* data)
  {
    parmsNewImage parms;
    parms.id = (int) (long) id;
    strncpy(parms.type,type,sizeof(parms.type)-1);
    parms.width = width;
    parms.height = height;

    size_t bytes = 0;
    if      (!strcasecmp(type,"RGB8"        )) bytes = width*height*3*sizeof(char);
    else if (!strcasecmp(type,"RGBA8"       )) bytes = width*height*4*sizeof(char);
    else if (!strcasecmp(type,"RGB_FLOAT32" )) bytes = width*height*3*sizeof(float);
    else if (!strcasecmp(type,"RGBA_FLOAT32")) bytes = width*height*4*sizeof(float);
    else throw std::runtime_error("unknown image type: "+std::string(type)); 

    COIBUFFER buffer;   
    COIRESULT result = COIBufferCreate(max(bytes,size_t(1)),COI_BUFFER_STREAMING_TO_SINK,0,data,1,&process,&buffer);
    if (result != COI_SUCCESS) throw std::runtime_error("COIBufferCreate failed: " + std::string(COIResultGetName(result)));
    
    COI_ACCESS_FLAGS flag = COI_SINK_READ;
    result = COIPipelineRunFunctionSync (pipeline, runNewImage, 1, &buffer, &flag, 0, NULL, &parms, sizeof(parms), NULL, 0);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));    

    result = COIBufferDestroy(buffer);
    if (result != COI_SUCCESS) throw std::runtime_error("COIBufferDestroy failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtNewTexture(Device::RTTexture id, const char* type)
  {
    parmsNewTexture parms;
    parms.id = (int) (long) id;
    strncpy(parms.type,type,sizeof(parms.type)-1);

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewTexture, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtNewMaterial(Device::RTMaterial id, const char* type)
  {
    parmsNewMaterial parms;
    parms.id = (int) (long) id;
    strncpy(parms.type,type,sizeof(parms.type)-1);

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewMaterial, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtNewShape(Device::RTShape id, const char* type)
  {
    parmsNewShape parms;
    parms.id = (int) (long) id;
    strncpy(parms.type,type,sizeof(parms.type)-1);

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewShape, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }
  
  void COIDevice::COIProcess::rtNewLight(Device::RTLight id, const char* type)
  {
    parmsNewLight parms;
    parms.id = (int) (long) id;
    strncpy(parms.type,type,sizeof(parms.type)-1);

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewLight, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtNewShapePrimitive(Device::RTPrimitive id, 
                                       Device::RTShape shape, 
                                       Device::RTMaterial material, 
                                       const float* transform)
  {
    parmsNewShapePrimitive parms;
    parms.id = (int) (long) id;
    parms.shape = (int) (long) shape;
    parms.material = (int) (long) material;
    if (transform) {
      for (size_t i=0; i<12; i++)
        parms.transform[i] = transform[i];
    } else {
      parms.transform[ 0] = 1.0f; parms.transform[ 1] = 0.0f; parms.transform[ 2] = 0.0f; 
      parms.transform[ 3] = 0.0f; parms.transform[ 4] = 1.0f; parms.transform[ 5] = 0.0f; 
      parms.transform[ 6] = 0.0f; parms.transform[ 7] = 0.0f; parms.transform[ 8] = 1.0f; 
      parms.transform[ 9] = 0.0f; parms.transform[10] = 0.0f; parms.transform[11] = 0.0f; 
    }

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewShapePrimitive, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }
  
  void COIDevice::COIProcess::rtNewLightPrimitive(Device::RTPrimitive id, 
                                       Device::RTLight light, 
                                       Device::RTMaterial material, 
                                       const float* transform)
  {
    parmsNewLightPrimitive parms;
    parms.id = (int) (long) id;
    parms.light = (int) (long) light;
    parms.material = (int) (long) material;
    if (transform) {
      for (size_t i=0; i<12; i++)
        parms.transform[i] = transform[i];
    } else {
      parms.transform[ 0] = 1.0f; parms.transform[ 1] = 0.0f; parms.transform[ 2] = 0.0f; 
      parms.transform[ 3] = 0.0f; parms.transform[ 4] = 1.0f; parms.transform[ 5] = 0.0f; 
      parms.transform[ 6] = 0.0f; parms.transform[ 7] = 0.0f; parms.transform[ 8] = 1.0f; 
      parms.transform[ 9] = 0.0f; parms.transform[10] = 0.0f; parms.transform[11] = 0.0f; 
    }

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewLightPrimitive, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtTransformPrimitive(Device::RTPrimitive id, Device::RTPrimitive primitive, const float* transform) 
  {
    parmsTransformPrimitive parms;
    parms.id = (int) (long) id;
    parms.primitive = (int) (long) primitive;
    if (transform) {
      for (size_t i=0; i<12; i++)
        parms.transform[i] = transform[i];
    } else {
      parms.transform[ 0] = 1.0f; parms.transform[ 1] = 0.0f; parms.transform[ 2] = 0.0f; 
      parms.transform[ 3] = 0.0f; parms.transform[ 4] = 1.0f; parms.transform[ 5] = 0.0f; 
      parms.transform[ 6] = 0.0f; parms.transform[ 7] = 0.0f; parms.transform[ 8] = 1.0f; 
      parms.transform[ 9] = 0.0f; parms.transform[10] = 0.0f; parms.transform[11] = 0.0f; 
    }

    COIRESULT result = COIPipelineRunFunction (pipeline, runTransformPrimitive, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtNewScene(Device::RTScene id, const char* type)
  {
    parmsNewScene parms;
    parms.id = (int) (long) id;
    strncpy(parms.type,type,sizeof(parms.type)-1);

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewScene, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetPrimitive(RTScene scene, size_t slot, RTPrimitive prim)
  {
    parmsSetPrimitive parms;
    parms.scene = (int) (long) scene;
    parms.slot = slot;
    parms.prim = (int) (long) prim;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetPrimitive, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtNewToneMapper(Device::RTToneMapper id, const char* type)
  {
    parmsNewToneMapper parms;
    parms.id = (int) (long) id;
    strncpy(parms.type,type,sizeof(parms.type)-1);

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewToneMapper, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtNewRenderer(Device::RTRenderer id, const char* type)
  {
    parmsNewRenderer parms;
    parms.id = (int) (long) id;
    strncpy(parms.type,type,sizeof(parms.type)-1);

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewRenderer, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtNewFrameBuffer(Device::RTFrameBuffer id, const char* type, size_t width, size_t height, size_t numBuffers)
  {
    SwapChain* swapchain = NULL;
    if      (!strcasecmp(type,"RGB_FLOAT32")) swapchain = new SwapChain(&process,width*height*3*sizeof(float),numBuffers);
    else if (!strcasecmp(type,"RGBA8"      )) swapchain = new SwapChain(&process,width*height*4,numBuffers);
    else if (!strcasecmp(type,"RGB8"       )) swapchain = new SwapChain(&process,width*height*3,numBuffers);
    else throw std::runtime_error("unknown framebuffer type: "+std::string(type));
    swapchains[(int)(long)id] = swapchain;

    parmsNewFrameBuffer parms;
    parms.id = (int) (long) id;
    strncpy(parms.type,type,sizeof(parms.type)-1);
    parms.width = width;
    parms.height = height;
    parms.depth = numBuffers;

    COIRESULT result = COIPipelineRunFunction (pipeline, runNewFrameBuffer, 
                                                   numBuffers, swapchain->buffer, swapchain->flag, 
                                                   0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    
    if (result != COI_SUCCESS) 
      throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void* COIDevice::COIProcess::rtMapFrameBuffer(RTFrameBuffer frameBuffer, int bufID)
  {
    Ref<SwapChain> swapchain = swapchains[(int)(long)frameBuffer];
    if (bufID < 0) bufID = swapchain->curBuffer;
    return swapchain->map(bufID);
  }

  void COIDevice::COIProcess::rtUnmapFrameBuffer(RTFrameBuffer frameBuffer, int bufID)
  {
    Ref<SwapChain> swapchain = swapchains[(int)(long)frameBuffer];
    if (bufID < 0) bufID = swapchain->curBuffer;
    return swapchain->unmap(bufID);
  }

  void COIDevice::COIProcess::rtSwapBuffers(Device::RTFrameBuffer frameBuffer) 
  {
    parmsSwapBuffers parms;
    parms.framebuffer = (int) (long) frameBuffer;
    
    COIRESULT result = COIPipelineRunFunction (pipeline, runSwapBuffers, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));

    Ref<SwapChain> swapchain = swapchains[(int)(long)frameBuffer];
    swapchain->swapBuffers();
  }
  
  void COIDevice::COIProcess::rtIncRef(Device::RTHandle handle)
  {
    parmsIncRef parms;
    parms.handle = (int) (long) handle;
    
    COIRESULT result = COIPipelineRunFunction (pipeline, runIncRef, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtDecRef(Device::RTHandle handle)
  {
    parmsDecRef parms;
    parms.handle = (int) (long) handle;
    
    COIRESULT result = COIPipelineRunFunction (pipeline, runDecRef, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  /*******************************************************************
                    setting of parameters
  *******************************************************************/

  void COIDevice::COIProcess::rtSetBool1(Device::RTHandle handle, const char* property, bool x)
  {
    parmsSetBoolN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = false;
    parms.z = false;
    parms.w = false;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetBool1, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }
  
  void COIDevice::COIProcess::rtSetBool2(Device::RTHandle handle, const char* property, bool x, bool y)
  {
    parmsSetBoolN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = y;
    parms.z = false;
    parms.w = false;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetBool2, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
 }

  void COIDevice::COIProcess::rtSetBool3(Device::RTHandle handle, const char* property, bool x, bool y, bool z)
  {
    parmsSetBoolN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = y;
    parms.z = z;
    parms.w = false;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetBool3, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetBool4(Device::RTHandle handle, const char* property, bool x, bool y, bool z, bool w)
  {
    parmsSetBoolN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = y;
    parms.z = z;
    parms.w = w;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetBool4, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetInt1(Device::RTHandle handle, const char* property, int x)
  {
    parmsSetIntN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = 0;
    parms.z = 0;
    parms.w = 0;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetInt1, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetInt2(Device::RTHandle handle, const char* property, int x, int y)
  {
    parmsSetIntN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = y;
    parms.z = 0;
    parms.w = 0;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetInt2, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetInt3(Device::RTHandle handle, const char* property, int x, int y, int z)
  {
    parmsSetIntN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = y;
    parms.z = z;
    parms.w = 0;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetInt3, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetInt4(Device::RTHandle handle, const char* property, int x, int y, int z, int w)
  {
    parmsSetIntN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = y;
    parms.z = z;
    parms.w = w;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetInt4, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetFloat1(Device::RTHandle handle, const char* property, float x)
  {
    parmsSetFloatN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = 0.0f;
    parms.z = 0.0f;
    parms.w = 0.0f;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetFloat1, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetFloat2(Device::RTHandle handle, const char* property, float x, float y)
  {
    parmsSetFloatN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = y;
    parms.z = 0.0f;
    parms.w = 0.0f;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetFloat2, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetFloat3(Device::RTHandle handle, const char* property, float x, float y, float z)
  {
    parmsSetFloatN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = y;
    parms.z = z;
    parms.w = 0.0f;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetFloat3, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetFloat4(Device::RTHandle handle, const char* property, float x, float y, float z, float w)
  {
    parmsSetFloatN parms;
    parms.handle = (int) (long) handle;
    parms.x = x;
    parms.y = y;
    parms.z = z;
    parms.w = w;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetFloat4, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetArray(Device::RTHandle handle, const char* property, const char* type, Device::RTData data, size_t size, size_t stride, size_t ofs)
  {
    parmsSetArray parms;
    parms.handle = (int) (long) handle;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    strncpy(parms.type,type,sizeof(parms.type)-1);
    parms.data = (int) (long) data;
    parms.size = size;
    parms.stride = stride;
    parms.ofs = ofs;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetArray, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetString(Device::RTHandle handle, const char* property, const char* str)
  {
    parmsSetString parms;
    parms.handle = (int) (long) handle;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    strncpy(parms.str,str,sizeof(parms.str)-1);

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetString, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetImage(Device::RTHandle handle, const char* property, Device::RTImage image)
  {
    parmsSetImage parms;
    parms.handle = (int) (long) handle;
    parms.image = (int) (long) image;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetImage, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetTexture(Device::RTHandle handle, const char* property, Device::RTTexture texture)
  {
    parmsSetTexture parms;
    parms.handle = (int) (long) handle;
    parms.texture = (int) (long) texture;
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetTexture, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtSetTransform(Device::RTHandle handle, const char* property, const float* transform)
  {
    parmsSetTransform parms;
    parms.handle = (int) (long) handle;
    for (size_t i=0; i<12; i++) parms.transform[i] = transform[i];
    strncpy(parms.property,property,sizeof(parms.property)-1);
    size_t zeros = sizeof(parms.property)-strlen(parms.property)-1;

    COIRESULT result = COIPipelineRunFunction (pipeline, runSetTransform, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms)-zeros, NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtClear(Device::RTHandle handle)
  {
    parmsClear parms;
    parms.handle = (int) (long) handle;

    COIRESULT result = COIPipelineRunFunction (pipeline, runClear, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  void COIDevice::COIProcess::rtCommit(Device::RTHandle handle)
  {
    parmsCommit parms;
    parms.handle = (int) (long) handle;

    COIRESULT result = COIPipelineRunFunction (pipeline, runCommit, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), NULL, 0, NULL);
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  /*******************************************************************
                          render calls
  *******************************************************************/

  void COIDevice::COIProcess::rtRenderFrame(Device::RTRenderer renderer, Device::RTCamera camera, Device::RTScene scene, 
                                            Device::RTToneMapper toneMapper, Device::RTFrameBuffer frameBuffer, int accumulate)
  {
    parmsRenderFrame parms;
    parms.renderer = (int) (long) renderer;
    parms.camera = (int) (long) camera;
    parms.scene = (int) (long) scene;
    parms.toneMapper = (int) (long) toneMapper;
    parms.frameBuffer = (int) (long) frameBuffer;
    parms.accumulate = accumulate;

    Ref<SwapChain> swapchain = swapchains[(int)(long)frameBuffer];

    COIRESULT result = COIPipelineRunFunction (pipeline, runRenderFrame, 
                                               1, &swapchain->buffer[swapchain->curBuffer], &swapchain->flag[swapchain->curBuffer], 
                                               0, NULL, &parms, sizeof(parms), NULL, 0, NULL);

    if (result != COI_SUCCESS) 
      throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));
  }

  bool COIDevice::COIProcess::rtPick(Device::RTCamera camera, float x, float y, Device::RTScene scene, float& px, float& py, float& pz) 
  {
    parmsPick parms;
    parms.camera = (int) (long) camera;
    parms.x = x;
    parms.y = y;
    parms.scene = (int) (long) scene;

    returnPick ret;
    COIRESULT result = COIPipelineRunFunctionSync (pipeline, runPick, 0, NULL, NULL, 0, NULL, &parms, sizeof(parms), &ret, sizeof(ret));
    if (result != COI_SUCCESS) throw std::runtime_error("COIPipelineRunFunction failed: "+std::string(COIResultGetName(result)));

    /*! return pick data */
    px = ret.x;
    py = ret.y;
    pz = ret.z;
    return ret.hit;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  COIDevice::COIDevice(const char* executable, size_t numThreads, size_t verbose)
    : nextHandle(1), serverID(0), serverCount(1) 
  {
    uint32_t engines = 0;
    COIEngineGetCount( COI_ISA_MIC, &engines );
    if (engines == 0) throw std::runtime_error("no Xeon Phi device found");

    /* initialize all devices */
    for (uint32_t i=0; i<engines; i++) 
      devices.push_back(new COIProcess(i,executable,numThreads,verbose));

    /* set server ID of devices */
    for (size_t i=0; i<devices.size(); i++) 
    {
      int id    = serverID*devices.size()+i;
      int count = serverCount*devices.size();
      devices[i]->rtSetInt1(NULL,"serverID",id);
      devices[i]->rtSetInt1(NULL,"serverCount",count);
    }

    /* dummy 0 handle */
    counters.push_back(0);
    buffers.push_back(NULL);
  }

  COIDevice::~COIDevice()
  {
    for (size_t i=0; i<devices.size(); i++) delete devices[i];
    devices.clear();
  }

  /*******************************************************************
                     handle ID allocations
  *******************************************************************/

  int COIDevice::allocHandle() 
  {
    if (pool.empty()) {
      pool.push_back(nextHandle++);
      counters.push_back(0);
      buffers.push_back(NULL);
    }
    int id = pool.back();
    counters[id] = 1;
    pool.pop_back();
    return id;
  }
  
  void COIDevice::incRef(int id) {
    Lock<MutexSys> lock(handleMutex);
    counters[id]++;
  }
  
  bool COIDevice::decRef(int id) 
  {
    Lock<MutexSys> lock(handleMutex);
    if (--counters[id] == 0) {
      pool.push_back((int)id);
      buffers[id] = null;
      for (size_t i=0; i<devices.size(); i++) 
        devices[i]->free((int)id);
      return true;
    }
    return false;
  }

  /*******************************************************************
                      creation of objects
  *******************************************************************/

  Device::RTCamera COIDevice::rtNewCamera(const char* type) 
  { 
    Device::RTCamera id = (Device::RTCamera) allocHandle();
    for (size_t i=0; i<devices.size(); i++) 
      devices[i]->rtNewCamera(id,type);
    return id;
  }

  Device::RTData COIDevice::rtNewData(const char* type, size_t bytes, const void* data) 
  { 
    Device::RTData id = (Device::RTData) allocHandle();

    for (size_t i=0; i<devices.size(); i++) 
      devices[i]->rtNewData(id,type,bytes,data);

    if (!strcasecmp(type,"immutable_managed")) 
      alignedFree((void*)data);

    return id;
  }

  Device::RTData COIDevice::rtNewDataFromFile(const char* type, const char* fileName, size_t offset, size_t bytes)
  { 
    if (strcasecmp(type,"immutable"))
      throw std::runtime_error("unknown data type: "+(std::string)type);
    
    /*! read data from file */
    FILE* file = fopen(fileName,"rb");
    if (!file) throw std::runtime_error("cannot open file "+(std::string)fileName);
    fseek(file,(long)offset,SEEK_SET);
    
    char* data = (char*) alignedMalloc(bytes);
    if (bytes != fread(data,1,sizeof(bytes),file))
      throw std::runtime_error("error filling data buffer from file");
    fclose(file);
    
    return rtNewData("immutable_managed", bytes, data);
  }

  Device::RTImage COIDevice::rtNewImage(const char* type, size_t width, size_t height, const void* data, const bool copy)
  { 
    Device::RTImage id = (Device::RTImage) allocHandle();

    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtNewImage(id,type,width,height,data);

    if (!copy) free((void*)data);

    return id;
  }

  Device::RTImage COIDevice::rtNewImageFromFile(const char* fileName)
  { 
    /*! load image locally */
    Ref<Image> image = loadImage(fileName);
    if (!image) throw std::runtime_error("cannot load image: "+std::string(fileName));
    else if (Ref<Image3c> cimg = image.dynamicCast<Image3c>())
      return rtNewImage("RGB8",cimg->width,cimg->height,cimg->steal_ptr(),false);
    else if (Ref<Image4c> cimg = image.dynamicCast<Image4c>())
      return rtNewImage("RGBA8",cimg->width,cimg->height,cimg->steal_ptr(),false);
    else if (Ref<Image3f> fimg = image.dynamicCast<Image3f>())
      return rtNewImage("RGB_FLOAT32",fimg->width,fimg->height,fimg->steal_ptr(),false);
    else if (Ref<Image4f> fimg = image.dynamicCast<Image4f>())
      return rtNewImage("RGBA_FLOAT32",fimg->width,fimg->height,fimg->steal_ptr(),false);
    else
      throw std::runtime_error("unknown image type");
  }

  Device::RTTexture COIDevice::rtNewTexture(const char* type)
  { 
    Device::RTTexture id = (Device::RTTexture) allocHandle();
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtNewTexture(id,type);
    return id;
  }

  Device::RTMaterial COIDevice::rtNewMaterial(const char* type)
  { 
    Device::RTMaterial id = (Device::RTMaterial) allocHandle();
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtNewMaterial(id,type);
    return id;
  }

  Device::RTShape COIDevice::rtNewShape(const char* type)
  { 
    Device::RTShape id = (Device::RTShape) allocHandle();
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtNewShape(id,type);
    return id;
  }

  Device::RTLight COIDevice::rtNewLight(const char* type)
  { 
    Device::RTLight id = (Device::RTLight) allocHandle();
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtNewLight(id,type);
    return id;
  }

  Device::RTPrimitive COIDevice::rtNewShapePrimitive(Device::RTShape shape, Device::RTMaterial material, const float* transform)
  { 
    Device::RTPrimitive id = (Device::RTPrimitive) allocHandle();
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtNewShapePrimitive(id,shape,material,transform);
    return id;
  }

  Device::RTPrimitive COIDevice::rtNewLightPrimitive(Device::RTLight light, Device::RTMaterial material, const float* transform)
  { 
    Device::RTPrimitive id = (Device::RTPrimitive) allocHandle();
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtNewLightPrimitive(id,light,material,transform);
    return id;
  }

  Device::RTPrimitive COIDevice::rtTransformPrimitive(Device::RTPrimitive primitive, const float* transform) 
  { 
    Device::RTPrimitive id = (Device::RTPrimitive) allocHandle();
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtTransformPrimitive(id,primitive,transform);
    return id;
  }

  Device::RTScene COIDevice::rtNewScene(const char* type)
  { 
    Device::RTScene id = (Device::RTScene) allocHandle();
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtNewScene(id,type);
    return id;
  }

  void COIDevice::rtSetPrimitive(RTScene scene, size_t slot, RTPrimitive prim)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetPrimitive(scene,slot,prim);
  }

  Device::RTToneMapper COIDevice::rtNewToneMapper(const char* type)
  { 
    Device::RTToneMapper id = (Device::RTToneMapper) allocHandle();
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtNewToneMapper(id,type);
    return id;
  }

  Device::RTRenderer COIDevice::rtNewRenderer(const char* type)
  { 
    Device::RTRenderer id = (Device::RTRenderer) allocHandle();
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtNewRenderer(id,type);
    return id;
  }

  Device::RTFrameBuffer COIDevice::rtNewFrameBuffer(const char* type, size_t width, size_t height, size_t numBuffers, void** ptrs)
  { 
    int id = allocHandle();
    Device::RTFrameBuffer hid = (Device::RTFrameBuffer) id;

    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtNewFrameBuffer(hid,type,width,height,numBuffers);

    if      (!strcasecmp(type,"RGB_FLOAT32")) buffers[id] = new SwapChain(type,width,height,numBuffers,ptrs,FrameBufferRGBFloat32::create);
    else if (!strcasecmp(type,"RGBA8"      )) buffers[id] = new SwapChain(type,width,height,numBuffers,ptrs,FrameBufferRGBA8     ::create);
    else if (!strcasecmp(type,"RGB8"       )) buffers[id] = new SwapChain(type,width,height,numBuffers,ptrs,FrameBufferRGB8      ::create);
    else throw std::runtime_error("unknown framebuffer type: "+std::string(type));

    return hid;
  }

  void* COIDevice::rtMapFrameBuffer(Device::RTFrameBuffer frameBuffer, int bufID) 
  { 
    Ref<SwapChain>& swapChain = buffers[(size_t)frameBuffer];
    if (bufID < 0) bufID = swapChain->id();

    if (devices.size() == 1) 
      return devices[0]->rtMapFrameBuffer(frameBuffer,bufID);

    /* map framebuffers of all devices */
    std::vector<char*> ptrs(devices.size());
    for (size_t i=0; i<devices.size(); i++)
      ptrs[i] = (char*) devices[i]->rtMapFrameBuffer(frameBuffer,bufID);
    
    /* merge images from different devices */
    Ref<FrameBuffer>& buffer = swapChain->buffer(bufID);
    char* dptr = (char*) buffer->getData();
    size_t stride = buffer->getStride();
    for (size_t y=0; y<buffer->getHeight(); y++) {
      size_t row = y/4, subrow = y%4;
      size_t devRow = row/devices.size();
      size_t devID  = row%devices.size();
      char* src = ptrs[devID]+stride*(4*devRow+subrow);
      char* dst = dptr+y*stride;
      memcpy(dst,src,stride);
    }

    /* unmap framebuffers of all devices */
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtUnmapFrameBuffer(frameBuffer,bufID);

    return dptr;
  }

  void COIDevice::rtUnmapFrameBuffer(Device::RTFrameBuffer frameBuffer, int bufID)
  { 
    Ref<SwapChain>& swapChain = buffers[(size_t)frameBuffer];
    if (bufID < 0) bufID = swapChain->id();
    
    if (devices.size() == 1) 
      devices[0]->rtUnmapFrameBuffer(frameBuffer,bufID);
  }

  void COIDevice::rtSwapBuffers(Device::RTFrameBuffer frameBuffer) 
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSwapBuffers(frameBuffer);

    Ref<SwapChain>& swapchain = buffers[(int)(long)frameBuffer];
    swapchain->swapBuffers();
  }
  
  void COIDevice::rtIncRef(Device::RTHandle handle)
  { 
    incRef((int)(size_t)handle);
  }

  void COIDevice::rtDecRef(Device::RTHandle handle)
  { 
    if (decRef((int)(size_t)handle)) {
      for (size_t i=0; i<devices.size(); i++)
        devices[i]->rtDecRef(handle);
    }
  }

  /*******************************************************************
                    setting of parameters
  *******************************************************************/

  void COIDevice::rtSetBool1(Device::RTHandle handle, const char* property, bool x)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetBool1(handle,property,x);
  }
  
  void COIDevice::rtSetBool2(Device::RTHandle handle, const char* property, bool x, bool y)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetBool2(handle,property,x,y);
 }

  void COIDevice::rtSetBool3(Device::RTHandle handle, const char* property, bool x, bool y, bool z)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetBool3(handle,property,x,y,z);
  }

  void COIDevice::rtSetBool4(Device::RTHandle handle, const char* property, bool x, bool y, bool z, bool w)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetBool4(handle,property,x,y,z,w);
  }

  void COIDevice::rtSetInt1(Device::RTHandle handle, const char* property, int x)
  { 
    if (!handle) {
      if (!strcmp(property,"serverID"   )) 
        serverID = x;
      else if (!strcmp(property,"serverCount")) {
        serverCount = x;
        for (size_t i=0; i<devices.size(); i++) {
          int id    = serverID*devices.size()+i;
          int count = serverCount*devices.size();
          devices[i]->rtSetInt1(NULL,"serverID",id);
          devices[i]->rtSetInt1(NULL,"serverCount",count);
        }
      }
      else
        for (size_t i=0; i<devices.size(); i++)
          devices[i]->rtSetInt1(handle,property,x);
    }
    else
      for (size_t i=0; i<devices.size(); i++)
        devices[i]->rtSetInt1(handle,property,x);
  }

  void COIDevice::rtSetInt2(Device::RTHandle handle, const char* property, int x, int y)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetInt2(handle,property,x,y);
  }

  void COIDevice::rtSetInt3(Device::RTHandle handle, const char* property, int x, int y, int z)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetInt3(handle,property,x,y,z);
  }

  void COIDevice::rtSetInt4(Device::RTHandle handle, const char* property, int x, int y, int z, int w)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetInt4(handle,property,x,y,z,w);
  }

  void COIDevice::rtSetFloat1(Device::RTHandle handle, const char* property, float x)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetFloat1(handle,property,x);
  }

  void COIDevice::rtSetFloat2(Device::RTHandle handle, const char* property, float x, float y)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetFloat2(handle,property,x,y);
  }

  void COIDevice::rtSetFloat3(Device::RTHandle handle, const char* property, float x, float y, float z)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetFloat3(handle,property,x,y,z);
  }

  void COIDevice::rtSetFloat4(Device::RTHandle handle, const char* property, float x, float y, float z, float w)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetFloat4(handle,property,x,y,z,w);
  }

  void COIDevice::rtSetArray(Device::RTHandle handle, const char* property, const char* type, Device::RTData data, size_t size, size_t stride, size_t ofs)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetArray(handle,property,type,data,size,stride,ofs);
  }

  void COIDevice::rtSetString(Device::RTHandle handle, const char* property, const char* str)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetString(handle,property,str);
  }

  void COIDevice::rtSetImage(Device::RTHandle handle, const char* property, Device::RTImage img)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetImage(handle,property,img);
  }

  void COIDevice::rtSetTexture(Device::RTHandle handle, const char* property, Device::RTTexture tex)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetTexture(handle,property,tex);
  }

  void COIDevice::rtSetTransform(Device::RTHandle handle, const char* property, const float* transform)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtSetTransform(handle,property,transform);
  }

  void COIDevice::rtClear(Device::RTHandle handle)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtClear(handle);
  }

  void COIDevice::rtCommit(Device::RTHandle handle)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtCommit(handle);
  }

  /*******************************************************************
                          render calls
  *******************************************************************/

  void COIDevice::rtRenderFrame(Device::RTRenderer renderer, Device::RTCamera camera, Device::RTScene scene, 
                                    Device::RTToneMapper toneMapper, Device::RTFrameBuffer frameBuffer, int accumulate)
  { 
    for (size_t i=0; i<devices.size(); i++)
      devices[i]->rtRenderFrame(renderer,camera,scene,toneMapper,frameBuffer,accumulate);
  }

  bool COIDevice::rtPick(Device::RTCamera camera, float x, float y, Device::RTScene scene, float& px, float& py, float& pz) 
  { 
    return devices[0]->rtPick(camera,x,y,scene,px,py,pz);
  }

  __dllexport Device* create(const char* parms, size_t numThreads, size_t verbose) 
  {
    return new COIDevice(parms,numThreads,verbose);
  }
}
