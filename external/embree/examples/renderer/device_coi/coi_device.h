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

#ifndef __EMBREE_COI_DEVICE_H__
#define __EMBREE_COI_DEVICE_H__

#include "device/device.h"
#include "coi_common.h"
#include "device_singleray/api/swapchain.h"

#include <source/COIProcess_source.h>
#include <source/COIEngine_source.h>
#include <source/COIPipeline_source.h>
#include <source/COIBuffer_source.h>
#include <source/COIEvent_source.h>

#include <map>

#define STREAM_BUFFER_SIZE (1024*1024)

namespace embree
{
  class COIDevice : public Device
  {
  public:
    
    class COIProcess : public RefCount
    {
    public:
      COIProcess (int cardID, const char* executable, size_t numThreads, size_t verbose);
      ~COIProcess ();
      void loadLibrary (const char* library);
      
      struct SwapChain : public RefCount
      {
        SwapChain (COIPROCESS* process, size_t bytes, size_t depth) 
        {
          curBuffer = 0;
          numBuffers = depth;
          numBytes = bytes;
          
          buffer = new COIBUFFER[numBuffers];
          flag = new COI_ACCESS_FLAGS[numBuffers];
          mapInst = new COIMAPINSTANCE[numBuffers];
          for (size_t i=0; i<numBuffers; i++) {
            //COIRESULT result = COIBufferCreate(max(bytes,size_t(1)),COI_BUFFER_STREAMING_TO_SOURCE,0,NULL,1,process,&buffer[i]);
            COIRESULT result = COIBufferCreate(max(bytes,size_t(1)),COI_BUFFER_NORMAL,COI_OPTIMIZE_SOURCE_READ|COI_OPTIMIZE_SINK_WRITE|COI_OPTIMIZE_HUGE_PAGE_SIZE,NULL,1,process,&buffer[i]);
            if (result != COI_SUCCESS) throw std::runtime_error("COIBufferCreate failed: " + std::string(COIResultGetName(result)));
          flag[i] = COI_SINK_WRITE_ENTIRE;
          }
        }
        
        ~SwapChain () 
        {
          for (size_t i=0; i<numBuffers; i++) {
            COIRESULT result = COIBufferDestroy(buffer[i]);
            if (result != COI_SUCCESS) throw std::runtime_error("COIBufferDestroy failed: "+std::string(COIResultGetName(result)));
          }
          delete[] buffer;
          delete[] flag;
          delete[] mapInst;
        }
        
        void swapBuffers () 
        {
          curBuffer++;
          if (curBuffer >= numBuffers) curBuffer = 0;
        }
        
        /* map the framebuffer */
        void* map (int bufID) 
        {
          void* ptr = NULL;
          COIRESULT result = COIBufferMap(buffer[bufID],0,numBytes,COI_MAP_READ_ONLY,0,NULL,NULL,&mapInst[bufID],&ptr);
          if (result != COI_SUCCESS)
            throw std::runtime_error("COIBufferMap failed: "+std::string(COIResultGetName(result)));
          return ptr;
        }
        
        /* unmap the framebuffer */
        void unmap (int bufID)
        {
          COIRESULT result = COIBufferUnmap(mapInst[bufID],0,NULL,NULL);
          if (result != COI_SUCCESS)
            throw std::runtime_error("COIBufferUnmap failed: "+std::string(COIResultGetName(result)));
        }
        
      public:
        size_t curBuffer;
        size_t numBuffers;
        size_t numBytes;
        COIBUFFER* buffer;
        COI_ACCESS_FLAGS* flag;
        COIMAPINSTANCE* mapInst;
      };
      
    public:
      
      /*******************************************************************
                         creation of objects
      *******************************************************************/
      
      void rtNewCamera(RTCamera id, const char* type);
      void rtNewData(RTData id, const char* type, size_t bytes, const void* data);
      void rtNewImage(RTImage id, const char* type, size_t width, size_t height, const void* data);
      void rtNewTexture(RTTexture id, const char* type);
      void rtNewMaterial(RTMaterial id, const char* type);
      void rtNewShape(RTShape id, const char* type);
      void rtNewLight(RTLight id, const char* type);
      void rtNewShapePrimitive(RTPrimitive id, RTShape shape, RTMaterial material, const float* transform = NULL);
      void rtNewLightPrimitive(RTPrimitive id, RTLight light, RTMaterial material, const float* transform = NULL);
      void rtTransformPrimitive(RTPrimitive id, RTPrimitive prim, const float* transform);
      void rtNewScene(RTScene id, const char* type);
      void rtSetPrimitive(RTScene scene, size_t slot, RTPrimitive prim);
      void rtNewToneMapper(RTToneMapper id, const char* type);
      void rtNewRenderer(RTRenderer id, const char* type);
      void rtNewFrameBuffer(RTFrameBuffer id, const char* type, size_t width, size_t height, size_t buffers);
      void* rtMapFrameBuffer(RTFrameBuffer frameBuffer, int bufID);
      void rtUnmapFrameBuffer(RTFrameBuffer frameBuffer, int bufID);
      void rtSwapBuffers(RTFrameBuffer frameBuffer);
      void rtIncRef(RTHandle handle);
      void rtDecRef(RTHandle handle);
      
      /*******************************************************************
                            setting of parameters
      *******************************************************************/
      
      void rtSetBool1(RTHandle handle, const char* property, bool x);
      void rtSetBool2(RTHandle handle, const char* property, bool x, bool y);
      void rtSetBool3(RTHandle handle, const char* property, bool x, bool y, bool z);
      void rtSetBool4(RTHandle handle, const char* property, bool x, bool y, bool z, bool w);
      void rtSetInt1(RTHandle handle, const char* property, int x);
      void rtSetInt2(RTHandle handle, const char* property, int x, int y);
      void rtSetInt3(RTHandle handle, const char* property, int x, int y, int z);
      void rtSetInt4(RTHandle handle, const char* property, int x, int y, int z, int w);
      void rtSetFloat1(RTHandle handle, const char* property, float x);
      void rtSetFloat2(RTHandle handle, const char* property, float x, float y);
      void rtSetFloat3(RTHandle handle, const char* property, float x, float y, float z);
      void rtSetFloat4(RTHandle handle, const char* property, float x, float y, float z, float w);
      void rtSetArray(RTHandle handle, const char* property, const char* type, RTData data, size_t size, size_t stride, size_t ofs);
      void rtSetString(RTHandle handle, const char* property, const char* str);
      void rtSetImage(RTHandle handle, const char* property, RTImage img);
      void rtSetTexture(RTHandle handle, const char* property, RTTexture tex);
      void rtSetTransform(RTHandle handle, const char* property, const float* transform);
      void rtClear(RTHandle handle);
      void rtCommit(RTHandle handle);
      
      /*******************************************************************
                            render calls
      *******************************************************************/
      
      void rtRenderFrame(RTRenderer renderer, RTCamera camera, RTScene scene, RTToneMapper toneMapper, RTFrameBuffer frameBuffer, int accumulate);
      bool rtPick(RTCamera camera, float x, float y, RTScene scene, float& px, float& py, float& pz);
    
      void free(int id);
      
    private:
      COIPROCESS  process;
      COIENGINE   engine;
      COIPIPELINE pipeline;
      COIBUFFER   stream;
      std::vector<COILIBRARY> libs;
      std::map<int, Ref<SwapChain> > swapchains;

    private:
      COIFUNCTION runNewCamera;
      COIFUNCTION runNewData;
      COIFUNCTION runNewImage;
      COIFUNCTION runNewTexture;
      COIFUNCTION runNewMaterial;
      COIFUNCTION runNewShape;
      COIFUNCTION runNewLight;
      COIFUNCTION runNewShapePrimitive;
      COIFUNCTION runNewLightPrimitive;
      COIFUNCTION runTransformPrimitive;
      COIFUNCTION runNewScene;
      COIFUNCTION runSetPrimitive;
      COIFUNCTION runNewToneMapper;
      COIFUNCTION runNewRenderer;
      COIFUNCTION runNewFrameBuffer;
      COIFUNCTION runSwapBuffers;
      COIFUNCTION runIncRef;
      COIFUNCTION runDecRef;
      COIFUNCTION runSetBool1;
      COIFUNCTION runSetBool2;
      COIFUNCTION runSetBool3;
      COIFUNCTION runSetBool4;
      COIFUNCTION runSetInt1;
      COIFUNCTION runSetInt2;
      COIFUNCTION runSetInt3;
      COIFUNCTION runSetInt4;
      COIFUNCTION runSetFloat1;
      COIFUNCTION runSetFloat2;
      COIFUNCTION runSetFloat3;
      COIFUNCTION runSetFloat4;
      COIFUNCTION runSetArray;
      COIFUNCTION runSetString;
      COIFUNCTION runSetImage;
      COIFUNCTION runSetTexture;
      COIFUNCTION runSetTransform;
      COIFUNCTION runClear;
      COIFUNCTION runCommit;
      COIFUNCTION runRenderFrame;
      COIFUNCTION runPick;
      COIFUNCTION runNewDataStart;
      COIFUNCTION runNewDataSet;
      COIFUNCTION runNewDataEnd;
    };
    
  public:
    COIDevice(const char* executable, size_t numThreads, size_t verbose);
    ~COIDevice();
  
    /*******************************************************************
                         creation of objects
    *******************************************************************/

    RTCamera rtNewCamera(const char* type);
    RTData rtNewData(const char* type, size_t bytes, const void* data);
    RTData rtNewDataFromFile(const char* type, const char* file, size_t offset, size_t bytes);
    RTImage rtNewImage(const char* type, size_t width, size_t height, const void* data, const bool copy);
    RTImage rtNewImageFromFile(const char* file);
    RTTexture rtNewTexture(const char* type);
    RTMaterial rtNewMaterial(const char* type);
    RTShape rtNewShape(const char* type);
    RTLight rtNewLight(const char* type);
    RTPrimitive rtNewShapePrimitive(RTShape shape, RTMaterial material, const float* transform = NULL);
    RTPrimitive rtNewLightPrimitive(RTLight light, RTMaterial material, const float* transform = NULL);
    RTPrimitive rtTransformPrimitive(RTPrimitive prim, const float* transform);
    RTScene rtNewScene(const char* type);
    void rtSetPrimitive(RTScene scene, size_t slot, RTPrimitive prim);
    RTToneMapper rtNewToneMapper(const char* type);
    RTRenderer rtNewRenderer(const char* type);
    RTFrameBuffer rtNewFrameBuffer(const char* type, size_t width, size_t height, size_t buffers, void** ptrs);
    void* rtMapFrameBuffer(RTFrameBuffer frameBuffer, int bufID);
    void rtUnmapFrameBuffer(RTFrameBuffer frameBuffer, int bufID);
    void rtSwapBuffers(RTFrameBuffer frameBuffer);
    void rtIncRef(RTHandle handle);
    void rtDecRef(RTHandle handle);
    
    /*******************************************************************
                            setting of parameters
    *******************************************************************/

    void rtSetBool1(RTHandle handle, const char* property, bool x);
    void rtSetBool2(RTHandle handle, const char* property, bool x, bool y);
    void rtSetBool3(RTHandle handle, const char* property, bool x, bool y, bool z);
    void rtSetBool4(RTHandle handle, const char* property, bool x, bool y, bool z, bool w);
    void rtSetInt1(RTHandle handle, const char* property, int x);
    void rtSetInt2(RTHandle handle, const char* property, int x, int y);
    void rtSetInt3(RTHandle handle, const char* property, int x, int y, int z);
    void rtSetInt4(RTHandle handle, const char* property, int x, int y, int z, int w);
    void rtSetFloat1(RTHandle handle, const char* property, float x);
    void rtSetFloat2(RTHandle handle, const char* property, float x, float y);
    void rtSetFloat3(RTHandle handle, const char* property, float x, float y, float z);
    void rtSetFloat4(RTHandle handle, const char* property, float x, float y, float z, float w);
    void rtSetArray(RTHandle handle, const char* property, const char* type, RTData data, size_t size, size_t stride, size_t ofs);
    void rtSetString(RTHandle handle, const char* property, const char* str);
    void rtSetImage(RTHandle handle, const char* property, RTImage img);
    void rtSetTexture(RTHandle handle, const char* property, RTTexture tex);
    void rtSetTransform(RTHandle handle, const char* property, const float* transform);
    void rtClear(RTHandle handle);
    void rtCommit(RTHandle handle);
    
    /*******************************************************************
                            render calls
    *******************************************************************/
    
    void rtRenderFrame(RTRenderer renderer, RTCamera camera, RTScene scene, RTToneMapper toneMapper, RTFrameBuffer frameBuffer, int accumulate);
    bool rtPick(RTCamera camera, float x, float y, RTScene scene, float& px, float& py, float& pz);

    /*! handle management */
  private:
    MutexSys handleMutex;
    int nextHandle;                         //!< next ID to take if pool is empty
    std::vector<int> pool;                  //!< pool of handle IDs
    std::vector<int> counters;              //!< reference counters of handles
    std::vector<Ref<SwapChain> > buffers;   //!< local framebuffer representation

    /*! allocate a new handle ID */
    int allocHandle();

    /*! increment local reference count */
    void incRef(int id);
    
    /*! decrement local reference count */
    bool decRef(int id);

    /* multi card rendering */
  private:
    int serverID;
    int serverCount;
     
  private:
    std::vector<COIProcess*> devices;
  };
}

#endif
