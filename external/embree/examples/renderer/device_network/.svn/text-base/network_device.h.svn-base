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

#ifndef __EMBREE_NETWORK_DEVICE_H__
#define __EMBREE_NETWORK_DEVICE_H__

#include "../device/device.h"
#include "../device_singleray/api/swapchain.h"
#include "network_common.h"

namespace embree
{
  class NetworkDevice : public Device
  {
  public:

    /*! construction */
    NetworkDevice();

    /*! construction */
    NetworkDevice(const std::vector<network::socket_t>& server);

    /*! destruction */
    ~NetworkDevice();

    void init ();

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
    RTFrameBuffer rtNewFrameBuffer(const char* type, size_t width, size_t height, size_t buffers, void** ptr);
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

    /*******************************************************************
                        network broadcast to servers
    *******************************************************************/

    void broadcast(bool value);
    void broadcast(char value);
    void broadcast(int value);
    void broadcast(float value);
    void broadcast(const std::string& value);
    void broadcast(const Vector3f& value);
    void broadcast(const void* data, size_t bytes);
    void flush();

  private:

    /*! read a single command from the socket */
    void receiveCommand(int id);

    /*! receive thread */
    bool terminate;
    static void _receive(void* data);
    void receive();
    Barrier receiveBarrier;

    /*! connection to rendering servers */
  protected:
    Atomic serverID;                         //!< to assign a unique server to each thread
    std::vector<thread_t> threads;           //!< each thread receives data of one rendering server
    std::vector<network::socket_t> servers;  //!< sockets of the rendering servers

    /*! handle management */
  private:
    MutexSys handleMutex;
    int nextHandle;                                      //!< next ID to take if pool is empty
    std::vector<int> pool;                               //!< pool of handle IDs
    std::vector<int> counters;                           //!< reference counters of handles
    std::vector<Ref<SwapChain> > buffers;  //!< local framebuffer representation

    /*! allocate a new handle ID */
    int allocHandle();

    /*! increment local reference count */
    void incRef(int id);
    
    /*! decrement local reference count */
    void decRef(int id);
    
    /*! return data for pick command */
  private:
    struct {
      Barrier barrier;                       //!< to wait for pick results
      bool hit;                              //!< if the ray hit something
      Vector3f pos;                             //!< location in space that got hit
    } pick;
  };
}

#endif
