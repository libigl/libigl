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

#ifndef __EMBREE_ISPC_RENDER_DEVICE_H__
#define __EMBREE_ISPC_RENDER_DEVICE_H__

#include "../default.h"
#include "device/device.h"
#include "sys/sync/mutex.h"

namespace embree
{
  class ISPCDevice : public Device
  {
  public:

    /*! construction */
    ISPCDevice(size_t numThreads, size_t verbose);

    /*! destruction */
    ~ISPCDevice();

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
    RTPrimitive rtNewShapePrimitive(RTShape shape, RTMaterial material, const float* transform);
    RTPrimitive rtNewLightPrimitive(RTLight light, RTMaterial material, const float* transform);
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

  private:
    MutexSys mutex;
  };
}

#endif

