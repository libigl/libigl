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

#include "coi_server.h"
#include "sys/sysinfo.h"

#include <sink/COIProcess_sink.h>
#include <sink/COIPipeline_sink.h>

#include <queue>
#include <vector>
#include <stdint.h>

extern "C" embree::Device* create(const char* parms, size_t numThreads, size_t verbose);
 
namespace embree 
{
  Device* g_device = NULL;
  int g_verbose = 0;

  /*! handle management =========================== */
  std::vector<Device::RTHandle> g_handles;
  std::vector<Ref<SwapChain> > g_swapChains;
  std::vector<int> g_numRefs;
  
  /*! return the handle associated with this ID */
  template<typename T> T get(size_t id) {
    if (id == 0) return (T) NULL;
    return((T) g_handles[id]);
  }

  /*! update ID to handle mapping */
  void set(int id, Device::RTHandle handle, SwapChain* swapChain = NULL) 
  {
    if (g_handles.size() <= (size_t) id) { 
      g_handles.resize(id + 1);  
      g_swapChains.resize(id + 1);
      g_numRefs.resize(id + 1); 
    }
    g_handles[id] = handle;
    g_swapChains[id] = swapChain;
    g_numRefs[id] = 1;
  }

  extern "C" void rtNewCamera(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewCamera* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewCamera(%s)\n", parms->id, parms->type);
      fflush(stdout);
    }
    set(parms->id, g_device->rtNewCamera(parms->type));
  }

  extern "C" void rtNewData(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewData* parms, uint16_t parmBytes, void* ret, uint16_t retBytess)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewData(%zu)\n", parms->id, parms->bytes);
      fflush(stdout);
    }
    set(parms->id, g_device->rtNewData("immutable",parms->bytes,buffers[0]));
  }

  char* g_data = NULL;

  extern "C" void rtNewDataStart(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewDataStart* parms, uint16_t parmBytes, void* ret, uint16_t retBytess)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewDataStart(%zu)\n", parms->id, parms->bytes);
      fflush(stdout);
    }
    g_data = (char*) alignedMalloc(parms->bytes);
  }

  extern "C" void rtNewDataSet(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewDataSet* parms, uint16_t parmBytes, void* ret, uint16_t retBytess)
  {
    if (g_verbose) {
      printf("  rtNewDataSet(%zu,%zu)\n", parms->offset, parms->bytes);
      fflush(stdout);
    }
    memcpy(g_data+parms->offset,buffers[0],parms->bytes);
  }

  extern "C" void rtNewDataEnd(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewDataStart* parms, uint16_t parmBytes, void* ret, uint16_t retBytess)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewDataEnd(%zu)\n", parms->id, parms->bytes);
      fflush(stdout);
    }
    set(parms->id, g_device->rtNewData("immutable_managed",parms->bytes,g_data));
    g_data = NULL;
  }

  extern "C" void rtNewImage(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewImage* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewImage(%s, %zu, %zu)\n", parms->id, parms->type, parms->width, parms->height);
      fflush(stdout);
    }
    set(parms->id, g_device->rtNewImage(parms->type,parms->width,parms->height,buffers[0],true));
  }

  extern "C" void rtNewTexture(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewTexture* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewTexture(%s)\n", parms->id, parms->type);
      fflush(stdout);
    }
    set(parms->id, g_device->rtNewTexture(parms->type));
  }

  extern "C" void rtNewMaterial(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewMaterial* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewMaterial(%s)\n", parms->id, parms->type);
      fflush(stdout);
    }
    set(parms->id, g_device->rtNewMaterial(parms->type));
  }

  extern "C" void rtNewShape(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewShape* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewShape(%s)\n", parms->id, parms->type);
      fflush(stdout);
    }
    set(parms->id, g_device->rtNewShape(parms->type));
  }

  extern "C" void rtNewLight(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewLight* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) printf("handle %06d = rtNewLight(%s)\n", parms->id, parms->type);
    set(parms->id, g_device->rtNewLight(parms->type));
  }

  extern "C" void rtNewShapePrimitive(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, 
                                      parmsNewShapePrimitive* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) { 
      printf("handle %06d = rtNewShapePrimitive(%06d, %06d, %s)\n", 
             parms->id, parms->shape, parms->material, std::stringOf(copyFromArray(parms->transform)).c_str());
      fflush(stdout);
    }
    set(parms->id, g_device->rtNewShapePrimitive(get<Device::RTShape>(parms->shape),get<Device::RTMaterial>(parms->material),parms->transform));
  }

  extern "C" void rtNewLightPrimitive(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, 
                                      parmsNewLightPrimitive* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewLightPrimitive(%06d, %06d, %s)\n", 
             parms->id, parms->light, parms->material, std::stringOf(copyFromArray(parms->transform)).c_str());
      fflush(stdout);
    }
    set(parms->id, g_device->rtNewLightPrimitive(get<Device::RTLight>(parms->light),get<Device::RTMaterial>(parms->material),parms->transform));
  }

  extern "C" void rtTransformPrimitive(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, 
                                       parmsTransformPrimitive* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("handle %06d = rtTransformPrimitive(%06d,%s)\n", 
             parms->id, parms->primitive, std::stringOf(copyFromArray(parms->transform)).c_str());
      fflush(stdout);
    }
    set(parms->id, g_device->rtTransformPrimitive(get<Device::RTPrimitive>(parms->primitive),parms->transform));
  }

  extern "C" void rtNewScene(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewScene* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewScene(%s)\n", parms->id, parms->type);
      fflush(stdout);
    }
    set(parms->id, g_device->rtNewScene(parms->type));
  }

  extern "C" void rtSetPrimitive(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetPrimitive* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetPrimitive(%06d, %06d, %06d)\n", parms->scene, parms->slot, parms->prim);
      fflush(stdout);
    }
    g_device->rtSetPrimitive(get<Device::RTScene>(parms->scene),parms->slot,get<Device::RTPrimitive>(parms->prim));
  }

  extern "C" void rtNewToneMapper(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewToneMapper* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewToneMapper(%s)\n", parms->id, parms->type);
      fflush(stdout);
    }
    set(parms->id, g_device->rtNewToneMapper(parms->type));
  }

  extern "C" void rtNewRenderer(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewRenderer* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) printf("handle %06d = rtNewRenderer(%s)\n", parms->id, parms->type);
    set(parms->id, g_device->rtNewRenderer(parms->type));
  }

  extern "C" void rtNewFrameBuffer(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsNewFrameBuffer* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("handle %06d = rtNewFrameBuffer(%s, %zu, %zu, %zu)\n", parms->id, parms->type, parms->width, parms->height, parms->depth);
      fflush(stdout);
    }                                         
    set(parms->id, g_device->rtNewFrameBuffer(parms->type, parms->width, parms->height, parms->depth, buffers));
  }

  extern "C" void rtSwapBuffers(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSwapBuffers* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSwapBuffers(%06d)\n", parms->framebuffer);
      fflush(stdout);
    }
    g_device->rtSwapBuffers(get<Device::RTFrameBuffer>(parms->framebuffer));
  }

  extern "C" void rtIncRef(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsIncRef* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtIncRef(%06d)\n", parms->handle);
      fflush(stdout);
    }
    g_device->rtIncRef(get<Device::RTHandle>(parms->handle));
  }

  extern "C" void rtDecRef(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsDecRef* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtDecRef(%06d)\n", parms->handle);
      fflush(stdout);
    }
    g_device->rtDecRef(get<Device::RTHandle>(parms->handle));
  }

  extern "C" void rtSetBool1(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetBoolN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetBool1(%06d, %s, [%d])\n", parms->handle, parms->property, parms->x);
      fflush(stdout);
    }
    g_device->rtSetBool1(get<Device::RTHandle>(parms->handle),parms->property,parms->x);
  }

  extern "C" void rtSetBool2(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetBoolN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetBool2(%06d, %s, [%d, %d])\n", parms->handle, parms->property, parms->x, parms->y);
      fflush(stdout);
    }
    g_device->rtSetBool2(get<Device::RTHandle>(parms->handle),parms->property,parms->x,parms->y);
  }

  extern "C" void rtSetBool3(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetBoolN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetBool3(%06d, %s, [%d, %d, %d])\n", parms->handle, parms->property, parms->x, parms->y, parms->z);
      fflush(stdout);
    }
    g_device->rtSetBool3(get<Device::RTHandle>(parms->handle),parms->property,parms->x,parms->y,parms->z);
  }

  extern "C" void rtSetBool4(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetBoolN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetBool4(%06d, %s, [%d, %d, %d, %d])\n", parms->handle, parms->property, parms->x, parms->y, parms->z, parms->w);
      fflush(stdout);
    }
    g_device->rtSetBool4(get<Device::RTHandle>(parms->handle),parms->property,parms->x,parms->y,parms->z,parms->w);
  }

  extern "C" void rtSetInt1(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetIntN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetInt1(%06d, %s, [%d])\n", parms->handle, parms->property, parms->x);
      fflush(stdout);
    }
    g_device->rtSetInt1(get<Device::RTHandle>(parms->handle),parms->property,parms->x);
  }

  extern "C" void rtSetInt2(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetIntN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) printf("rtSetInt2(%06d, %s, [%d, %d])\n", parms->handle, parms->property, parms->x, parms->y);
    g_device->rtSetInt2(get<Device::RTHandle>(parms->handle),parms->property,parms->x,parms->y);
  }

  extern "C" void rtSetInt3(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetIntN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetInt3(%06d, %s, [%d, %d, %d])\n", parms->handle, parms->property, parms->x, parms->y, parms->z);
      fflush(stdout);
    }
    g_device->rtSetInt3(get<Device::RTHandle>(parms->handle),parms->property,parms->x,parms->y,parms->z);
  }

  extern "C" void rtSetInt4(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetIntN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetInt4(%06d, %s, [%d, %d, %d, %d])\n", parms->handle, parms->property, parms->x, parms->y, parms->z, parms->w);
      fflush(stdout);
    }
    g_device->rtSetInt4(get<Device::RTHandle>(parms->handle),parms->property,parms->x,parms->y,parms->z,parms->w);
  }

  extern "C" void rtSetFloat1(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetFloatN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetFloat1(%06d, %s, [%f])\n", parms->handle, parms->property, parms->x);
      fflush(stdout);
    }
    g_device->rtSetFloat1(get<Device::RTHandle>(parms->handle),parms->property,parms->x);
  }

  extern "C" void rtSetFloat2(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetFloatN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetFloat2(%06d, %s, [%f, %f])\n", parms->handle, parms->property, parms->x, parms->y);
      fflush(stdout);
    }
    g_device->rtSetFloat2(get<Device::RTHandle>(parms->handle),parms->property,parms->x,parms->y);
  }

  extern "C" void rtSetFloat3(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetFloatN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetFloat3(%06d, %s, [%f, %f, %f])\n", parms->handle, parms->property, parms->x, parms->y, parms->z);
      fflush(stdout);
    }
    g_device->rtSetFloat3(get<Device::RTHandle>(parms->handle),parms->property,parms->x,parms->y,parms->z);
  }

  extern "C" void rtSetFloat4(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetFloatN* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetFloat4(%06d, %s, [%f, %f, %f, %f])\n", parms->handle, parms->property, parms->x, parms->y, parms->z, parms->w);
      fflush(stdout);
    }
    g_device->rtSetFloat4(get<Device::RTHandle>(parms->handle),parms->property,parms->x,parms->y,parms->z,parms->w);
  }

  extern "C" void rtSetArray(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetArray* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetArray(%06d, %s, %s, %06d, %zu, %zu, %zu)\n", parms->handle, parms->property, parms->type, parms->data, parms->size, parms->stride, parms->ofs);
      fflush(stdout);
    }
    g_device->rtSetArray(get<Device::RTHandle>(parms->handle),parms->property,parms->type,get<Device::RTData>(parms->data),parms->size,parms->stride,parms->ofs);
  }

  extern "C" void rtSetString(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetString* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetString(%06d, %s, %s)\n", parms->handle, parms->property, parms->str);
      fflush(stdout);
    }
    g_device->rtSetString(get<Device::RTHandle>(parms->handle),parms->property,parms->str);
  }

  extern "C" void rtSetImage(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetImage* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {\
      printf("rtSetImage(%06d, %s, %06d)\n", parms->handle, parms->property, parms->image);
      fflush(stdout);
    }
    g_device->rtSetImage(get<Device::RTHandle>(parms->handle),parms->property,get<Device::RTImage>(parms->image));
  }

  extern "C" void rtSetTexture(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetTexture* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetTexture(%06d, %s, %06d)\n", parms->handle, parms->property, parms->texture);
      fflush(stdout);
    }
    g_device->rtSetTexture(get<Device::RTHandle>(parms->handle),parms->property,get<Device::RTTexture>(parms->texture));
  }

  extern "C" void rtSetTransform(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsSetTransform* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtSetTransform(%06d, %s, %s)\n", parms->handle, parms->property, std::stringOf(copyFromArray(parms->transform)).c_str());
      fflush(stdout);
    }
    g_device->rtSetTransform(get<Device::RTHandle>(parms->handle),parms->property,parms->transform);
  }

  extern "C" void rtClear(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsClear* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtClear(%06d)\n", parms->handle);
      fflush(stdout);
    }
    g_device->rtClear(get<Device::RTHandle>(parms->handle));
  }

  extern "C" void rtCommit(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsCommit* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtCommit(%06d)\n", parms->handle);
      fflush(stdout);
    }
    g_device->rtCommit(get<Device::RTHandle>(parms->handle));
  }

  extern "C" void rtRenderFrame(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsRenderFrame* parms, uint16_t parmBytes, void* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtRenderFrame(%06d, %06d, %06d, %06d, %06d, %d)\n", 
             parms->renderer, parms->camera, parms->scene, parms->toneMapper, parms->frameBuffer, parms->accumulate);
      fflush(stdout);
    }

    g_device->rtRenderFrame(get<Device::RTRenderer>(parms->renderer),
                            get<Device::RTCamera>(parms->camera),
                            get<Device::RTScene>(parms->scene),
                            get<Device::RTToneMapper>(parms->toneMapper),
                            get<Device::RTFrameBuffer>(parms->frameBuffer),
                            parms->accumulate);
  }

  extern "C" void rtPick(uint32_t numBuffers, void** buffers, uint64_t* bufferBytes, parmsPick* parms, uint16_t parmBytes, returnPick* ret, uint16_t retBytes)
  {
    if (g_verbose) {
      printf("rtPick([%f, %f], %06d, %06d)\n", parms->x, parms->y, parms->camera, parms->scene);
      fflush(stdout);
    }

    float px, py, pz;
    bool isHit = g_device->rtPick(get<Device::RTCamera>(parms->camera), parms->x, parms->y, get<Device::RTScene>(parms->scene), px, py, pz);
    if (g_verbose) {
      printf("  returning %s, [%f, %f, %f]\n", (isHit ? "TRUE" : "FALSE"), px, py, pz);
      fflush(stdout);
    }

    ret->x = px;
    ret->y = py;
    ret->z = pz;
    ret->hit = isHit;
  }

  void main(int argc, char **argv) 
  {
    /* get number of threads to use */
    size_t numThreads = 0;
    if (argc >= 2) numThreads = atoi(argv[1]);
    if (argc >= 3) g_verbose = atoi(argv[2]);
    
#if defined(__MIC__)
    if (numThreads == 0) numThreads = getNumberOfLogicalThreads()-4;
#else
    if (numThreads == 0) numThreads = getNumberOfLogicalThreads();
#endif

    /* enable wait to attach with debugger */
#if 0
    std::cout << "waiting " << std::flush;
    for (int i=0; i<20; i++) {
      sleep(1);
      std::cout << "." << std::flush;
    }
    std::cout << " [DONE]" << std::endl;
#endif
    
    /*! create device */
    g_device = create("",numThreads,g_verbose);
  }
}

int main(int argc, char **argv) try 
{
  embree::main(argc, argv);

  // Functions enqueued on the sink side will not start executing until
  // you call COIPipelineStartExecutingRunFunctions(). This call is to
  // synchronize any initialization required on the sink side
  COIRESULT result = COIPipelineStartExecutingRunFunctions();
  if (result != COI_SUCCESS) 
    throw std::runtime_error("COIPipelineStartExecutingRunFunctions failed: "+std::string(COIResultGetName(result)));
  
  // This call will wait until COIProcessDestroy() gets called on the source
  // side. If COIProcessDestroy is called without force flag set, this call
  // will make sure all the functions enqueued are executed and does all
  // clean up required to exit gracefully.
  COIProcessWaitForShutdown();
  return 0;
  
 } catch (const std::exception &e) {
  std::cout << "main(): Error: " << e.what() << std::endl;
  return(1);

 } catch (...) {
  std::cout << "main(): Error: unknown exception caught." << std::endl;
  return(1);
 }
