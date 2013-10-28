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

#include "singleray_device.h"
#include "image/image.h"
#include "sys/taskscheduler.h"

/* include general stuff */
#include "api/handle.h"
#include "api/data.h"
#include "api/parms.h"
#include "sys/sync/mutex.h"

/* include all scenes */
#include "api/instance.h"
#include "api/scene.h"
#include "api/scene_flat.h"
#include "api/scene_instancing.h"

/* include all cameras */
#include "cameras/pinholecamera.h"
#include "cameras/depthoffieldcamera.h"

/* include all lights */
#include "lights/ambientlight.h"
#include "lights/pointlight.h"
#include "lights/spotlight.h"
#include "lights/directionallight.h"
#include "lights/distantlight.h"
#include "lights/hdrilight.h"
#include "lights/trianglelight.h"

/* include all materials */
#include "materials/matte.h"
#include "materials/plastic.h"
#include "materials/dielectric.h"
#include "materials/thindielectric.h"
#include "materials/mirror.h"
#include "materials/metal.h"
#include "materials/brushedmetal.h"
#include "materials/metallicpaint.h"
#include "materials/matte_textured.h"
#include "materials/obj.h"
#include "materials/velvet.h"

/* include all shapes */
#include "shapes/triangle.h"
#include "shapes/trianglemesh.h"
#include "shapes/sphere.h"
#include "shapes/disk.h"

/* include all textures */
#include "textures/nearestneighbor.h"

/* include all tonemappers */
#include "tonemappers/defaulttonemapper.h"

/* include all renderers */
#include "renderers/debugrenderer.h"
#include "renderers/integratorrenderer.h"

/* include ray tracing core interface */
#include "embree/include/embree.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#define strcasecmp lstrcmpiA
#pragma warning(disable:4297) // function assumed not to throw an exception but does
#endif

#define RT_COMMAND_HEADER Lock<MutexSys> lock(mutex); g_time++

namespace embree
{
  /*******************************************************************
                  creation of device
  *******************************************************************/

  __dllexport Device* create(const char* parms, size_t numThreads, size_t verbose) 
  {
    return new SingleRayDevice(numThreads,verbose);
  }

  int g_serverCount = 1;
  int g_serverID = 0;
  size_t g_time = 0;

  /*******************************************************************
                  type definitions
  *******************************************************************/

  /* camera handle */
  struct _RTCamera : public _RTHandle { };

  /* data handle */
  struct _RTData : public _RTHandle { };

  /* image handle */
  struct _RTImage : public _RTHandle { };

  /* texture handle */
  struct _RTTexture : public _RTHandle { };

  /* material handle */
  struct _RTMaterial : public _RTHandle { };

  /* shape handle */
  struct _RTShape : public _RTHandle { };

  /* light handle */
  struct _RTLight : public _RTHandle { };

  /* primitive handle */
  struct _RTPrimitive : public _RTHandle { };

  /* scene handle */
  struct _RTScene : public _RTHandle { };

  /* renderer handle */
  struct _RTRenderer : public _RTHandle { };

  /* framebuffer handle */
  struct _RTFrameBuffer : public _RTHandle { };


  /*******************************************************************
                  construction
  *******************************************************************/
  
  SingleRayDevice::SingleRayDevice(size_t numThreads, size_t verbose)
  {
    rtcInit();
    rtcSetVerbose(verbose);
    rtcStartThreads(numThreads);
  }

  SingleRayDevice::~SingleRayDevice() 
  {
    rtcStopThreads();
    rtcExit();
  }

  /*******************************************************************
                    creation of handles
  *******************************************************************/

  Device::RTCamera SingleRayDevice::rtNewCamera(const char* type)
  {
    RT_COMMAND_HEADER;
    if      (!strcasecmp(type,"pinhole")) return (Device::RTCamera) new ConstructorHandle<PinHoleCamera,Camera>;
    else if (!strcasecmp(type,"depthoffield")) return (Device::RTCamera) new ConstructorHandle<DepthOfFieldCamera,Camera>;
    else throw std::runtime_error("unknown camera type: "+std::string(type));
  }

  Device::RTData SingleRayDevice::rtNewData(const char* type, size_t bytes, const void* data)
  {
    RT_COMMAND_HEADER;
    if (!strcasecmp(type,"immutable")) 
      return (Device::RTData) new ConstHandle<Data>(new Data(bytes,data,true));
    else if (!strcasecmp(type,"immutable_managed")) 
      return (Device::RTData) new ConstHandle<Data>(new Data(bytes,data,false));
    else
      throw std::runtime_error("unknown data buffer type: "+std::string(type));
  }

  Device::RTData SingleRayDevice::rtNewDataFromFile(const char* type, const char* fileName, size_t offset, size_t bytes)
  {
    RT_COMMAND_HEADER;

    /*! we always load locally */
    if (!strncmp(fileName,"server:",7)) 
      fileName += 7;

    if (!strcasecmp(type,"immutable")) 
    {
      FILE* file = fopen(fileName,"rb");
      if (!file) throw std::runtime_error("cannot open file "+(std::string)fileName);
      
      Data* data = new Data(bytes);
      fseek(file,(long)offset,SEEK_SET);
      if (bytes != fread(data->map(),1,bytes,file))
        throw std::runtime_error("error filling data buffer from file");
      fclose(file);

      return (Device::RTData) new ConstHandle<Data>(data);
    }
    else
      throw std::runtime_error("unknown data buffer type: "+std::string(type));
  }

  Device::RTImage SingleRayDevice::rtNewImage(const char* type, size_t width, size_t height, const void* data, const bool copy)
  {
    RT_COMMAND_HEADER;
    if      (!strcasecmp(type,"RGB8"        )) return (Device::RTImage) new ConstHandle<Image>(new Image3c(width,height,(Col3c*)data,copy));
    else if (!strcasecmp(type,"RGBA8"       )) return (Device::RTImage) new ConstHandle<Image>(new Image4c(width,height,(Col4c*)data,copy));
    else if (!strcasecmp(type,"RGB_FLOAT32" )) return (Device::RTImage) new ConstHandle<Image>(new Image3f(width,height,(Col3f*)data,copy));
    else if (!strcasecmp(type,"RGBA_FLOAT32")) return (Device::RTImage) new ConstHandle<Image>(new Image4f(width,height,(Col4f*)data,copy));
    else throw std::runtime_error("unknown image type: "+std::string(type));
  }

  Device::RTImage SingleRayDevice::rtNewImageFromFile(const char* file)
  {
    RT_COMMAND_HEADER;
#if defined(__MIC__)      
      throw std::runtime_error("rtNewImageFromFile not supported on MIC");
#else
    if (!strncmp(file,"server:",7)) file += 7;
    Ref<Image> image = loadImage(file);
    if (image) 
      return (Device::RTImage) new ConstHandle<Image>(image);
    else
      return (Device::RTImage) new ConstHandle<Image>(new Image3c(1,1,Col3c(255,255,255)));
#endif
  }

  Device::RTTexture SingleRayDevice::rtNewTexture(const char* type) {
    RT_COMMAND_HEADER;
    if (!strcasecmp(type,"nearest")) return (Device::RTTexture) new ConstructorHandle<NearestNeighbor,Texture>;
    else if (!strcasecmp(type,"image")) return (Device::RTTexture) new ConstructorHandle<NearestNeighbor,Texture>;
    else throw std::runtime_error("unsupported texture type: "+std::string(type));
  }

  Device::RTMaterial SingleRayDevice::rtNewMaterial(const char* type)
  {
    RT_COMMAND_HEADER;
    if      (!strcasecmp(type,"Matte")         ) return (Device::RTMaterial) new ConstructorHandle<Matte,Material>;
    else if (!strcasecmp(type,"Plastic")       ) return (Device::RTMaterial) new ConstructorHandle<Plastic,Material>;
    else if (!strcasecmp(type,"Dielectric")    ) return (Device::RTMaterial) new ConstructorHandle<Dielectric,Material>;
    else if (!strcasecmp(type,"Glass")         ) return (Device::RTMaterial) new ConstructorHandle<Dielectric,Material>;
    else if (!strcasecmp(type,"ThinDielectric")) return (Device::RTMaterial) new ConstructorHandle<ThinDielectric,Material>;
    else if (!strcasecmp(type,"ThinGlass")     ) return (Device::RTMaterial) new ConstructorHandle<ThinDielectric,Material>;
    else if (!strcasecmp(type,"Mirror")        ) return (Device::RTMaterial) new ConstructorHandle<Mirror,Material>;
    else if (!strcasecmp(type,"Metal")         ) return (Device::RTMaterial) new ConstructorHandle<Metal,Material>;
    else if (!strcasecmp(type,"BrushedMetal")  ) return (Device::RTMaterial) new ConstructorHandle<BrushedMetal,Material>;
    else if (!strcasecmp(type,"MetallicPaint") ) return (Device::RTMaterial) new ConstructorHandle<MetallicPaint,Material>;
    else if (!strcasecmp(type,"MatteTextured") ) return (Device::RTMaterial) new ConstructorHandle<MatteTextured,Material>;
    else if (!strcasecmp(type,"Obj")           ) return (Device::RTMaterial) new ConstructorHandle<Obj,Material>;
    else if (!strcasecmp(type,"Velvet")        ) return (Device::RTMaterial) new ConstructorHandle<Velvet,Material>;
    else throw std::runtime_error("unknown material type: "+std::string(type));
  }

  Device::RTShape SingleRayDevice::rtNewShape(const char* type) {
    RT_COMMAND_HEADER;
    if      (!strcasecmp(type,"trianglemesh")) return (Device::RTShape) new CreateHandle<TriangleMesh,Shape>;
    else if (!strcasecmp(type,"triangle")    ) return (Device::RTShape) new ConstructorHandle<Triangle,Shape>;
    else if (!strcasecmp(type,"sphere")      ) return (Device::RTShape) new ConstructorHandle<Sphere,Shape>;
    else if (!strcasecmp(type,"disk")        ) return (Device::RTShape) new ConstructorHandle<Disk,Shape>;
    else throw std::runtime_error("unknown shape type: "+std::string(type));
  }

  Device::RTLight SingleRayDevice::rtNewLight(const char* type)
  {
    RT_COMMAND_HEADER;
    if      (!strcasecmp(type,"ambientlight"    )) return (Device::RTLight) new ConstructorHandle<AmbientLight,Light>;
    else if (!strcasecmp(type,"pointlight"      )) return (Device::RTLight) new ConstructorHandle<PointLight,Light>;
    else if (!strcasecmp(type,"spotlight"       )) return (Device::RTLight) new ConstructorHandle<SpotLight,Light>;
    else if (!strcasecmp(type,"directionallight")) return (Device::RTLight) new ConstructorHandle<DirectionalLight,Light>;
    else if (!strcasecmp(type,"distantlight"    )) return (Device::RTLight) new ConstructorHandle<DistantLight,Light>;
    else if (!strcasecmp(type,"hdrilight"       )) return (Device::RTLight) new ConstructorHandle<HDRILight,Light>;
    else if (!strcasecmp(type,"trianglelight"   )) return (Device::RTLight) new ConstructorHandle<TriangleLight,Light>;
    else throw std::runtime_error("unknown light type: "+std::string(type));
  }

  Device::RTPrimitive SingleRayDevice::rtNewShapePrimitive(Device::RTShape shape_i, 
                                                           Device::RTMaterial material_i, 
                                                           const float* transform)
  {
    RT_COMMAND_HEADER;
    Ref<InstanceHandle<Shape> > shape = castHandle<InstanceHandle<Shape>    >(shape_i   ,"shape"   );
    Ref<InstanceHandle<Material> > material = castHandle<InstanceHandle<Material> >(material_i,"material");
    AffineSpace3f space = transform ? copyFromArray(transform) : AffineSpace3f(one);
    return (Device::RTPrimitive) new PrimitiveHandle(shape,material,space);
  }

  Device::RTPrimitive SingleRayDevice::rtNewLightPrimitive(Device::RTLight light_i, 
                                                        Device::RTMaterial material_i, 
                                                        const float* transform)
  {
    RT_COMMAND_HEADER;
    Ref<InstanceHandle<Light> > light = castHandle<InstanceHandle<Light> >(light_i,"light");
    Ref<InstanceHandle<Material> > material = NULL;
    if (material_i) material = castHandle<InstanceHandle<Material> >(material_i,"material");
    AffineSpace3f space = transform ? copyFromArray(transform) : AffineSpace3f(one);
    return (Device::RTPrimitive) new PrimitiveHandle(light,material,space);
  }

  Device::RTPrimitive SingleRayDevice::rtTransformPrimitive(Device::RTPrimitive primitive, const float* transform) 
  {
    RT_COMMAND_HEADER;
    Ref<PrimitiveHandle> prim = dynamic_cast<PrimitiveHandle*>((_RTHandle*)primitive);
    AffineSpace3f space = transform ? copyFromArray(transform) : AffineSpace3f(one);
    return (Device::RTPrimitive) new PrimitiveHandle(space, prim);
  }
  
  Device::RTScene SingleRayDevice::rtNewScene(const char* type) 
  {
    RT_COMMAND_HEADER;
    if      (!strcmp(type,"default" )) return (Device::RTScene) new BackendSceneFlat::Handle;
    else if (!strcmp(type,"flat"    )) return (Device::RTScene) new BackendSceneFlat::Handle;
    else if (!strcmp(type,"twolevel")) return (Device::RTScene) new BackendSceneInstancing::Handle;
    else throw std::runtime_error("unknown scene type: "+std::string(type));
  }
     
  void SingleRayDevice::rtSetPrimitive(RTScene hscene, size_t slot, RTPrimitive hprim) 
  {
    RT_COMMAND_HEADER;
    Ref<BackendScene::Handle> scene = castHandle<BackendScene::Handle>(hscene,"scene");
    if (hprim == NULL) { scene->setPrimitive(slot,NULL); return; }
    Ref<PrimitiveHandle> prim = dynamic_cast<PrimitiveHandle*>((_RTHandle*)hprim);
    scene->setPrimitive(slot,prim);
  }

  Device::RTToneMapper SingleRayDevice::rtNewToneMapper(const char* type)
  {
    RT_COMMAND_HEADER;
    if      (!strcasecmp(type,"default")) return (Device::RTToneMapper) new ConstructorHandle<DefaultToneMapper,ToneMapper>;
    else throw std::runtime_error("unknown tonemapper type: "+std::string(type));
  }

  Device::RTRenderer SingleRayDevice::rtNewRenderer(const char* type)
  {
    RT_COMMAND_HEADER;
    if (!strcasecmp(type,"debug"     )) return (Device::RTRenderer) new ConstructorHandle<DebugRenderer,Renderer>;
    if (!strcasecmp(type,"pathtracer")) {
      ConstructorHandle<IntegratorRenderer,Renderer>* handle = new ConstructorHandle<IntegratorRenderer,Renderer>;
      handle->set("integrator",Variant("pathtracer"));
      return (Device::RTRenderer) handle;
    }
    else throw std::runtime_error("unknown renderer type: " + std::string(type));
  }

  Device::RTFrameBuffer SingleRayDevice::rtNewFrameBuffer(const char* type, size_t width, size_t height, size_t buffers, void** ptrs) 
  {
    RT_COMMAND_HEADER;
    Ref<SwapChain> swapchain = null;
    if      (!strcasecmp(type,"RGB_FLOAT32")) swapchain = new SwapChain(type,width,height,buffers,ptrs,FrameBufferRGBFloat32::create);
    else if (!strcasecmp(type,"RGBA8"      )) swapchain = new SwapChain(type,width,height,buffers,ptrs,FrameBufferRGBA8     ::create);
    else if (!strcasecmp(type,"RGB8"       )) swapchain = new SwapChain(type,width,height,buffers,ptrs,FrameBufferRGB8      ::create);
    else throw std::runtime_error("unknown framebuffer type: "+std::string(type));
    return (Device::RTFrameBuffer) new ConstHandle<SwapChain>(swapchain);
  }

  void* SingleRayDevice::rtMapFrameBuffer(Device::RTFrameBuffer frameBuffer_i, int bufID) 
  {
    RT_COMMAND_HEADER;
    Ref<ConstHandle<SwapChain> > frameBuffer = castHandle<ConstHandle<SwapChain> >(frameBuffer_i,"framebuffer");
    Ref<SwapChain> instance = frameBuffer->getInstance();
    if (bufID < 0) bufID = instance->id();
    instance->buffer(bufID)->wait();
    return instance->buffer(bufID)->getData();
  }

  void SingleRayDevice::rtUnmapFrameBuffer(Device::RTFrameBuffer frameBuffer_i, int bufID) 
  {
    RT_COMMAND_HEADER;
    castHandle<ConstHandle<SwapChain> >(frameBuffer_i,"framebuffer");
  }

  void SingleRayDevice::rtSwapBuffers(Device::RTFrameBuffer frameBuffer_i) 
  {
    RT_COMMAND_HEADER;
    Ref<SwapChain> frameBuffer = castHandle<ConstHandle<SwapChain> >(frameBuffer_i,"framebuffer")->getInstance();
    frameBuffer->swapBuffers();
  }

  void SingleRayDevice::rtIncRef(Device::RTHandle handle) {
    ((_RTHandle*)handle)->refInc();
  }

  void SingleRayDevice::rtDecRef(Device::RTHandle handle) {
    ((_RTHandle*)handle)->refDec();
  }

  /*******************************************************************
                  setting of parameters
  *******************************************************************/

  void SingleRayDevice::rtSetBool1(Device::RTHandle handle, const char* property, bool x) {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(x));
  }

  void SingleRayDevice::rtSetBool2(Device::RTHandle handle, const char* property, bool x, bool y)  {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(x,y));
  }

  void SingleRayDevice::rtSetBool3(Device::RTHandle handle, const char* property, bool x, bool y, bool z) {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(x,y,z));
  }

  void SingleRayDevice::rtSetBool4(Device::RTHandle handle, const char* property, bool x, bool y, bool z, bool w)  {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(x,y,z,w));
  }

  void SingleRayDevice::rtSetInt1(Device::RTHandle handle, const char* property, int x)  {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) {
      if      (!strcmp(property,"serverID"   )) g_serverID = x;
      else if (!strcmp(property,"serverCount")) g_serverCount = x;
      return;
    }
    ((_RTHandle*)handle)->set(property,Variant(x));
  }

  void SingleRayDevice::rtSetInt2(Device::RTHandle handle, const char* property, int x, int y)  {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(x,y));
  }

  void SingleRayDevice::rtSetInt3(Device::RTHandle handle, const char* property, int x, int y, int z)  {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(x,y,z));
  }

  void SingleRayDevice::rtSetInt4(Device::RTHandle handle, const char* property, int x, int y, int z, int w)  {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(x,y,z,w));
  }

  void SingleRayDevice::rtSetFloat1(Device::RTHandle handle, const char* property, float x)  {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(x));
  }

  void SingleRayDevice::rtSetFloat2(Device::RTHandle handle, const char* property, float x, float y)  {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(x,y));
  }

  void SingleRayDevice::rtSetFloat3(Device::RTHandle handle, const char* property, float x, float y, float z)  {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(x,y,z));
  }

  void SingleRayDevice::rtSetFloat4(Device::RTHandle handle, const char* property, float x, float y, float z, float w)  {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(x,y,z,w));
  }

  void SingleRayDevice::rtSetArray(Device::RTHandle handle_i, const char* property, const char* type, Device::RTData data_i, size_t size, size_t stride, size_t ofs)
  {
    RT_COMMAND_HEADER;
    Ref<ConstHandle<Data> > data    = castHandle<ConstHandle<Data> >(data_i,"data");
    _RTHandle* handle = (_RTHandle*)handle_i;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    if      (!strcasecmp(type,"bool1" )) handle->set(property,Variant(data->getInstance(),Variant::BOOL1 ,size,stride == size_t(-1) ? 1*sizeof(bool ) : stride, ofs));
    else if (!strcasecmp(type,"bool2" )) handle->set(property,Variant(data->getInstance(),Variant::BOOL2 ,size,stride == size_t(-1) ? 2*sizeof(bool ) : stride, ofs));
    else if (!strcasecmp(type,"bool3" )) handle->set(property,Variant(data->getInstance(),Variant::BOOL3 ,size,stride == size_t(-1) ? 3*sizeof(bool ) : stride, ofs));
    else if (!strcasecmp(type,"bool4" )) handle->set(property,Variant(data->getInstance(),Variant::BOOL4 ,size,stride == size_t(-1) ? 4*sizeof(bool ) : stride, ofs)); 
    else if (!strcasecmp(type,"int1"  )) handle->set(property,Variant(data->getInstance(),Variant::INT1  ,size,stride == size_t(-1) ? 1*sizeof(int  ) : stride, ofs));
    else if (!strcasecmp(type,"int2"  )) handle->set(property,Variant(data->getInstance(),Variant::INT2  ,size,stride == size_t(-1) ? 2*sizeof(int  ) : stride, ofs));
    else if (!strcasecmp(type,"int3"  )) handle->set(property,Variant(data->getInstance(),Variant::INT3  ,size,stride == size_t(-1) ? 3*sizeof(int  ) : stride, ofs));
    else if (!strcasecmp(type,"int4"  )) handle->set(property,Variant(data->getInstance(),Variant::INT4  ,size,stride == size_t(-1) ? 4*sizeof(int  ) : stride, ofs));
    else if (!strcasecmp(type,"float1")) handle->set(property,Variant(data->getInstance(),Variant::FLOAT1,size,stride == size_t(-1) ? 1*sizeof(float) : stride, ofs));
    else if (!strcasecmp(type,"float2")) handle->set(property,Variant(data->getInstance(),Variant::FLOAT2,size,stride == size_t(-1) ? 2*sizeof(float) : stride, ofs));
    else if (!strcasecmp(type,"float3")) handle->set(property,Variant(data->getInstance(),Variant::FLOAT3,size,stride == size_t(-1) ? 3*sizeof(float) : stride, ofs)); 
    else if (!strcasecmp(type,"float4")) handle->set(property,Variant(data->getInstance(),Variant::FLOAT4,size,stride == size_t(-1) ? 4*sizeof(float) : stride, ofs));
    else throw std::runtime_error("unknown array type: "+std::string(type));
  }

  void SingleRayDevice::rtSetString(Device::RTHandle handle, const char* property, const char* str) {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(str));
  }

  void SingleRayDevice::rtSetImage(Device::RTHandle handle, const char* property, Device::RTImage img) {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    if (ConstHandle<Image>* image = dynamic_cast<ConstHandle<Image>*>((_RTHandle*)img)) {
      if (!image->getInstance()) throw std::runtime_error("invalid image value");
      ((_RTHandle*)handle)->set(property,Variant(image->getInstance()));
    } 
    else throw std::runtime_error("invalid image handle");
  }

  void SingleRayDevice::rtSetTexture(Device::RTHandle handle, const char* property, Device::RTTexture tex) {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    Ref<InstanceHandle<Texture> > texture = castHandle<InstanceHandle<Texture> >(tex,"texture");
    ((_RTHandle*)handle)->set(property,Variant(texture->getInstance()));
  }

  void SingleRayDevice::rtSetTransform(Device::RTHandle handle, const char* property, const float* transform)
  {
    RT_COMMAND_HEADER;
    if (!property) throw std::runtime_error("invalid property");
    if (!handle  ) return;
    ((_RTHandle*)handle)->set(property,Variant(copyFromArray(transform)));
  }

  void SingleRayDevice::rtClear(Device::RTHandle handle) {
    RT_COMMAND_HEADER;
    if (!handle) throw std::runtime_error("invalid handle");
    ((_RTHandle*)handle)->clear();
  }

  void SingleRayDevice::rtCommit(Device::RTHandle handle) {
    RT_COMMAND_HEADER;
    if (!handle) throw std::runtime_error("invalid handle");
    ((_RTHandle*)handle)->create();
  }

  /*******************************************************************
                            render call
  *******************************************************************/

  void SingleRayDevice::rtRenderFrame(Device::RTRenderer renderer_i, Device::RTCamera camera_i,
                                      Device::RTScene scene_i, Device::RTToneMapper toneMapper_i, 
                                      Device::RTFrameBuffer frameBuffer_i, int accumulate)
  {
    RT_COMMAND_HEADER;

    /* extract objects from handles */
    Ref<InstanceHandle<Renderer> >      renderer    = castHandle<InstanceHandle<Renderer     > >(renderer_i   ,"renderer"   );
    Ref<InstanceHandle<Camera> >        camera      = castHandle<InstanceHandle<Camera       > >(camera_i     ,"camera"     );
    Ref<BackendScene::Handle >          scene       = castHandle<BackendScene::Handle>          (scene_i      ,"scene"      );
    Ref<InstanceHandle<ToneMapper> >    toneMapper  = castHandle<InstanceHandle<ToneMapper   > >(toneMapper_i ,"tonemapper" );
    Ref<ConstHandle<SwapChain> > frameBuffer = castHandle<ConstHandle<SwapChain> >(frameBuffer_i,"framebuffer");

    /* render the frame */
    renderer->getInstance()->renderFrame(camera->getInstance(),scene->getInstance(),toneMapper->getInstance(),frameBuffer->getInstance(),accumulate);
  }

  bool SingleRayDevice::rtPick(Device::RTCamera camera_i, float x, float y, Device::RTScene scene_i, float& px, float& py, float& pz)
  {
    RT_COMMAND_HEADER;

    /* extract objects from handles */
    Ref<InstanceHandle<Camera> > camera = castHandle<InstanceHandle<Camera> >(camera_i,"camera");
    Ref<BackendScene::Handle >   scene  = castHandle<BackendScene::Handle>   (scene_i ,"scene" );

    /* trace ray */
    Ray ray; camera->getInstance()->ray(Vec2f(x,y), Vec2f(0.5f, 0.5f), ray);
    scene->getInstance()->intersector->intersect(ray);
    Vector3f p = ray.org + ray.tfar * ray.dir; 
    px = p.x; py = p.y; pz = p.z; 
    return (bool)ray;
  }
}
