// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "render_device.h"
#include "image/image.h"
#include "sys/taskscheduler_standard.h"

/* include general stuff */
#include "api/handle.h"
#include "api/data.h"
#include "api/parms.h"
#include "api/scene.h"
#include "sys/sync/mutex.h"

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
#include "materials/metallicpaint.h"
#include "materials/matte_textured.h"
#include "materials/obj.h"
#include "materials/velvet.h"

/* include all shapes */
#include "shapes/triangle.h"
#include "shapes/trianglemesh.h"
#include "shapes/sphere.h"

/* include all textures */
#include "textures/nearestneighbor.h"

/* include all tonemappers */
#include "tonemappers/defaulttonemapper.h"

/* include all renderers */
#include "renderers/debugrenderer.h"
#include "renderers/integratorrenderer.h"

/* include ray tracing core interface */
#include "rtcore/common/accel.h"
#include "rtcore/common/intersector.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#define strcasecmp lstrcmpiA
#pragma warning(disable:4297) // function assumed not to throw an exception but does
#endif

namespace embree
{
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

  /*! Primitive Handle */
  class PrimitiveHandle : public _RTHandle {
  public:

    /*! Constructs shape primitive. */
    PrimitiveHandle (const Ref<Shape>& shape, const Ref<Material>& material, const AffineSpace3f& transform)
      : shape(shape), material(material), transform(transform) {}

    /*! Constructs light primitive. */
    PrimitiveHandle (const Ref<Light>& light, const AffineSpace3f& transform)
      : light(light), transform(transform) {}

    /*! Constructs new primitive. */
    PrimitiveHandle (const Ref<Shape>& shape, const Ref<Light>& light, const Ref<Material>& material, const AffineSpace3f& transform)
      : shape(shape), light(light), material(material), transform(transform) {}

    /*! Creation is not allowed. */
    void create() { throw std::runtime_error("cannot modify constant handle"); }

    /*! Setting parameters is not allowed. */
    void set(const std::string& property, const Variant& data) { throw std::runtime_error("cannot modify constant handle"); }
  public:
    Ref<Shape> shape;       //!< Shape in case of a shape primitive
    Ref<Light> light;       //!< Light in case of a light primitive
    Ref<Material> material; //!< Material of shape primitive
    AffineSpace3f transform;  //!< Transformation of primitive
  };

  /*******************************************************************
                  construction
  *******************************************************************/
  
  RenderDevice::RenderDevice(size_t numThreads) 
  {
    if (numThreads != 0) {
      if (scheduler) delete scheduler;
      scheduler = new TaskSchedulerStandard(numThreads); // FIXME: this is only safe when creating a single device
    }
  }

  RenderDevice::~RenderDevice() {
    rtcFreeMemory();
  }

  /*******************************************************************
                    creation of handles
  *******************************************************************/

  Device::RTCamera RenderDevice::rtNewCamera(const char* type)
  {
    Lock<MutexSys> lock(mutex);
    if      (!strcasecmp(type,"pinhole")) return (Device::RTCamera) new NormalHandle<PinholeCamera,Camera>;
    else if (!strcasecmp(type,"depthoffield")) return (Device::RTCamera) new NormalHandle<DepthOfFieldCamera,Camera>;
    else throw std::runtime_error("unknown camera type: "+std::string(type));
  }

  Device::RTData RenderDevice::rtNewData(const char* type, size_t bytes, const void* data)
  {
    if (!strcasecmp(type,"immutable")) 
      return (Device::RTData) new ConstHandle<Data>(new Data(bytes,data,true));
    else if (!strcasecmp(type,"immutable_managed")) 
      return (Device::RTData) new ConstHandle<Data>(new Data(bytes,data,false));
    else
      throw std::runtime_error("unknown data buffer type: "+std::string(type));
  }

  Device::RTData RenderDevice::rtNewDataFromFile(const char* type, const char* fileName, size_t offset, size_t bytes)
  {
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

  Device::RTImage RenderDevice::rtNewImage(const char* type, size_t width, size_t height, const void* data)
  {
    Lock<MutexSys> lock(mutex);
    if (!strcasecmp(type,"RGB8")) {
      return (Device::RTImage) new ConstHandle<Image>(new Image3c(width,height,(const Col3c*)data));
    }
    else if (!strcasecmp(type,"RGB_FLOAT32")) {
      return (Device::RTImage) new ConstHandle<Image>(new Image3f(width,height,(const Col3f*)data));
    }
    else throw std::runtime_error("unknown image type: "+std::string(type));
  }

  Device::RTImage RenderDevice::rtNewImageFromFile(const char* file)
  {
    Lock<MutexSys> lock(mutex);
    if (!strncmp(file,"server:",7)) file += 7;
    Ref<Image> image = loadImage(file);
    if (image) 
      return (Device::RTImage) new ConstHandle<Image>(image);
    else
      return (Device::RTImage) new ConstHandle<Image>(new Image3c(1,1,Col3c(255,255,255)));
  }

  Device::RTTexture RenderDevice::rtNewTexture(const char* type) {
    Lock<MutexSys> lock(mutex);
    if (!strcasecmp(type,"nearest")) return (Device::RTTexture) new NormalHandle<NearestNeighbor,Texture>;
    else if (!strcasecmp(type,"image")) return (Device::RTTexture) new NormalHandle<NearestNeighbor,Texture>;
    else throw std::runtime_error("unsupported texture type: "+std::string(type));
  }

  Device::RTMaterial RenderDevice::rtNewMaterial(const char* type)
  {
    Lock<MutexSys> lock(mutex);
    if      (!strcasecmp(type,"Matte")         ) return (Device::RTMaterial) new NormalHandle<Matte,Material>;
    else if (!strcasecmp(type,"Plastic")       ) return (Device::RTMaterial) new NormalHandle<Plastic,Material>;
    else if (!strcasecmp(type,"Dielectric")    ) return (Device::RTMaterial) new NormalHandle<Dielectric,Material>;
    else if (!strcasecmp(type,"Glass")         ) return (Device::RTMaterial) new NormalHandle<Dielectric,Material>;
    else if (!strcasecmp(type,"ThinDielectric")) return (Device::RTMaterial) new NormalHandle<ThinDielectric,Material>;
    else if (!strcasecmp(type,"ThinGlass")     ) return (Device::RTMaterial) new NormalHandle<ThinDielectric,Material>;
    else if (!strcasecmp(type,"Mirror")        ) return (Device::RTMaterial) new NormalHandle<Mirror,Material>;
    else if (!strcasecmp(type,"Metal")         ) return (Device::RTMaterial) new NormalHandle<Metal,Material>;
    else if (!strcasecmp(type,"MetallicPaint") ) return (Device::RTMaterial) new NormalHandle<MetallicPaint,Material>;
    else if (!strcasecmp(type,"MatteTextured") ) return (Device::RTMaterial) new NormalHandle<MatteTextured,Material>;
    else if (!strcasecmp(type,"Obj")           ) return (Device::RTMaterial) new NormalHandle<Obj,Material>;
    else if (!strcasecmp(type,"Velvet")        ) return (Device::RTMaterial) new NormalHandle<Velvet,Material>;
    else throw std::runtime_error("unknown material type: "+std::string(type));
  }

  Device::RTShape RenderDevice::rtNewShape(const char* type) {
    Lock<MutexSys> lock(mutex);
    if      (!strcasecmp(type,"trianglemesh")) return (Device::RTShape) new NormalHandle2<TriangleMesh::Handle,Shape>;
    else if (!strcasecmp(type,"triangle")    ) return (Device::RTShape) new NormalHandle<Triangle,Shape>;
    else if (!strcasecmp(type,"sphere")      ) return (Device::RTShape) new NormalHandle<Sphere,Shape>;
    else throw std::runtime_error("unknown shape type: "+std::string(type));
  }

  Device::RTLight RenderDevice::rtNewLight(const char* type)
  {
    Lock<MutexSys> lock(mutex);
    if      (!strcasecmp(type,"ambientlight"    )) return (Device::RTLight) new NormalHandle<AmbientLight,Light>;
    else if (!strcasecmp(type,"pointlight"      )) return (Device::RTLight) new NormalHandle<PointLight,Light>;
    else if (!strcasecmp(type,"spotlight"       )) return (Device::RTLight) new NormalHandle<SpotLight,Light>;
    else if (!strcasecmp(type,"directionallight")) return (Device::RTLight) new NormalHandle<DirectionalLight,Light>;
    else if (!strcasecmp(type,"distantlight"    )) return (Device::RTLight) new NormalHandle<DistantLight,Light>;
    else if (!strcasecmp(type,"hdrilight"       )) return (Device::RTLight) new NormalHandle<HDRILight,Light>;
    else if (!strcasecmp(type,"trianglelight"   )) return (Device::RTLight) new NormalHandle<TriangleLight,Light>;
    else throw std::runtime_error("unknown light type: "+std::string(type));
  }

  Device::RTPrimitive RenderDevice::rtNewShapePrimitive(Device::RTShape shape_i, Device::RTMaterial material_i, const float* transform)
  {
    Lock<MutexSys> lock(mutex);
    BaseHandle<Shape>*    shape    = castHandle<BaseHandle<Shape>    >(shape_i   ,"shape"   );
    BaseHandle<Material>* material = castHandle<BaseHandle<Material> >(material_i,"material");
    AffineSpace3f space = transform ? copyFromArray(transform) : AffineSpace3f(one);
    return (Device::RTPrimitive) new PrimitiveHandle(shape->instance,material->instance,space);
  }

  Device::RTPrimitive RenderDevice::rtNewLightPrimitive(Device::RTLight light_i, const float* transform)
  {
    Lock<MutexSys> lock(mutex);
    BaseHandle<Light>* light = castHandle<BaseHandle<Light> >(light_i,"light");
    AffineSpace3f space = transform ? copyFromArray(transform) : AffineSpace3f(one);
    return (Device::RTPrimitive) new PrimitiveHandle(light->instance,space);
  }

  Device::RTPrimitive RenderDevice::rtTransformPrimitive(Device::RTPrimitive primitive, const float* transform) 
  {
    Lock<MutexSys> lock(mutex);
    PrimitiveHandle* prim = dynamic_cast<PrimitiveHandle*>((_RTHandle*)primitive);
    AffineSpace3f space = transform ? copyFromArray(transform) : AffineSpace3f(one);
    return (Device::RTPrimitive) new PrimitiveHandle(prim->shape,prim->light,prim->material,space * prim->transform);
  }

  void calculateSize(Device::RTPrimitive primitive, size_t& numTriangles, size_t& numVertices)
  {
    PrimitiveHandle* prim = dynamic_cast<PrimitiveHandle*>((_RTHandle*)primitive);
    if (!prim) throw std::runtime_error("invalid primitive");
    
    if (prim->shape) {
      numVertices  += prim->shape->numVertices();
      numTriangles += prim->shape->numTriangles();
    }
    else if (prim->light)
    {
      if (Ref<Shape> shape = prim->light->shape()) {
        numVertices  += shape->numVertices();
        numTriangles += shape->numTriangles();
      }
    }
    else throw std::runtime_error("invalid primitive");
  }
  
  BBox3f extractTriangles(Ref<BackendScene>& scene, Device::RTPrimitive primitive, 
                          size_t& instOfs, BuildTriangle* triangles, size_t& numTriangles, BuildVertex* vertices, size_t& numVertices)
  {
    BBox3f bounds = empty;
    PrimitiveHandle* prim = dynamic_cast<PrimitiveHandle*>((_RTHandle*)primitive);
    if (!prim) throw std::runtime_error("invalid primitive");

    /* extract geometry */
    if (prim->shape)
    {
      Ref<Shape> shape = prim->shape->transform(prim->transform);
      size_t id = scene->add(new Instance(instOfs++,shape,prim->material,null));
      bounds.grow(shape->extract(id,triangles,numTriangles,vertices,numVertices));
    }
    
    /* extract lights */
    else if (prim->light)
    {
      Ref<Light> light = prim->light->transform(prim->transform);
      scene->add(light);
      if (Ref<Shape> shape = light->shape()) 
      {
        size_t id = scene->add(new Instance(instOfs++,shape,null,light.dynamicCast<AreaLight>()));
        bounds.grow(shape->extract(id,triangles,numTriangles,vertices,numVertices));
      }
    }
    else throw std::runtime_error("invalid primitive");

    return bounds;
  }
  
  Device::RTScene RenderDevice::rtNewScene(const char* type, Device::RTPrimitive* prims, size_t size)
  {
    Lock<MutexSys> lock(mutex);
    Ref<BackendScene> scene = new BackendScene;
    
    /* count number of vertices and triangles */
    size_t numAllocatedTriangles = 0;
    size_t numAllocatedVertices = 0;
    for (size_t i=0; i<size; i++)
      calculateSize(prims[i],numAllocatedTriangles,numAllocatedVertices);

    /* allocate triangle buffers */
    size_t numMeshes = 0;
    size_t numTriangles = 0;
    size_t numVertices = 0;
    BuildTriangle* triangles = (BuildTriangle*) rtcMalloc(numAllocatedTriangles*sizeof(BuildTriangle));
    BuildVertex*   vertices  = (BuildVertex*  ) rtcMalloc(numAllocatedVertices *sizeof(BuildVertex));

    /* extract all primitives */
    BBox3f bounds = empty;
    for (size_t i=0; i<size; i++) {
      bounds.grow(extractTriangles(scene,prims[i],numMeshes,triangles,numTriangles,vertices,numVertices));
      if (numTriangles > numAllocatedTriangles) throw std::runtime_error("internal error");
      if (numVertices  > numAllocatedVertices ) throw std::runtime_error("internal error");
    }

    /* build acceleration structure */
    size_t split = strcspn(type, " \t\r");
    std::string accelType(type,split);
    std::string triType(type+split+1);
    scene->setAccel(rtcCreateAccel(accelType.c_str(), triType.c_str(), triangles, numTriangles, vertices, numVertices, bounds));

    return (Device::RTScene) new ConstHandle<BackendScene>(scene);
  }

  Device::RTToneMapper RenderDevice::rtNewToneMapper(const char* type)
  {
    if      (!strcasecmp(type,"default")) return (Device::RTToneMapper) new NormalHandle<DefaultToneMapper,ToneMapper>;
    else throw std::runtime_error("unknown tonemapper type: "+std::string(type));
  }

  Device::RTRenderer RenderDevice::rtNewRenderer(const char* type)
  {
    Lock<MutexSys> lock(mutex);
    if (!strcasecmp(type,"debug"     )) return (Device::RTRenderer) new NormalHandle<DebugRenderer,Renderer>;
    if (!strcasecmp(type,"pathtracer")) {
      NormalHandle<IntegratorRenderer,Renderer>* handle = new NormalHandle<IntegratorRenderer,Renderer>;
      handle->set("integrator",Variant("pathtracer"));
      return (Device::RTRenderer) handle;
    }
    else throw std::runtime_error("unknown renderer type: " + std::string(type));
  }

  Device::RTFrameBuffer RenderDevice::rtNewFrameBuffer(const char* type, size_t width, size_t height, size_t buffers) {
    Lock<MutexSys> lock(mutex);
    if (!strcasecmp(type,"RGB_FLOAT32")) return (Device::RTFrameBuffer) new ConstHandle<SwapChain<> >(new SwapChain<>(width,height,buffers));
    else throw std::runtime_error("unknown framebuffer type: "+std::string(type));
  }

  void RenderDevice::rtGetFrameBufferSize(Device::RTFrameBuffer frameBuffer_i, size_t& width, size_t& height) 
  {
    ConstHandle<SwapChain<> >* frameBuffer = castHandle<ConstHandle<SwapChain<> > >(frameBuffer_i,"framebuffer");
    width = frameBuffer->instance->getWidth();
    height = frameBuffer->instance->getHeight();
  }

  void* RenderDevice::rtMapFrameBuffer(Device::RTFrameBuffer frameBuffer_i) {
    ConstHandle<SwapChain<> >* frameBuffer = castHandle<ConstHandle<SwapChain<> > >(frameBuffer_i,"framebuffer");
    return &frameBuffer->instance->buffer()->get(0,0);
  }

  void RenderDevice::rtUnmapFrameBuffer(Device::RTFrameBuffer frameBuffer_i) {
    Lock<MutexSys> lock(mutex);
    castHandle<ConstHandle<SwapChain<> > >(frameBuffer_i,"framebuffer");
  }

  void RenderDevice::rtSwapBuffers(Device::RTFrameBuffer frameBuffer_i) 
  {
    Lock<MutexSys> lock(mutex);
    Ref<SwapChain<> > frameBuffer = castHandle<ConstHandle<SwapChain<> > >(frameBuffer_i,"framebuffer")->instance;
    frameBuffer->swapBuffers();
  }

  void RenderDevice::rtIncRef(Device::RTHandle handle) {
    ((_RTHandle*)handle)->incRef();
  }

  void RenderDevice::rtDecRef(Device::RTHandle handle) {
    ((_RTHandle*)handle)->decRef();
  }

  /*******************************************************************
                  setting of parameters
  *******************************************************************/

  void RenderDevice::rtSetBool1(Device::RTHandle handle, const char* property, bool x) {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x));
  }

  void RenderDevice::rtSetBool2(Device::RTHandle handle, const char* property, bool x, bool y)  {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x,y));
  }

  void RenderDevice::rtSetBool3(Device::RTHandle handle, const char* property, bool x, bool y, bool z) {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x,y,z));
  }

  void RenderDevice::rtSetBool4(Device::RTHandle handle, const char* property, bool x, bool y, bool z, bool w)  {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x,y,z,w));
  }

  void RenderDevice::rtSetInt1(Device::RTHandle handle, const char* property, int x)  {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x));
  }

  void RenderDevice::rtSetInt2(Device::RTHandle handle, const char* property, int x, int y)  {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x,y));
  }

  void RenderDevice::rtSetInt3(Device::RTHandle handle, const char* property, int x, int y, int z)  {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x,y,z));
  }

  void RenderDevice::rtSetInt4(Device::RTHandle handle, const char* property, int x, int y, int z, int w)  {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x,y,z,w));
  }

  void RenderDevice::rtSetFloat1(Device::RTHandle handle, const char* property, float x)  {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x));
  }

  void RenderDevice::rtSetFloat2(Device::RTHandle handle, const char* property, float x, float y)  {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x,y));
  }

  void RenderDevice::rtSetFloat3(Device::RTHandle handle, const char* property, float x, float y, float z)  {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x,y,z));
  }

  void RenderDevice::rtSetFloat4(Device::RTHandle handle, const char* property, float x, float y, float z, float w)  {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(x,y,z,w));
  }

  void RenderDevice::rtSetArray(Device::RTHandle handle_i, const char* property, const char* type, Device::RTData data_i, size_t size, size_t stride, size_t ofs)
  {
    Lock<MutexSys> lock(mutex);
    ConstHandle<Data>* data    = castHandle<ConstHandle<Data> >(data_i,"data");
    _RTHandle* handle = (_RTHandle*)handle_i;
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    if      (!strcasecmp(type,"bool1" )) handle->set(property,Variant(data->instance,Variant::BOOL1 ,size,stride == size_t(-1) ? 1*sizeof(bool ) : stride, ofs));
    else if (!strcasecmp(type,"bool2" )) handle->set(property,Variant(data->instance,Variant::BOOL2 ,size,stride == size_t(-1) ? 2*sizeof(bool ) : stride, ofs));
    else if (!strcasecmp(type,"bool3" )) handle->set(property,Variant(data->instance,Variant::BOOL3 ,size,stride == size_t(-1) ? 3*sizeof(bool ) : stride, ofs));
    else if (!strcasecmp(type,"bool4" )) handle->set(property,Variant(data->instance,Variant::BOOL4 ,size,stride == size_t(-1) ? 4*sizeof(bool ) : stride, ofs)); 
    else if (!strcasecmp(type,"int1"  )) handle->set(property,Variant(data->instance,Variant::INT1  ,size,stride == size_t(-1) ? 1*sizeof(int  ) : stride, ofs));
    else if (!strcasecmp(type,"int2"  )) handle->set(property,Variant(data->instance,Variant::INT2  ,size,stride == size_t(-1) ? 2*sizeof(int  ) : stride, ofs));
    else if (!strcasecmp(type,"int3"  )) handle->set(property,Variant(data->instance,Variant::INT3  ,size,stride == size_t(-1) ? 3*sizeof(int  ) : stride, ofs));
    else if (!strcasecmp(type,"int4"  )) handle->set(property,Variant(data->instance,Variant::INT4  ,size,stride == size_t(-1) ? 4*sizeof(int  ) : stride, ofs));
    else if (!strcasecmp(type,"float1")) handle->set(property,Variant(data->instance,Variant::FLOAT1,size,stride == size_t(-1) ? 1*sizeof(float) : stride, ofs));
    else if (!strcasecmp(type,"float2")) handle->set(property,Variant(data->instance,Variant::FLOAT2,size,stride == size_t(-1) ? 2*sizeof(float) : stride, ofs));
    else if (!strcasecmp(type,"float3")) handle->set(property,Variant(data->instance,Variant::FLOAT3,size,stride == size_t(-1) ? 3*sizeof(float) : stride, ofs)); 
    else if (!strcasecmp(type,"float4")) handle->set(property,Variant(data->instance,Variant::FLOAT4,size,stride == size_t(-1) ? 4*sizeof(float) : stride, ofs));
    else throw std::runtime_error("unknown array type: "+std::string(type));
  }

  void RenderDevice::rtSetString(Device::RTHandle handle, const char* property, const char* str) {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(str));
  }

  void RenderDevice::rtSetImage(Device::RTHandle handle, const char* property, Device::RTImage img) {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    if (ConstHandle<Image>* image = dynamic_cast<ConstHandle<Image>*>((_RTHandle*)img)) {
      if (!image->instance) throw std::runtime_error("invalid image value");
      ((_RTHandle*)handle)->set(property,Variant(image->instance));
    } 
    else throw std::runtime_error("invalid image handle");
  }

  void RenderDevice::rtSetTexture(Device::RTHandle handle, const char* property, Device::RTTexture tex) {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    BaseHandle<Texture>* texture = castHandle<BaseHandle<Texture> >(tex,"texture");
    ((_RTHandle*)handle)->set(property,Variant(texture->instance));
  }

  void RenderDevice::rtSetTransform(Device::RTHandle handle, const char* property, const float* transform)
  {
    Lock<MutexSys> lock(mutex);
    if (!handle  ) throw std::runtime_error("invalid handle"  );
    if (!property) throw std::runtime_error("invalid property");
    ((_RTHandle*)handle)->set(property,Variant(copyFromArray(transform)));
  }

  void RenderDevice::rtClear(Device::RTHandle handle) {
    Lock<MutexSys> lock(mutex);
    if (!handle) throw std::runtime_error("invalid handle");
    ((_RTHandle*)handle)->clear();
  }

  void RenderDevice::rtCommit(Device::RTHandle handle) {
    Lock<MutexSys> lock(mutex);
    if (!handle) throw std::runtime_error("invalid handle");
    ((_RTHandle*)handle)->create();
  }

  /*******************************************************************
                            render call
  *******************************************************************/

  void RenderDevice::rtRenderFrame(Device::RTRenderer renderer_i, Device::RTCamera camera_i,
                                   Device::RTScene scene_i, Device::RTToneMapper toneMapper_i, 
                                   Device::RTFrameBuffer frameBuffer_i, int accumulate)
  {
    Lock<MutexSys> lock(mutex);

    /* extract objects from handles */
    BaseHandle<Renderer>*      renderer    = castHandle<BaseHandle<Renderer     > >(renderer_i   ,"renderer"   );
    BaseHandle<Camera>*        camera      = castHandle<BaseHandle<Camera       > >(camera_i     ,"camera"     );
    ConstHandle<BackendScene>* scene       = castHandle<ConstHandle<BackendScene> >(scene_i      ,"scene"      );
    BaseHandle<ToneMapper>*    toneMapper  = castHandle<BaseHandle<ToneMapper   > >(toneMapper_i ,"tonemapper" );
    ConstHandle<SwapChain<> >* frameBuffer = castHandle<ConstHandle<SwapChain<> > >(frameBuffer_i,"framebuffer");
    
    /* render the frame */
    renderer->instance->renderFrame(camera->instance,scene->instance,toneMapper->instance,frameBuffer->instance,accumulate);
  }

  bool RenderDevice::rtPick(Device::RTCamera camera_i, float x, float y, Device::RTScene scene_i, float& px, float& py, float& pz)
  {
    Lock<MutexSys> lock(mutex);

    /* extract objects from handles */
    BaseHandle<Camera>*        camera = castHandle<BaseHandle<Camera> >(camera_i,"camera");
    ConstHandle<BackendScene>* scene  = castHandle<ConstHandle<BackendScene> >(scene_i,"scene");

    /* trace ray */
    Ray ray; camera->instance->ray(Vec2f(x,y), Vec2f(0.5f, 0.5f), ray);
    Hit hit; scene->instance->intersector->intersect(ray,hit);
    Vec3f p = ray.org + hit.t * ray.dir; 
    px = p.x; py = p.y; pz = p.z; 
    return (bool)hit;
  }
}
