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

#ifndef __EMBREE_NATIVE_DEVICE_H__
#define __EMBREE_NATIVE_DEVICE_H__

#include <stddef.h>

/*! \file api.h This file implements the interface to the renderer
 *  backend. The library has to get initialized by calling rtInit() at
 *  the beginning and rtExit() at the end of your application. The API
 *  is completely functional, meaning that objects can NOT be
 *  modified. Handles are references to objects, and objects are
 *  internally reference counted, thus destroyed when no longer
 *  needed. Calling rtDelete does only delete the handle, not the
 *  object the handle references. The first parameter of every
 *  rtNewXXX function is a type that selects the object to create. For
 *  an overview of the parameters supported by a specific object see
 *  the header file of the object implementation. Using the rtSetXXX
 *  functions buffers the set values inside the handle. A subsequent
 *  call to rtCommit will set the handle reference to a new object
 *  with the changed parameters. The original object is not changed by
 *  this process. The semantics of modifying an object A used by
 *  another object B can only be achieved by creating A' and a new B'
 *  that uses A'. The RTImage, RTFramebuffer, and RTScene handles are
 *  constant, thus the rtSetXXX and rtCommit function cannot be used
 *  for them. */

namespace embree
{
  /*******************************************************************
                           Embree Device
  *******************************************************************/
  typedef int light_mask_t;

  class Device
  {
  public:
    virtual ~Device() {}

    /*! Create new Embree device. */
    static Device* rtCreateDevice(const char* type = "default", size_t numThreads = 0, size_t verbose = 0);

    /*******************************************************************
                        type definitions
    *******************************************************************/
    
    /*! Generic handle. */
    typedef struct __RTHandle {}* RTHandle;
    
    /*! Transformation handle. */
    typedef struct __RTTransform : public __RTHandle { }* RTTransform;
    
    /*! Camera handle. */
    typedef struct __RTCamera : public __RTHandle { }* RTCamera;
    
    /*! Data handle (constant). */
    typedef struct __RTData : public __RTHandle { }* RTData;
    
    /*! Image handle (constant). */
    typedef struct __RTImage : public __RTHandle { }* RTImage;
    
    /*! Texture handle. */
    typedef struct __RTTexture : public __RTHandle { }* RTTexture;
    
    /*! Material handle. */
    typedef struct __RTMaterial : public __RTHandle { }* RTMaterial;
    
    /*! Shape handle. */
    typedef struct __RTShape : public __RTHandle { }* RTShape;
    
    /*! Light handle. */
    typedef struct __RTLight : public __RTHandle { }* RTLight;
    
    /*! Primitive handle. */
    typedef struct __RTPrimitive : public __RTHandle { }* RTPrimitive;
    
    /*! Scene handle (constant). */
    typedef struct __RTScene : public __RTHandle { }* RTScene;

    /*! Tonemapper handle. */
    typedef struct __RTToneMapper : public __RTHandle { }* RTToneMapper;

    /*! Renderer handle. */
    typedef struct __RTRenderer : public __RTHandle { }* RTRenderer;
    
    /*! Framebuffer handle (constant). */
    typedef struct __RTFrameBuffer : public __RTHandle { }* RTFrameBuffer;
    
    /*! Ray structure. Describes ray layout for ray queries. */
    struct RTRay {
      struct vec3 { float x,y,z; };
      vec3 org;   //!< origin of ray
      float near; //!< start of ray segment
      vec3 dir;   //!< direction of ray
      float far;  //!< end of ray segment
    };
    
    /*! Hit structure. Describes hit layout for ray queries. */
    struct RTHit {
      int prim;    //!< ID of hit primitive
      float dist;  //!< distance of hit
    };
    
    /*******************************************************************
                         creation of objects
    *******************************************************************/

    /*! Creates a new camera. \param type is the type of camera to
     *  create (e.g. "pinhole"). \returns camera handle. */
    virtual RTCamera rtNewCamera(const char* type) = 0;

    /*! Creates a new data object. \param type is the type of the data
     *  buffer, can be "immutable" (constant buffer and data gets copied
     *  on construction) or "immutable_managed" (constant buffer but
     *  provided data will be deleted by the system using free). \param
     *  bytes is the number of bytes for the data buffer \param data is
     *  the data to fill the buffer with. \returns data buffer handle */
    virtual RTData rtNewData(const char* type, size_t bytes, const void* data) = 0;

    /*! Creates a new data object and initializes its content from a
     *  file. \param type is the type of the data buffer, can only be
     *  "immutable" (constant buffer). \param file is the name of the
     *  file to open. If the filename starts with "server:" the data
     *  is loaded directly on the rendering servers in network
     *  mode. \param offset is the location on the file to start
     *  reading \param bytes is the number of bytes for the data
     *  buffer. \returns data buffer handle */
    virtual RTData rtNewDataFromFile(const char* type, const char* file, size_t offset, size_t bytes) = 0;

    /*! Creates a new image. The data gets directly copied by this
     *  function. \param type is the type of image to create
     *  (e.g. "RGB8", "RGB_FLOAT32"). \param width is the width of the
     *  image \param height is the height of the image \param data is a
     *  pointer to the image data. \returns image handle */
    virtual RTImage rtNewImage(const char* type, size_t width, size_t height, const void* data, const bool copy = true) = 0;

    /*! Creates a new image from a file. \param file is the file to
     *  read the image from. If the filename starts with "server:" the
     *  image is loaded directly on rendering servers in network
     *  mode. \returns image handle */
    virtual RTImage rtNewImageFromFile(const char* file) = 0;

    /*! Creates a texture. \param type is the type of texture to create
     *  (e.g. "image" mapped). \returns texture handle */
    virtual RTTexture rtNewTexture(const char* type) = 0;

    /*! Creates a new material. \param type is the type of material to
     *  create (e.g. "Matte", "Plastic", "Dielectric", "ThinDielectric",
     *  "Mirror", "Metal", "MetallicPaint", "MatteTextured",
     *  "Obj"). \returns material handle */
    virtual RTMaterial rtNewMaterial(const char* type) = 0;

      /*! Creates a new shape. \param type is the type of shape to create
       *  (e.g. "trianglemesh", "triangle", "sphere") \returns shape
       *  handle */
    virtual RTShape rtNewShape(const char* type) = 0;

    /*! Creates a new light source. \param type is the type of shape to
     *  create (e.g. "ambientlight", "pointlight", "spotlight",
     *  "directionallight", "distantlight", "hdrilight",
     *  "trianglelight"). \returns light handle */
    virtual RTLight rtNewLight(const char* type) = 0;

    /*! Creates a new shape primitive. \param shape is the shape to
     *  instantiate \param material is the material to attach to the
     *  shape \param transform is an optional pointer to a
     *  transformation to transform the shape \returns primitive
     *  handle */
    virtual RTPrimitive rtNewShapePrimitive(RTShape shape, RTMaterial material, const float* transform = NULL) = 0;

    /*! Creates a new light primitive. \param light is the light to
     *  instantiate \param transform is an optional pointer to a
     *  transformation to transform the shape \returns primitive
     *  handle */
    virtual RTPrimitive rtNewLightPrimitive(RTLight light, RTMaterial material, const float* transform = NULL) = 0;

    /*! Transforms a primitive. \param prim is the primitive to
     *  transform \param transform is a pointer to a
     *  transformation \returns primitive
     *  handle */
    virtual RTPrimitive rtTransformPrimitive(RTPrimitive prim, const float* transform) = 0;

    /*! Creates a new scene. \returns scene handle */
    virtual RTScene rtNewScene(const char* type) = 0; 

    /*! Adds or deletes a primitive to/from the scene. Primitives are
        deleted by adding NULL to a primitive slot. */
    virtual void rtSetPrimitive(RTScene scene, size_t slot, RTPrimitive prim) = 0;

    /*! Creates a new tonemapper. \returns tonemapper handle. */
    virtual RTToneMapper rtNewToneMapper(const char* type) = 0;

    /*! Creates a new renderer. \param type is the type of renderer to
     *  create (e.g. "debug", "pathtracer"). \returns renderer handle */
    virtual RTRenderer rtNewRenderer(const char* type) = 0;

    /*! Creates a new framebuffer. \param type is the type of
     *  framebuffer to create (e.g. "RGB_FLOAT32"). \param width is
     *  the width of the framebuffer in pixels \param height is the
     *  height of the framebuffer in pixels \param iteration is the
     *  current iteration for this frame buffer \param buffers are the
     *  number of buffers for multi-buffering \returns framebuffer
     *  handle */
    virtual RTFrameBuffer rtNewFrameBuffer(const char* type, size_t width, size_t height, size_t buffers = 1, void** ptrs = NULL) = 0;

    /*! Map the framebuffer data. \param frameBuffer is the framebuffer
     *  to map \returns pointer to framebuffer data */
    virtual void* rtMapFrameBuffer(RTFrameBuffer frameBuffer, int bufID = -1) = 0;

    /*! Unmap the framebuffer data. \param frameBuffer is the
     *  framebuffer to unmap. */
    virtual void rtUnmapFrameBuffer(RTFrameBuffer frameBuffer, int bufID = -1) = 0;

    /*! Switch to the next buffer in the chain. */
    virtual void rtSwapBuffers(RTFrameBuffer frameBuffer) = 0;

    /*! Increases the reference counter of the handle. */
    virtual void rtIncRef(RTHandle handle) = 0;

    /*! Decreases the reference counter of the handle. If zero is reached 
        the handle is destroyed, but not necessarily the referenced object. */
    virtual void rtDecRef(RTHandle handle) = 0;
    
    /*******************************************************************
                            setting of parameters
    *******************************************************************/

    /*! Sets a boolean parameter of the handle. */
    virtual void rtSetBool1(RTHandle handle, const char* property, bool x) = 0;

    /*! Sets a bool2 parameter of the handle. */
    virtual void rtSetBool2(RTHandle handle, const char* property, bool x, bool y) = 0;

    /*! Sets a bool3 parameter of the handle. */
    virtual void rtSetBool3(RTHandle handle, const char* property, bool x, bool y, bool z) = 0;

    /*! Sets a bool4 parameter of the handle. */
    virtual void rtSetBool4(RTHandle handle, const char* property, bool x, bool y, bool z, bool w) = 0;

    /*! Sets an integer parameter of the handle. */
    virtual void rtSetInt1(RTHandle handle, const char* property, int x) = 0;

    /*! Sets an int2 parameter of the handle. */
    virtual void rtSetInt2(RTHandle handle, const char* property, int x, int y) = 0;

    /*! Sets an int3 parameter of the handle. */
    virtual void rtSetInt3(RTHandle handle, const char* property, int x, int y, int z) = 0;

    /*! Sets an int4 parameter of the handle. */
    virtual void rtSetInt4(RTHandle handle, const char* property, int x, int y, int z, int w) = 0;

    /*! Sets a float parameter of the handle. */
    virtual void rtSetFloat1(RTHandle handle, const char* property, float x) = 0;

    /*! Sets a float2 parameter of the handle. */
    virtual void rtSetFloat2(RTHandle handle, const char* property, float x, float y) = 0;

    /*! Sets a float3 parameter of the handle. */
    virtual void rtSetFloat3(RTHandle handle, const char* property, float x, float y, float z) = 0;

    /*! Sets a float4 parameter of the handle. */
    virtual void rtSetFloat4(RTHandle handle, const char* property, float x, float y, float z, float w) = 0;

    /*! Sets an typed array parameter of the handle. The data is copied when calling rtCommit. */
    virtual void rtSetArray(RTHandle handle, const char* property, const char* type, RTData data, size_t size, size_t stride, size_t ofs) = 0;

    /*! Sets a string parameter of the handle. */
    virtual void rtSetString(RTHandle handle, const char* property, const char* str) = 0;

    /*! Sets an image parameter of the handle. */
    virtual void rtSetImage(RTHandle handle, const char* property, RTImage img) = 0;

    /*! Sets a texture parameter of the handle. */
    virtual void rtSetTexture(RTHandle handle, const char* property, RTTexture tex) = 0;

    /*! Sets a transformation of the handle. */
    virtual void rtSetTransform(RTHandle handle, const char* property, const float* transform) = 0;

    /*! Clear all parameters cached in the handle to reduce memory
     *  consumption. For instance, mesh handles no longer reference data
     *  buffers passed at creation time. */
    virtual void rtClear(RTHandle handle) = 0;

    /*! Commits all changes by setting the reference of the handle to a
     *  new object with specified parameters. */
    virtual void rtCommit(RTHandle handle) = 0;

    /*******************************************************************
                            render calls
    *******************************************************************/
    
    /*! Renders a frame. \param renderer is the renderer to use \param
     *  camera is the camera to use \param scene is the scene to
     *  render \param tonemapper is the tonemapper to use \parm
     *  frameBuffer is the framebuffer to render into */
    virtual void rtRenderFrame(RTRenderer renderer, RTCamera camera, RTScene scene, RTToneMapper tonemapper, RTFrameBuffer frameBuffer, int accumulate) = 0;

    /*! Pick a 3D point. \returns true if a point was picked, false otherwise
     *  \parm x is the x coordinate [0:1] in the image plane
     *  \parm y is the y coordinate [0:1] in the image plane
     *  \parm p is the world space position of the picked point, if any
     *  \parm camera is the camera to use \param scene is the scene for picking */
    virtual bool rtPick(RTCamera camera, float x, float y, RTScene scene, float& px, float& py, float& pz) = 0;
  };
}

#endif
