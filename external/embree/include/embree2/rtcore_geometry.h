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

#ifndef __RTCORE_GEOMETRY_H__
#define __RTCORE_GEOMETRY_H__

/*! \ingroup embree_kernel_api */
/*! \{ */

/*! invalid geometry ID */
#define RTC_INVALID_GEOMETRY_ID ((unsigned)-1)

/*! \brief Specifies the type of buffers when mapping buffers */
enum RTCBufferType {
  RTC_INDEX_BUFFER    = 0x01000000,
  RTC_VERTEX_BUFFER   = 0x02000000,
  RTC_VERTEX_BUFFER0  = 0x02000000,
  RTC_VERTEX_BUFFER1  = 0x02000001,
};

/*! \brief Supported types of matrix layout for functions involving matrices */
enum RTCMatrixType {
  RTC_MATRIX_ROW_MAJOR = 0,
  RTC_MATRIX_COLUMN_MAJOR = 1,
  RTC_MATRIX_COLUMN_MAJOR_ALIGNED16 = 2,
};

/*! \brief Supported geometry flags to specify handling in dynamic scenes. */
enum RTCGeometryFlags 
{
  RTC_GEOMETRY_STATIC     = 0,    //!< specifies static geometry that will change rarely
  RTC_GEOMETRY_DEFORMABLE = 1,    //!< specifies dynamic geometry with deformable motion (BVH refit possible)
  RTC_GEOMETRY_DYNAMIC    = 2,    //!< specifies dynamic geometry with arbitrary motion (BVH refit not possible)
};

/*! Intersection filter function for single rays. */
typedef void (*RTCFilterFunc)(void* ptr,           /*!< pointer to user data */
                              RTCRay& ray          /*!< intersection to filter */);

/*! Intersection filter function for ray packets of size 4. */
typedef void (*RTCFilterFunc4)(const void* valid,  /*!< pointer to valid mask */
                               void* ptr,          /*!< pointer to user data */
                               RTCRay4& ray        /*!< intersection to filter */);

/*! Intersection filter function for ray packets of size 8. */
typedef void (*RTCFilterFunc8)(const void* valid,  /*!< pointer to valid mask */
                               void* ptr,          /*!< pointer to user data */
                               RTCRay8& ray        /*!< intersection to filter */);

/*! Intersection filter function for ray packets of size 16. */
typedef void (*RTCFilterFunc16)(const void* valid, /*!< pointer to valid mask */
                                void* ptr,         /*!< pointer to user data */
                                RTCRay16& ray      /*!< intersection to filter */);

/*! \brief Creates a new scene instance. 

  A scene instance contains a reference to a scene to instantiate and
  the transformation to instantiate the scene with. An implementation
  will typically transform the ray with the inverse of the provided
  transformation and continue traversing the ray through the provided
  scene. If any geometry is hit, the instance ID (instID) member of
  the ray will get set to the geometry ID of the instance. */
RTCORE_API unsigned rtcNewInstance (RTCScene target,                  //!< the scene the instance belongs to
                                    RTCScene source                   //!< the scene to instantiate
  );

/*! \brief Sets transformation of the instance */
RTCORE_API void rtcSetTransform (RTCScene scene,                          //!< scene handle
                                 unsigned geomID,                         //!< ID of geometry
                                 RTCMatrixType layout,                    //!< layout of transformation matrix
                                 const float* xfm                         //!< transformation matrix
                                 );

/*! \brief Creates a new triangle mesh. The number of triangles
  (numTriangles), number of vertices (numVertices), and number of time
  steps (1 for normal meshes, and 2 for linear motion blur), have to
  get specified. The triangle indices can be set be mapping and
  writing to the index buffer (RTC_INDEX_BUFFER) and the triangle
  vertices can be set by mapping and writing into the vertex buffer
  (RTC_VERTEX_BUFFER). In case of linear motion blur, two vertex
  buffers have to get filled (RTC_VERTEX_BUFFER0, RTC_VERTEX_BUFFER1),
  one for each time step. The index buffer has the default layout of
  three 32 bit integer indices for each triangle. An index points to
  the ith vertex. The vertex buffer stores single precision x,y,z
  floating point coordinates aligned to 16 bytes. The value of the 4th
  float used for alignment can be arbitrary. */
RTCORE_API unsigned rtcNewTriangleMesh (RTCScene scene,                    //!< the scene the mesh belongs to
                                        RTCGeometryFlags flags,            //!< geometry flags
                                        size_t numTriangles,               //!< number of triangles
                                        size_t numVertices,                //!< number of vertices
                                        size_t numTimeSteps = 1            //!< number of motion blur time steps
  );

/*! \brief Sets 32 bit ray mask. */
RTCORE_API void rtcSetMask (RTCScene scene, unsigned geomID, int mask);

/*! \brief Maps specified buffer. This function can be used to set index and
 *  vertex buffers of geometries. */
RTCORE_API void* rtcMapBuffer(RTCScene scene, unsigned geomID, RTCBufferType type);

/*! \brief Unmaps specified buffer. 

  A buffer has to be unmapped before the rtcEnable, rtcDisable,
  rtcUpdate, or rtcDeleteGeometry calls are executed. */
RTCORE_API void rtcUnmapBuffer(RTCScene scene, unsigned geomID, RTCBufferType type);

/*! \brief Shares a data buffer between the application and
 *  Embree. The passed buffer is used by Embree to store index and
 *  vertex data. It has to remain valid as long as the mesh exists,
 *  and the user is responsible to free the data when the mesh gets
 *  deleted. One can optionally speficy a byte offset and byte stride
 *  of the elements stored inside the buffer. The addresses
 *  ptr+offset+i*stride have to be aligned to 4 bytes on Xeon CPUs and
 *  16 bytes on Xeon Phi accelerators. For vertex buffers, the 4 bytes
 *  after the z-coordinate of the last vertex have to be readable memory,
 *  thus padding is required for some layouts. If this function is not
 *  called, Embree will allocate and manage buffers of the default
 *  layout. */
RTCORE_API void rtcSetBuffer(RTCScene scene, unsigned geomID, RTCBufferType type, 
                             void* ptr, size_t offset = 0, size_t stride = 16);

/*! \brief Enable geometry. Enabled geometry can be hit by a ray. */
RTCORE_API void rtcEnable (RTCScene scene, unsigned geomID);

/*! \brief Update geometry. 

  This function has to get called, each time the user modifies some
  geometry for dynamic scenes. The function does not have to get
  called after initializing some geometry for the first time. */
RTCORE_API void rtcUpdate (RTCScene scene, unsigned geomID);

/*! \brief Disable geometry. 

  Disabled geometry is not hit by any ray. Disabling and enabling
  geometry gives higher performance than deleting and recreating
  geometry. */
RTCORE_API void rtcDisable (RTCScene scene, unsigned geomID);

/*! \brief Sets the intersection filter function for single rays. */
RTCORE_API void rtcSetIntersectionFilterFunction (RTCScene scene, unsigned geomID, RTCFilterFunc func);

/*! \brief Sets the intersection filter function for ray packets of size 4. */
RTCORE_API void rtcSetIntersectionFilterFunction4 (RTCScene scene, unsigned geomID, RTCFilterFunc4 func);

/*! \brief Sets the intersection filter function for ray packets of size 8. */
RTCORE_API void rtcSetIntersectionFilterFunction8 (RTCScene scene, unsigned geomID, RTCFilterFunc8 func);

/*! \brief Sets the intersection filter function for ray packets of size 16. */
RTCORE_API void rtcSetIntersectionFilterFunction16 (RTCScene scene, unsigned geomID, RTCFilterFunc16 func);

/*! \brief Sets the occlusion filter function for single rays. */
RTCORE_API void rtcSetOcclusionFilterFunction (RTCScene scene, unsigned geomID, RTCFilterFunc func);

/*! \brief Sets the occlusion filter function for ray packets of size 4. */
RTCORE_API void rtcSetOcclusionFilterFunction4 (RTCScene scene, unsigned geomID, RTCFilterFunc4 func);

/*! \brief Sets the occlusion filter function for ray packets of size 8. */
RTCORE_API void rtcSetOcclusionFilterFunction8 (RTCScene scene, unsigned geomID, RTCFilterFunc8 func);

/*! \brief Sets the occlusion filter function for ray packets of size 16. */
RTCORE_API void rtcSetOcclusionFilterFunction16 (RTCScene scene, unsigned geomID, RTCFilterFunc16 func);

/*! \brief Deletes the geometry. */
RTCORE_API void rtcDeleteGeometry (RTCScene scene, unsigned geomID);

/*! @} */

#endif
