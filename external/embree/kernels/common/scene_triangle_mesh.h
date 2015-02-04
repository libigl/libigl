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

#pragma once

#include "common/geometry.h"
#include "common/buffer.h"

namespace embree
{
  /*! Triangle Mesh */
  struct TriangleMesh : public Geometry
  {
    static const GeometryTy geom_type = TRIANGLE_MESH;

    struct Triangle {
      unsigned int v[3];
    };
    
  public:
    TriangleMesh (Scene* parent, RTCGeometryFlags flags, size_t numTriangles, size_t numVertices, size_t numTimeSteps); 
  
    void write(std::ofstream& file);

    /* geometry interface */
  public:
    void enabling();
    void disabling();
    void setMask (unsigned mask);
    void setBuffer(RTCBufferType type, void* ptr, size_t offset, size_t stride);
    void* map(RTCBufferType type);
    void unmap(RTCBufferType type);
    void setUserData (void* ptr, bool ispc);
    void immutable ();
    bool verify ();

  public:

    /*! returns number of triangles */
    __forceinline size_t size() const {
      return numTriangles;
    }
    
    /*! returns i'th triangle*/
    __forceinline const Triangle& triangle(size_t i) const {
      assert(i < numTriangles);
      return triangles[i];
    }

    /*! returns i'th vertex of j'th timestep */
    __forceinline const Vec3fa vertex(size_t i, size_t j = 0) const 
    {
      assert(i < numVertices);
      assert(j < numTimeSteps);
      return vertices[j][i];
    }

    /*! returns i'th vertex of j'th timestep */
    __forceinline const char* vertexPtr(size_t i, size_t j = 0) const 
    {
      assert(i < numVertices);
      assert(j < numTimeSteps);
      return vertices[j].getPtr(i);
    }
    
    /*! returns the stride in bytes of the triangle buffer */
    __forceinline size_t getTriangleBufferStride() const {
      return triangles.getBufferStride();
    }

    /*! returns the stride in butes of the vertex buffer */
    __forceinline size_t getVertexBufferStride() const {
      return vertices[0].getBufferStride();
    }

    /*! check if the i'th primitive is valid */
    __forceinline bool valid(size_t i, BBox3fa* bbox = NULL) const 
    {
      const Triangle& tri = triangle(i);
      if (tri.v[0] >= numVertices) return false;
      if (tri.v[1] >= numVertices) return false;
      if (tri.v[2] >= numVertices) return false;

      for (size_t j=0; j<numTimeSteps; j++) {
	const Vec3fa v0 = vertex(tri.v[0],j);
	const Vec3fa v1 = vertex(tri.v[1],j);
	const Vec3fa v2 = vertex(tri.v[2],j);
	if (!inFloatRange(v0) || !inFloatRange(v1) || !inFloatRange(v2))
	  return false;
      }

      if (bbox) {
	const Vec3fa v0 = vertex(tri.v[0]);
	const Vec3fa v1 = vertex(tri.v[1]);
	const Vec3fa v2 = vertex(tri.v[2]);
	*bbox = BBox3fa(min(v0,v1,v2),max(v0,v1,v2));
      }
      return true;
    }

    /*! calculates the bounds of the i'th triangle */
    __forceinline BBox3fa bounds(size_t i) const 
    {
      const Triangle& tri = triangle(i);
      const Vec3fa v0 = vertex(tri.v[0]);
      const Vec3fa v1 = vertex(tri.v[1]);
      const Vec3fa v2 = vertex(tri.v[2]);
      return BBox3fa(min(v0,v1,v2),max(v0,v1,v2));
    }

#if defined(__MIC__)


    template<unsigned int HINT=0>
      __forceinline mic3f getTriangleVertices(const Triangle &tri,const size_t dim=0) const 
      {
	assert( tri.v[0] < numVertices );
	assert( tri.v[1] < numVertices );
	assert( tri.v[2] < numVertices );

#if !defined(RTCORE_BUFFER_STRIDE)
	
	const float *__restrict__ const vptr0 = (float*) vertices[dim].getPtr(tri.v[0]);
	const float *__restrict__ const vptr1 = (float*) vertices[dim].getPtr(tri.v[1]);
	const float *__restrict__ const vptr2 = (float*) vertices[dim].getPtr(tri.v[2]);

	const mic_f v0 = broadcast4to16f(vptr0); 
	const mic_f v1 = broadcast4to16f(vptr1); 
	const mic_f v2 = broadcast4to16f(vptr2); 
	return mic3f(v0,v1,v2);
#else
	const mic_i stride = vertices[dim].getBufferStride();

	const mic_i offset0_64 = mul_uint64(stride,mic_i(tri.v[0]));
	const mic_i offset1_64 = mul_uint64(stride,mic_i(tri.v[1]));
	const mic_i offset2_64 = mul_uint64(stride,mic_i(tri.v[2]));

	const char  *__restrict__ const base  = vertices[dim].getPtr();
	const size_t off0 = offset0_64.uint64(0);
	const size_t off1 = offset1_64.uint64(0);
	const size_t off2 = offset2_64.uint64(0);

	const float *__restrict__ const vptr0_64 = (float*)(base + off0);
	const float *__restrict__ const vptr1_64 = (float*)(base + off1);
	const float *__restrict__ const vptr2_64 = (float*)(base + off2);

	if (HINT)
	{
	  prefetch<HINT>(vptr1_64);
	  prefetch<HINT>(vptr2_64);
	}

	assert( vptr0_64 == (float*)vertexPtr(tri.v[0],dim) );
	assert( vptr1_64 == (float*)vertexPtr(tri.v[1],dim) );
	assert( vptr2_64 == (float*)vertexPtr(tri.v[2],dim) );
	
	const mic_m m_3f = 0x7;
	const mic_f v0 = permute<0,0,0,0>(uload16f(m_3f,vptr0_64));
	const mic_f v1 = permute<0,0,0,0>(uload16f(m_3f,vptr1_64));
	const mic_f v2 = permute<0,0,0,0>(uload16f(m_3f,vptr2_64));
	 //FIXME: there should be no need to zero the last component

	return mic3f(select(0x7777,v0,mic_f::zero()),select(0x7777,v1,mic_f::zero()),select(0x7777,v2,mic_f::zero()));
#endif	
	
	
      }
    
#endif
    
  public:
    unsigned int mask;                //!< for masking out geometry
    unsigned int numTimeSteps;        //!< number of time steps (1 or 2)
    
    BufferT<Triangle> triangles;      //!< array of triangles
    size_t numTriangles;              //!< number of triangles
    
    BufferT<Vec3fa> vertices[2];    //!< vertex array
    size_t numVertices;               //!< number of vertices
  };

  __forceinline std::ostream &operator<<(std::ostream &o, const TriangleMesh::Triangle &t)
  {
    o << "tri " << t.v[0] << " " << t.v[1] << " " << t.v[2] << std::endl;
    return o;
  } 

}
