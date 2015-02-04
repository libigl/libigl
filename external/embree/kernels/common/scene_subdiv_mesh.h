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
#include "algorithms/pmap.h"
#include "algorithms/pset.h"

namespace embree
{
  class SubdivMesh : public Geometry
  {
    ALIGNED_CLASS;
  public:

    static const size_t MAX_VALENCE = 5;                //!< maximal number of vertices of a patch
    static const size_t MAX_RING_FACE_VALENCE = 32;     //!< maximal number of faces per ring
    static const size_t MAX_RING_EDGE_VALENCE = 2*32;   //!< maximal number of edges per ring
    static const GeometryTy geom_type = SUBDIV_MESH;

    struct Edge 
    {
      /*! edge constructor */
      __forceinline Edge(const uint32 v0, const uint32 v1)
	: v0(v0), v1(v1) {}

      /*! create an 64 bit identifier that is unique for the not oriented edge */
      __forceinline operator uint64() const       
      {
	uint32 p0 = v0, p1 = v1;
	if (p0<p1) std::swap(p0,p1);
	return (((uint64)p0) << 32) | (uint64)p1;
      }

    public:
      uint32 v0,v1;    //!< start and end vertex of the edge
    };

    class __aligned(32) HalfEdge
    {
      friend class SubdivMesh;
    public:

      HalfEdge () 
        : vtx_index(-1), next_half_edge_ofs(0), prev_half_edge_ofs(0), opposite_half_edge_ofs(0), edge_crease_weight(0), 
          vertex_crease_weight(0), edge_level(0), align(0) 
	{
	  assert(sizeof(HalfEdge) == 32);
	}

      __forceinline bool hasOpposite() const { return opposite_half_edge_ofs != 0; }
      __forceinline void setOpposite(HalfEdge* opposite) { opposite_half_edge_ofs = opposite-this; }

      __forceinline       HalfEdge* next()       { assert( next_half_edge_ofs != 0 ); return &this[next_half_edge_ofs]; }
      __forceinline const HalfEdge* next() const { assert( next_half_edge_ofs != 0 ); return &this[next_half_edge_ofs]; }

      __forceinline       HalfEdge* prev()       { assert( prev_half_edge_ofs != 0 ); return &this[prev_half_edge_ofs]; }
      __forceinline const HalfEdge* prev() const { assert( prev_half_edge_ofs != 0 ); return &this[prev_half_edge_ofs]; }

      __forceinline       HalfEdge* opposite()       { assert( opposite_half_edge_ofs != 0 ); return &this[opposite_half_edge_ofs]; }
      __forceinline const HalfEdge* opposite() const { assert( opposite_half_edge_ofs != 0 ); return &this[opposite_half_edge_ofs]; }

      __forceinline       HalfEdge* rotate()       { return opposite()->next(); }
      __forceinline const HalfEdge* rotate() const { return opposite()->next(); }

      __forceinline unsigned int getStartVertexIndex() const { return vtx_index; }
      __forceinline unsigned int getEndVertexIndex  () const { return next()->vtx_index; }

      /*! tests if the start vertex of the edge is regular */
      __forceinline bool isRegularVertex() const 
      {
	const HalfEdge* p = this;
	if (p->hasCreases()) return false;
	
        if (!p->hasOpposite()) return false;
        if ((p = p->rotate()) == this) return false;

        if (!p->hasOpposite()) return false;
        if ((p = p->rotate()) == this) return false;

        if (!p->hasOpposite()) return false;
        if ((p = p->rotate()) == this) return false;

        if (!p->hasOpposite()) return false;
        if ((p = p->rotate()) != this) return false;

        return true;
      }

      /*! tests if the face is a regular b-spline face */
      __forceinline bool isRegularFace() const 
      {
	const HalfEdge* p = this;

        if (!p->isRegularVertex()) return false;
        if ((p = p->next()) == this) return false;

        if (!p->isRegularVertex()) return false;
        if ((p = p->next()) == this) return false;

        if (!p->isRegularVertex()) return false;
        if ((p = p->next()) == this) return false;
        
        if (!p->isRegularVertex()) return false;
        if ((p = p->next()) != this) return false;

	return true;
      }

      /*! tests if the face is a no regular b-spline face (0=regular, 1=irregular quad, 2=non-quad) */
      __forceinline int noRegularFace() const 
      {
	int ret = 0;
	const HalfEdge* p = this;

        if (!p->isRegularVertex()) ret = 1;
        if ((p = p->next()) == this) return 2;

        if (!p->isRegularVertex()) ret = 1;
        if ((p = p->next()) == this) return 2;

        if (!p->isRegularVertex()) ret = 1;
        if ((p = p->next()) == this) return 2;
        
        if (!p->isRegularVertex()) ret = 1;
        if ((p = p->next()) != this) return 2;

	return ret;
      }

      /*! tests if the face can be diced (using bspline or gregory patch) */
      __forceinline bool isGregoryFace() const 
      {
	const HalfEdge* p = this;

        if (p->hasCreases()) return false;
        if ((p = p->next()) == this) return false;

        if (p->hasCreases()) return false;
        if ((p = p->next()) == this) return false;

        if (p->hasCreases()) return false;
        if ((p = p->next()) == this) return false;
        
        if (p->hasCreases()) return false;
        if ((p = p->next()) != this) return false;

	return true;
      }

      /*! calculates conservative bounds of a catmull clark subdivision face */
      __forceinline BBox3fa bounds(const BufferT<Vec3fa>& vertices) const
      {
        BBox3fa bounds = this->get1RingBounds(vertices);
        for (const HalfEdge* p=this->next(); p!=this; p=p->next())
          bounds.extend(p->get1RingBounds(vertices));
        return bounds;
      }

      /*! tests if this is a valid patch */
      __forceinline bool valid(const BufferT<Vec3fa>& vertices) const
      {
        size_t N = 1;
        if (!this->validRing(vertices)) return false;
        for (const HalfEdge* p=this->next(); p!=this; p=p->next(), N++) {
          if (!p->validRing(vertices)) return false;
        }
        return N >= 3 && N <= MAX_VALENCE;
      }

      /*! stream output */
      friend __forceinline std::ostream &operator<<(std::ostream &o, const SubdivMesh::HalfEdge &h)
      {
        return o << "{ " << 
          "vertex = " << h.vtx_index << ", " << //" -> " << h.next()->vtx_index << ", " << 
          "prev = " << h.prev_half_edge_ofs << ", " << 
          "next = " << h.next_half_edge_ofs << ", " << 
          "opposite = " << h.opposite_half_edge_ofs << ", " << 
          //"edge_crease = " << h.edge_crease_weight << ", " << 
          //"vertex_crease = " << h.vertex_crease_weight << ", " << 
          //"edge_level = " << h.edge_level << 
          " }";
      } 

    private:

      /*! calculates the bounds of the face associated with the half-edge */
      __forceinline BBox3fa getFaceBounds(const BufferT<Vec3fa>& vertices) const 
      {
        BBox3fa b = vertices[getStartVertexIndex()];
        for (const HalfEdge* p = next(); p!=this; p=p->next()) {
          b.extend(vertices[p->getStartVertexIndex()]);
        }
        return b;
      }

      /*! calculates the bounds of the 1-ring associated with the vertex of the half-edge */
      __forceinline BBox3fa get1RingBounds(const BufferT<Vec3fa>& vertices) const 
      {
        BBox3fa bounds = empty;
        const HalfEdge* p = this;
        do 
        {
          /* calculate bounds of current face */
          bounds.extend(p->getFaceBounds(vertices));
          p = p->prev();
          
          /* continue with next face */
          if (likely(p->hasOpposite())) 
            p = p->opposite();
          
          /* if there is no opposite go the long way to the other side of the border */
          else {
            p = this;
            while (p->hasOpposite()) 
              p = p->opposite()->next();
          }
          
        } while (p != this); 
        
        return bounds;
      }

      /*! tests if this is a valid face */
      __forceinline bool validFace(const BufferT<Vec3fa>& vertices, size_t& N) const 
      {
        const Vec3fa v = vertices[getStartVertexIndex()];
        if (!inFloatRange(v)) return false;
        size_t n = 1;
        for (const HalfEdge* p = next(); p!=this; p=p->next(), n++) {
          const Vec3fa v = vertices[p->getStartVertexIndex()];
          if (!inFloatRange(v)) return false;
        }
        N += n-2;
        return n >= 3 && n <= MAX_VALENCE;
      }

      /*! tests if this is a valid ring */
      __forceinline bool validRing(const BufferT<Vec3fa>& vertices) const 
      {
        size_t faceValence = 0;
        size_t edgeValence = 0;

        BBox3fa bounds = empty;
        const HalfEdge* p = this;
        do 
        {
          /* calculate bounds of current face */
          if (!p->validFace(vertices,edgeValence)) 
            return false;
          
          faceValence++;
          p = p->prev();
          
          /* continue with next face */
          if (likely(p->hasOpposite())) 
            p = p->opposite();
          
          /* if there is no opposite go the long way to the other side of the border */
          else {
            faceValence++;
            edgeValence++;
            p = this;
            while (p->hasOpposite()) 
              p = p->opposite()->next();
          }
          
        } while (p != this); 
        
        return faceValence <= MAX_RING_FACE_VALENCE && edgeValence <= MAX_RING_EDGE_VALENCE;
      }

      /*! tests if the edge has creases */
      __forceinline bool hasCreases() const {
	return max(edge_crease_weight,vertex_crease_weight) != 0.0f;
      }

    private:
      unsigned int vtx_index;         //!< index of edge start vertex
      int next_half_edge_ofs;         //!< relative offset to next half edge of face
      int prev_half_edge_ofs;         //!< relative offset to previous half edge of face
      int opposite_half_edge_ofs;     //!< relative offset to opposite half edge

    public:
      float edge_crease_weight;       //!< crease weight attached to edge
      float vertex_crease_weight;     //!< crease weight attached to start vertex
      float edge_level;               //!< subdivision factor for edge
      float align;                    //!< aligns the structure to 32 bytes
    };

    struct KeyHalfEdge 
    {
      KeyHalfEdge() {}
      
      KeyHalfEdge (uint64 key, HalfEdge* edge) 
      : key(key), edge(edge) {}
      
      __forceinline operator uint64() const { 
	return key; 
      }

      friend __forceinline bool operator<(const KeyHalfEdge& e0, const KeyHalfEdge& e1) {
        return e0.key < e1.key;
      }
      
    public:
      uint64 key;
      HalfEdge* edge;
    };

  public:
    SubdivMesh(Scene* parent, RTCGeometryFlags flags, size_t numFaces, size_t numEdges, size_t numVertices, 
               size_t numCreases, size_t numCorners, size_t numHoles, size_t numTimeSteps);
    ~SubdivMesh();

    void enabling();
    void disabling();
    void setMask (unsigned mask);
    void setBuffer(RTCBufferType type, void* ptr, size_t offset, size_t stride);
    void* map(RTCBufferType type);
    void unmap(RTCBufferType type);
    void update ();
    void updateBuffer (RTCBufferType type);
    void setUserData (void* ptr, bool ispc);
    void immutable ();
    bool verify ();
    void setDisplacementFunction (RTCDisplacementFunc func, RTCBounds* bounds);

  public:

    /*! return the number of faces */
    size_t size() const { 
      return numFaces; 
    };

    /*! calculates the bounds of the i'th subdivision patch */
    __forceinline BBox3fa bounds(size_t i) const {
      return halfEdges[faceStartEdge[i]].bounds(vertices[0]);
    }

    /*! check if the i'th primitive is valid */
    __forceinline bool valid(size_t i, BBox3fa* bbox = NULL) const {
      if (bbox) *bbox = bounds(i);
      return getHalfEdge(i)->valid(vertices[0]) && !holeSet.lookup(i);
    }

    /*! initializes the half edge data structure */
    void initializeHalfEdgeStructures ();

  private:

    /*! recalculates the half edges */
    void calculateHalfEdges();

    /*! updates half edges when recalculation is not necessary */
    void updateHalfEdges();

  public:
    /*! returns the start half edge for some face */
    __forceinline const HalfEdge* getHalfEdge ( const size_t f ) const { 
      return &halfEdges[faceStartEdge[f]]; 
    }    

    /*! returns the vertex buffer for some time step */
    __forceinline const BufferT<Vec3fa>& getVertexBuffer( const size_t t = 0 ) const {
      return vertices[t];
    }
    
    //__forceinline const Vec3fa*   getVertexPositionPtr( const unsigned t = 0 ) const { return (Vec3fa*)vertices[t].getPtr(); } // FIXME: this function should never get used, always pass BufferT<Vec3fa> object

  public:
    unsigned int mask;                //!< for masking out geometry
    unsigned int numTimeSteps;        //!< number of time steps (1 or 2)  

    RTCDisplacementFunc displFunc;    //!< displacement function
    BBox3fa             displBounds;  //!< bounds for maximal displacement 

  private:
    size_t numFaces;                  //!< number of faces
    size_t numEdges;                  //!< number of edges
    size_t numVertices;               //!< number of vertices
    size_t numHalfEdges;              //!< number of half edges used by faces

    /*! all buffers in this section are provided by the application */
  private:
    
    /*! buffer containing the number of vertices for each face */
    BufferT<int> faceVertices;

    /*! indices of the vertices composing each face */
    BufferT<unsigned> vertexIndices;

    /*! vertex buffer (one buffer for each time step) */
    BufferT<Vec3fa> vertices[2];

    /*! edge crease buffer containing edges (pairs of vertices) that carry edge crease weights */
    BufferT<Edge> edge_creases;

    /*! edge crease weights for each edge of the edge_creases buffer */
    BufferT<float> edge_crease_weights;

    /*! vertex crease buffer containing all vertices that carry vertex crease weights */
    BufferT<unsigned> vertex_creases;

    /*! vertex crease weights for each vertex of the vertex_creases buffer */
    BufferT<float> vertex_crease_weights;

    /*! subdivision level for each half edge of the vertexIndices buffer */
    BufferT<float> levels;

    /*! buffer that marks specific faces as holes */
    BufferT<unsigned> holes;

    /*! all data in this section if generated by initializeHalfEdgeStructures function */
  private:

    /*! fast lookup table to find the first half edge for some face */
    std::vector<uint32> faceStartEdge;

    /*! Half edge structure. */
    std::vector<HalfEdge> halfEdges;

    /*! set with all holes */
    pset<uint32> holeSet;

    /*! the following data is only required during construction of the
     *  half edge structure and can be cleared for static scenes */
  private:

    /*! two arrays used to sort the half edges */
    std::vector<KeyHalfEdge> halfEdges0;
    std::vector<KeyHalfEdge> halfEdges1;

    /*! map with all vertex creases */
    pmap<uint32,float> vertexCreaseMap;

    /*! map with all edge creases */
    pmap<uint64,float> edgeCreaseMap;

    /*! flag whether only the edge levels have changed and the mesh has no creases,
     *  allows for simple bvh update instead of full rebuild in cached mode */
    bool levelUpdate;

  public:
    /* check for simple edge level update */
    __forceinline bool checkLevelUpdate() { return levelUpdate; }

  };
};
