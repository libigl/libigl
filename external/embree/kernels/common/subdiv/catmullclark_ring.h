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
#include "common/scene_subdiv_mesh.h"

namespace embree
{
  struct __aligned(64) FinalQuad {
    Vec3fa vtx[4];
  };
  
  struct __aligned(64) CatmullClark1Ring
  {
    static const size_t MAX_FACE_VALENCE = SubdivMesh::MAX_RING_FACE_VALENCE;
    static const size_t MAX_EDGE_VALENCE = SubdivMesh::MAX_RING_EDGE_VALENCE;
    
    Vec3fa ring [MAX_EDGE_VALENCE];
    float crease_weight[MAX_FACE_VALENCE];
    
    int border_index;
    Vec3fa vtx;
    unsigned int face_valence;
    unsigned int edge_valence;
    float vertex_crease_weight;
    float vertex_level; // maximal level of all adjacent edges
    float edge_level; // level of first edge
    bool noForcedSubdivision; // varying edge crease weight stitching fix

  public:
    CatmullClark1Ring () {}

    __forceinline bool hasBorder() const {
      return border_index != -1;
    }
    
    __forceinline const Vec3fa& front(size_t i) const {
      assert(edge_valence>i);
      return ring[i];
    }
    
    __forceinline const Vec3fa& back(size_t i) const {
      assert(i>0 && edge_valence>=i);
      return ring[edge_valence-i];
    }
    
    __forceinline bool has_last_face() const {
      return border_index != edge_valence-2;
    }
    
    __forceinline bool has_second_face() const {
      return (border_index == -1) || (border_index >= 4);
    }
    
    __forceinline Vec3fa regular_border_vertex_4() const 
    {
      assert(border_index != -1);
      assert(face_valence == 2 || face_valence == 3);
      if (face_valence == 3 && border_index == 4) return ring[4];
      else return 2.0f*vtx-ring[0];
    }

    __forceinline Vec3fa regular_border_vertex_5() const 
    {
      assert(border_index != -1);
      assert(face_valence == 2 || face_valence == 3);
      if (face_valence == 2) return 2.0f*vtx-ring[1];
      else if (border_index == 2) return 2.0f*ring[4]-ring[5];
      else { assert(border_index == 4); return 2.0f*ring[4]-ring[3]; }
    }

    __forceinline Vec3fa regular_border_vertex_6() const 
    {
      assert(border_index != -1);
      assert(face_valence == 2 || face_valence == 3);
      if (face_valence == 3 && border_index == 2) return ring[4];
      else return 2.0f*vtx-ring[2];
    }

    __forceinline BBox3fa bounds() const
    {
      BBox3fa bounds ( vtx );
      for (size_t i = 0; i<edge_valence ; i++)
	bounds.extend( ring[i] );
      return bounds;
    }
    
    __forceinline void init(const SubdivMesh::HalfEdge* const h, const BufferT<Vec3fa>& vertices) 
    {
      noForcedSubdivision = true;
      border_index = -1;
      vtx = (Vec3fa_t)vertices[ h->getStartVertexIndex() ];
      vertex_crease_weight = h->vertex_crease_weight;
      
      SubdivMesh::HalfEdge* p = (SubdivMesh::HalfEdge*) h;
      
      size_t i=0;
      crease_weight[i/2] = p->edge_crease_weight;
      edge_level = vertex_level = p->edge_level;
      if (!p->hasOpposite()) crease_weight[i/2] = inf;
      
      do // FIXME: assertions!?
      {
        /* store first two vertices of face */
        p = p->next();
        ring[i++] = (Vec3fa_t) vertices[ p->getStartVertexIndex() ];
        p = p->next();
        ring[i++] = (Vec3fa_t) vertices[ p->getStartVertexIndex() ];
        p = p->next();
        crease_weight[i/2] = p->edge_crease_weight;
	vertex_level = max(vertex_level,p->edge_level);
	
        /* continue with next face */
        if (likely(p->hasOpposite())) 
          p = p->opposite();
        
        /* if there is no opposite go the long way to the other side of the border */
        else
        {
          /*! mark first border edge and store dummy vertex for face between the two border edges */
          border_index = i;
          crease_weight[i/2] = inf; 
          ring[i++] = (Vec3fa_t) vertices[ p->getStartVertexIndex() ];
          ring[i++] = vtx; // dummy vertex
          crease_weight[i/2] = inf;
	  
          /*! goto other side of border */
          p = (SubdivMesh::HalfEdge*) h;
          while (p->hasOpposite()) 
            p = p->opposite()->next();
        }
	
      } while (p != h); 
      
      edge_valence = i;
      face_valence = i >> 1;
    }
    
    
    __forceinline void subdivide(CatmullClark1Ring& dest) const
    {
      dest.noForcedSubdivision = true;
      dest.edge_level = 0.5f*edge_level;
      dest.vertex_level = 0.5f*vertex_level;
      dest.face_valence         = face_valence;
      dest.edge_valence         = edge_valence;
      dest.border_index = border_index;
      dest.vertex_crease_weight   = max(0.0f,vertex_crease_weight-1.0f);
      
      /* calculate face points */
      Vec3fa_t S = Vec3fa_t(0.0f);
      for (size_t i=0; i<face_valence-1; i++) {
        S += dest.ring[2*i+1] = ((vtx + ring[2*i]) + (ring[2*i+1] + ring[2*i+2])) * 0.25f;
      }
      S += dest.ring[edge_valence-1] = ((vtx + ring[edge_valence-2]) + (ring[edge_valence-1] + ring[0])) * 0.25f;
      
      /* calculate new edge points */
      size_t num_creases = 0;
      size_t crease_id[MAX_FACE_VALENCE];
      Vec3fa_t C = Vec3fa_t(0.0f);
      for (size_t i=1; i<face_valence; i++)
      {
        const Vec3fa_t v = vtx + ring[2*i];
        const Vec3fa_t f = dest.ring[2*i-1] + dest.ring[2*i+1];
        S += ring[2*i];
        dest.crease_weight[i] = max(crease_weight[i]-1.0f,0.0f);
        //dest.crease_weight[i] = crease_weight[i] < 1.0f ? 0.0f : 0.5f*crease_weight[i];
	dest.noForcedSubdivision &= crease_weight[i] == 0.0f || crease_weight[i] > 1.0f;
        
        /* fast path for regular edge points */
        if (likely(crease_weight[i] <= 0.0f)) {
          dest.ring[2*i] = (v+f) * 0.25f;
        }
        
        /* slower path for hard edge rule */
        else {
          C += ring[2*i]; crease_id[num_creases++] = i;
          dest.ring[2*i] = v*0.5f;
	  
          /* even slower path for blended edge rule */
          if (unlikely(crease_weight[i] < 1.0f)) {
            const float w0 = crease_weight[i], w1 = 1.0f-w0;
            dest.ring[2*i] = w1*((v+f)*0.25f) + w0*(v*0.5f);
          }
        }
      }
      {
        const size_t i=0;
        const Vec3fa_t v = vtx + ring[2*i];
        const Vec3fa_t f = dest.ring[edge_valence-1] + dest.ring[2*i+1];
        S += ring[2*i];
        dest.crease_weight[i] = max(crease_weight[i]-1.0f,0.0f);
        //dest.crease_weight[i] = crease_weight[i] < 1.0f ? 0.0f : 0.5f*crease_weight[i];
	dest.noForcedSubdivision &= crease_weight[i] == 0.0f || crease_weight[i] > 1.0f;

        /* fast path for regular edge points */
        if (likely(crease_weight[i] <= 0.0f)) {
          dest.ring[2*i] = (v+f) * 0.25f;
        }
        
        /* slower path for hard edge rule */
        else {
          C += ring[2*i]; crease_id[num_creases++] = i;
          dest.ring[2*i] = v*0.5f;
	  
          /* even slower path for blended edge rule */
          if (unlikely(crease_weight[i] < 1.0f)) {
            const float w0 = crease_weight[i], w1 = 1.0f-w0;
            dest.ring[2*i] = w1*((v+f)*0.25f) + w0*(v*0.5f);
          }
        }
      }
      
      /* compute new vertex using smooth rule */
      const float inv_face_valence = 1.0f / (float)face_valence;
      const Vec3fa_t v_smooth = (Vec3fa_t)(S*inv_face_valence + (float(face_valence)-2.0f)*vtx)*inv_face_valence;
      dest.vtx = v_smooth;
      
      /* compute new vertex using vertex_crease_weight rule */
      if (unlikely(vertex_crease_weight > 0.0f)) 
      {
        if (vertex_crease_weight >= 1.0f) {
          dest.vtx = vtx;
        } else {
          const float t0 = vertex_crease_weight, t1 = 1.0f-t0;
          dest.vtx = t0*vtx + t1*v_smooth;
        }
        return;
      }
      
      if (likely(num_creases <= 1))
        return;
      
      /* compute new vertex using crease rule */
      if (likely(num_creases == 2)) {
        const Vec3fa_t v_sharp = (Vec3fa_t)(C + 6.0f * vtx) * (1.0f / 8.0f);
        const float crease_weight0 = crease_weight[crease_id[0]];
        const float crease_weight1 = crease_weight[crease_id[1]];
        dest.vtx = v_sharp;
        dest.crease_weight[crease_id[0]] = max(0.25f*(3.0f*crease_weight0 + crease_weight1)-1.0f,0.0f);
        dest.crease_weight[crease_id[1]] = max(0.25f*(3.0f*crease_weight1 + crease_weight0)-1.0f,0.0f);
	dest.noForcedSubdivision = (dest.crease_weight[crease_id[0]] != 0.0f) || (dest.crease_weight[crease_id[1]] != 0.0f);
        //dest.crease_weight[crease_id[0]] = max(0.5f*(crease_weight0 + crease_weight1)-1.0f,0.0f);
        //dest.crease_weight[crease_id[1]] = max(0.5f*(crease_weight1 + crease_weight0)-1.0f,0.0f);
        const float t0 = 0.5f*(crease_weight0+crease_weight1), t1 = 1.0f-t0;
        //dest.crease_weight[crease_id[0]] = t0 < 1.0f ? 0.0f : 0.5f*t0;
        //dest.crease_weight[crease_id[1]] = t0 < 1.0f ? 0.0f : 0.5f*t0;
        if (unlikely(t0 < 1.0f)) {
          dest.vtx = t0*v_sharp + t1*v_smooth;
        }
      }
      
      /* compute new vertex using corner rule */
      else {
        dest.vtx = vtx;
      }
    }
    
    __forceinline bool isRegular() const 
    {
      if (border_index == -1) {
	if (face_valence == 4) return true;
      } else {
	if (face_valence < 4) return true;
      }
      return false;
    }

     /* returns true if the vertex can be part of a dicable B-Spline patch or is a final Quad */
    __forceinline bool isRegularOrFinal(const size_t depth) const 
    {
      if (depth < 10)
      //if (vertex_level > 1.0f) 
      {
	if (border_index == -1) 
	{
	  if (face_valence != 4)
	    return false;
	  if (vertex_crease_weight > 0.0f) 
	    return false;
	} 
	else {
	  if (face_valence == 2 && vertex_crease_weight > 1E5); // FIXME: use inf
	  else if (face_valence == 3 && vertex_crease_weight == 0.0f);
	  else return false;
	}

	for (size_t i=1; i<face_valence; i++)
	  if (crease_weight[i] > 0.0f && (2*i != border_index) && (2*(i-1) != border_index)) 
	    return false;

	if (crease_weight[0] > 0.0f && (2*(face_valence-1) != border_index)) 
	  return false;
	
	if (!noForcedSubdivision)
	  return false;
      }
      return true;
    }

    /* returns true if the vertex can be part of a dicable B-Spline patch or is a final Quad */
    __forceinline bool isRegularOrFinal2(const size_t depth) const 
    {
      //if (depth < 2)
      if (vertex_level > 1.0f) 
      {
	if (border_index == -1) 
	{
	  if (face_valence != 4)
	    return false;
	  if (vertex_crease_weight > 0.0f) 
	    return false;
	} 
	else {
	  if (face_valence == 2 && vertex_crease_weight > 1E5); // FIXME: use inf
	  else if (face_valence == 3 && vertex_crease_weight == 0.0f);
	  else return false;
	}

	for (size_t i=1; i<face_valence; i++)
	  if (crease_weight[i] > 0.0f && (2*i != border_index) && (2*(i-1) != border_index)) 
	    return false;

	if (crease_weight[0] > 0.0f && (2*(face_valence-1) != border_index)) 
	  return false;
	
	if (!noForcedSubdivision)
	  return false;
      }
      return true;
    }

    /* returns true if the vertex can be part of a dicable gregory patch (using gregory patches) */
    __forceinline bool isGregoryOrFinal(const size_t depth) const 
    {
      if (depth < 10)
      //if (vertex_level > 1.0f) 
      {
	if (vertex_crease_weight > 0.0f) 
	  return false;
	
	for (size_t i=1; i<face_valence; i++) 
	  if (crease_weight[i] > 0.0f && (2*i != border_index) && (2*(i-1) != border_index)) 
	    return false;
	  
       	if (crease_weight[0] > 0.0f && (2*(face_valence-1) != border_index)) 
	  return false;

	if (!noForcedSubdivision)
	  return false;
      }
      return true;
    }
    
    /* computes the limit vertex */
    __forceinline Vec3fa getLimitVertex() const
    {
      /* border vertex rule */
      if (unlikely(border_index != -1))
      {
	if (unlikely(std::isinf(vertex_crease_weight)))
	  return vtx;
	
	const unsigned int second_border_index = border_index+2 >= edge_valence ? 0 : border_index+2;
	return (4.0f * vtx + ring[border_index] + ring[second_border_index]) * 1.0f/6.0f;
      }
      
      Vec3fa_t F( 0.0f );
      Vec3fa_t E( 0.0f );
      
      for (size_t i=0; i<face_valence; i++) {
        F += ring[2*i+1];
        E += ring[2*i];
      }
      
      const float n = (float)face_valence;
      return (Vec3fa_t)(n*n*vtx+4*E+F) / ((n+5.0f)*n);      
    }
    
    /* gets limit tangent in the direction of egde vtx -> ring[0] */
    __forceinline Vec3fa getLimitTangent() const 
    {
      /* border vertex rule */
      if (unlikely(border_index != -1))
      {
	if (unlikely(std::isinf(vertex_crease_weight)))
	  return ring[0] - vtx;
	
	if (border_index != edge_valence-2 && face_valence != 2) {
	  return ring[0] - vtx; 
	}
	else
	{
	  const unsigned int second_border_index = border_index+2 >= edge_valence ? 0 : border_index+2;
	  return (ring[second_border_index] - ring[border_index]) * 0.5f;
	}
      }
      
      Vec3fa_t alpha( 0.0f );
      Vec3fa_t beta ( 0.0f );
      
      const float n = (float)face_valence;
      //const float delta = 1.0f / sqrtf(4.0f + cos(M_PI/n)*cos(M_PI/n));
      //const float c0 = 2.0f/n * delta;
      //const float c1 = 1.0f/n * (1.0f - delta*cosf(M_PI/n));
      const float c0 = 1.0f/n * 1.0f / sqrtf(4.0f + cosf(M_PI/n)*cosf(M_PI/n));  
      const float c1 = (1.0f/n + cosf(M_PI/n) * c0); // FIXME: plus or minus
      for (size_t i=0; i<face_valence; i++)
      {
	const float a = c1 * cosf(2.0f*M_PI*i/n);
	const float b = c0 * cosf((2.0f*M_PI*i+M_PI)/n); // FIXME: factor of 2 missing?
	alpha +=  a * ring[2*i];
	beta  +=  b * ring[2*i+1];
      }
      return alpha +  beta;
    }
    
    /* gets limit tangent in the direction of egde vtx -> ring[edge_valence-2] */
    __forceinline Vec3fa getSecondLimitTangent() const 
    {
      /* border vertex rule */
      if (unlikely(border_index != -1))
      {
        if (unlikely(std::isinf(vertex_crease_weight)))
          return ring[2] - vtx;
        
        //if (border_index == 0 && face_valence != 2) {
        if (border_index == edge_valence-2 && face_valence != 2) {
          return ring[2] - vtx;
        }
        else {
          const unsigned int second_border_index = border_index+2 >= edge_valence ? 0 : border_index+2;
          return (ring[border_index] - ring[second_border_index]) * 0.5f;
        }
      }
      
      Vec3fa_t alpha( 0.0f );
      Vec3fa_t beta ( 0.0f );
      const float n = (float)face_valence;
      //const float delta = 1.0f / sqrtf(4.0f + cos(M_PI/n)*cos(M_PI/n));
      //const float c0 = 2.0f/n * delta;
      //const float c1 = 1.0f/n * (1.0f - delta*cosf(M_PI/n));
      const float c0 = 1.0f/n * 1.0f / sqrtf(4.0f + cosf(M_PI/n)*cosf(M_PI/n));  
      const float c1 = (1.0f/n + cosf(M_PI/n) * c0);
      for (size_t i=0; i<face_valence; i++)
      {
	const float a = c1 * cosf(2.0f*M_PI*(float(i)-1.0f)/n);
	const float b = c0 * cosf((2.0f*M_PI*(float(i)-1.0f)+M_PI)/n);
	alpha += a * ring[2*i];
	beta  += b * ring[2*i+1];
      }
      return alpha +  beta;      
    }

    /* gets surface normal */
    const Vec3fa getNormal() const  {
      return cross(getSecondLimitTangent(),getLimitTangent());
    }
    
    /* returns center of the n-th quad in the 1-ring */
    __forceinline Vec3fa getQuadCenter(const size_t index) const
    {
      const Vec3fa_t &p0 = vtx;
      const Vec3fa_t &p1 = ring[2*index+0];
      const Vec3fa_t &p2 = ring[2*index+1];
      const Vec3fa_t &p3 = index == face_valence-1 ? ring[0] : ring[2*index+2];
      const Vec3fa p = (p0+p1+p2+p3) * 0.25f;
      return p;
    }
    
    /* returns center of the n-th edge in the 1-ring */
    __forceinline Vec3fa getEdgeCenter(const size_t index) const {
      return (vtx + ring[index*2]) * 0.5f;
    }
    
    friend __forceinline std::ostream &operator<<(std::ostream &o, const CatmullClark1Ring &c)
    {
      o << "vtx " << c.vtx << " size = " << c.edge_valence << ", " << 
	"hard_edge = " << c.border_index << ", face_valence " << c.face_valence << 
	", edge_level = " << c.edge_level << ", vertex_level = " << c.vertex_level << ", ring: " << std::endl;
      
      for (size_t i=0; i<c.edge_valence; i++) {
        o << i << " -> " << c.ring[i];
        if (i % 2 == 0) o << " crease = " << c.crease_weight[i/2];
        o << std::endl;
      }
      return o;
    } 
  };
  
  struct __aligned(64) GeneralCatmullClark1Ring
  {
    static const size_t MAX_FACE_VALENCE = SubdivMesh::MAX_RING_FACE_VALENCE;
    static const size_t MAX_EDGE_VALENCE = SubdivMesh::MAX_RING_EDGE_VALENCE;
    
    Vec3fa vtx;
    Vec3fa ring[MAX_EDGE_VALENCE]; 
    int face_size[MAX_FACE_VALENCE];       // number of vertices-2 of nth face in ring
    float crease_weight[MAX_FACE_VALENCE]; // FIXME: move into 4th component of ring entries
    unsigned int face_valence;
    unsigned int edge_valence;
    int border_face;
    float vertex_crease_weight;
    float vertex_level;                      //!< maximal level of adjacent edges
    float edge_level; // level of first edge
    bool only_quads;

    GeneralCatmullClark1Ring() {}
    
    __forceinline bool has_last_face() const {
      return border_face != face_valence-1;
    }
    
    __forceinline bool has_second_face() const {
      return (border_face == -1) || (border_face >= 2);
    }
    
    __forceinline void init(const SubdivMesh::HalfEdge* const h, const BufferT<Vec3fa>& vertices)
    {
      only_quads = true;
      border_face = -1;
      vtx = (Vec3fa_t)vertices[ h->getStartVertexIndex() ];
      vertex_crease_weight = h->vertex_crease_weight;
      SubdivMesh::HalfEdge* p = (SubdivMesh::HalfEdge*) h;
      
      size_t e=0, f=0;
      crease_weight[f] = p->edge_crease_weight;
      edge_level = vertex_level = p->edge_level;
      if (!p->hasOpposite()) crease_weight[f] = inf;
      
      do 
      {
	/* store first N-2 vertices of face */
	size_t vn = 0;
	SubdivMesh::HalfEdge* p_prev = p->prev();
        for (SubdivMesh::HalfEdge* v = p->next(); v!=p_prev; v=v->next()) {
          assert(e < MAX_EDGE_VALENCE);
          ring[e++] = (Vec3fa_t) vertices[ v->getStartVertexIndex() ];
	  vn++;
	}
	assert(f < MAX_FACE_VALENCE);
	face_size[f] = vn;
	only_quads &= (vn == 2);
	p = p_prev;
	crease_weight[++f] = p->edge_crease_weight;
	vertex_level = max(vertex_level,p->edge_level);
	
        /* continue with next face */
        if (likely(p->hasOpposite())) 
          p = p->opposite();
        
        /* if there is no opposite go the long way to the other side of the border */
        else
        {
          /*! mark first border edge and store dummy vertex for face between the two border edges */
	  assert(f < MAX_FACE_VALENCE);
          border_face = f;
	  face_size[f] = 2;
          crease_weight[f] = inf; 
	  assert(e < MAX_EDGE_VALENCE);
          ring[e++] = (Vec3fa_t) vertices[ p->getStartVertexIndex() ];
	  assert(e < MAX_EDGE_VALENCE);
          ring[e++] = vtx; // dummy vertex
          crease_weight[++f] = inf;
	  
          /*! goto other side of border */
          p = (SubdivMesh::HalfEdge*) h;
          while (p->hasOpposite()) 
            p = p->opposite()->next();
        }
	
      } while (p != h); 
      
      edge_valence = e;
      face_valence = f;
    }
    
    __forceinline void subdivide(CatmullClark1Ring& dest) const
    {
      dest.noForcedSubdivision = true;
      dest.edge_level = 0.5f*edge_level;
      dest.vertex_level = 0.5f*vertex_level;
      dest.face_valence = face_valence;
      dest.edge_valence = 2*face_valence;
      dest.border_index = border_face == -1 ? -1 : 2*border_face; // FIXME:
      dest.vertex_crease_weight   = max(0.0f,vertex_crease_weight-1.0f);
      assert(dest.face_valence <= CatmullClark1Ring::MAX_FACE_VALENCE);

      /* calculate face points */
      Vec3fa_t S = Vec3fa_t(0.0f);
      for (size_t f=0, v=0; f<face_valence; v+=face_size[f++]) {
        Vec3fa_t F = vtx;
        for (size_t k=v; k<=v+face_size[f]; k++) F += ring[k%edge_valence]; // FIXME: optimize
        S += dest.ring[2*f+1] = F/float(face_size[f]+2);
      }
      
      /* calculate new edge points */
      size_t num_creases = 0;
      size_t crease_id[MAX_FACE_VALENCE];
      Vec3fa_t C = Vec3fa_t(0.0f);
      for (size_t i=0, j=0; i<face_valence; j+=face_size[i++])
      {
        const Vec3fa_t v = vtx + ring[j];
        Vec3fa_t f = dest.ring[2*i+1];
        if (i == 0) f += dest.ring[dest.edge_valence-1]; 
        else        f += dest.ring[2*i-1];
        S += ring[j];
        dest.crease_weight[i] = max(crease_weight[i]-1.0f,0.0f);
        
        /* fast path for regular edge points */
        if (likely(crease_weight[i] <= 0.0f)) {
          dest.ring[2*i] = (v+f) * 0.25f;
        }
        
        /* slower path for hard edge rule */
        else {
          C += ring[j]; crease_id[num_creases++] = i;
          dest.ring[2*i] = v*0.5f;
	  
          /* even slower path for blended edge rule */
          if (unlikely(crease_weight[i] < 1.0f)) {
            const float w0 = crease_weight[i], w1 = 1.0f-w0;
            dest.ring[2*i] = w1*((v+f)*0.25f) + w0*(v*0.5f);
          }
        }
      }
      
      /* compute new vertex using smooth rule */
      const float inv_face_valence = 1.0f / (float)face_valence;
      const Vec3fa_t v_smooth = (Vec3fa_t)(S*inv_face_valence + (float(face_valence)-2.0f)*vtx)*inv_face_valence;
      dest.vtx = v_smooth;
      
      /* compute new vertex using vertex_crease_weight rule */
      if (unlikely(vertex_crease_weight > 0.0f)) 
      {
        if (vertex_crease_weight >= 1.0f) {
          dest.vtx = vtx;
        } else {
          const float t0 = vertex_crease_weight, t1 = 1.0f-t0;
          dest.vtx = t0*vtx + t1*v_smooth;
        }
        return;
      }
      
      if (likely(num_creases <= 1))
        return;
      
      /* compute new vertex using crease rule */
      if (likely(num_creases == 2)) {
        const Vec3fa_t v_sharp = (Vec3fa_t)(C + 6.0f * vtx) * (1.0f / 8.0f);
        const float crease_weight0 = crease_weight[crease_id[0]];
        const float crease_weight1 = crease_weight[crease_id[1]];
        dest.vtx = v_sharp;
        dest.crease_weight[crease_id[0]] = max(0.25f*(3.0f*crease_weight0 + crease_weight1)-1.0f,0.0f);
        dest.crease_weight[crease_id[1]] = max(0.25f*(3.0f*crease_weight1 + crease_weight0)-1.0f,0.0f);
        //dest.crease_weight[crease_id[0]] = max(0.5f*(crease_weight0 + crease_weight1)-1.0f,0.0f);
        //dest.crease_weight[crease_id[1]] = max(0.5f*(crease_weight1 + crease_weight0)-1.0f,0.0f);
        const float t0 = 0.5f*(crease_weight0+crease_weight1), t1 = 1.0f-t0;
        //dest.crease_weight[crease_id[0]] = t0 < 1.0f ? 0.0f : 0.5f*t0;
        //dest.crease_weight[crease_id[1]] = t0 < 1.0f ? 0.0f : 0.5f*t0;
        if (unlikely(t0 < 1.0f)) {
          dest.vtx = t0*v_sharp + t1*v_smooth;
        }
      }
      
      /* compute new vertex using corner rule */
      else {
        dest.vtx = vtx;
      }
    }

    void convert(CatmullClark1Ring& dst) const
    {
      assert(only_quads);
      assert(std::all_of(&face_size[0],&face_size[face_valence],[](int i) { return i == 2; }));
      
      dst.edge_level = edge_level;
      dst.vertex_level = vertex_level;
      dst.vtx = vtx;
      dst.face_valence = face_valence;
      dst.edge_valence = 2*face_valence;
      dst.border_index = border_face == -1 ? -1 : 2*border_face;
      for (size_t i=0; i<face_valence; i++) 
	dst.crease_weight[i] = crease_weight[i];
      dst.vertex_crease_weight = vertex_crease_weight;
      dst.noForcedSubdivision = true;
      for (size_t i=0; i<edge_valence; i++) dst.ring[i] = ring[i];
    }
    
    friend __forceinline std::ostream &operator<<(std::ostream &o, const GeneralCatmullClark1Ring &c)
    {
      o << "vtx " << c.vtx << " size = " << c.edge_valence << ", border_face = " << c.border_face << ", " << " face_valence = " << c.face_valence << 
	", edge_level = " << c.edge_level << ", vertex_level = " << c.vertex_level << ", ring: " << std::endl;
      for (size_t v=0, f=0; f<c.face_valence; v+=c.face_size[f++]) {
        for (size_t i=v; i<v+c.face_size[f]; i++) {
          o << i << " -> " << c.ring[i];
          if (i == v) o << " crease = " << c.crease_weight[f];
          o << std::endl;
        }
      }
      return o;
    } 
  };
}
