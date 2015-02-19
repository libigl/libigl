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

#include "catmullclark_ring.h"

namespace embree
{
  class __aligned(64) CatmullClarkPatch
  {
  public:
    CatmullClark1Ring ring[4];

    __forceinline CatmullClarkPatch () {}

    __forceinline CatmullClarkPatch (const SubdivMesh::HalfEdge* first_half_edge, const BufferT<Vec3fa>& vertices) 
    {
      for (size_t i=0; i<4; i++)
        ring[i].init(first_half_edge+i,vertices);
    }

    __forceinline Vec3fa normal(const float uu, const float vv) const // FIXME: remove
    {
      const int iu = (int) uu;
      const int iv = (int) vv;
      const int index = iv > 0 ? 3-iu : iu;
      const Vec3fa tu = ring[index].getLimitTangent();
      const Vec3fa tv = ring[index].getSecondLimitTangent();
      return cross(tv,tu);
    }   

    __forceinline Vec3fa eval(const float uu, const float vv) const
    {
      const int iu = (int) uu;
      const int iv = (int) vv;
      const int index = iv > 0 ? 3-iu : iu;
      return ring[index].getLimitVertex();
    }   

    __forceinline Vec3fa getLimitVertex(const size_t index) const {
      return ring[index].getLimitVertex();
    }

    __forceinline Vec3fa getLimitTangent(const size_t index) const {
      return ring[index].getLimitTangent();
    }

    __forceinline Vec3fa getSecondLimitTangent(const size_t index) const {
      return ring[index].getSecondLimitTangent();
    }

    __forceinline BBox3fa bounds() const
    {
      BBox3fa bounds (ring[0].bounds());
      for (size_t i=1; i<4; i++)
	bounds.extend(ring[i].bounds());
      return bounds;
    }

    /*! returns true if the patch is a B-spline patch */
    __forceinline bool isRegular() const {
      return ring[0].isRegular() && ring[1].isRegular() && ring[2].isRegular() && ring[3].isRegular();
    }

    /*! returns true if the patch is a B-spline patch or final Quad */
    __forceinline bool isRegularOrFinal(const size_t depth) const {
      return ring[0].isRegularOrFinal(depth) && ring[1].isRegularOrFinal(depth) && ring[2].isRegularOrFinal(depth) && ring[3].isRegularOrFinal(depth);
    }

    /*! returns true if the patch is a B-spline patch or final Quad */
    __forceinline bool isRegularOrFinal2(const size_t depth) const {
      return ring[0].isRegularOrFinal2(depth) && ring[1].isRegularOrFinal2(depth) && ring[2].isRegularOrFinal2(depth) && ring[3].isRegularOrFinal2(depth);
    }

    /*! returns true of the patch is a Gregory patch */
    __forceinline bool isGregoryOrFinal(const size_t depth) const {
      return ring[0].isGregoryOrFinal(depth) && ring[1].isGregoryOrFinal(depth) && ring[2].isGregoryOrFinal(depth) && ring[3].isGregoryOrFinal(depth);
    }

    static __forceinline void init_regular(const CatmullClark1Ring& p0,
					   const CatmullClark1Ring& p1,
					   CatmullClark1Ring& dest0,
					   CatmullClark1Ring& dest1) 
    {
      dest1.vertex_level = dest0.vertex_level = p0.edge_level;
      dest1.face_valence = dest0.face_valence = 4;
      dest1.edge_valence = dest0.edge_valence = 8;
      dest1.border_index = dest0.border_index = -1;
      dest1.vtx = dest0.vtx = (Vec3fa_t)p0.ring[0];
      dest1.vertex_crease_weight = dest0.vertex_crease_weight = 0.0f;
      dest1.noForcedSubdivision = dest0.noForcedSubdivision = true;

      dest1.ring[2] = dest0.ring[0] = (Vec3fa_t)p0.ring[1];
      dest1.ring[1] = dest0.ring[7] = (Vec3fa_t)p1.ring[0];
      dest1.ring[0] = dest0.ring[6] = (Vec3fa_t)p1.vtx;
      dest1.ring[7] = dest0.ring[5] = (Vec3fa_t)p1.ring[4];
      dest1.ring[6] = dest0.ring[4] = (Vec3fa_t)p0.ring[p0.edge_valence-1];
      dest1.ring[5] = dest0.ring[3] = (Vec3fa_t)p0.ring[p0.edge_valence-2];
      dest1.ring[4] = dest0.ring[2] = (Vec3fa_t)p0.vtx;
      dest1.ring[3] = dest0.ring[1] = (Vec3fa_t)p0.ring[2];

      dest1.crease_weight[1] = dest0.crease_weight[0] = 0.0f;
      dest1.crease_weight[0] = dest0.crease_weight[3] = p1.crease_weight[1];
      dest1.crease_weight[3] = dest0.crease_weight[2] = 0.0f;
      dest1.crease_weight[2] = dest0.crease_weight[1] = p0.crease_weight[0];
    }


    static __forceinline void init_border(const CatmullClark1Ring &p0,
                                          const CatmullClark1Ring &p1,
                                          CatmullClark1Ring &dest0,
                                          CatmullClark1Ring &dest1) 
    {
      dest1.vertex_level = dest0.vertex_level = p0.edge_level;
      dest1.face_valence = dest0.face_valence = 3;
      dest1.edge_valence = dest0.edge_valence = 6;
      dest0.border_index = 2;
      dest1.border_index = 4;
      dest1.vtx  = dest0.vtx = (Vec3fa_t)p0.ring[0];
      dest1.vertex_crease_weight = dest0.vertex_crease_weight = 0.0f;
      dest1.noForcedSubdivision = dest0.noForcedSubdivision = true;

      dest1.ring[2] = dest0.ring[0] = (Vec3fa_t)p0.ring[1];
      dest1.ring[1] = dest0.ring[5] = (Vec3fa_t)p1.ring[0];
      dest1.ring[0] = dest0.ring[4] = (Vec3fa_t)p1.vtx;
      dest1.ring[5] = dest0.ring[3] = (Vec3fa_t)p0.ring[p0.border_index+1]; // dummy
      dest1.ring[4] = dest0.ring[2] = (Vec3fa_t)p0.vtx;
      dest1.ring[3] = dest0.ring[1] = (Vec3fa_t)p0.ring[2];

      dest1.crease_weight[1] = dest0.crease_weight[0] = 0.0f;
      dest1.crease_weight[0] = dest0.crease_weight[2] = p1.crease_weight[1];
      dest1.crease_weight[2] = dest0.crease_weight[1] = p0.crease_weight[0];
    }

    static __forceinline void init_regular(const Vec3fa_t &center, const Vec3fa_t center_ring[8], const size_t offset, CatmullClark1Ring &dest)
    {
      dest.vertex_level = 0.0f;
      dest.face_valence = 4;
      dest.edge_valence = 8;
      dest.border_index = -1;
      dest.vtx     = (Vec3fa_t)center;
      dest.vertex_crease_weight = 0.0f;
      dest.noForcedSubdivision = true;
      for (size_t i=0; i<8; i++) 
	dest.ring[i] = (Vec3fa_t)center_ring[(offset+i)%8];
      for (size_t i=0; i<4; i++) 
        dest.crease_weight[i] = 0.0f;
    }

    void subdivide(CatmullClarkPatch patch[4]) const
    {
      ring[0].subdivide(patch[0].ring[0]);
      ring[1].subdivide(patch[1].ring[1]);
      ring[2].subdivide(patch[2].ring[2]);
      ring[3].subdivide(patch[3].ring[3]);

      patch[0].ring[0].edge_level = 0.5f*ring[0].edge_level;
      patch[0].ring[1].edge_level = 0.25f*(ring[1].edge_level+ring[3].edge_level);
      patch[0].ring[2].edge_level = 0.25f*(ring[0].edge_level+ring[2].edge_level);
      patch[0].ring[3].edge_level = 0.5f*ring[3].edge_level;

      patch[1].ring[0].edge_level = 0.5f*ring[0].edge_level;
      patch[1].ring[1].edge_level = 0.5f*ring[1].edge_level;
      patch[1].ring[2].edge_level = 0.25f*(ring[0].edge_level+ring[2].edge_level);
      patch[1].ring[3].edge_level = 0.25f*(ring[1].edge_level+ring[3].edge_level);

      patch[2].ring[0].edge_level = 0.25f*(ring[0].edge_level+ring[2].edge_level);
      patch[2].ring[1].edge_level = 0.5f*ring[1].edge_level;
      patch[2].ring[2].edge_level = 0.5f*ring[2].edge_level;
      patch[2].ring[3].edge_level = 0.25f*(ring[1].edge_level+ring[3].edge_level);

      patch[3].ring[0].edge_level = 0.25f*(ring[0].edge_level+ring[2].edge_level);
      patch[3].ring[1].edge_level = 0.25f*(ring[1].edge_level+ring[3].edge_level);
      patch[3].ring[2].edge_level = 0.5f*ring[2].edge_level;
      patch[3].ring[3].edge_level = 0.5f*ring[3].edge_level;

      if (likely(ring[0].has_last_face()))
        init_regular(patch[0].ring[0],patch[1].ring[1],patch[0].ring[1],patch[1].ring[0]);
      else
        init_border(patch[0].ring[0],patch[1].ring[1],patch[0].ring[1],patch[1].ring[0]);

      if (likely(ring[1].has_last_face()))
        init_regular(patch[1].ring[1],patch[2].ring[2],patch[1].ring[2],patch[2].ring[1]);
      else
        init_border(patch[1].ring[1],patch[2].ring[2],patch[1].ring[2],patch[2].ring[1]);

      if (likely(ring[2].has_last_face()))
        init_regular(patch[2].ring[2],patch[3].ring[3],patch[2].ring[3],patch[3].ring[2]);
      else
        init_border(patch[2].ring[2],patch[3].ring[3],patch[2].ring[3],patch[3].ring[2]);

      if (likely(ring[3].has_last_face()))
        init_regular(patch[3].ring[3],patch[0].ring[0],patch[3].ring[0],patch[0].ring[3]);
      else
        init_border(patch[3].ring[3],patch[0].ring[0],patch[3].ring[0],patch[0].ring[3]);

      Vec3fa_t center = (ring[0].vtx + ring[1].vtx + ring[2].vtx + ring[3].vtx) * 0.25f;
      Vec3fa_t center_ring[8];

      // counter-clockwise
      center_ring[0] = (Vec3fa_t)patch[3].ring[3].ring[0];
      center_ring[7] = (Vec3fa_t)patch[3].ring[3].vtx;
      center_ring[6] = (Vec3fa_t)patch[2].ring[2].ring[0];
      center_ring[5] = (Vec3fa_t)patch[2].ring[2].vtx;
      center_ring[4] = (Vec3fa_t)patch[1].ring[1].ring[0];
      center_ring[3] = (Vec3fa_t)patch[1].ring[1].vtx;
      center_ring[2] = (Vec3fa_t)patch[0].ring[0].ring[0];
      center_ring[1] = (Vec3fa_t)patch[0].ring[0].vtx;

      init_regular(center,center_ring,0,patch[0].ring[2]);
      init_regular(center,center_ring,2,patch[1].ring[3]);
      init_regular(center,center_ring,4,patch[2].ring[0]);
      init_regular(center,center_ring,6,patch[3].ring[1]);
    }

    __forceinline void init( FinalQuad& quad ) const
    {
      quad.vtx[0] = (Vec3fa_t)ring[0].vtx;
      quad.vtx[1] = (Vec3fa_t)ring[1].vtx;
      quad.vtx[2] = (Vec3fa_t)ring[2].vtx;
      quad.vtx[3] = (Vec3fa_t)ring[3].vtx;
    };

    friend __forceinline std::ostream &operator<<(std::ostream &o, const CatmullClarkPatch &p)
    {
      o << "CatmullClarkPatch { " << std::endl;
      for (size_t i=0; i<4; i++)
	o << "ring" << i << ": " << p.ring[i] << std::endl;
      o << "}" << std::endl;
      return o;
    }
  };

  class __aligned(64) GeneralCatmullClarkPatch
  {
  public:
    static const size_t SIZE = SubdivMesh::MAX_VALENCE;
    GeneralCatmullClark1Ring ring[SIZE];
    size_t N;

    __forceinline GeneralCatmullClarkPatch () 
      : N(0) {}

    __forceinline size_t size() const { 
      return N; 
    }

    __forceinline bool isQuadPatch() const {
      return (N == 4) && ring[0].only_quads && ring[1].only_quads && ring[2].only_quads && ring[3].only_quads;
    }

    __forceinline GeneralCatmullClarkPatch (const SubdivMesh::HalfEdge* h, const BufferT<Vec3fa>& vertices) 
    {
      size_t i = 0;
      const SubdivMesh::HalfEdge* edge = h; 
      do {
	ring[i].init(edge,vertices);
        edge = edge->next();
        i++;
      } while ((edge != h) && (i < SIZE));
      N = i;
    }

    static __forceinline void init_regular(const CatmullClark1Ring& p0,
					   const CatmullClark1Ring& p1,
					   CatmullClark1Ring& dest0,
					   CatmullClark1Ring& dest1) 
    {
      dest1.vertex_level = dest0.vertex_level = p0.edge_level;
      dest1.face_valence = dest0.face_valence = 4;
      dest1.edge_valence = dest0.edge_valence = 8;
      dest1.border_index = dest0.border_index = -1;
      dest1.vtx = dest0.vtx = (Vec3fa_t)p0.ring[0];
      dest1.vertex_crease_weight = dest0.vertex_crease_weight = 0.0f;
      dest1.noForcedSubdivision = dest0.noForcedSubdivision = true;

      dest1.ring[2] = dest0.ring[0] = (Vec3fa_t)p0.ring[1];
      dest1.ring[1] = dest0.ring[7] = (Vec3fa_t)p1.ring[0];
      dest1.ring[0] = dest0.ring[6] = (Vec3fa_t)p1.vtx;
      dest1.ring[7] = dest0.ring[5] = (Vec3fa_t)p1.ring[4];
      dest1.ring[6] = dest0.ring[4] = (Vec3fa_t)p0.ring[p0.edge_valence-1];
      dest1.ring[5] = dest0.ring[3] = (Vec3fa_t)p0.ring[p0.edge_valence-2];
      dest1.ring[4] = dest0.ring[2] = (Vec3fa_t)p0.vtx;
      dest1.ring[3] = dest0.ring[1] = (Vec3fa_t)p0.ring[2];

      dest1.crease_weight[1] = dest0.crease_weight[0] = 0.0f;
      dest1.crease_weight[0] = dest0.crease_weight[3] = p1.crease_weight[1];
      dest1.crease_weight[3] = dest0.crease_weight[2] = 0.0f;
      dest1.crease_weight[2] = dest0.crease_weight[1] = p0.crease_weight[0];
    }


    static __forceinline void init_border(const CatmullClark1Ring &p0,
                                          const CatmullClark1Ring &p1,
                                          CatmullClark1Ring &dest0,
                                          CatmullClark1Ring &dest1) 
    {
      dest1.vertex_level = dest0.vertex_level = p0.edge_level;
      dest1.face_valence = dest0.face_valence = 3;
      dest1.edge_valence = dest0.edge_valence = 6;
      dest0.border_index = 2;
      dest1.border_index = 4;
      dest1.vtx  = dest0.vtx = (Vec3fa_t)p0.ring[0];
      dest1.vertex_crease_weight = dest0.vertex_crease_weight = 0.0f;
      dest1.noForcedSubdivision = dest0.noForcedSubdivision = true;

      dest1.ring[2] = dest0.ring[0] = (Vec3fa_t)p0.ring[1];
      dest1.ring[1] = dest0.ring[5] = (Vec3fa_t)p1.ring[0];
      dest1.ring[0] = dest0.ring[4] = (Vec3fa_t)p1.vtx;
      dest1.ring[5] = dest0.ring[3] = (Vec3fa_t)p0.ring[p0.border_index+1]; // dummy
      dest1.ring[4] = dest0.ring[2] = (Vec3fa_t)p0.vtx;
      dest1.ring[3] = dest0.ring[1] = (Vec3fa_t)p0.ring[2];

      dest1.crease_weight[1] = dest0.crease_weight[0] = 0.0f;
      dest1.crease_weight[0] = dest0.crease_weight[2] = p1.crease_weight[1];
      dest1.crease_weight[2] = dest0.crease_weight[1] = p0.crease_weight[0];
    }

    static __forceinline void init_regular(const Vec3fa_t &center, const Vec3fa_t center_ring[2*SIZE], const float vertex_level, const size_t N, const size_t offset, CatmullClark1Ring &dest)
    {
      assert(N<CatmullClark1Ring::MAX_FACE_VALENCE);
      dest.vertex_level = vertex_level;
      dest.face_valence = N;
      dest.edge_valence = 2*N;
      dest.border_index = -1;
      dest.vtx     = (Vec3fa_t)center;
      dest.vertex_crease_weight = 0.0f;
      dest.noForcedSubdivision = true;
      for (size_t i=0; i<2*N; i++) 
	dest.ring[i] = (Vec3fa_t)center_ring[(2*N+offset+i-1)%(2*N)];
      for (size_t i=0; i<N; i++) 
        dest.crease_weight[i] = 0.0f;
    }
 
    void subdivide(CatmullClarkPatch patch[SIZE], size_t& N_o) const
    {
      N_o = N;

      for (size_t i=0; i<N; i++) {
        size_t ip1 = (i+1)%N; // FIXME: %
        ring[i].subdivide(patch[i].ring[0]);
        patch[i]  .ring[0].edge_level = 0.5f*ring[i].edge_level;
        patch[ip1].ring[3].edge_level = 0.5f*ring[i].edge_level;
      }

      Vec3fa_t center = Vec3fa_t(0.0f);
      Vec3fa_t center_ring[2*SIZE];
      float center_vertex_level = 2.0f; // guarantees that irregular vertices get always isolated also for non-quads

      for (size_t i=0; i<N; i++)
      {
        size_t ip1 = (i+1)%N; // FIXME: %
        size_t im1 = (i+N-1)%N; // FIXME: %
        if (likely(ring[i].has_last_face())) init_regular(patch[i].ring[0],patch[ip1].ring[0],patch[i].ring[1],patch[ip1].ring[3]); 
        else                                 init_border (patch[i].ring[0],patch[ip1].ring[0],patch[i].ring[1],patch[ip1].ring[3]);

	float level = 0.25f*(ring[im1].edge_level+ring[ip1].edge_level);
        patch[i].ring[1].edge_level = patch[ip1].ring[2].edge_level = level;
	center_vertex_level = max(center_vertex_level,level);

        center += ring[i].vtx;
        center_ring[2*i+0] = (Vec3fa_t)patch[i].ring[0].vtx;
        center_ring[2*i+1] = (Vec3fa_t)patch[i].ring[0].ring[0];
      }
      center /= float(N);

      for (size_t i=0; i<N; i++) {
        init_regular(center,center_ring,center_vertex_level,N,2*i,patch[i].ring[2]);
      }
    }

    void init(CatmullClarkPatch& patch) const
    {
      assert(size() == 4);
      ring[0].convert(patch.ring[0]);
      ring[1].convert(patch.ring[1]);
      ring[2].convert(patch.ring[2]);
      ring[3].convert(patch.ring[3]);
    }

    friend __forceinline std::ostream &operator<<(std::ostream &o, const GeneralCatmullClarkPatch &p)
    {
      o << "GeneralCatmullClarkPatch { " << std::endl;
      for (size_t i=0; i<p.N; i++)
	o << "ring" << i << ": " << p.ring[i] << std::endl;
      o << "}" << std::endl;
      return o;
    }
  };
}


