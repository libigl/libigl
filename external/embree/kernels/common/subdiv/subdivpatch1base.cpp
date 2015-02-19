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

#include "common/scene_subdiv_mesh.h"
#include "subdivpatch1base.h"

namespace embree
{

    // if (grid_size_simd_blocks == 1)
    //   {
    //     mic_m m_active = 0xffff;
    //     for (unsigned int i=grid_u_res-1;i<16;i+=grid_u_res)
    //       m_active ^= (unsigned int)1 << i;
    //     m_active &= ((unsigned int)1 << (grid_u_res * (grid_v_res-1)))-1;
    //     grid_mask = m_active;
    //   }

  void SubdivPatch1Base::updateEdgeLevels(const float edge_level[4],const SubdivMesh *const mesh)
  {
    /* init discrete edge tessellation levels and grid resolution */

    assert( edge_level[0] >= 0.0f );
    assert( edge_level[1] >= 0.0f );
    assert( edge_level[2] >= 0.0f );
    assert( edge_level[3] >= 0.0f );
      
    level[0] = max(ceilf(edge_level[0]),1.0f);
    level[1] = max(ceilf(edge_level[1]),1.0f);
    level[2] = max(ceilf(edge_level[2]),1.0f);
    level[3] = max(ceilf(edge_level[3]),1.0f);

    grid_u_res = max(level[0],level[2])+1; // n segments -> n+1 points
    grid_v_res = max(level[1],level[3])+1;

#if defined(__MIC__)
    grid_size_simd_blocks        = ((grid_u_res*grid_v_res+15)&(-16)) / 16;
    grid_subtree_size_64b_blocks = 5; // single leaf with u,v,x,y,z      
#else
    /* 8-wide SIMD is default on Xeon */
    grid_size_simd_blocks        = ((grid_u_res*grid_v_res+7)&(-8)) / 8;
    grid_subtree_size_64b_blocks = (sizeof(Quad2x2)+63) / 64; // single Quad2x2

#endif
    /* need stiching? */

    flags &= ~TRANSITION_PATCH;

    const unsigned int int_edge_points0 = (unsigned int)level[0] + 1;
    const unsigned int int_edge_points1 = (unsigned int)level[1] + 1;
    const unsigned int int_edge_points2 = (unsigned int)level[2] + 1;
    const unsigned int int_edge_points3 = (unsigned int)level[3] + 1;
      
    if (int_edge_points0 < (unsigned int)grid_u_res ||
	int_edge_points2 < (unsigned int)grid_u_res ||
	int_edge_points1 < (unsigned int)grid_v_res ||
	int_edge_points3 < (unsigned int)grid_v_res)
      flags |= TRANSITION_PATCH;

    /* tessellate into grid blocks for larger grid resolutions, generate bvh4 subtree over grid blocks*/

#if defined(__MIC__)
    const size_t leafBlocks = 4;
#else
    const size_t leafBlocks = (sizeof(Quad2x2)+63) / 64;
#endif
    grid_subtree_size_64b_blocks = getSubTreeSize64bBlocks( leafBlocks ); // u,v,x,y,z 

    /* has displacements? */
    flags &= ~HAS_DISPLACEMENT;
    if (mesh->displFunc != NULL)
      flags |= HAS_DISPLACEMENT;

  }


  /*! Construction from vertices and IDs. */
  SubdivPatch1Base::SubdivPatch1Base (const CatmullClarkPatch& ipatch,
                                      const unsigned int gID,
                                      const unsigned int pID,
                                      const SubdivMesh *const mesh,
                                      const Vec2f uv[4],
                                      const float edge_level[4]) 
    : geom(gID),
      prim(pID),  
      flags(0)
  {
    assert(sizeof(SubdivPatch1Base) == 5 * 64);

    for (size_t i=0;i<4;i++)
      {
        /* need to reverse input here */
        u[i] = (unsigned short)(uv[i].y * 65535.0f);
        v[i] = (unsigned short)(uv[i].x * 65535.0f);
      }


    updateEdgeLevels(edge_level,mesh);
     

    /* determine whether patch is regular or not */

    if (ipatch.isRegularOrFinal(0) && mesh->displFunc == NULL)
      {
        flags |= REGULAR_PATCH;
        patch.init( ipatch );
      }
    else
      {
        GregoryPatch gpatch; 
        gpatch.init( ipatch ); 
        gpatch.exportDenseConrolPoints( patch.v );
      }
#if 0
    DBG_PRINT( grid_u_res );
    DBG_PRINT( grid_v_res );
    DBG_PRINT( grid_size_16wide_blocks );
    DBG_PRINT( grid_mask );
    DBG_PRINT( grid_subtree_size_64b_blocks );
#endif
  }

  BBox3fa SubdivPatch1Base::bounds(const SubdivMesh* const mesh) const
  {
#if FORCE_TESSELLATION_BOUNDS == 1

#if !defined(_MSC_VER) || defined(__INTEL_COMPILER)
    __aligned(64) float u_array[(grid_size_simd_blocks + 1) * 16]; // +16 for unaligned access
    __aligned(64) float v_array[(grid_size_simd_blocks + 1) * 16]; // +16 for unaligned access

#else
    const size_t array_elements = (grid_size_simd_blocks + 1) * 8;

    float *const ptr = (float*)_malloca(2 *array_elements * sizeof(float) + 64);
    float *const uv_arrays = (float*)ALIGN_PTR(ptr, 64);
    float *const u_array = &uv_arrays[0];
    float *const v_array = &uv_arrays[grid_size_simd_blocks * 8];

#endif

    const unsigned int real_grid_size = grid_u_res*grid_v_res;
    gridUVTessellator(level, grid_u_res, grid_v_res, u_array, v_array);

    if (unlikely(needsStiching()))
      stichUVGrid(level, grid_u_res, grid_v_res, u_array, v_array);

    // FIXME: remove
#if defined(__MIC__)
    for (size_t i = real_grid_size; i<grid_size_simd_blocks * 16; i++)
#else
      for (size_t i = real_grid_size; i<grid_size_simd_blocks * 8; i++)
#endif
        {
          u_array[i] = 1.0f;
          v_array[i] = 1.0f;
        }

    BBox3fa b(empty);
    assert(grid_size_simd_blocks >= 1);

#if defined(__MIC__)
    for (size_t i = 0; i<grid_size_simd_blocks; i++)
      {
        const mic_f u = load16f(&u_array[i * 16]);
        const mic_f v = load16f(&v_array[i * 16]);

        mic3f vtx = eval16(u, v);


        /* eval displacement function */
        if (unlikely(mesh->displFunc != NULL))
          {
            mic3f normal = normal16(u, v);
            normal = normalize(normal);

            const Vec2f uv0 = getUV(0);
            const Vec2f uv1 = getUV(1);
            const Vec2f uv2 = getUV(2);
            const Vec2f uv3 = getUV(3);

            const mic_f patch_uu = bilinear_interpolate(uv0.x, uv1.x, uv2.x, uv3.x, u, v);
            const mic_f patch_vv = bilinear_interpolate(uv0.y, uv1.y, uv2.y, uv3.y, u, v);

            mesh->displFunc(mesh->userPtr,
                            geom,
                            prim,
                            (const float*)&patch_uu,
                            (const float*)&patch_vv,
                            (const float*)&normal.x,
                            (const float*)&normal.y,
                            (const float*)&normal.z,
                            (float*)&vtx.x,
                            (float*)&vtx.y,
                            (float*)&vtx.z,
                            16);

          }

        /* extend bounding box */
        b.extend(getBBox3fa(vtx));
      }

#else

#if !defined(__AVX__)

    for (size_t i = 0; i<grid_size_simd_blocks * 2; i++)
      {
        ssef u = load4f(&u_array[i * 4]);
        ssef v = load4f(&v_array[i * 4]);

        sse3f vtx = eval4(u, v);

        /* eval displacement function */
        if (unlikely(mesh->displFunc != NULL))
          {
            const Vec2f uv0 = getUV(0);
            const Vec2f uv1 = getUV(1);
            const Vec2f uv2 = getUV(2);
            const Vec2f uv3 = getUV(3);

            const ssef patch_uu = bilinear_interpolate(uv0.x, uv1.x, uv2.x, uv3.x, u, v);
            const ssef patch_vv = bilinear_interpolate(uv0.y, uv1.y, uv2.y, uv3.y, u, v);

            sse3f nor = normal4(u, v);

            nor = normalize_safe(nor);

            mesh->displFunc(mesh->userPtr,
                            geom,
                            prim,
                            (const float*)&patch_uu,
                            (const float*)&patch_vv,
                            (const float*)&nor.x,
                            (const float*)&nor.y,
                            (const float*)&nor.z,
                            (float*)&vtx.x,
                            (float*)&vtx.y,
                            (float*)&vtx.z,
                            4);
          }
        b.extend(getBBox3fa(vtx));
      }

#else

    for (size_t i = 0; i<grid_size_simd_blocks; i++)
      {
        avxf u = load8f(&u_array[i * 8]);
        avxf v = load8f(&v_array[i * 8]);

        avx3f vtx = eval8(u, v);

        /* eval displacement function */
        if (unlikely(mesh->displFunc != NULL))
          {
            const Vec2f uv0 = getUV(0);
            const Vec2f uv1 = getUV(1);
            const Vec2f uv2 = getUV(2);
            const Vec2f uv3 = getUV(3);

            const avxf patch_uu = bilinear_interpolate(uv0.x, uv1.x, uv2.x, uv3.x, u, v);
            const avxf patch_vv = bilinear_interpolate(uv0.y, uv1.y, uv2.y, uv3.y, u, v);

            avx3f nor = normal8(u, v);

            nor = normalize_safe(nor);

            mesh->displFunc(mesh->userPtr,
                            geom,
                            prim,
                            (const float*)&patch_uu,
                            (const float*)&patch_vv,
                            (const float*)&nor.x,
                            (const float*)&nor.y,
                            (const float*)&nor.z,
                            (float*)&vtx.x,
                            (float*)&vtx.y,
                            (float*)&vtx.z,
                            8);
          }
        b.extend(getBBox3fa(vtx));
      }
#endif

#endif
    b.lower.a = 0.0f;
    b.upper.a = 0.0f;

#if defined(DEBUG)
    isfinite(b.lower.x);
    isfinite(b.lower.y);
    isfinite(b.lower.z);

    isfinite(b.upper.x);
    isfinite(b.upper.y);
    isfinite(b.upper.z);
#endif

#else
    BBox3fa b = patch.bounds();
    if (unlikely(isGregoryPatch()))
      {
        b.extend(GregoryPatch::extract_f_m_Vec3fa(patch.v, 0));
        b.extend(GregoryPatch::extract_f_m_Vec3fa(patch.v, 1));
        b.extend(GregoryPatch::extract_f_m_Vec3fa(patch.v, 2));
        b.extend(GregoryPatch::extract_f_m_Vec3fa(patch.v, 3));
      }
#endif

#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    _freea(ptr);
#endif
    return b;
  }


}
