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

#include "bvh4i_intersector16_subdiv.h"
#include "bvh4i_leaf_intersector.h"
#include "geometry/subdivpatch1.h"
#include "common/subdiv/tessellation_cache.h"

#define TIMER(x) 

namespace embree
{
  namespace isa
  {
    __thread TessellationCache *tess_cache = NULL;

    struct __aligned(64) Quad4x4 {
      mic3f vtx;
      unsigned short uu[16];
      unsigned short vv[16];

      __forceinline void prefetchData() const
      {
	prefetch<PFHINT_NT>(&vtx.x);
	prefetch<PFHINT_NT>(&vtx.y);
	prefetch<PFHINT_NT>(&vtx.z);
	prefetch<PFHINT_NT>(uu);
      }
      
      __forceinline mic_f getU() const
      {
	return load16f_uint16(uu);
      }

      __forceinline mic_f getV() const
      {
	return load16f_uint16(vv);
      }

      __forceinline void set(const mic3f &v3, const mic_f &u, const mic_f &v)
      {
	vtx = v3;
	store16f_uint16(uu,u*65535.0f);
	store16f_uint16(vv,v*65535.0f);
      }

    };

    __forceinline void createSubPatchBVH4iLeaf(BVH4i::NodeRef &ref,
					       const unsigned int patchIndex) 
    {
      *(volatile unsigned int*)&ref = (patchIndex << BVH4i::encodingBits) | BVH4i::leaf_mask;
    }


    BBox3fa createSubTree(BVH4i::NodeRef &curNode,
			  mic_f *const lazymem,
			  const SubdivPatch1 &patch,
			  const float *const grid_u_array,
			  const float *const grid_v_array,
			  const GridRange &range,
			  unsigned int &localCounter,
			  const SubdivMesh* const geom)
    {
      if (range.hasLeafSize())
	{
	  const unsigned int u_start = range.u_start;
	  const unsigned int u_end   = range.u_end;
	  const unsigned int v_start = range.v_start;
	  const unsigned int v_end   = range.v_end;

	  const unsigned int u_size = u_end-u_start+1;
	  const unsigned int v_size = v_end-v_start+1;
	  const unsigned int currentIndex = localCounter;
	  localCounter += 4; // x,y,z + 16bit u,v


	  mic_f uu = mic_f::inf();
	  mic_f vv = mic_f::inf();

	  for (unsigned int v=v_start;v<=v_end;v++)
	    for (unsigned int u=u_start;u<=u_end;u++)
	      {
		const unsigned int local_v = v - v_start;
		const unsigned int local_u = u - u_start;
		uu[4 * local_v + local_u] = grid_u_array[ v * patch.grid_u_res + u ];
		vv[4 * local_v + local_u] = grid_v_array[ v * patch.grid_u_res + u ];
	      }

	  /* set invalid grid u,v value to border elements */

	  for (unsigned int x=u_size-1;x<4;x++)
	    for (unsigned int y=0;y<4;y++)
	      {
		uu[4 * y + x] = uu[4 * y + u_size-1];
		vv[4 * y + x] = vv[4 * y + u_size-1];
	      }

	  for (unsigned int y=v_size-1;y<4;y++)
	    for (unsigned int x=0;x<4;x++)
	      {
		uu[4 * y + x] = uu[4 * (v_size-1) + x];
		vv[4 * y + x] = vv[4 * (v_size-1) + x];
	      }
	  
	  mic3f vtx = patch.eval16(uu,vv);

	  const Vec2f uv0 = patch.getUV(0);
	  const Vec2f uv1 = patch.getUV(1);
	  const Vec2f uv2 = patch.getUV(2);
	  const Vec2f uv3 = patch.getUV(3);

	  if (unlikely(geom->displFunc != NULL))
	    {
	      mic3f normal      = patch.normal16(uu,vv);
	      normal = normalize(normal);

	      const mic_f patch_uu = bilinear_interpolate(uv0.x,uv1.x,uv2.x,uv3.x,uu,vv);
	      const mic_f patch_vv = bilinear_interpolate(uv0.y,uv1.y,uv2.y,uv3.y,uu,vv);

	      geom->displFunc(geom->userPtr,
			      patch.geom,
			      patch.prim,
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
	  const BBox3fa leafGridBounds = getBBox3fa(vtx);

	  Quad4x4 *quad4x4 = (Quad4x4*)&lazymem[currentIndex];
	  quad4x4->set(vtx,uu,vv);

	  createSubPatchBVH4iLeaf( curNode, currentIndex);	  
	  return leafGridBounds;
	}
      /* allocate new bvh4i node */
      const size_t num64BytesBlocksPerNode = 2;
      const size_t currentIndex = localCounter;
      localCounter += num64BytesBlocksPerNode;

      createBVH4iNode<2>(curNode,currentIndex);

      BVH4i::Node &node = *(BVH4i::Node*)curNode.node((BVH4i::Node*)lazymem);

      node.setInvalid();
      GridRange r[4];
      
      const unsigned int children = range.splitIntoSubRanges(r);
      
      /* create four subtrees */
      BBox3fa bounds( empty );
      for (unsigned int i=0;i<children;i++)
	{
	  BBox3fa bounds_subtree = createSubTree( node.child(i), 
						  lazymem, 
						  patch, 
						  grid_u_array,
						  grid_v_array,
						  r[i],
						  localCounter,
						  geom);
	  node.setBounds(i, bounds_subtree);
	  bounds.extend( bounds_subtree );
	}
      return bounds;      
    }

    BVH4i::NodeRef initLocalLazySubdivTree(const SubdivPatch1 &patch,
					   unsigned int currentIndex,
					   mic_f *lazymem,
					   const SubdivMesh* const geom)
    {

      TIMER(double msec = 0.0);
      TIMER(msec = getSeconds());

      __aligned(64) float u_array[(patch.grid_size_simd_blocks+1)*16]; // for unaligned access
      __aligned(64) float v_array[(patch.grid_size_simd_blocks+1)*16];

      gridUVTessellator(patch.level,
			patch.grid_u_res,
			patch.grid_v_res,
			u_array,
			v_array);

      /* stich different tessellation levels in u/v grid */
      if (patch.needsStiching())
	stichUVGrid(patch.level,patch.grid_u_res,patch.grid_v_res,u_array,v_array);


      BVH4i::NodeRef subtree_root = 0;
      const unsigned int oldIndex = currentIndex;

      BBox3fa bounds = createSubTree( subtree_root,
				      lazymem,
				      patch,
				      u_array,
				      v_array,
				      GridRange(0,patch.grid_u_res-1,0,patch.grid_v_res-1),
				      currentIndex,
				      geom);

      assert(currentIndex - oldIndex == patch.grid_subtree_size_64b_blocks);
      TIMER(msec = getSeconds()-msec);    
      return subtree_root;
    }



    static AtomicMutex mtx;

    void createTessellationCache()
    {
      TessellationCache *cache = (TessellationCache *)_mm_malloc(sizeof(TessellationCache),64);
      assert( (size_t)cache % 64 == 0 );
      cache->init();	
#if DEBUG
      mtx.lock();
      std::cout << "Enabling tessellation cache with " << cache->allocated64ByteBlocks() << " blocks = " << cache->allocated64ByteBlocks()*64 << " bytes as default size" << std::endl;
      mtx.unlock();
#endif
      tess_cache = cache;
    }

    __forceinline BVH4i::NodeRef lookUpTessellationCache(TessellationCache *local_cache,
							 const unsigned int patchIndex,
							 const unsigned int commitCounter,
							 const SubdivPatch1* const patches,
							 Scene *const scene)
    {
      TessellationCache::InputTagType tag = (TessellationCache::InputTagType)patchIndex;

      BVH4i::NodeRef subtree_root = local_cache->lookup(tag,commitCounter);
      if (unlikely(subtree_root == BVH4i::invalidNode))
	{
	  const SubdivPatch1& subdiv_patch = patches[patchIndex];
		  
	  subdiv_patch.prefetchData();

	  const SubdivMesh* const geom = (SubdivMesh*)scene->get(subdiv_patch.geom); // FIXME: test flag first

	  const unsigned int blocks = subdiv_patch.grid_subtree_size_64b_blocks;

	  TessellationCache::CacheTag &t = local_cache->request(tag,commitCounter,blocks);		      
	  mic_f *local_mem = (mic_f*)local_cache->getPtr(); 

	  unsigned int currentIndex = t.getRootRef();

	  subtree_root = initLocalLazySubdivTree(subdiv_patch,currentIndex,local_mem,geom);		      
	  assert( subtree_root != BVH4i::invalidNode);

	  local_cache->updateRootRef(t,subtree_root);
	}
      return subtree_root;
    }



    static unsigned int BVH4I_LEAF_MASK = BVH4i::leaf_mask; // needed due to compiler efficiency bug
    static unsigned int M_LANE_7777 = 0x7777;               // needed due to compiler efficiency bug

    // ============================================================================================
    // ============================================================================================
    // ============================================================================================


    void BVH4iIntersector16Subdiv::intersect(mic_i* valid_i, BVH4i* bvh, Ray16& ray16)
    {
      /* near and node stack */
      __aligned(64) float   stack_dist[4*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[4*BVH4i::maxDepth+1];

      /* query per thread tessellation cache */
      TessellationCache *local_cache = NULL;
      if (!tess_cache)
	createTessellationCache();
      local_cache = tess_cache;


      /* setup */
      const mic_m m_valid    = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir16     = rcp_safe(ray16.dir);
      const mic_f inf        = mic_f(pos_inf);
      const mic_f zero       = mic_f::zero();

      store16f(stack_dist,inf);
      ray16.primID = select(m_valid,mic_i(-1),ray16.primID);
      ray16.geomID = select(m_valid,mic_i(-1),ray16.geomID);

      Scene *const scene                         = (Scene*)bvh->geometry;
      const Node      * __restrict__ const nodes = (Node     *)bvh->nodePtr();
      Triangle1 * __restrict__ const accel       = (Triangle1*)bvh->triPtr();
      const unsigned int commitCounter           = scene->commitCounter;

      stack_node[0] = BVH4i::invalidNode;
      long rayIndex = -1;
      while((rayIndex = bitscan64(rayIndex,toInt(m_valid))) != BITSCAN_NO_BIT_SET_64)	    
        {
	  stack_node[1] = bvh->root;
	  size_t sindex = 2;

	  const mic_f org_xyz      = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
	  const mic_f dir_xyz      = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);
	  const mic_f rdir_xyz     = loadAOS4to16f(rayIndex,rdir16.x,rdir16.y,rdir16.z);
	  const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
	  const mic_f min_dist_xyz = broadcast1to16f(&ray16.tnear[rayIndex]);
	  mic_f       max_dist_xyz = broadcast1to16f(&ray16.tfar[rayIndex]);

	  const unsigned int leaf_mask = BVH4I_LEAF_MASK;

	  while (1)
	    {

	      NodeRef curNode = stack_node[sindex-1];
	      sindex--;

	      traverse_single_intersect<false>(curNode,
					       sindex,
					       rdir_xyz,
					       org_rdir_xyz,
					       min_dist_xyz,
					       max_dist_xyz,
					       stack_node,
					       stack_dist,
					       nodes,
					       leaf_mask);
		   

	      /* return if stack is empty */
	      if (unlikely(curNode == BVH4i::invalidNode)) break;

	      STAT3(normal.trav_leaves,1,1,1);
	      STAT3(normal.trav_prims,1,1,1);

	      //////////////////////////////////////////////////////////////////////////////////////////////////

	      const unsigned int patchIndex = curNode.offsetIndex();

	      BVH4i::NodeRef subtree_root = lookUpTessellationCache(local_cache,
								    patchIndex,
								    commitCounter,
								    (SubdivPatch1*)accel,
								    scene);
	      mic_f     * const __restrict__ lazymem     = (mic_f*)local_cache->getPtr(); /* lazymem could change to realloc */

	      // -------------------------------------
	      // -------------------------------------
	      // -------------------------------------

	      float   * __restrict__ const sub_stack_dist = &stack_dist[sindex];
	      NodeRef * __restrict__ const sub_stack_node = &stack_node[sindex];
	      sub_stack_node[0] = BVH4i::invalidNode;
	      sub_stack_node[1] = subtree_root;
	      ustore16f(sub_stack_dist,inf);
	      size_t sub_sindex = 2;

	      while (1)
		{
		  curNode = sub_stack_node[sub_sindex-1];
		  sub_sindex--;

		  traverse_single_intersect<false>(curNode,
						   sub_sindex,
						   rdir_xyz,
						   org_rdir_xyz,
						   min_dist_xyz,
						   max_dist_xyz,
						   sub_stack_node,
						   sub_stack_dist,
						   (BVH4i::Node*)lazymem,
						   leaf_mask);
		 		   

		  /* return if stack is empty */
		  if (unlikely(curNode == BVH4i::invalidNode)) break;

		  const unsigned int uvIndex = curNode.offsetIndex();
		  
		  const mic_m m_active = 0x777;
		  const Quad4x4 *__restrict__ const quad4x4 = (Quad4x4*)&lazymem[uvIndex];
		  quad4x4->prefetchData();
		  const mic_f uu = quad4x4->getU();
		  const mic_f vv = quad4x4->getV();

		  intersect1_quad16(rayIndex, 
				    dir_xyz,
				    org_xyz,
				    ray16,
				    quad4x4->vtx,
				    uu,
				    vv,
				    4,
				    m_active,
				    patchIndex);
		}

	      // -------------------------------------
	      // -------------------------------------
	      // -------------------------------------

	      compactStack(stack_node,stack_dist,sindex,mic_f(ray16.tfar[rayIndex]));

	      // ------------------------
	    }
	}


      /* update primID/geomID and compute normals */
      mic_m m_hit = (ray16.primID != -1) & m_valid;
      rayIndex = -1;
      while((rayIndex = bitscan64(rayIndex,toInt(m_hit))) != BITSCAN_NO_BIT_SET_64)	    
        {
	  const SubdivPatch1& subdiv_patch = ((SubdivPatch1*)accel)[ray16.primID[rayIndex]];
	  ray16.primID[rayIndex] = subdiv_patch.prim;
	  ray16.geomID[rayIndex] = subdiv_patch.geom;
#if defined(RTCORE_RETURN_SUBDIV_NORMAL)

	  if (unlikely(!subdiv_patch.hasDisplacement()))
	    {
	      const Vec3fa normal    = subdiv_patch.normal(ray16.u[rayIndex],ray16.v[rayIndex]);
	      ray16.Ng.x[rayIndex]   = normal.x;
	      ray16.Ng.y[rayIndex]   = normal.y;
	      ray16.Ng.z[rayIndex]   = normal.z;
	    }
#endif

	  const Vec2f uv0 = subdiv_patch.getUV(0);
	  const Vec2f uv1 = subdiv_patch.getUV(1);
	  const Vec2f uv2 = subdiv_patch.getUV(2);
	  const Vec2f uv3 = subdiv_patch.getUV(3);
	  
	  const float patch_u = bilinear_interpolate(uv0.x,uv1.x,uv2.x,uv3.x,ray16.u[rayIndex],ray16.v[rayIndex]);
	  const float patch_v = bilinear_interpolate(uv0.y,uv1.y,uv2.y,uv3.y,ray16.u[rayIndex],ray16.v[rayIndex]);
	  ray16.u[rayIndex] = patch_u;
	  ray16.v[rayIndex] = patch_v;
	}

    }

    void BVH4iIntersector16Subdiv::occluded(mic_i* valid_i, BVH4i* bvh, Ray16& ray16)
    {
      /* near and node stack */
      __aligned(64) NodeRef stack_node[4*BVH4i::maxDepth+1];

      /* query per thread tessellation cache */
      TessellationCache *local_cache = NULL;
      if (!tess_cache)
	createTessellationCache();
      local_cache = tess_cache;

      /* setup */
      const mic_m m_valid = *(mic_i*)valid_i != mic_i(0);
      const mic3f rdir16  = rcp_safe(ray16.dir);
      mic_m terminated    = !m_valid;
      const mic_f inf     = mic_f(pos_inf);
      const mic_f zero    = mic_f::zero();

      Scene *const scene                   = (Scene*)bvh->geometry;
      const Node      * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();
      const unsigned int commitCounter           = scene->commitCounter;

      stack_node[0] = BVH4i::invalidNode;
      ray16.primID = select(m_valid,mic_i(-1),ray16.primID);
      ray16.geomID = select(m_valid,mic_i(-1),ray16.geomID);

      long rayIndex = -1;
      while((rayIndex = bitscan64(rayIndex,toInt(m_valid))) != BITSCAN_NO_BIT_SET_64)	    
        {
	  stack_node[1] = bvh->root;
	  size_t sindex = 2;

	  const mic_f org_xyz      = loadAOS4to16f(rayIndex,ray16.org.x,ray16.org.y,ray16.org.z);
	  const mic_f dir_xyz      = loadAOS4to16f(rayIndex,ray16.dir.x,ray16.dir.y,ray16.dir.z);
	  const mic_f rdir_xyz     = loadAOS4to16f(rayIndex,rdir16.x,rdir16.y,rdir16.z);
	  const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
	  const mic_f min_dist_xyz = broadcast1to16f(&ray16.tnear[rayIndex]);
	  const mic_f max_dist_xyz = broadcast1to16f(&ray16.tfar[rayIndex]);
	  const mic_i v_invalidNode(BVH4i::invalidNode);
	  const unsigned int leaf_mask = BVH4I_LEAF_MASK;

	  while (1)
	    {
	      NodeRef curNode = stack_node[sindex-1];
	      sindex--;

	      traverse_single_occluded< false >(curNode,
						sindex,
						rdir_xyz,
						org_rdir_xyz,
						min_dist_xyz,
						max_dist_xyz,
						stack_node,
						nodes,
						leaf_mask);

	      /* return if stack is empty */
	      if (unlikely(curNode == BVH4i::invalidNode)) break;

	      STAT3(shadow.trav_leaves,1,1,1);
	      STAT3(shadow.trav_prims,1,1,1);

	      //////////////////////////////////////////////////////////////////////////////////////////////////

	      

	      const unsigned int patchIndex = curNode.offsetIndex();
	      BVH4i::NodeRef subtree_root = lookUpTessellationCache(local_cache,
								    patchIndex,
								    commitCounter,
								    (SubdivPatch1*)accel,
								    scene);
	      mic_f     * const __restrict__ lazymem     = (mic_f*)local_cache->getPtr(); /* lazymem could change to realloc */

		{
		    
		  // -------------------------------------
		  // -------------------------------------
		  // -------------------------------------

		  __aligned(64) NodeRef sub_stack_node[64];
		  sub_stack_node[0] = BVH4i::invalidNode;
		  sub_stack_node[1] = subtree_root;
		  size_t sub_sindex = 2;

		  ///////////////////////////////////////////////////////////////////////////////////////////////////////
		  while (1)
		    {
		      curNode = sub_stack_node[sub_sindex-1];
		      sub_sindex--;

		      traverse_single_occluded<false>(curNode,
						      sub_sindex,
						      rdir_xyz,
						      org_rdir_xyz,
						      min_dist_xyz,
						      max_dist_xyz,
						      sub_stack_node,
						      (BVH4i::Node*)lazymem,
						      leaf_mask);
		 		   

		      /* return if stack is empty */
		      if (unlikely(curNode == BVH4i::invalidNode)) break;

		      const unsigned int uvIndex = curNode.offsetIndex();
		  
		      const mic_m m_active = 0x777;
		      const Quad4x4 *__restrict__ const quad4x4 = (Quad4x4*)&lazymem[uvIndex];
		      quad4x4->prefetchData();
		      const mic_f uu = quad4x4->getU();
		      const mic_f vv = quad4x4->getV();
		  
		      if (unlikely(occluded1_quad16(rayIndex, 
						    dir_xyz,
						    org_xyz,
						    ray16,
						    quad4x4->vtx,
						    uu,
						    vv,
						    4,
						    m_active,
						    patchIndex)))
			{
			  terminated |= (mic_m)((unsigned int)1 << rayIndex);
			  break;
			}
		    }
		}

		if (unlikely(terminated & (mic_m)((unsigned int)1 << rayIndex))) break;
	      //////////////////////////////////////////////////////////////////////////////////////////////////

	    }

	  if (unlikely(all(toMask(terminated)))) break;
	}


      store16i(m_valid & toMask(terminated),&ray16.geomID,0);

    }


    void BVH4iIntersector1Subdiv::intersect(BVH4i* bvh, Ray& ray)
    {
      /* near and node stack */
      __aligned(64) float   stack_dist[3*BVH4i::maxDepth+1];
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      /* query per thread tessellation cache */
      TessellationCache *local_cache = NULL;
      if (!tess_cache)
	createTessellationCache();
      local_cache = tess_cache;

      /* setup */
      const mic3f rdir16     = rcp_safe(mic3f(mic_f(ray.dir.x),mic_f(ray.dir.y),mic_f(ray.dir.z)));
      const mic_f inf        = mic_f(pos_inf);
      const mic_f zero       = mic_f::zero();

      store16f(stack_dist,inf);

      const Node      * __restrict__ nodes = (Node    *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();
      Scene *const scene                   = (Scene*)bvh->geometry;
      const unsigned int commitCounter     = scene->commitCounter;

      stack_node[0] = BVH4i::invalidNode;      
      stack_node[1] = bvh->root;

      size_t sindex = 2;

      const mic_f org_xyz      = loadAOS4to16f(ray.org.x,ray.org.y,ray.org.z);
      const mic_f dir_xyz      = loadAOS4to16f(ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f rdir_xyz     = loadAOS4to16f(rdir16.x[0],rdir16.y[0],rdir16.z[0]);
      const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
      const mic_f min_dist_xyz = broadcast1to16f(&ray.tnear);
      mic_f       max_dist_xyz = broadcast1to16f(&ray.tfar);
	  
      const unsigned int leaf_mask = BVH4I_LEAF_MASK;
	  
      while (1)
	{
	  NodeRef curNode = stack_node[sindex-1];
	  sindex--;

	  traverse_single_intersect<false>(curNode,
					   sindex,
					   rdir_xyz,
					   org_rdir_xyz,
					   min_dist_xyz,
					   max_dist_xyz,
					   stack_node,
					   stack_dist,
					   nodes,
					   leaf_mask);            		    

	  /* return if stack is empty */
	  if (unlikely(curNode == BVH4i::invalidNode)) break;

	  STAT3(normal.trav_leaves,1,1,1);
	  STAT3(normal.trav_prims,1,1,1);


	  //////////////////////////////////////////////////////////////////////////////////////////////////

	  const unsigned int patchIndex = curNode.offsetIndex();

	  BVH4i::NodeRef subtree_root = lookUpTessellationCache(local_cache,
								patchIndex,
								commitCounter,
								(SubdivPatch1*)accel,
								scene);
	  mic_f     * const __restrict__ lazymem     = (mic_f*)local_cache->getPtr(); /* lazymem could change to realloc */

	  // -------------------------------------
	  // -------------------------------------
	  // -------------------------------------

	  float   * __restrict__ const sub_stack_dist = &stack_dist[sindex];
	  NodeRef * __restrict__ const sub_stack_node = &stack_node[sindex];
	  sub_stack_node[0] = BVH4i::invalidNode;
	  sub_stack_node[1] = subtree_root;
	  ustore16f(sub_stack_dist,inf);
	  size_t sub_sindex = 2;

	  while (1)
	    {
	      curNode = sub_stack_node[sub_sindex-1];
	      sub_sindex--;

	      traverse_single_intersect<false>(curNode,
					       sub_sindex,
					       rdir_xyz,
					       org_rdir_xyz,
					       min_dist_xyz,
					       max_dist_xyz,
					       sub_stack_node,
					       sub_stack_dist,
					       (BVH4i::Node*)lazymem,
					       leaf_mask);
		 		   

	      /* return if stack is empty */
	      if (unlikely(curNode == BVH4i::invalidNode)) break;

	      const unsigned int uvIndex = curNode.offsetIndex();
		  
	      const mic_m m_active = 0x777;
	      const Quad4x4 *__restrict__ const quad4x4 = (Quad4x4*)&lazymem[uvIndex];
	      quad4x4->prefetchData();
	      const mic_f uu = quad4x4->getU();
	      const mic_f vv = quad4x4->getV();

	      intersect1_quad16(dir_xyz,
				org_xyz,
				ray,
				quad4x4->vtx,
				uu,
				vv,
				4,
				m_active,
				patchIndex);
	    }



	  compactStack(stack_node,stack_dist,sindex,max_dist_xyz);
	}
    }

    void BVH4iIntersector1Subdiv::occluded(BVH4i* bvh, Ray& ray)
    {
      /* near and node stack */
      __aligned(64) NodeRef stack_node[3*BVH4i::maxDepth+1];

      /* query per thread tessellation cache */
      TessellationCache *local_cache = NULL;
      if (!tess_cache)
	createTessellationCache();
      local_cache = tess_cache;

      /* setup */
      const mic3f rdir16      = rcp_safe(mic3f(ray.dir.x,ray.dir.y,ray.dir.z));
      const mic_f inf         = mic_f(pos_inf);
      const mic_f zero        = mic_f::zero();

      const Node      * __restrict__ nodes = (Node     *)bvh->nodePtr();
      const Triangle1 * __restrict__ accel = (Triangle1*)bvh->triPtr();
      Scene *const scene                   = (Scene*)bvh->geometry;
      const unsigned int commitCounter     = scene->commitCounter;

      stack_node[0] = BVH4i::invalidNode;
      stack_node[1] = bvh->root;
      size_t sindex = 2;

      const mic_f org_xyz      = loadAOS4to16f(ray.org.x,ray.org.y,ray.org.z);
      const mic_f dir_xyz      = loadAOS4to16f(ray.dir.x,ray.dir.y,ray.dir.z);
      const mic_f rdir_xyz     = loadAOS4to16f(rdir16.x[0],rdir16.y[0],rdir16.z[0]);
      const mic_f org_rdir_xyz = org_xyz * rdir_xyz;
      const mic_f min_dist_xyz = broadcast1to16f(&ray.tnear);
      const mic_f max_dist_xyz = broadcast1to16f(&ray.tfar);

      const unsigned int leaf_mask = BVH4I_LEAF_MASK;
	  
      while (1)
	{
	  NodeRef curNode = stack_node[sindex-1];
	  sindex--;
            
	  
	  traverse_single_occluded< false >(curNode,
					    sindex,
					    rdir_xyz,
					    org_rdir_xyz,
					    min_dist_xyz,
					    max_dist_xyz,
					    stack_node,
					    nodes,
					    leaf_mask);	    

	  /* return if stack is empty */
	  if (unlikely(curNode == BVH4i::invalidNode)) break;


	      STAT3(shadow.trav_leaves,1,1,1);
	      STAT3(shadow.trav_prims,1,1,1);

	      //////////////////////////////////////////////////////////////////////////////////////////////////

	      

	      const unsigned int patchIndex = curNode.offsetIndex();
	      BVH4i::NodeRef subtree_root = lookUpTessellationCache(local_cache,
								    patchIndex,
								    commitCounter,
								    (SubdivPatch1*)accel,
								    scene);
	      mic_f     * const __restrict__ lazymem     = (mic_f*)local_cache->getPtr(); /* lazymem could change to realloc */

		    
	      // -------------------------------------
	      // -------------------------------------
	      // -------------------------------------

	      __aligned(64) NodeRef sub_stack_node[64];
	      sub_stack_node[0] = BVH4i::invalidNode;
	      sub_stack_node[1] = subtree_root;
	      size_t sub_sindex = 2;

	      ///////////////////////////////////////////////////////////////////////////////////////////////////////
	      while (1)
		{
		  curNode = sub_stack_node[sub_sindex-1];
		  sub_sindex--;

		  traverse_single_occluded<false>(curNode,
						  sub_sindex,
						  rdir_xyz,
						  org_rdir_xyz,
						  min_dist_xyz,
						  max_dist_xyz,
						  sub_stack_node,
						  (BVH4i::Node*)lazymem,
						  leaf_mask);
		 		   

		  /* return if stack is empty */
		  if (unlikely(curNode == BVH4i::invalidNode)) break;

		  const unsigned int uvIndex = curNode.offsetIndex();
		  
		  const mic_m m_active = 0x777;
		  const Quad4x4 *__restrict__ const quad4x4 = (Quad4x4*)&lazymem[uvIndex];
		  quad4x4->prefetchData();
		  const mic_f uu = quad4x4->getU();
		  const mic_f vv = quad4x4->getV();
		  
		  if (unlikely(occluded1_quad16(dir_xyz,
						org_xyz,
						ray,
						quad4x4->vtx,
						uu,
						vv,
						4,
						m_active,
						patchIndex)))
		    {
		      ray.geomID = 0;
		      return;
		    }
		}

	  //////////////////////////////////////////////////////////////////////////////////////////////////

	}
    }

    // ----------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------------------



    typedef BVH4iIntersector16Subdiv SubdivIntersector16SingleMoellerFilter;
    typedef BVH4iIntersector16Subdiv SubdivIntersector16SingleMoellerNoFilter;

    DEFINE_INTERSECTOR16   (BVH4iSubdivMeshIntersector16        , SubdivIntersector16SingleMoellerFilter);
    DEFINE_INTERSECTOR16   (BVH4iSubdivMeshIntersector16NoFilter, SubdivIntersector16SingleMoellerNoFilter);

    typedef BVH4iIntersector1Subdiv SubdivMeshIntersector1MoellerFilter;
    typedef BVH4iIntersector1Subdiv SubdivMeshIntersector1MoellerNoFilter;

    DEFINE_INTERSECTOR1    (BVH4iSubdivMeshIntersector1        , SubdivMeshIntersector1MoellerFilter);
    DEFINE_INTERSECTOR1    (BVH4iSubdivMeshIntersector1NoFilter, SubdivMeshIntersector1MoellerNoFilter);

  }
}
