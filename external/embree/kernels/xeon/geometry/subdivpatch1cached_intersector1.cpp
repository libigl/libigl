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

#include "subdivpatch1cached_intersector1.h"
#include "xeon/bvh4/bvh4.h"
#include "xeon/bvh4/bvh4_intersector1.h"

#define TIMER(x) 

namespace embree
{
  namespace isa
  {  
    
    __thread TessellationCache *SubdivPatch1CachedIntersector1::thread_cache = NULL;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BVH4::NodeRef SubdivPatch1CachedIntersector1::buildSubdivPatchTree(const SubdivPatch1Cached &patch,
                                                                       void *const lazymem,
                                                                       const SubdivMesh* const geom)
    {
        
      TIMER(double msec = 0.0);
      TIMER(msec = getSeconds());
        
      assert( patch.grid_size_simd_blocks >= 1 );
#if !defined(_MSC_VER) || defined(__INTEL_COMPILER)
      __aligned(64) float grid_x[(patch.grid_size_simd_blocks+1)*8]; 
      __aligned(64) float grid_y[(patch.grid_size_simd_blocks+1)*8];
      __aligned(64) float grid_z[(patch.grid_size_simd_blocks+1)*8]; 
        
      __aligned(64) float grid_u[(patch.grid_size_simd_blocks+1)*8]; 
      __aligned(64) float grid_v[(patch.grid_size_simd_blocks+1)*8];
     
#else
      const size_t array_elements = (patch.grid_size_simd_blocks + 1) * 8;
      float *const ptr = (float*)_malloca(5 * array_elements * sizeof(float) + 64);
      float *const grid_arrays = (float*)ALIGN_PTR(ptr,64);

      float *grid_x = &grid_arrays[array_elements * 0];
      float *grid_y = &grid_arrays[array_elements * 1];
      float *grid_z = &grid_arrays[array_elements * 2];
      float *grid_u = &grid_arrays[array_elements * 3];
      float *grid_v = &grid_arrays[array_elements * 4];

        
#endif   
      evalGrid(patch,grid_x,grid_y,grid_z,grid_u,grid_v,geom);
        
      BVH4::NodeRef subtree_root = BVH4::encodeNode( (BVH4::Node*)lazymem);
      unsigned int currentIndex = 0;
      BBox3fa bounds = createSubTree( subtree_root,
				      (float*)lazymem,
				      patch,
				      grid_x,
				      grid_y,
				      grid_z,
				      grid_u,
				      grid_v,
				      GridRange(0,patch.grid_u_res-1,0,patch.grid_v_res-1),
				      currentIndex,
				      geom);
        
      assert(currentIndex == patch.grid_subtree_size_64b_blocks);

      TIMER(msec = getSeconds()-msec);            

#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
      _freea(ptr);
#endif
      return subtree_root;
    }
    
    
    BBox3fa SubdivPatch1CachedIntersector1::createSubTree(BVH4::NodeRef &curNode,
                                                          float *const lazymem,
                                                          const SubdivPatch1Cached &patch,
                                                          const float *const grid_x_array,
                                                          const float *const grid_y_array,
                                                          const float *const grid_z_array,
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
        
	  assert(u_size >= 1);
	  assert(v_size >= 1);
        
	  assert(u_size*v_size <= 9);
        
	  const unsigned int currentIndex = localCounter;
	  localCounter +=  (sizeof(Quad2x2)+63) / 64; 
        
	  Quad2x2 *qquad = (Quad2x2*)&lazymem[currentIndex*16];
                
	  ssef leaf_x_array[3];
	  ssef leaf_y_array[3];
	  ssef leaf_z_array[3];
	  ssef leaf_u_array[3];
	  ssef leaf_v_array[3];
        
	  for (unsigned int v=v_start;v<=v_end;v++)
	    {
	      const size_t offset = v * patch.grid_u_res + u_start;
	      const unsigned int local_v = v - v_start;
	      leaf_x_array[local_v] = loadu4f(&grid_x_array[ offset ]);
	      leaf_y_array[local_v] = loadu4f(&grid_y_array[ offset ]);
	      leaf_z_array[local_v] = loadu4f(&grid_z_array[ offset ]);
	      leaf_u_array[local_v] = loadu4f(&grid_u_array[ offset ]);
	      leaf_v_array[local_v] = loadu4f(&grid_v_array[ offset ]);            
	    }
        
	  /* set invalid grid u,v value to border elements */
	  for (unsigned int x=u_size-1;x<3;x++)
	    for (unsigned int y=0;y<3;y++)
	      {
		leaf_x_array[y][x] = leaf_x_array[y][u_size-1];
		leaf_y_array[y][x] = leaf_y_array[y][u_size-1];
		leaf_z_array[y][x] = leaf_z_array[y][u_size-1];
		leaf_u_array[y][x] = leaf_u_array[y][u_size-1];
		leaf_v_array[y][x] = leaf_v_array[y][u_size-1];
	      }
        
	  for (unsigned int y=v_size-1;y<3;y++)
	    for (unsigned int x=0;x<3;x++)
	      {
		leaf_x_array[y][x] = leaf_x_array[v_size-1][x];
		leaf_y_array[y][x] = leaf_y_array[v_size-1][x];
		leaf_z_array[y][x] = leaf_z_array[v_size-1][x];
		leaf_u_array[y][x] = leaf_u_array[v_size-1][x];
		leaf_v_array[y][x] = leaf_v_array[v_size-1][x];
	      }
        
                
	  qquad->init( leaf_x_array, 
		       leaf_y_array, 
		       leaf_z_array, 
		       leaf_u_array, 
		       leaf_v_array);
        
#if 0
	  DBG_PRINT("LEAF");
	  DBG_PRINT(u_start);
	  DBG_PRINT(v_start);
	  DBG_PRINT(u_end);
	  DBG_PRINT(v_end);
        
	  for (unsigned int y=0;y<3;y++)
	    for (unsigned int x=0;x<3;x++)
	      std::cout << y << " " << x 
			<< " ->  x = " << leaf_x_array[y][x] << " y = " << leaf_v_array[y][x] << " z = " << leaf_z_array[y][x]
			<< "   u = " << leaf_u_array[y][x] << " v = " << leaf_v_array[y][x] << std::endl;
        
	  DBG_PRINT( *qquad );
        
#endif          
        
	  BBox3fa bounds = qquad->bounds();
	  curNode = BVH4::encodeLeaf(qquad,2);
        
	  return bounds;
	}
      
      
      /* allocate new bvh4 node */
      const size_t currentIndex = localCounter;
      
      /* 128 bytes == 2 x 64 bytes cachelines */
      localCounter += 2; 
      
      BVH4::Node *node = (BVH4::Node *)&lazymem[currentIndex*16];
      
      curNode = BVH4::encodeNode( node );
      
      node->clear();
      
      GridRange r[4];
      
      const unsigned int children = range.splitIntoSubRanges(r);
      
      /* create four subtrees */
      BBox3fa bounds( empty );
      
      for (unsigned int i=0;i<children;i++)
	{
	  BBox3fa bounds_subtree = createSubTree( node->child(i), 
						  lazymem, 
						  patch, 
						  grid_x_array,
						  grid_y_array,
						  grid_z_array,
						  grid_u_array,
						  grid_v_array,
						  r[i],						  
						  localCounter,
						  geom);
	  node->set(i, bounds_subtree);
	  bounds.extend( bounds_subtree );
	}
      
      return bounds;
    }
    
    void SubdivPatch1CachedIntersector1::createTessellationCache()
    {
      TessellationCache *cache = (TessellationCache *)_mm_malloc(sizeof(TessellationCache),64);
      assert( (size_t)cache % 64 == 0 );
      cache->init();	
#if defined(DEBUG) && 0
      static AtomicMutex mtx;
      mtx.lock();
      std::cout << "Enabling tessellation cache with " << cache->allocated64ByteBlocks() << " blocks = " << cache->allocated64ByteBlocks()*64 << " bytes as default size" << std::endl;
      mtx.unlock();
#endif
      thread_cache = cache;
    }
    
  };
}
