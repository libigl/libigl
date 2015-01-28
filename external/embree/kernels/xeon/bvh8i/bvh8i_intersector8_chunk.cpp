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

#include "bvh8i_intersector8_chunk.h"
#include "geometry/triangle8_intersector8_moeller.h"


#define DBG(x) 

#define START_SINDEX 5

namespace embree
{
  namespace isa
  {    
    
    template<typename TriangleIntersector8>    
    void BVH8iIntersector8Chunk<TriangleIntersector8>::intersect(avxb* valid_i, BVH8i* bvh, Ray8& ray)
    {
#if defined(__AVX__)
     /* load node and primitive array */
      const Node     * __restrict__ nodes = (Node    *)bvh->nodePtr();
      const Triangle * __restrict__ accel = (Triangle*)bvh->triPtr();
      
      /* load ray */
      const avxb valid0 = *valid_i;
      const avx3f rdir = rcp_safe(ray.dir);
      const avx3f org_rdir = ray.org * rdir;
      avxf ray_tnear = select(valid0,ray.tnear,pos_inf);
      avxf ray_tfar  = select(valid0,ray.tfar ,neg_inf);
      const avxf inf = avxf(pos_inf);
      
      /* allocate stack and push root node */
      avxf    stack_near[3*BVH4i::maxDepth+1];
      NodeRef stack_node[3*BVH4i::maxDepth+1];
      stack_node[0] = BVH4i::invalidNode;
      stack_near[0] = inf;
      stack_node[1] = bvh->root;
      stack_near[1] = ray_tnear; 
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      avxf*    __restrict__ sptr_near = stack_near + 2;
      
      while (1)
      {
        /* pop next node from stack */
        sptr_node--;
        sptr_near--;
        NodeRef curNode = *sptr_node;
        if (unlikely(curNode == BVH4i::invalidNode)) 
          break;
        
        /* cull node if behind closest hit point */
        avxf curDist = *sptr_near;
        if (unlikely(none(ray_tfar > curDist))) 
          continue;
        
        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(curNode.isLeaf()))
            break;
          
          const avxb valid_node = ray_tfar > curDist;
          STAT3(normal.trav_nodes,1,popcnt(valid_node),8);
          const Node* __restrict__ const node = (BVH8i::Node*)curNode.node(nodes);
          
          /* pop of next node */
          sptr_node--;
          sptr_near--;
          curNode = *sptr_node; // FIXME: this trick creates issues with stack depth
          curDist = *sptr_near;
          
          for (unsigned i=0; i<8; i++)
          {
            const NodeRef child = node->children[i];
            if (unlikely(child == BVH4i::emptyNode)) break;
            
#if defined(__AVX2__)
            const avxf lclipMinX = msub(node->lower_x[i],rdir.x,org_rdir.x);
            const avxf lclipMinY = msub(node->lower_y[i],rdir.y,org_rdir.y);
            const avxf lclipMinZ = msub(node->lower_z[i],rdir.z,org_rdir.z);
            const avxf lclipMaxX = msub(node->upper_x[i],rdir.x,org_rdir.x);
            const avxf lclipMaxY = msub(node->upper_y[i],rdir.y,org_rdir.y);
            const avxf lclipMaxZ = msub(node->upper_z[i],rdir.z,org_rdir.z);
            const avxf lnearP = maxi(maxi(mini(lclipMinX, lclipMaxX), mini(lclipMinY, lclipMaxY)), mini(lclipMinZ, lclipMaxZ));
            const avxf lfarP  = mini(mini(maxi(lclipMinX, lclipMaxX), maxi(lclipMinY, lclipMaxY)), maxi(lclipMinZ, lclipMaxZ));
            const avxb lhit   = maxi(lnearP,ray_tnear) <= mini(lfarP,ray_tfar);      
#else
            const avxf lclipMinX = node->lower_x[i] * rdir.x - org_rdir.x;
            const avxf lclipMinY = node->lower_y[i] * rdir.y - org_rdir.y;
            const avxf lclipMinZ = node->lower_z[i] * rdir.z - org_rdir.z;
            const avxf lclipMaxX = node->upper_x[i] * rdir.x - org_rdir.x;
            const avxf lclipMaxY = node->upper_y[i] * rdir.y - org_rdir.y;
            const avxf lclipMaxZ = node->upper_z[i] * rdir.z - org_rdir.z;
            const avxf lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
            const avxf lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
            const avxb lhit   = max(lnearP,ray_tnear) <= min(lfarP,ray_tfar);      
#endif
            
            /* if we hit the child we choose to continue with that child if it 
               is closer than the current next child, or we push it onto the stack */
            if (likely(any(lhit)))
            {
              const avxf childDist = select(lhit,lnearP,inf);
              const NodeRef child = node->children[i];
              
              /* push cur node onto stack and continue with hit child */
              if (any(childDist < curDist))
              {
                *sptr_node = curNode;
                *sptr_near = curDist; 
		sptr_node++;
		sptr_near++;

                curDist = childDist;
                curNode = child;
              }
              
              /* push hit child onto stack*/
              else {
                *sptr_node = child;
                *sptr_near = childDist; 
		sptr_node++;
		sptr_near++;

              }
              assert(sptr_node - stack_node < BVH4i::maxDepth);
            }	      
          }
        }
        
        /* return if stack is empty */
        if (unlikely(curNode == BVH4i::invalidNode)) 
          break;
        
        /* intersect leaf */
        const avxb valid_leaf = ray_tfar > curDist;
        STAT3(normal.trav_leaves,1,popcnt(valid_leaf),8);
        size_t items; const Triangle* tri  = (Triangle*) curNode.leaf(accel, items);
        TriangleIntersector8::intersect(valid_leaf,ray,tri,items,bvh->geometry);
        ray_tfar = select(valid_leaf,ray.tfar,ray_tfar);
      }
      AVX_ZERO_UPPER();
#endif       
    }
    
     template<typename TriangleIntersector8>
    void BVH8iIntersector8Chunk<TriangleIntersector8>::occluded(avxb* valid_i, BVH8i* bvh, Ray8& ray)
    {
#if defined(__AVX__)
      /* load node and primitive array */
      const Node     * __restrict__ nodes = (Node    *)bvh->nodePtr();
      const Triangle * __restrict__ accel = (Triangle*)bvh->triPtr();
      
      /* load ray */
      const avxb valid = *valid_i;
      avxb terminated = !valid;
      const avx3f rdir = rcp_safe(ray.dir);
      const avx3f org_rdir = ray.org * rdir;
      avxf ray_tnear = select(valid,ray.tnear,pos_inf);
      avxf ray_tfar  = select(valid,ray.tfar ,neg_inf);
      const avxf inf = avxf(pos_inf);
      
      /* allocate stack and push root node */
      avxf    stack_near[3*BVH4i::maxDepth+1];
      NodeRef stack_node[3*BVH4i::maxDepth+1];
      stack_node[0] = BVH4i::invalidNode;
      stack_near[0] = inf;
      stack_node[1] = bvh->root;
      stack_near[1] = ray_tnear; 
      NodeRef* __restrict__ sptr_node = stack_node + 2;
      avxf*    __restrict__ sptr_near = stack_near + 2;
      
      while (1)
      {
        /* pop next node from stack */
        sptr_node--;
        sptr_near--;
        NodeRef curNode = *sptr_node;
        if (unlikely(curNode == BVH4i::invalidNode)) 
          break;
        
        /* cull node if behind closest hit point */
        avxf curDist = *sptr_near;
        if (unlikely(none(ray_tfar > curDist))) 
          continue;
        
        while (1)
        {
          /* test if this is a leaf node */
          if (unlikely(curNode.isLeaf()))
            break;
          
          const avxb valid_node = ray_tfar > curDist;
          STAT3(shadow.trav_nodes,1,popcnt(valid_node),8);
          const Node* __restrict__ const node = (Node*)curNode.node(nodes);
          
          /* pop of next node */
          sptr_node--;
          sptr_near--;
          curNode = *sptr_node; // FIXME: this trick creates issues with stack depth
          curDist = *sptr_near;
          
          for (unsigned i=0; i<8; i++)
          {
            const NodeRef child = node->children[i];
            if (unlikely(child == BVH4i::emptyNode)) break;
            
#if defined(__AVX2__)
            const avxf lclipMinX = msub(node->lower_x[i],rdir.x,org_rdir.x);
            const avxf lclipMinY = msub(node->lower_y[i],rdir.y,org_rdir.y);
            const avxf lclipMinZ = msub(node->lower_z[i],rdir.z,org_rdir.z);
            const avxf lclipMaxX = msub(node->upper_x[i],rdir.x,org_rdir.x);
            const avxf lclipMaxY = msub(node->upper_y[i],rdir.y,org_rdir.y);
            const avxf lclipMaxZ = msub(node->upper_z[i],rdir.z,org_rdir.z);
            const avxf lnearP = maxi(maxi(mini(lclipMinX, lclipMaxX), mini(lclipMinY, lclipMaxY)), mini(lclipMinZ, lclipMaxZ));
            const avxf lfarP  = mini(mini(maxi(lclipMinX, lclipMaxX), maxi(lclipMinY, lclipMaxY)), maxi(lclipMinZ, lclipMaxZ));
            const avxb lhit   = maxi(lnearP,ray_tnear) <= mini(lfarP,ray_tfar);      
#else
            const avxf lclipMinX = node->lower_x[i] * rdir.x - org_rdir.x;
            const avxf lclipMinY = node->lower_y[i] * rdir.y - org_rdir.y;
            const avxf lclipMinZ = node->lower_z[i] * rdir.z - org_rdir.z;
            const avxf lclipMaxX = node->upper_x[i] * rdir.x - org_rdir.x;
            const avxf lclipMaxY = node->upper_y[i] * rdir.y - org_rdir.y;
            const avxf lclipMaxZ = node->upper_z[i] * rdir.z - org_rdir.z;
            const avxf lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
            const avxf lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
            const avxb lhit   = max(lnearP,ray_tnear) <= min(lfarP,ray_tfar);      
#endif
            
            /* if we hit the child we choose to continue with that child if it 
               is closer than the current next child, or we push it onto the stack */
            if (likely(any(lhit)))
            {
              const avxf childDist = select(lhit,lnearP,inf);
              sptr_node++;
              sptr_near++;
              
              /* push cur node onto stack and continue with hit child */
              if (any(childDist < curDist))
              {
                *(sptr_node-1) = curNode;
                *(sptr_near-1) = curDist; 
                curDist = childDist;
                curNode = child;
              }
              
              /* push hit child onto stack*/
              else {
                *(sptr_node-1) = child;
                *(sptr_near-1) = childDist; 
              }
              assert(sptr_node - stack_node < BVH4i::maxDepth);
            }	      
          }
        }
        
        /* return if stack is empty */
        if (unlikely(curNode == BVH4i::invalidNode)) 
          break;
        
        /* intersect leaf */
        const avxb valid_leaf = ray_tfar > curDist;
        STAT3(shadow.trav_leaves,1,popcnt(valid_leaf),8);
        size_t items; const Triangle* tri  = (Triangle*) curNode.leaf(accel, items);
        terminated |= TriangleIntersector8::occluded(!terminated,ray,tri,items,bvh->geometry);
        if (all(terminated)) break;
        ray_tfar = select(terminated,neg_inf,ray_tfar);
      }
      store8i(valid & terminated,&ray.geomID,0);
      AVX_ZERO_UPPER();
#endif      
    }
    
    DEFINE_INTERSECTOR8(BVH8iTriangle8Intersector8ChunkMoeller,BVH8iIntersector8Chunk<Triangle8Intersector8MoellerTrumbore>);

  }
}  


//       DBG(PING);
//       avxf stack_near[QBVH_MAX_STACK_DEPTH];
//       unsigned int stack_node[QBVH_MAX_STACK_DEPTH];
      
//       avxb m_active = *valid_i;
      
//       const avx3f &origin    = ray.org;
//       const avx3f &direction = ray.dir;
//       const avxf   max_dist  = ray.tfar;
      
//       const BVH8i::BVH8iNode *const __restrict__ bvh8  = (BVH8i::BVH8iNode*)bvh->qbvh;
      
//       const Triangle1 *const __restrict__ accel = (Triangle1*)bvh->accel;
      
//       const avxf orgx = origin.x;
//       const avxf orgy = origin.y;
//       const avxf orgz = origin.z;  
      
//       const avxf ulp(6E-8f);
//       const avxf direction_x = select(direction.x == 0,ulp,direction.x);
//       const avxf direction_y = select(direction.y == 0,ulp,direction.y);
//       const avxf direction_z = select(direction.z == 0,ulp,direction.z);
      
//       const avxf rdirx = rcp(direction_x);
//       const avxf rdiry = rcp(direction_y);
//       const avxf rdirz = rcp(direction_z);
      
//       const avxf org_rdirx = orgx * rdirx;
//       const avxf org_rdiry = orgy * rdiry;
//       const avxf org_rdirz = orgz * rdirz;
      
//       const avxf inf  = avxf(pos_inf);
      
//       /* --  dummy node at index 0, removes branch in inner-loop -- */
      
//       stack_node[0] = QBVH_TERMINAL_TOKEN;
//       stack_near[0] = inf;
      
//       unsigned int* __restrict__ sptr_node = stack_node + 1;
//       avxf       * __restrict__ sptr_near = stack_near + 1;
      
//       *sptr_node++ = bvh8[0].min_d[0];
//       *sptr_near++ = ray.tnear; 
      
//       avxf max_distance   = select(m_active,max_dist,avxf(zero));
      
//       while (1)
//       {
// 	sptr_node--;
// 	sptr_near--;
        
// 	unsigned int curNode = *sptr_node;
        
// 	if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) break;
        
// 	avxf minDist = *sptr_near;
// 	const avxb m_dist = max_distance > minDist;
        
// 	if (unlikely( none(m_dist) )) continue;
        
// 	assert(curNode != QBVH_TERMINAL_TOKEN);
        
// #if defined(USE_AVXI_OR_MIN_MAX_CMP)
// 	const avxi min_distance = cast(ray.tnear);
// #else
// 	const avxf min_distance = ray.tnear;
// #endif
        
// 	while(1)
//         {
//           DBG(DBG_PRINT(curNode));
          
//           if (unlikely(qbvhLeaf(curNode) != 0)) break;
          
//           assert (curNode != QBVH_TERMINAL_TOKEN);
          
//           sptr_node--;
//           sptr_near--;
//           STAT3(normal.trav_nodes,1,popcnt(m_active),8);
//           const BVH8i::BVH8iNode* __restrict__ const bptr = BVH8i::bvh8ChildPtrNoMask(bvh8,curNode);
          
//           curNode = *sptr_node;
//           minDist = *sptr_near;
          
//           DBG(
//             DBG_PRINT(rdirx);
//             DBG_PRINT(rdiry);
//             DBG_PRINT(rdirz);
            
//             DBG_PRINT(org_rdirx);
//             DBG_PRINT(org_rdiry);
//             DBG_PRINT(org_rdirz);
            
//             DBG_PRINT(*qptr);
//             );
          
// #pragma unroll(8)
//           for (unsigned int i=0;i<8;i++)
//           {
//             const unsigned int dd = bptr->min_d[i];
//             if (unlikely(dd == QBVH_LEAF_MASK)) { break; }
            
//             const avxf lclipMinX = bptr->min_x[i] * rdirx - org_rdirx;
//             const avxf lclipMinY = bptr->min_y[i] * rdiry - org_rdiry;
//             const avxf lclipMinZ = bptr->min_z[i] * rdirz - org_rdirz;
//             const avxf lclipMaxX = bptr->max_x[i] * rdirx - org_rdirx;
//             const avxf lclipMaxY = bptr->max_y[i] * rdiry - org_rdiry;
//             const avxf lclipMaxZ = bptr->max_z[i] * rdirz - org_rdirz;
            
// #if !defined(USE_AVXI_OR_MIN_MAX_CMP)
//             const avxf lnearP = max(max(min(lclipMinX, lclipMaxX), min(lclipMinY, lclipMaxY)), min(lclipMinZ, lclipMaxZ));
//             const avxf lfarP  = min(min(max(lclipMinX, lclipMaxX), max(lclipMinY, lclipMaxY)), max(lclipMinZ, lclipMaxZ));
//             const avxb lhit   = max(lnearP,min_distance) <= min(lfarP,max_distance);      
// #else
//             const avxi iclipMinX = cast(lclipMinX);
//             const avxi iclipMaxX = cast(lclipMaxX);
//             const avxi iclipMinY = cast(lclipMinY);
//             const avxi iclipMaxY = cast(lclipMaxY);
//             const avxi iclipMinZ = cast(lclipMinZ);
//             const avxi iclipMaxZ = cast(lclipMaxZ);
            
//             const avxi lnearP = max(max(min(iclipMinX, iclipMaxX), min(iclipMinY, iclipMaxY)), min(iclipMinZ, iclipMaxZ));
//             const avxi lfarP  = min(min(max(iclipMinX, iclipMaxX), max(iclipMinY, iclipMaxY)), max(iclipMinZ, iclipMaxZ));      
//             //const avxb lhit   = cast(max(lnearP,min_distance)) <= cast(min(lfarP,cast(max_distance)));
//             const avxb lhit   = max(lnearP,min_distance) <= min(lfarP,cast(max_distance));
            
// #endif
//             //const avxb lhit   = max(lnearP,ray.tnear) <= min(lfarP,ray.t);
            
//             DBG(
//               DBG_PRINT(i);
//               DBG_PRINT(lhit);
//               );
            
//             if (likely(any(lhit)))
//             {
// #if defined(USE_AVXI_OR_MIN_MAX_CMP)
//               const avxf boxDist = select(lhit,cast(lnearP),inf);
// #else
//               const avxf boxDist = select(lhit,lnearP,inf);
// #endif
              
              
//               sptr_node++;
//               sptr_near++;
              
              
//               if (any(boxDist < minDist))
//               {
//                 *(sptr_node-1) = curNode;
//                 *(sptr_near-1) = minDist; 
                
//                 minDist = boxDist;
//                 curNode = dd;
                
//               }
//               else
//               {
//                 *(sptr_node-1) = dd;
//                 *(sptr_near-1) = boxDist; 
//               }
              
//               assert(sptr_node - stack_node < QBVH_MAX_STACK_DEPTH);
//             }	      
//           }
          
//         }
        
// 	if (unlikely(curNode == QBVH_TERMINAL_TOKEN)) {break; }
        
// 	assert (curNode != QBVH_TERMINAL_TOKEN);
        
// 	unsigned int itemOffset    = qbvhItemOffset(curNode);
        
// 	const Triangle1  *__restrict__ tptr    = (Triangle1*)((char*)accel + itemOffset);
// 	const unsigned int items = qbvhItems(curNode);
//         STAT3(normal.trav_leaves,1,popcnt(m_active),8);
        
// 	assert(items > 0);
// 	assert(items <= 4);
        
// 	//prefetch<PFHINT_L1>(tptr +  1); 
// 	//prefetch<PFHINT_L1>(tptr +  2); 
// 	//prefetch<PFHINT_L1>(tptr +  3); 
        
// 	const avxf zero_ = avxf(zero);
        
// 	for(unsigned int i=0;i<items;i++,tptr++) 
// 	{
// 	  DBG(DBG_PRINT(*tptr));
          
//           STAT3(normal.trav_prims,1,popcnt(m_active),8);
//           const avxf v0 = broadcast4f(&tptr->v0);
//           const avxf v1 = broadcast4f(&tptr->v1);
//           const avxf v2 = broadcast4f(&tptr->v2);
//           const avxf e1 = v0-v1;
//           const avxf e2 = v2-v0;
          
//           const avxf v0_x = shuffle<0>(v0);
//           const avxf v0_y = shuffle<1>(v0);
//           const avxf v0_z = shuffle<2>(v0);
          
          
//           const avx3f org =  avx3f(avxf(v0_x) - orgx,avxf(v0_y) - orgy,avxf(v0_z) - orgz);
          
//           //const vec3f normal = cross(e1,e2);
          
//           const avxf nz = avxf(tptr->normal.z);
//           const avxf nx = avxf(tptr->normal.x);
//           const avxf ny = avxf(tptr->normal.y);
          
//           const avxf den = direction_z * nz + direction_y * ny + direction_x * nx;
          
//           DBG(DBG_PRINT(den));
          
//           const avxf sgnDet = signmsk(den);
//           const avxf absDet = abs(den);
          
//           //const avxf rcp_den = rcp(den);
          
          
//           const avxf odz = direction_x*org.y - org.x*direction_y;
//           const avxf odx = direction_y*org.z - org.y*direction_z;
//           const avxf ody = direction_z*org.x - org.z*direction_x;	  
          
//           const avxf e1_x = shuffle<0>(e1);
//           const avxf e1_y = shuffle<1>(e1);
//           const avxf e1_z = shuffle<2>(e1);
          
//           const avxf vv = odz*e1_z + ody*e1_y + odx*e1_x;
          
//           DBG(DBG_PRINT(vv));
          
//           const avxf e2_x = shuffle<0>(e2);
//           const avxf e2_y = shuffle<1>(e2);
//           const avxf e2_z = shuffle<2>(e2);
          
//           const avxf uu = odz*e2_z + ody*e2_y + odx*e2_x;
          
//           DBG(DBG_PRINT(uu));
          
//           const avxf u = uu^sgnDet; //uu * rcp_den;
//           const avxf v = vv^sgnDet; //vv * rcp_den;
          
          
//           const avxb valid_u = m_active & (u >= zero_);
//           const avxb valid_v = valid_u  & (v >= zero_);
          
          
//           const avxb m_aperture = valid_v & (u+v <= absDet); 
          
          
//           if ( unlikely(none(m_aperture)) ) continue;
          
//           const avxf rcp_den = rcp(den);
//           const avxf nom = org.z       * nz +       org.y * ny +       org.x * nx;
          
//           const avxf t = rcp_den*nom;
//           const avxb m_t_min = m_aperture & (t >= ray.tnear);
//           const avxb m_final  = m_t_min & (t < max_distance);
          
          
	  
//           if ( unlikely(none(m_final)) ) continue;
          
          
//           const avxf geomID(*(float*)&tptr->geomID);
//           const avxf primID(*(float*)&tptr->primID);
          
//           max_distance = select(m_final,t,max_distance);
          
//           store8f(m_final,(float*)&ray.geomID,geomID);
//           store8f(m_final,(float*)&ray.primID,primID);
	  
//           store8f(m_final,(float*)&ray.Ng.x,nx);
//           store8f(m_final,(float*)&ray.Ng.y,ny);
//           store8f(m_final,(float*)&ray.Ng.z,nz);
          
//           store8f(m_final,(float*)&ray.u,u*rcp_den);
//           store8f(m_final,(float*)&ray.v,v*rcp_den);
          
//           DBG(DBG_PRINT(u));
//           DBG(DBG_PRINT(v));
//           DBG(DBG_PRINT(nx));
//           DBG(DBG_PRINT(ny));
//           DBG(DBG_PRINT(nz));
//           DBG(DBG_PRINT(geomID));
//           DBG(DBG_PRINT(primID));
//           DBG(DBG_PRINT(t));
//           DBG(DBG_PRINT(ray));
	  
// 	}
        
// 	// ---------------------------------------------
// 	// ---------------------------------------------
        
//       }
      
//       DBG(DBG_PRINT(ray));
      
//       DBG(exit(0));
      
//       ray.tfar        = select(m_active,max_distance,ray.tfar);
  
  // while(1)
  //   {
  
  //     const size_t posN = bitscan(i_lr_hit);
  //     DBG(DBG_PRINT(posN));
  
  //     const unsigned int distN   = ((unsigned int*)&near8)[posN];
  //     const unsigned int nodeN   = bptr->min_d[posN];		  
  
  //     if (curNode_dist <= distN)
  //       {
  // 	stack[sindex].node = nodeN;
  // 	*(unsigned int*)&stack[sindex].dist = distN;
  //       }
  //     else
  //       {
  // 	stack[sindex].node = curNode;
  // 	*(unsigned int*)&stack[sindex].dist = curNode_dist;
  // 	curNode_dist = distN;
  // 	curNode = nodeN;
  //       }
  
  //     sindex++;
  //     assert(sindex < BVH8_MAX_STACK_DEPTH);
  
  //     i_lr_hit &= i_lr_hit-1;
  //     if (likely(i_lr_hit == 0)) break;
  
  //   }
  
