// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "bvh4mb.h"
#include "../triangle/triangles.h"
#include "bvh4mb/bvh4mb_intersector.h"

namespace embree
{
  BVH4MB::BVH4MB (const TriangleType& trity, const std::string& intTy, const Vec3fa* vertices, size_t numVertices, bool freeVertices) 
  : Accel(intTy), trity(trity), maxLeafTris(maxLeafBlocks*trity.blockSize), root(NULL), vertices(NULL), numVertices(0), freeVertices(freeVertices)
  {
    if (trity.needVertices) {
      this->vertices = vertices;
      this->numVertices = numVertices;
    }
  }

  BVH4MB::~BVH4MB () 
  {
    if (freeVertices && vertices) 
      alignedFree(vertices); vertices  = NULL;
  }

  Ref<RefCount> BVH4MB::query(const char* name)
  {
    if (trity.name == "triangle4i") {
      if (!strcmp(name,Intersector::name)) {
        if (intTy == "default" ) return new BVH4MBIntersector<Triangle4iIntersectorMoellerTrumboreMB>(this);
        if (intTy == "fast"    ) return new BVH4MBIntersector<Triangle4iIntersectorMoellerTrumboreMB>(this);
        if (intTy == "accurate") return new BVH4MBIntersector<Triangle4iIntersectorPlueckerMB>(this);
        if (intTy == "moeller" ) return new BVH4MBIntersector<Triangle4iIntersectorMoellerTrumboreMB>(this);
        if (intTy == "pluecker") return new BVH4MBIntersector<Triangle4iIntersectorPlueckerMB>(this);
        throw std::runtime_error("unknown triangle intersector \""+intTy+"\" for triangle4i");
      }
      throw std::runtime_error("unknown triangle intersector interface \""+std::string(name)+"\"");
    }
    throw std::runtime_error("unknown BVH4MB triangle type \""+std::string(trity.name)+"\"");
    return null;
  }

  size_t BVH4MB::rotate(Base* nodeID, size_t depth)
  {
    /*! nothing to rotate if we reached a leaf node. */
    if (nodeID->isLeaf()) return 0;
    Node* parent = nodeID->node();

    /*! rotate all children first */
    ssei cdepth;
    for (size_t c=0; c<4; c++)
      cdepth[c] = (int)rotate(parent->child[c],depth+1);

    /* compute current area of all children */
    ssef sizeX = parent->upper_x-parent->lower_x;
    ssef sizeY = parent->upper_y-parent->lower_y;
    ssef sizeZ = parent->upper_z-parent->lower_z;
    ssef childArea = sizeX*(sizeY + sizeZ) + sizeY*sizeZ;

    /*! transpose node bounds */
    ssef plower0,plower1,plower2,plower3; transpose(parent->lower_x,parent->lower_y,parent->lower_z,ssef(zero),plower0,plower1,plower2,plower3);
    ssef pupper0,pupper1,pupper2,pupper3; transpose(parent->upper_x,parent->upper_y,parent->upper_z,ssef(zero),pupper0,pupper1,pupper2,pupper3);
    BBox<ssef> other0(plower0,pupper0), other1(plower1,pupper1), other2(plower2,pupper2), other3(plower3,pupper3);

    /*! Find best rotation. We pick a target child of a first child,
      and swap this with an other child. We perform the best such
      swap. */
    float bestCost = pos_inf;
    int bestChild = -1, bestTarget = -1, bestOther = -1;
    for (size_t c=0; c<4; c++)
    {
      /*! ignore leaf nodes as we cannot descent into */
      if (parent->child[c]->isLeaf()) continue;
      Node* child = parent->child[c]->node();

      /*! transpose child bounds */
      ssef clower0,clower1,clower2,clower3; transpose(child->lower_x,child->lower_y,child->lower_z,ssef(zero),clower0,clower1,clower2,clower3);
      ssef cupper0,cupper1,cupper2,cupper3; transpose(child->upper_x,child->upper_y,child->upper_z,ssef(zero),cupper0,cupper1,cupper2,cupper3);
      BBox<ssef> target0(clower0,cupper0), target1(clower1,cupper1), target2(clower2,cupper2), target3(clower3,cupper3);

      /*! put other0 at each target position */
      float cost00 = halfArea3f(merge(other0 ,target1,target2,target3));
      float cost01 = halfArea3f(merge(target0,other0 ,target2,target3));
      float cost02 = halfArea3f(merge(target0,target1,other0 ,target3));
      float cost03 = halfArea3f(merge(target0,target1,target2,other0 ));
      ssef cost0 = ssef(cost00,cost01,cost02,cost03);
      ssef min0 = vreduce_min(cost0);
      int pos0 = (int)__bsf(movemask(min0 == cost0));

      /*! put other1 at each target position */
      float cost10 = halfArea3f(merge(other1 ,target1,target2,target3));
      float cost11 = halfArea3f(merge(target0,other1 ,target2,target3));
      float cost12 = halfArea3f(merge(target0,target1,other1 ,target3));
      float cost13 = halfArea3f(merge(target0,target1,target2,other1 ));
      ssef cost1 = ssef(cost10,cost11,cost12,cost13);
      ssef min1 = vreduce_min(cost1);
      int pos1 = (int)__bsf(movemask(min1 == cost1));

      /*! put other2 at each target position */
      float cost20 = halfArea3f(merge(other2 ,target1,target2,target3));
      float cost21 = halfArea3f(merge(target0,other2 ,target2,target3));
      float cost22 = halfArea3f(merge(target0,target1,other2 ,target3));
      float cost23 = halfArea3f(merge(target0,target1,target2,other2 ));
      ssef cost2 = ssef(cost20,cost21,cost22,cost23);
      ssef min2 = vreduce_min(cost2);
      int pos2 = (int)__bsf(movemask(min2 == cost2));

      /*! put other3 at each target position */
      float cost30 = halfArea3f(merge(other3 ,target1,target2,target3));
      float cost31 = halfArea3f(merge(target0,other3 ,target2,target3));
      float cost32 = halfArea3f(merge(target0,target1,other3 ,target3));
      float cost33 = halfArea3f(merge(target0,target1,target2,other3 ));
      ssef cost3 = ssef(cost30,cost31,cost32,cost33);
      ssef min3 = vreduce_min(cost3);
      int pos3 = (int)__bsf(movemask(min3 == cost3));

      /*! find best other child */
      ssef otherCost = ssef(extract<0>(min0),extract<0>(min1),extract<0>(min2),extract<0>(min3));
      int pos[4] = { pos0,pos1,pos2,pos3 };
      sseb valid = ssei(int(depth+1))+cdepth <= ssei(maxDepth); // only select swaps that fulfill depth constraints
      if (none(valid)) continue;
      
      size_t n = select_min(valid,otherCost);
      float cost = otherCost[n]-childArea[c]; //< increasing the original child bound is bad, decreasing good

      /*! accept a swap when it reduces cost and is not swapping a node with itself */
      if (cost < bestCost && n != c) {
        bestCost = cost;
        bestChild = (int)c;
        bestOther = (int)n;
        bestTarget = pos[n];
      }
    }

    /*! if we did not find a swap that improves the SAH then do nothing */
    if (bestCost >= 0) return 1+reduce_max(cdepth);

    /*! perform the best found tree rotation */
    Node* child = parent->child[bestChild]->node();
    swap(parent,bestOther,child,bestTarget);
    parent->lower_x[bestChild] = reduce_min(child->lower_x);
    parent->lower_y[bestChild] = reduce_min(child->lower_y);
    parent->lower_z[bestChild] = reduce_min(child->lower_z);
    parent->upper_x[bestChild] = reduce_max(child->upper_x);
    parent->upper_y[bestChild] = reduce_max(child->upper_y);
    parent->upper_z[bestChild] = reduce_max(child->upper_z);
    parent->lower_dx[bestChild] = reduce_min(child->lower_dx);
    parent->lower_dy[bestChild] = reduce_min(child->lower_dy);
    parent->lower_dz[bestChild] = reduce_min(child->lower_dz);
    parent->upper_dx[bestChild] = reduce_max(child->upper_dx);
    parent->upper_dy[bestChild] = reduce_max(child->upper_dy);
    parent->upper_dz[bestChild] = reduce_max(child->upper_dz);

    /*! This returned depth is conservative as the child that was
     *  pulled up in the tree could have been on the critical path. */
    cdepth[bestOther]++; // bestOther was pushed down one level
    return 1+reduce_max(cdepth); 
  }

  float BVH4MB::sort(Base* node, int maxDepth)
  {
    /*! sum up triangle areas */
    if (node->isLeaf()) 
    {
      float A = 0.0f;
      size_t num; char* tri = node->leaf(num);
      for (size_t i=0; i<num; i++)
        A += trity.area(tri+i*trity.bytes,vertices);
      return A;
    }
    /*! sort node children based on approximated opacity */
    else 
    {
      /*! calculate half of node surface area */
      Node* n = node->node();
      const ssef dx = n->upper_x-n->lower_x;
      const ssef dy = n->upper_y-n->lower_y;
      const ssef dz = n->upper_z-n->lower_z;
      const ssef area = dx*(dy+dz)+dy*dz;
      if (maxDepth < 0) return reduce_add(area);
      const ssef rcpBoundArea = rcp(area);
      
      /*! recurse into children first */
      float A = 0.0f;
      float opacity[4];
      for (size_t c=0; c<4; c++) {
        float dA = sort(n->child[c],maxDepth-1);
        opacity[c] = dA*rcpBoundArea[c];
        A += dA;
      }
      
      /*! sort children */
      for (size_t i=3; i>0; i--) {
        for (size_t j=0; j<i; j++) {
          if (opacity[j+0] > opacity[j+1]) {
            std::swap(opacity[j+0],opacity[j+1]);
            swap(n,j+0,n,j+1);
          }
        }
      }

      return A;
    }
  }

  std::pair<BBox3f,BBox3f> BVH4MB::refit(Base* node)
  {
    /*! merge bounds of triangles for both time steps */
    if (node->isLeaf()) 
    {
      BBox3f bounds0 = empty;
      BBox3f bounds1 = empty;
      size_t num; char* tri = node->leaf(num);
      for (size_t i=0; i<num; i++) {
        std::pair<BBox3f,BBox3f> bounds = trity.bounds(tri+i*trity.bytes,vertices);
        bounds0.grow(bounds.first);
        bounds1.grow(bounds.second);
      }
      return std::pair<BBox3f,BBox3f>(bounds0,bounds1);
    }
    /*! set and propagate merged bounds for both time steps */
    else
    {
      Node* n = node->node();
      if (!n->hasBounds()) {
        for (size_t i=0; i<4; i++) {
          std::pair<BBox3f,BBox3f> bounds = refit(n->child[i]);
          n->set(i,bounds.first,bounds.second);
        }
      }
      BBox3f bounds0 = merge(n->bounds0(0),n->bounds0(1),n->bounds0(2),n->bounds0(3));
      BBox3f bounds1 = merge(n->bounds1(0),n->bounds1(1),n->bounds1(2),n->bounds1(3));
      return std::pair<BBox3f,BBox3f>(bounds0,bounds1);
    }
  }

  float BVH4MB::statistics(Base* node, float ap, size_t& depth)
  {
    if (node->isNode())
    {
      numNodes++;
      depth = 0;
      size_t cdepth = 0;
      const Node* n = node->node();
      const ssef dx = n->upper_x-n->lower_x;
      const ssef dy = n->upper_y-n->lower_y;
      const ssef dz = n->upper_z-n->lower_z;
      const ssef ac = 2.0f*(dx*(dy+dz)+dy*dz);
      const float sah0 = statistics(n->child[0],ac[0],cdepth); depth=max(depth,cdepth);
      const float sah1 = statistics(n->child[1],ac[1],cdepth); depth=max(depth,cdepth);
      const float sah2 = statistics(n->child[2],ac[2],cdepth); depth=max(depth,cdepth);
      const float sah3 = statistics(n->child[3],ac[3],cdepth); depth=max(depth,cdepth);
      depth++;
      return ap*travCost + sah0 + sah1 + sah2 + sah3;
    }
    else
    {
      depth = 0;
      size_t num; char* tri = node->leaf(num);
      if (!num) return 0.0f;

      numLeaves++;
      numPrimBlocks += num;
      for (size_t i=0; i<num; i++)
        numPrims += trity.size(tri+i*trity.bytes);
      return trity.intCost * ap * num;
    }
  }

  void BVH4MB::print(std::ostream& cout)
  {
    /* calculate statistics */
    numNodes = numLeaves = numPrimBlocks = numPrims = depth = 0;
    bvhSAH = statistics(root,0.0f,depth);

    /* output statistics */
    std::ios::fmtflags flags = std::cout.flags();
    size_t bytesNodes = numNodes     *sizeof(Node);
    size_t bytesTris  = numPrimBlocks*trity.bytes;
    size_t bytesVertices = numVertices*sizeof(Vec3f);
    size_t bytesTotal = bytesNodes+bytesTris+bytesVertices;
    cout.setf(std::ios::scientific, std::ios::floatfield);
    cout.precision(2);
    cout << "sah = " << bvhSAH << std::endl;
    cout.setf(std::ios::fixed, std::ios::floatfield);
    cout.precision(1);
    cout << "depth = " << depth << std::endl;
    cout << "size = " << bytesTotal/1E6 << " MB" << std::endl;
    cout << "nodes = "  << numNodes << " "
         << "(" << bytesNodes/1E6  << " MB) "
         << "(" << 100.0*double(bytesNodes)/double(bytesTotal) << "% of total) "
         << "(" << 100.0*(numNodes-1+numLeaves)/(4.0*numNodes) << "% used)" 
         << std::endl;
    cout << "leaves = " << numLeaves << " "
         << "(" << bytesTris/1E6  << " MB) "
         << "(" << 100.0*double(bytesTris)/double(bytesTotal) << "% of total) "
         << "(" << 100.0*numPrims/(trity.blockSize*numPrimBlocks) << "% used)" 
         << std::endl;
    cout << "vertices = " << numVertices << " "
         << "(" << bytesVertices/1E6 << " MB) " 
         << "(" << 100.0*double(bytesVertices)/double(bytesTotal) << "% of total) "
         << "(" << 100.0*12.0f/float(sizeof(Vec3f)) << "% used)" 
         << std::endl;
    cout.setf(flags);
  }
}
