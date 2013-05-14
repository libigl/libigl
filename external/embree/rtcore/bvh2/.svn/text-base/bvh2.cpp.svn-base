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

#include "bvh2.h"
#include "triangle/triangles.h"
#include "bvh2/bvh2_intersector.h"

namespace embree
{
  BVH2::BVH2 (const TriangleType& trity, const std::string& intTy, const Vec3fa* vertices, size_t numVertices, bool freeVertices) 
  : Accel(intTy), trity(trity), maxLeafTris(maxLeafBlocks*trity.blockSize), root(NULL), vertices(NULL), numVertices(0), freeVertices(freeVertices)
  {
    if (trity.needVertices) {
      this->vertices = vertices;
      this->numVertices = numVertices;
    }
  }

  BVH2::~BVH2 () 
  {
    if (freeVertices && vertices) 
      alignedFree(vertices); vertices  = NULL;
  }

  Ref<RefCount> BVH2::query(const char* name) 
  {
    if (trity.name == "triangle1i") {
      if (!strcmp(name,Intersector::name)) {
        if (intTy == "default" ) return new BVH2Intersector<Triangle1iIntersectorMoellerTrumbore>(this);
        if (intTy == "fast"    ) return new BVH2Intersector<Triangle1iIntersectorMoellerTrumbore>(this);
        if (intTy == "accurate") return new BVH2Intersector<Triangle1iIntersectorPluecker>(this);
        if (intTy == "moeller" ) return new BVH2Intersector<Triangle1iIntersectorMoellerTrumbore>(this);
        if (intTy == "pluecker") return new BVH2Intersector<Triangle1iIntersectorPluecker>(this);
        throw std::runtime_error("unknown triangle intersector \""+intTy+"\" for triangle1i");
      }
      throw std::runtime_error("unknown triangle intersector interface \""+std::string(name)+"\"");
    }
    else if (trity.name == "triangle4i") {
      if (!strcmp(name,Intersector::name)) {
        if (intTy == "default" ) return new BVH2Intersector<Triangle4iIntersectorMoellerTrumbore>(this);
        if (intTy == "fast"    ) return new BVH2Intersector<Triangle4iIntersectorMoellerTrumbore>(this);
        if (intTy == "accurate") return new BVH2Intersector<Triangle4iIntersectorPluecker>(this);
        if (intTy == "moeller" ) return new BVH2Intersector<Triangle4iIntersectorMoellerTrumbore>(this);
        if (intTy == "pluecker") return new BVH2Intersector<Triangle4iIntersectorPluecker>(this);
        throw std::runtime_error("unknown triangle intersector \""+intTy+"\" for triangle4i");
      }
      throw std::runtime_error("unknown triangle intersector interface \""+std::string(name)+"\"");
    }
    else if (trity.name == "triangle1v") {
      if (!strcmp(name,Intersector::name)) {
        if (intTy == "default" ) return new BVH2Intersector<Triangle1vIntersectorMoellerTrumbore>(this);
        if (intTy == "fast"    ) return new BVH2Intersector<Triangle1vIntersectorMoellerTrumbore>(this);
        if (intTy == "accurate") return new BVH2Intersector<Triangle1vIntersectorPluecker>(this);
        if (intTy == "moeller" ) return new BVH2Intersector<Triangle1vIntersectorMoellerTrumbore>(this);
        if (intTy == "pluecker") return new BVH2Intersector<Triangle1vIntersectorPluecker>(this);
        throw std::runtime_error("unknown triangle intersector \""+intTy+"\" for triangle1v");
      }
      throw std::runtime_error("unknown triangle intersector interface \""+std::string(name)+"\"");
    }
    else if (trity.name == "triangle4v") {
      if (!strcmp(name,Intersector::name)) {
        if (intTy == "default" ) return new BVH2Intersector<Triangle4vIntersectorMoellerTrumbore>(this);
        if (intTy == "fast"    ) return new BVH2Intersector<Triangle4vIntersectorMoellerTrumbore>(this);
        if (intTy == "accurate") return new BVH2Intersector<Triangle4vIntersectorPluecker>(this);
        if (intTy == "moeller" ) return new BVH2Intersector<Triangle4vIntersectorMoellerTrumbore>(this);
        if (intTy == "pluecker") return new BVH2Intersector<Triangle4vIntersectorPluecker>(this);
        throw std::runtime_error("unknown triangle intersector \""+intTy+"\" for triangle4v");
      }
      throw std::runtime_error("unknown triangle intersector interface \""+std::string(name)+"\"");
    }
    else if (trity.name == "triangle1") {
      if (!strcmp(name,Intersector::name)) {
        if (intTy == "default" ) return new BVH2Intersector<Triangle1IntersectorMoellerTrumbore>(this);
        if (intTy == "fast"    ) return new BVH2Intersector<Triangle1IntersectorMoellerTrumbore>(this);
        if (intTy == "moeller" ) return new BVH2Intersector<Triangle1IntersectorMoellerTrumbore>(this);
        throw std::runtime_error("unknown triangle intersector \""+intTy+"\" for triangle1");
      }
      throw std::runtime_error("unknown triangle intersector interface \""+std::string(name)+"\"");
    }
    else if (trity.name == "triangle4") {
      if (!strcmp(name,Intersector::name)) {
        if (intTy == "default" ) return new BVH2Intersector<Triangle4IntersectorMoellerTrumbore>(this);
        if (intTy == "fast"    ) return new BVH2Intersector<Triangle4IntersectorMoellerTrumbore>(this);
        if (intTy == "moeller" ) return new BVH2Intersector<Triangle4IntersectorMoellerTrumbore>(this);
        throw std::runtime_error("unknown triangle intersector \""+intTy+"\" for triangle4");
      }
      throw std::runtime_error("unknown triangle intersector interface \""+std::string(name)+"\"");
    }
    else if (trity.name == "triangle8") {
      if (!strcmp(name,Intersector::name)) {
        if (intTy == "default") return new BVH2Intersector<Triangle8IntersectorMoellerTrumbore>(this);
        if (intTy == "fast"   ) return new BVH2Intersector<Triangle8IntersectorMoellerTrumbore>(this);
        if (intTy == "moeller") return new BVH2Intersector<Triangle8IntersectorMoellerTrumbore>(this);
        throw std::runtime_error("unknown triangle intersector \""+intTy+"\" for triangle8");
      }
      throw std::runtime_error("unknown triangle intersector interface \""+std::string(name)+"\"");
    }
    throw std::runtime_error("unknown BVH2 triangle type \""+std::string(trity.name)+"\"");
    return null;
  }

  size_t BVH2::rotate(Base* node, size_t depth)
  {
    /*! nothing to rotate if we reached a leaf node. */
    if (node->isBarrier()) return maxLocalDepth;
    if (node->isLeaf()) return 0;
    Node* parent = node->node();

    /*! rotate all children first */
    size_t cdepth[2];
    for (size_t c=0; c<2; c++)
      cdepth[c] = rotate(parent->child[c],depth+1);

    /*! compute current area of all children */
    ssef sizeX = shuffle<2,3,0,1>(parent->lower_upper_x)-parent->lower_upper_x;
    ssef sizeY = shuffle<2,3,0,1>(parent->lower_upper_y)-parent->lower_upper_y;
    ssef sizeZ = shuffle<2,3,0,1>(parent->lower_upper_z)-parent->lower_upper_z;
    ssef childArea = sizeX*(sizeY + sizeZ) + sizeY*sizeZ;

    /*! transpose node bounds */
    ssef plower0,plower1,pupper0,pupper1;
    transpose(parent->lower_upper_x,parent->lower_upper_y,parent->lower_upper_z,ssef(zero),plower0,plower1,pupper0,pupper1);
    BBox<ssef> others[2] = { BBox<ssef>(plower0,pupper0), BBox<ssef>(plower1,pupper1) };

    /*! Find best rotation. We pick a target child of a first child,
      and swap this with an other child. We perform the best such
      swap. */
    float bestCost = pos_inf;
    int bestChild = -1, bestTarget = -1, bestOther = -1;
    for (size_t c=0; c<2; c++)
    {
      /*! ignore leaf nodes as we cannot descent into */
      if (parent->child[c]->isBarrier()) continue;
      if (parent->child[c]->isLeaf()) continue;

      /*! only select swaps that fulfill depth constraints */
      bool valid = depth+1+cdepth[1-c] <= maxDepth; 
      if (!valid) continue;

      Node* child = parent->child[c]->node();
      BBox<ssef>& other = others[1-c];

      /*! transpose child bounds */
      ssef clower0,clower1,cupper0,cupper1;
      transpose(child->lower_upper_x,child->lower_upper_y,child->lower_upper_z,ssef(zero),clower0,clower1,cupper0,cupper1);
      BBox<ssef> target0(clower0,cupper0), target1(clower1,cupper1);
      
      /*! compute cost for both possible swaps */
      float cost0 = halfArea3f(merge(other ,target1))-childArea[c];
      float cost1 = halfArea3f(merge(target0,other ))-childArea[c];
      
      if (min(cost0,cost1) < bestCost)
      {
        bestChild = (int)c;
        bestOther = (int)(1-c);
        if (cost0 < cost1) {
          bestCost = cost0;
          bestTarget = 0;
        } else {
          bestCost = cost0;
          bestTarget = 1;
        }
      }
    }

    /*! if we did not find a swap that improves the SAH then do nothing */
    if (bestCost >= 0) return 1+max(cdepth[0],cdepth[1]);

    /*! perform the best found tree rotation */
    Node* child = parent->child[bestChild]->node();
    swap(parent,bestOther,child,bestTarget);
    parent->lower_upper_x[bestChild+0] = min(child->lower_upper_x[0],child->lower_upper_x[1]);
    parent->lower_upper_y[bestChild+0] = min(child->lower_upper_y[0],child->lower_upper_y[1]);
    parent->lower_upper_z[bestChild+0] = min(child->lower_upper_z[0],child->lower_upper_z[1]);
    parent->lower_upper_x[bestChild+2] = max(child->lower_upper_x[2],child->lower_upper_x[3]);
    parent->lower_upper_y[bestChild+2] = max(child->lower_upper_y[2],child->lower_upper_y[3]);
    parent->lower_upper_z[bestChild+2] = max(child->lower_upper_z[2],child->lower_upper_z[3]);

    /*! This returned depth is conservative. */
    cdepth[bestOther]++; // bestOther was pushed down one level
    return 1+max(cdepth[0],cdepth[1]); 
  }

  float BVH2::sort(Base* node, int maxDepth)
  {
    if (node->isBarrier()) 
      return 0.0f;

    if (node->isLeaf()) 
    {
      float A = 0.0f;
      size_t num; char* tri = node->leaf(num);
      for (size_t i=0; i<num; i++)
        A += trity.area(tri+i*trity.bytes,vertices);
      return A;
    }
    else 
    {
      Node* n = node->node();
      const ssef dx = shuffle<2,3,0,1>(n->lower_upper_x)-n->lower_upper_x;
      const ssef dy = shuffle<2,3,0,1>(n->lower_upper_y)-n->lower_upper_y;
      const ssef dz = shuffle<2,3,0,1>(n->lower_upper_z)-n->lower_upper_z;
      const ssef area = dx*(dy+dz)+dy*dz;
      if (maxDepth < 0) return area[0]+area[1];
      const ssef rcpBoundArea = rcp(area);

      /*! recurse into children first */
      float A = 0.0f;
      float opacity[2];
      for (size_t c=0; c<2; c++) {
        float dA = sort(n->child[c],maxDepth-1);
        opacity[c] = dA*rcpBoundArea[c];
        A += dA;
      }
      
      /*! sort children */
      if (opacity[0] > opacity[1]) swap(n,0,n,1);
      
      return A;
    }
  }

  void BVH2::clearBarrier(Base*& node)
  {
    if (node->isBarrier()) {
      node = node->clearBarrier();
      return;
    }
    else if (node->isLeaf()) 
      return;
    else {
      Node* n = node->node();
      for (size_t c=0; c<2; c++)
        clearBarrier(n->child[c]);
    }
  }

  float BVH2::statistics(Base* node, float ap, size_t& depth)
  {
    if (node->isNode())
    {
      numNodes++;
      depth = 0;
      size_t cdepth = 0;
      const Node* n = node->node();
      const float a0 = area(n->bounds(0));
      const float a1 = area(n->bounds(1));
      const float sah0 = statistics(n->child[0],a0,cdepth); depth=max(depth,cdepth);
      const float sah1 = statistics(n->child[1],a1,cdepth); depth=max(depth,cdepth);
      depth++;
      return ap*travCost + sah0 + sah1;
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

  void BVH2::print(std::ostream& cout)
  {
    /* calculate statistics */
    numNodes = numLeaves = numPrimBlocks = numPrims = depth = 0;
    bvhSAH = statistics(root,0.0f,depth);
    assert(depth <= BVH2::maxDepth);

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
         << "(" << 100.0*(numNodes-1+numLeaves)/(2.0*numNodes) << "% used)" 
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
