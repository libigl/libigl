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

#include "heuristic_binning.h"

namespace embree
{
  template<int logBlockSize>
  const size_t HeuristicBinning<logBlockSize>::maxBins;

  template<int logBlockSize>
  void HeuristicBinning<logBlockSize>::bin(const PrimRef* prims, size_t num)
  {
    if (num == 0) return;
    
    size_t i; for (i=0; i<num-1; i+=2)
    {
      /*! map even and odd primitive to bin */
      const BBox3f prim0 = prims[i+0].bounds(); const Vector3i bin0 = mapping.bin(prim0); const Vector3f center0 = Vector3f(center2(prim0));
      const BBox3f prim1 = prims[i+1].bounds(); const Vector3i bin1 = mapping.bin(prim1); const Vector3f center1 = Vector3f(center2(prim1));
      
      /*! increase bounds for bins for even primitive */
      const int b00 = bin0.x; counts[b00][0]++; geomBounds[b00][0].grow(prim0); centBounds[b00][0].grow(center0);
      const int b01 = bin0.y; counts[b01][1]++; geomBounds[b01][1].grow(prim0); centBounds[b01][1].grow(center0);
      const int b02 = bin0.z; counts[b02][2]++; geomBounds[b02][2].grow(prim0); centBounds[b02][2].grow(center0);
      
      /*! increase bounds of bins for odd primitive */
      const int b10 = bin1.x; counts[b10][0]++; geomBounds[b10][0].grow(prim1); centBounds[b10][0].grow(center1);
      const int b11 = bin1.y; counts[b11][1]++; geomBounds[b11][1].grow(prim1); centBounds[b11][1].grow(center1);
      const int b12 = bin1.z; counts[b12][2]++; geomBounds[b12][2].grow(prim1); centBounds[b12][2].grow(center1);
    }
    
    /*! for uneven number of primitives */
    if (i < num)
    {
      /*! map primitive to bin */
      const BBox3f prim0 = prims[i].bounds(); const Vector3i bin0 = mapping.bin(prim0); const Vector3f center0 = Vector3f(center2(prim0));
      
      /*! increase bounds of bins */
      const int b00 = bin0.x; counts[b00][0]++; geomBounds[b00][0].grow(prim0); centBounds[b00][0].grow(center0);
      const int b01 = bin0.y; counts[b01][1]++; geomBounds[b01][1].grow(prim0); centBounds[b01][1].grow(center0);
      const int b02 = bin0.z; counts[b02][2]++; geomBounds[b02][2].grow(prim0); centBounds[b02][2].grow(center0);
    }
  }
  
  template<int logBlockSize>
  void HeuristicBinning<logBlockSize>::best(Split& split)
  {
    Vector3f rAreas [maxBins];      //!< area of bounds of primitives on the right
    Vector3f rCounts[maxBins];      //!< blocks of primitives on the right
    
    /* sweep from right to left and compute parallel prefix of merged bounds */
    assert(mapping.size() > 0);
    Vector3i count = 0; BBox3f bx = empty; BBox3f by = empty; BBox3f bz = empty;
    for (size_t i=mapping.size()-1; i>0; i--)
    {
      count += counts[i];
      rCounts[i] = Vector3f(blocks(count));
      bx = merge(bx,geomBounds[i][0]); rAreas[i][0] = halfArea(bx);
      by = merge(by,geomBounds[i][1]); rAreas[i][1] = halfArea(by);
      bz = merge(bz,geomBounds[i][2]); rAreas[i][2] = halfArea(bz);
    }
    
    /* sweep from left to right and compute SAH */
    Vector3i ii = 1; Vector3f bestSAH = pos_inf; Vector3i bestPos = 0;
    count = 0; bx = empty; by = empty; bz = empty;
    for (size_t i=1; i<mapping.size(); i++, ii+=1)
    {
      count += counts[i-1];
      bx = merge(bx,geomBounds[i-1][0]); float Ax = halfArea(bx);
      by = merge(by,geomBounds[i-1][1]); float Ay = halfArea(by);
      bz = merge(bz,geomBounds[i-1][2]); float Az = halfArea(bz);
      const Vector3f lCount = Vector3f(blocks(count));
      const Vector3f lArea = Vector3f(Ax,Ay,Az);
      const Vector3f sah = lArea*lCount + rAreas[i]*rCounts[i];
      bestPos = select(lt_mask(sah,bestSAH),ii ,bestPos);
      bestSAH = select(lt_mask(sah,bestSAH),sah,bestSAH);
    }

    /* find best dimension */
    for (int i=0; i<3; i++) 
    {
      if (unlikely(pinfo.centBounds.lower[i] >= pinfo.centBounds.upper[i])) 
        continue;
      
      if (bestSAH[i] < split.cost && bestPos[i] != 0) {
        split.dim = i;
        split.pos = bestPos[i];
        split.cost = bestSAH[i];
      }
    }

    split.mapping = mapping;
    if (split.pos == 0) return;
    
    /* calculate geometry info from binning data */
    size_t pos = split.pos, dim = split.dim;
    size_t numLeft = 0, numRight = 0;
    BBox3f lcentBounds = empty, rcentBounds = empty;
    BBox3f lgeomBounds = empty, rgeomBounds = empty;
    for (size_t i=0; i<pos; i++) {
      numLeft += counts[i][dim];
      lcentBounds.grow(centBounds[i][dim]);
      lgeomBounds.grow(geomBounds[i][dim]);
    }
    for (size_t i=pos; i<mapping.size(); i++) {
      numRight += counts[i][dim];
      rcentBounds.grow(centBounds[i][dim]);
      rgeomBounds.grow(geomBounds[i][dim]);
    }
    assert(numLeft + numRight == pinfo.size());
    new (&split.linfo) PrimInfo(numLeft ,lgeomBounds,lcentBounds);
    new (&split.rinfo) PrimInfo(numRight,rgeomBounds,rcentBounds);
  }
  
  template<int logBlockSize>
  void HeuristicBinning<logBlockSize>::reduce(const HeuristicBinning binners[], size_t num, HeuristicBinning& binner_o)
  {
    binner_o = binners[0];
    for (size_t tid=1; tid<num; tid++) {
      const HeuristicBinning& binner = binners[tid];
      for (size_t bin=0; bin<binner.mapping.size(); bin++) 
      {
        for (size_t dim=0; dim<3; dim++) {
          binner_o.counts    [bin][dim] += binner.counts[bin][dim];
          binner_o.geomBounds[bin][dim].grow(binner.geomBounds[bin][dim]);
          binner_o.centBounds[bin][dim].grow(binner.centBounds[bin][dim]);
        }
      }
    }
  }
  
  /*! explicit template instantiations */
  template class HeuristicBinning<0>;
  template class HeuristicBinning<2>;
  template class HeuristicBinning<3>;
}
