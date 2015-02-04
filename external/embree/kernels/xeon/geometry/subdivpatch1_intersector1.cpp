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

#include "subdivpatch1_intersector1.h"
#include "bicubic_bezier_patch.h"

namespace embree
{
  namespace isa
  {  
    void SubdivPatch1Intersector1::subdivide_intersect1_bspline(const Precalculations& pre,
                                                                Ray& ray,
                                                                const BSplinePatch &patch,
                                                                const unsigned int geomID,
                                                                const unsigned int primID,
                                                                const Vec2f &s,
                                                                const Vec2f &t,
                                                                const unsigned int subdiv_level)
    {
      if (subdiv_level == 0)
      {
        Vec3fa vtx[4];
        vtx[0] = patch.eval(s[0],t[0]);
        
        vtx[1] = patch.eval(s[1],t[0]);
        vtx[2] = patch.eval(s[1],t[1]);
        vtx[3] = patch.eval(s[0],t[1]);
        intersectTri(vtx[0],
                     vtx[1],
                     vtx[2],
                     ray,
                     geomID,
                     primID,NULL); 
        
        intersectTri(vtx[2],
                     vtx[3],
                     vtx[0],
                     ray,
                     geomID,
                     primID,NULL); 
      }
      else
      {
        const float mid_s = 0.5f * (s[0]+s[1]);
        const float mid_t = 0.5f * (t[0]+t[1]);
        Vec2f s_left(s[0],mid_s);
        Vec2f s_right(mid_s,s[1]);
        Vec2f t_left(t[0],mid_t);
        Vec2f t_right(mid_t,t[1]);
        subdivide_intersect1_bspline(pre,ray,patch,geomID,primID,s_left ,t_left,subdiv_level - 1);
        subdivide_intersect1_bspline(pre,ray,patch,geomID,primID,s_right,t_left,subdiv_level - 1);
        subdivide_intersect1_bspline(pre,ray,patch,geomID,primID,s_right,t_right,subdiv_level - 1);
        subdivide_intersect1_bspline(pre,ray,patch,geomID,primID,s_left ,t_right,subdiv_level - 1);
      }
    }
    
    void SubdivPatch1Intersector1::subdivide_intersect1(const Precalculations& pre,
                                                        Ray& ray,
                                                        const CatmullClarkPatch &patch,
                                                        const unsigned int geomID,
                                                        const unsigned int primID,
                                                        const unsigned int subdiv_level)
    {
      if (subdiv_level == 0)
      {
        __aligned(64) FinalQuad finalQuad;
        patch.init( finalQuad );
        
        intersectTri(finalQuad.vtx[0],
                     finalQuad.vtx[1],
                     finalQuad.vtx[2],
                     ray,
                     geomID,
                     primID,NULL); 
        
        intersectTri(finalQuad.vtx[2],
                     finalQuad.vtx[3],
                     finalQuad.vtx[0],
                     ray,
                     geomID,
                     primID,NULL); 
      }
      else
      {
        CatmullClarkPatch subpatches[4];
        patch.subdivide(subpatches);
        for (size_t i=0;i<4;i++)
          if (intersectBounds(pre,ray,subpatches[i].bounds()))
            subdivide_intersect1(pre,
                                 ray,
                                 subpatches[i],
                                 geomID,
                                 primID,
                                 subdiv_level - 1);	    
      }   
    }
    
    bool SubdivPatch1Intersector1::subdivide_occluded1(const Precalculations& pre,
                                                       Ray& ray,
                                                       const CatmullClarkPatch &patch,
                                                       const unsigned int geomID,
                                                       const unsigned int primID,						     
                                                       const unsigned int subdiv_level)
    {
      if (subdiv_level == 0)
      {
        __aligned(64) FinalQuad finalQuad;
        patch.init( finalQuad );
        
        if (occludedTri(finalQuad.vtx[0],
                        finalQuad.vtx[1],
                        finalQuad.vtx[2],
                        ray,
                        geomID,
                        primID,NULL)) return true; 
        
        if (occludedTri(finalQuad.vtx[2],
                        finalQuad.vtx[3],
                        finalQuad.vtx[0],
                        ray,
                        geomID,
                        primID,NULL)) return false;
      }
      else
      {
        CatmullClarkPatch subpatches[4];
        patch.subdivide(subpatches);
        for (size_t i=0;i<4;i++)
          if (intersectBounds(pre,ray,subpatches[i].bounds()))
            if (subdivide_occluded1(pre,
                                    ray,
                                    subpatches[i],
                                    geomID,
                                    primID,
                                    subdiv_level - 1)) return true;
      }   
      return false;
    }
    
    
    
    void SubdivPatch1Intersector1::subdivide_intersect1(const Precalculations& pre,
                                                        Ray& ray,
                                                        const BSplinePatch &patch,
                                                        const unsigned int geomID,
                                                        const unsigned int primID,
                                                        const unsigned int subdiv_level)
    {
      if (subdiv_level == 0)
      {
        __aligned(64) FinalQuad finalQuad;
        patch.init( finalQuad );
        
        intersectTri(finalQuad.vtx[0],
                     finalQuad.vtx[1],
                     finalQuad.vtx[2],
                     ray,
                     geomID,
                     primID,NULL); 
        
        intersectTri(finalQuad.vtx[2],
                     finalQuad.vtx[3],
                     finalQuad.vtx[0],
                     ray,
                     geomID,
                     primID,NULL); 
      }
      else
      {
        BSplinePatch subpatches[4];
        patch.subdivide(subpatches);
        for (size_t i=0;i<4;i++)
          if (intersectBounds(pre,ray,subpatches[i].bounds()))
            subdivide_intersect1(pre,
                                 ray,
                                 subpatches[i],
                                 geomID,
                                 primID,
                                 subdiv_level - 1);	    
        
      } 
    }
    
    bool SubdivPatch1Intersector1::subdivide_occluded1(const Precalculations& pre,
                                                       Ray& ray,
                                                       const BSplinePatch &patch,
                                                       const unsigned int geomID,
                                                       const unsigned int primID,
                                                       const unsigned int subdiv_level)
    {
      if (subdiv_level == 0)
      {
        __aligned(64) FinalQuad finalQuad;
        patch.init( finalQuad );
        
        if (occludedTri(finalQuad.vtx[0],
                        finalQuad.vtx[1],
                        finalQuad.vtx[2],
                        ray,
                        geomID,
                        primID,NULL)) return true; 
        
        if (occludedTri(finalQuad.vtx[2],
                        finalQuad.vtx[3],
                        finalQuad.vtx[0],
                        ray,
                        geomID,
                        primID,NULL)) return false;
      }
      else
      {
        BSplinePatch subpatches[4];
        patch.subdivide(subpatches);
        for (size_t i=0;i<4;i++)
          if (intersectBounds(pre,ray,subpatches[i].bounds()))
            if (subdivide_occluded1(pre,
                                    ray,
                                    subpatches[i],
                                    geomID,
                                    primID,
                                    subdiv_level - 1)) return true;
      }   
      return false;
    }
  };
}
