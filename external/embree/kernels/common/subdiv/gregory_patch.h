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

#include "catmullclark_patch.h"

namespace embree
{  
  class __aligned(64) GregoryPatch : public BSplinePatchT<Vec3fa> 
  {
  public:
    
    Vec3fa f[2][2]; // need 16 + 4 = 20 control points
    
    GregoryPatch() {
      memset(this,0,sizeof(GregoryPatch));
    }
    
    GregoryPatch(const Vec3fa matrix[4][4],
		 const Vec3fa f_m[2][2]) 
    {
      for (size_t y=0;y<4;y++)
	for (size_t x=0;x<4;x++)
	  v[y][x] = (Vec3fa_t)matrix[y][x];
      
      for (size_t y=0;y<2;y++)
	for (size_t x=0;x<2;x++)
	  f[y][x] = (Vec3fa_t)f_m[y][x];
    }
    
    Vec3fa& p0() { return v[0][0]; }
    Vec3fa& p1() { return v[0][3]; }
    Vec3fa& p2() { return v[3][3]; }
    Vec3fa& p3() { return v[3][0]; }
    
    Vec3fa& e0_p() { return v[0][1]; }
    Vec3fa& e0_m() { return v[1][0]; }
    Vec3fa& e1_p() { return v[1][3]; }
    Vec3fa& e1_m() { return v[0][2]; }
    Vec3fa& e2_p() { return v[3][2]; }
    Vec3fa& e2_m() { return v[2][3]; }
    Vec3fa& e3_p() { return v[2][0]; }
    Vec3fa& e3_m() { return v[3][1]; }
    
    Vec3fa& f0_p() { return v[1][1]; }
    Vec3fa& f1_p() { return v[1][2]; }
    Vec3fa& f2_p() { return v[2][2]; }
    Vec3fa& f3_p() { return v[2][1]; }
    Vec3fa& f0_m() { return f[0][0]; }
    Vec3fa& f1_m() { return f[0][1]; }
    Vec3fa& f2_m() { return f[1][1]; }
    Vec3fa& f3_m() { return f[1][0]; }
    
    const Vec3fa& p0() const { return v[0][0]; }
    const Vec3fa& p1() const { return v[0][3]; }
    const Vec3fa& p2() const { return v[3][3]; }
    const Vec3fa& p3() const { return v[3][0]; }
    
    const Vec3fa& e0_p() const { return v[0][1]; }
    const Vec3fa& e0_m() const { return v[1][0]; }
    const Vec3fa& e1_p() const { return v[1][3]; }
    const Vec3fa& e1_m() const { return v[0][2]; }
    const Vec3fa& e2_p() const { return v[3][2]; }
    const Vec3fa& e2_m() const { return v[2][3]; }
    const Vec3fa& e3_p() const { return v[2][0]; }
    const Vec3fa& e3_m() const { return v[3][1]; }
    
    const Vec3fa& f0_p() const { return v[1][1]; }
    const Vec3fa& f1_p() const { return v[1][2]; }
    const Vec3fa& f2_p() const { return v[2][2]; }
    const Vec3fa& f3_p() const { return v[2][1]; }
    const Vec3fa& f0_m() const { return f[0][0]; }
    const Vec3fa& f1_m() const { return f[0][1]; }
    const Vec3fa& f2_m() const { return f[1][1]; }
    const Vec3fa& f3_m() const { return f[1][0]; }
    
    
    Vec3fa initCornerVertex(const CatmullClarkPatch &irreg_patch, const size_t index)
    {
      return irreg_patch.ring[index].getLimitVertex();
    }
    
    
    Vec3fa initPositiveEdgeVertex(const CatmullClarkPatch &irreg_patch, const size_t index, const Vec3fa &p_vtx)
    {
      const Vec3fa tangent = irreg_patch.ring[index].getLimitTangent();
      return 1.0f/3.0f * tangent + p_vtx;
    }
    
    Vec3fa initNegativeEdgeVertex(const CatmullClarkPatch &irreg_patch, const size_t index, const Vec3fa &p_vtx)
    {
      const Vec3fa tangent = irreg_patch.ring[index].getSecondLimitTangent();
      return 1.0f/3.0f * tangent + p_vtx;
    }
    
    
    void initFaceVertex(const CatmullClarkPatch &irreg_patch,
			const size_t index,
			const Vec3fa &p_vtx,
			const Vec3fa &e0_p_vtx,
			const Vec3fa &e1_m_vtx,
			const unsigned int face_valence_p1,
			const Vec3fa &e0_m_vtx,
			const Vec3fa &e3_p_vtx,
			const unsigned int face_valence_p3,
			Vec3fa &f_p_vtx,
			Vec3fa &f_m_vtx)
    {
      const unsigned int face_valence         = irreg_patch.ring[index].face_valence;
      const unsigned int edge_valence         = irreg_patch.ring[index].edge_valence;
      const unsigned int border_index = irreg_patch.ring[index].border_index;
      
      const Vec3fa &vtx     = irreg_patch.ring[index].vtx;
      const Vec3fa e_i      = irreg_patch.ring[index].getEdgeCenter( 0 );
      const Vec3fa c_i_m_1  = irreg_patch.ring[index].getQuadCenter( 0 );
      const Vec3fa e_i_m_1  = irreg_patch.ring[index].getEdgeCenter( 1 );
      
      Vec3fa c_i, e_i_p_1;
      if (unlikely(border_index == edge_valence-2))
      {
        /* mirror quad center and edge mid-point */
        c_i     = c_i_m_1 + 2 * (e_i - c_i_m_1);
        e_i_p_1 = e_i_m_1 + 2 * (vtx - e_i_m_1);
      }
      else
      {
        c_i     = irreg_patch.ring[index].getQuadCenter( face_valence-1 );
        e_i_p_1 = irreg_patch.ring[index].getEdgeCenter( face_valence-1 );
      }
      
      Vec3fa c_i_m_2, e_i_m_2;
      if (unlikely(border_index == 2 || face_valence == 2))
      {
        /* mirror quad center and edge mid-point */
        c_i_m_2  = c_i_m_1 + 2 * (e_i_m_1 - c_i_m_1);
        e_i_m_2  = e_i + 2 * (vtx - e_i);	  
      }
      else
      {
        c_i_m_2  = irreg_patch.ring[index].getQuadCenter( 1 );
        e_i_m_2  = irreg_patch.ring[index].getEdgeCenter( 2 );
      }      
      
      const float d = 3.0f;
      const float c     = cosf(2.0*M_PI/(float)face_valence);
      const float c_e_p = cosf(2.0*M_PI/(float)face_valence_p1);
      const float c_e_m = cosf(2.0*M_PI/(float)face_valence_p3);
      
      const Vec3fa r_e_p = 1.0f/3.0f * (e_i_m_1 - e_i_p_1) + 2.0f/3.0f * (c_i_m_1 - c_i);
      
      f_p_vtx =  1.0f / d * (c_e_p * p_vtx + (d - 2.0f*c - c_e_p) * e0_p_vtx + 2.0f*c* e1_m_vtx + r_e_p);
      
      const Vec3fa r_e_m = 1.0f/3.0f * (e_i - e_i_m_2) + 2.0f/3.0f * (c_i_m_1 - c_i_m_2);
      
      f_m_vtx = 1.0f / d * (c_e_m * p_vtx + (d - 2.0f*c - c_e_m) * e0_m_vtx + 2.0f*c* e3_p_vtx + r_e_m);      
    }
    
    void init(const CatmullClarkPatch& patch)
    {
      p0() = initCornerVertex(patch,0);
      p1() = initCornerVertex(patch,1);
      p2() = initCornerVertex(patch,2);
      p3() = initCornerVertex(patch,3);
      
      e0_p() = initPositiveEdgeVertex(patch,0, p0());
      e1_p() = initPositiveEdgeVertex(patch,1, p1());
      e2_p() = initPositiveEdgeVertex(patch,2, p2());
      e3_p() = initPositiveEdgeVertex(patch,3, p3());
      
      e0_m() = initNegativeEdgeVertex(patch,0, p0());
      e1_m() = initNegativeEdgeVertex(patch,1, p1());
      e2_m() = initNegativeEdgeVertex(patch,2, p2());
      e3_m() = initNegativeEdgeVertex(patch,3, p3());
      
      const unsigned int face_valence_p0 = patch.ring[0].face_valence;
      const unsigned int face_valence_p1 = patch.ring[1].face_valence;
      const unsigned int face_valence_p2 = patch.ring[2].face_valence;
      const unsigned int face_valence_p3 = patch.ring[3].face_valence;
      
      initFaceVertex(patch,0,p0(),e0_p(),e1_m(),face_valence_p1,e0_m(),e3_p(),face_valence_p3,f0_p(),f0_m() );
      initFaceVertex(patch,1,p1(),e1_p(),e2_m(),face_valence_p2,e1_m(),e0_p(),face_valence_p0,f1_p(),f1_m() );
      initFaceVertex(patch,2,p2(),e2_p(),e3_m(),face_valence_p3,e2_m(),e1_p(),face_valence_p1,f2_p(),f2_m() );
      initFaceVertex(patch,3,p3(),e3_p(),e0_m(),face_valence_p0,e3_m(),e2_p(),face_valence_p3,f3_p(),f3_m() );
    }
    
    void init(const GeneralCatmullClarkPatch& patch)
    {
      assert(patch.size() == 4);
      CatmullClarkPatch qpatch; patch.init(qpatch);
      init(qpatch);
    }
    
    __forceinline void exportConrolPoints( Vec3fa matrix[4][4], Vec3fa f_m[2][2] ) const
    {
      for (size_t y=0;y<4;y++)
	for (size_t x=0;x<4;x++)
	  matrix[y][x] = (Vec3fa_t)v[y][x];
      
      for (size_t y=0;y<2;y++)
	for (size_t x=0;x<2;x++)
	  f_m[y][x] = (Vec3fa_t)f[y][x];
    }
    
    __forceinline void exportDenseConrolPoints( Vec3fa matrix[4][4] ) const //store all f_m into 4th component of Vec3fa matrix
    {
      for (size_t y=0;y<4;y++)
	for (size_t x=0;x<4;x++)
	  matrix[y][x] = (Vec3fa_t)v[y][x];
      
      matrix[0][0].w = f[0][0].x;
      matrix[0][1].w = f[0][0].y;
      matrix[0][2].w = f[0][0].z;
      matrix[0][3].w = 0.0f;
      
      matrix[1][0].w = f[0][1].x;
      matrix[1][1].w = f[0][1].y;
      matrix[1][2].w = f[0][1].z;
      matrix[1][3].w = 0.0f;
      
      matrix[2][0].w = f[1][1].x;
      matrix[2][1].w = f[1][1].y;
      matrix[2][2].w = f[1][1].z;
      matrix[2][3].w = 0.0f;
      
      matrix[3][0].w = f[1][0].x;
      matrix[3][1].w = f[1][0].y;
      matrix[3][2].w = f[1][0].z;
      matrix[3][3].w = 0.0f;
    }
    
    template<class T, class S>
      static __forceinline T deCasteljau_t(const S &uu, const T &v0, const T &v1, const T &v2, const T &v3)
    {
      const S one_minus_uu = 1.0f - uu;      
      const T v0_1 = one_minus_uu * v0   + uu * v1;
      const T v1_1 = one_minus_uu * v1   + uu * v2;
      const T v2_1 = one_minus_uu * v2   + uu * v3;      
      const T v0_2 = one_minus_uu * v0_1 + uu * v1_1;
      const T v1_2 = one_minus_uu * v1_1 + uu * v2_1;      
      const T v0_3 = one_minus_uu * v0_2 + uu * v1_2;
      return v0_3;
    }
    
    template<class T, class S>
      static __forceinline T deCasteljau_tangent_t(const S &uu, const T &v0, const T &v1, const T &v2, const T &v3)
    {
      const S one_minus_uu = 1.0f - uu;      
      const T v0_1         = one_minus_uu * v0   + uu * v1;
      const T v1_1         = one_minus_uu * v1   + uu * v2;
      const T v2_1         = one_minus_uu * v2   + uu * v3;      
      const T v0_2         = one_minus_uu * v0_1 + uu * v1_1;
      const T v1_2         = one_minus_uu * v1_1 + uu * v2_1;      
      return v1_2 - v0_2;      
    }
    
    static __forceinline void computeInnerVertices(const Vec3fa matrix[4][4],
						   const Vec3fa f_m[2][2],
						   const float uu,
						   const float vv,
						   Vec3fa_t &matrix_11,
						   Vec3fa_t &matrix_12,
						   Vec3fa_t &matrix_22,
						   Vec3fa_t &matrix_21)
    {
      if (unlikely(uu == 0.0f || uu == 1.0f || vv == 0.0f || vv == 1.0f)) 
      {
	matrix_11 = matrix[1][1];
	matrix_12 = matrix[1][2];
	matrix_22 = matrix[2][2];
	matrix_21 = matrix[2][1];	 
      }
      else
      {
	const Vec3fa_t f0_p = matrix[1][1];
	const Vec3fa_t f1_p = matrix[1][2];
	const Vec3fa_t f2_p = matrix[2][2];
	const Vec3fa_t f3_p = matrix[2][1];
        
	const Vec3fa_t f0_m = f_m[0][0];
	const Vec3fa_t f1_m = f_m[0][1];
	const Vec3fa_t f2_m = f_m[1][1];
	const Vec3fa_t f3_m = f_m[1][0];
        
	const Vec3fa_t F0 = (      uu  * f0_p +       vv  * f0_m) * 1.0f/(uu+vv);
	const Vec3fa_t F1 = ((1.0f-uu) * f1_m +       vv  * f1_p) * 1.0f/(1.0f-uu+vv);
	const Vec3fa_t F2 = ((1.0f-uu) * f2_p + (1.0f-vv) * f2_m) * 1.0f/(2.0f-uu-vv);
	const Vec3fa_t F3 = (      uu  * f3_m + (1.0f-vv) * f3_p) * 1.0f/(1.0f+uu-vv);
        
	matrix_11 = F0;
	matrix_12 = F1;
	matrix_22 = F2;
	matrix_21 = F3;     
      }
    } 
    
    static __forceinline Vec3fa normal(const Vec3fa matrix[4][4],
				       const Vec3fa f_m[2][2],
				       const float uu,
				       const float vv) 
    {
      
      Vec3fa_t matrix_11, matrix_12, matrix_22, matrix_21;
      computeInnerVertices(matrix,f_m,uu,vv,matrix_11, matrix_12, matrix_22, matrix_21);
      
      /* tangentU */
      const Vec3fa_t col0 = deCasteljau_t(vv, (Vec3fa_t)matrix[0][0], (Vec3fa_t)matrix[1][0], (Vec3fa_t)matrix[2][0], (Vec3fa_t)matrix[3][0]);
      const Vec3fa_t col1 = deCasteljau_t(vv, (Vec3fa_t)matrix[0][1], (Vec3fa_t)matrix_11   , (Vec3fa_t)matrix_21   , (Vec3fa_t)matrix[3][1]);
      const Vec3fa_t col2 = deCasteljau_t(vv, (Vec3fa_t)matrix[0][2], (Vec3fa_t)matrix_12   , (Vec3fa_t)matrix_22   , (Vec3fa_t)matrix[3][2]);
      const Vec3fa_t col3 = deCasteljau_t(vv, (Vec3fa_t)matrix[0][3], (Vec3fa_t)matrix[1][3], (Vec3fa_t)matrix[2][3], (Vec3fa_t)matrix[3][3]);
      
      const Vec3fa_t tangentU = deCasteljau_tangent_t(uu, col0, col1, col2, col3);
      
      /* tangentV */
      const Vec3fa_t row0 = deCasteljau_t(uu, (Vec3fa_t)matrix[0][0], (Vec3fa_t)matrix[0][1], (Vec3fa_t)matrix[0][2], (Vec3fa_t)matrix[0][3]);
      const Vec3fa_t row1 = deCasteljau_t(uu, (Vec3fa_t)matrix[1][0], (Vec3fa_t)matrix_11   , (Vec3fa_t)matrix_12   , (Vec3fa_t)matrix[1][3]);
      const Vec3fa_t row2 = deCasteljau_t(uu, (Vec3fa_t)matrix[2][0], (Vec3fa_t)matrix_21   , (Vec3fa_t)matrix_22   , (Vec3fa_t)matrix[2][3]);
      const Vec3fa_t row3 = deCasteljau_t(uu, (Vec3fa_t)matrix[3][0], (Vec3fa_t)matrix[3][1], (Vec3fa_t)matrix[3][2], (Vec3fa_t)matrix[3][3]);
      
      const Vec3fa_t tangentV = deCasteljau_tangent_t(vv, row0, row1, row2, row3);
      
      /* normal = tangentU x tangentV */
      const Vec3fa_t n = cross(tangentV,tangentU);
      
      return n;     
    }
    
    __forceinline Vec3fa normal( const float uu, const float vv) const
    {
      return normal(v,f,uu,vv);
    }
    
    static __forceinline float extract_f_m(const Vec3fa matrix[4][4],
					   const size_t y,
					   const size_t x)
    {
      return matrix[y][x].w;
    }
    
    static __forceinline Vec3fa extract_f_m_Vec3fa(const Vec3fa matrix[4][4],
						   const size_t n)
    {
      return Vec3fa( extract_f_m(matrix,n,0), extract_f_m(matrix,n,1), extract_f_m(matrix,n,2) );
    }
    
    
    
    static __forceinline Vec3fa eval(const Vec3fa matrix[4][4],
                                     const float &uu,
                                     const float &vv) 
    {
      Vec3fa f[2][2];
      f[0][0] = extract_f_m_Vec3fa( matrix, 0 );
      f[0][1] = extract_f_m_Vec3fa( matrix, 1 );
      f[1][1] = extract_f_m_Vec3fa( matrix, 2 );
      f[1][0] = extract_f_m_Vec3fa( matrix, 3 );
      
      Vec3fa_t v_11, v_12, v_22, v_21;
      computeInnerVertices(matrix,f,uu,vv,v_11, v_12, v_22, v_21);
      
      const float one_minus_uu = 1.0f - uu;
      const float one_minus_vv = 1.0f - vv;      
      
      const float B0_u = one_minus_uu * one_minus_uu * one_minus_uu;
      const float B0_v = one_minus_vv * one_minus_vv * one_minus_vv;
      const float B1_u = 3.0f * one_minus_uu * one_minus_uu * uu;
      const float B1_v = 3.0f * one_minus_vv * one_minus_vv * vv;
      const float B2_u = 3.0f * one_minus_uu * uu * uu;
      const float B2_v = 3.0f * one_minus_vv * vv * vv;
      const float B3_u = uu * uu * uu;
      const float B3_v = vv * vv * vv;
      
      const Vec3fa_t res = 
	(B0_u * matrix[0][0] + B1_u * matrix[0][1] + B2_u * matrix[0][2] + B3_u * matrix[0][3]) * B0_v + 
	(B0_u * matrix[1][0] + B1_u * v_11    + B2_u * v_12    + B3_u * matrix[1][3]) * B1_v + 
	(B0_u * matrix[2][0] + B1_u * v_21    + B2_u * v_22    + B3_u * matrix[2][3]) * B2_v + 
	(B0_u * matrix[3][0] + B1_u * matrix[3][1] + B2_u * matrix[3][2] + B3_u * matrix[3][3]) * B3_v; 
      
      return res;
    }
    
    template<class M, class T>
      static __forceinline Vec3<T> eval_t(const Vec3fa matrix[4][4],
                                          const T &uu,
                                          const T &vv) 
    {
      const M m_border = (uu == 0.0f) | (uu == 1.0f) | (vv == 0.0f) | (vv == 1.0f);
      
      const Vec3<T> f0_p = Vec3<T>(matrix[1][1].x,matrix[1][1].y,matrix[1][1].z);
      const Vec3<T> f1_p = Vec3<T>(matrix[1][2].x,matrix[1][2].y,matrix[1][2].z);
      const Vec3<T> f2_p = Vec3<T>(matrix[2][2].x,matrix[2][2].y,matrix[2][2].z);
      const Vec3<T> f3_p = Vec3<T>(matrix[2][1].x,matrix[2][1].y,matrix[2][1].z);
      
      const Vec3<T> f0_m = Vec3<T>( extract_f_m(matrix,0,0), extract_f_m(matrix,0,1), extract_f_m(matrix,0,2) );
      const Vec3<T> f1_m = Vec3<T>( extract_f_m(matrix,1,0), extract_f_m(matrix,1,1), extract_f_m(matrix,1,2) );
      const Vec3<T> f2_m = Vec3<T>( extract_f_m(matrix,2,0), extract_f_m(matrix,2,1), extract_f_m(matrix,2,2) );
      const Vec3<T> f3_m = Vec3<T>( extract_f_m(matrix,3,0), extract_f_m(matrix,3,1), extract_f_m(matrix,3,2) );
      
      const T one_minus_uu = T(1.0f) - uu;
      const T one_minus_vv = T(1.0f) - vv;      
      
#if 1
      const T inv0 = rcp(uu+vv);
      const T inv1 = rcp(one_minus_uu+vv);
      const T inv2 = rcp(one_minus_uu+one_minus_vv);
      const T inv3 = rcp(uu+one_minus_vv);
#else
      const T inv0 = 1.0f/(uu+vv);
      const T inv1 = 1.0f/(one_minus_uu+vv);
      const T inv2 = 1.0f/(one_minus_uu+one_minus_vv);
      const T inv3 = 1.0f/(uu+one_minus_vv);
#endif
      
      const Vec3<T> f0_i = (          uu * f0_p +           vv * f0_m) * inv0;
      const Vec3<T> f1_i = (one_minus_uu * f1_m +           vv * f1_p) * inv1;
      const Vec3<T> f2_i = (one_minus_uu * f2_p + one_minus_vv * f2_m) * inv2;
      const Vec3<T> f3_i = (          uu * f3_m + one_minus_vv * f3_p) * inv3;
      
      const Vec3<T> F0( select(m_border,f0_p.x,f0_i.x), select(m_border,f0_p.y,f0_i.y), select(m_border,f0_p.z,f0_i.z) );
      const Vec3<T> F1( select(m_border,f1_p.x,f1_i.x), select(m_border,f1_p.y,f1_i.y), select(m_border,f1_p.z,f1_i.z) );
      const Vec3<T> F2( select(m_border,f2_p.x,f2_i.x), select(m_border,f2_p.y,f2_i.y), select(m_border,f2_p.z,f2_i.z) );
      const Vec3<T> F3( select(m_border,f3_p.x,f3_i.x), select(m_border,f3_p.y,f3_i.y), select(m_border,f3_p.z,f3_i.z) );
      
      const T B0_u = one_minus_uu * one_minus_uu * one_minus_uu;
      const T B0_v = one_minus_vv * one_minus_vv * one_minus_vv;
      const T B1_u = 3.0f * one_minus_uu * one_minus_uu * uu;
      const T B1_v = 3.0f * one_minus_vv * one_minus_vv * vv;
      const T B2_u = 3.0f * one_minus_uu * uu * uu;
      const T B2_v = 3.0f * one_minus_vv * vv * vv;
      const T B3_u = uu * uu * uu;
      const T B3_v = vv * vv * vv;
      
      const T x = 
	(B0_u * matrix[0][0].x + B1_u * matrix[0][1].x + B2_u * matrix[0][2].x + B3_u * matrix[0][3].x) * B0_v + 
	(B0_u * matrix[1][0].x + B1_u * F0.x           + B2_u * F1.x           + B3_u * matrix[1][3].x) * B1_v + 
	(B0_u * matrix[2][0].x + B1_u * F3.x           + B2_u * F2.x           + B3_u * matrix[2][3].x) * B2_v + 
	(B0_u * matrix[3][0].x + B1_u * matrix[3][1].x + B2_u * matrix[3][2].x + B3_u * matrix[3][3].x) * B3_v; 
      
      const T y = 
	(B0_u * matrix[0][0].y + B1_u * matrix[0][1].y + B2_u * matrix[0][2].y + B3_u * matrix[0][3].y) * B0_v + 
	(B0_u * matrix[1][0].y + B1_u * F0.y           + B2_u * F1.y           + B3_u * matrix[1][3].y) * B1_v + 
	(B0_u * matrix[2][0].y + B1_u * F3.y           + B2_u * F2.y           + B3_u * matrix[2][3].y) * B2_v + 
	(B0_u * matrix[3][0].y + B1_u * matrix[3][1].y + B2_u * matrix[3][2].y + B3_u * matrix[3][3].y) * B3_v; 
      
      const T z = 
	(B0_u * matrix[0][0].z + B1_u * matrix[0][1].z + B2_u * matrix[0][2].z + B3_u * matrix[0][3].z) * B0_v + 
	(B0_u * matrix[1][0].z + B1_u * F0.z           + B2_u * F1.z           + B3_u * matrix[1][3].z) * B1_v + 
	(B0_u * matrix[2][0].z + B1_u * F3.z           + B2_u * F2.z           + B3_u * matrix[2][3].z) * B2_v + 
	(B0_u * matrix[3][0].z + B1_u * matrix[3][1].z + B2_u * matrix[3][2].z + B3_u * matrix[3][3].z) * B3_v; 
      
      
      return Vec3<T>(x,y,z);
    }
    
    
    template<class M, class T>
      static __forceinline Vec3<T> normal_t(const Vec3fa matrix[4][4],
                                            const T &uu,
                                            const T &vv) 
    {
      const M m_border = (uu == 0.0f) | (uu == 1.0f) | (vv == 0.0f) | (vv == 1.0f);
      
      const Vec3<T> f0_p = Vec3<T>(matrix[1][1].x,matrix[1][1].y,matrix[1][1].z);
      const Vec3<T> f1_p = Vec3<T>(matrix[1][2].x,matrix[1][2].y,matrix[1][2].z);
      const Vec3<T> f2_p = Vec3<T>(matrix[2][2].x,matrix[2][2].y,matrix[2][2].z);
      const Vec3<T> f3_p = Vec3<T>(matrix[2][1].x,matrix[2][1].y,matrix[2][1].z);
      
      const Vec3<T> f0_m = Vec3<T>( extract_f_m(matrix,0,0), extract_f_m(matrix,0,1), extract_f_m(matrix,0,2) );
      const Vec3<T> f1_m = Vec3<T>( extract_f_m(matrix,1,0), extract_f_m(matrix,1,1), extract_f_m(matrix,1,2) );
      const Vec3<T> f2_m = Vec3<T>( extract_f_m(matrix,2,0), extract_f_m(matrix,2,1), extract_f_m(matrix,2,2) );
      const Vec3<T> f3_m = Vec3<T>( extract_f_m(matrix,3,0), extract_f_m(matrix,3,1), extract_f_m(matrix,3,2) );
      
      const T one_minus_uu = T(1.0f) - uu;
      const T one_minus_vv = T(1.0f) - vv;      
      
      const T inv0 = rcp(uu+vv);
      const T inv1 = rcp(one_minus_uu+vv);
      const T inv2 = rcp(one_minus_uu+one_minus_vv);
      const T inv3 = rcp(uu+one_minus_vv);
      
      const Vec3<T> f0_i = (          uu * f0_p +           vv * f0_m) * inv0;
      const Vec3<T> f1_i = (one_minus_uu * f1_m +           vv * f1_p) * inv1;
      const Vec3<T> f2_i = (one_minus_uu * f2_p + one_minus_vv * f2_m) * inv2;
      const Vec3<T> f3_i = (          uu * f3_m + one_minus_vv * f3_p) * inv3;

#if 1
      const M m_border0 = (uu == 0.0f) & (vv == 0.0f);
      const M m_border1 = (uu == 1.0f) & (vv == 0.0f);
      const M m_border2 = (uu == 1.0f) & (vv == 1.0f);
      const M m_border3 = (uu == 0.0f) & (vv == 1.0f);
      
      const Vec3<T> matrix_11( select(m_border0,f0_p.x,f0_i.x), select(m_border0,f0_p.y,f0_i.y), select(m_border0,f0_p.z,f0_i.z) );
      const Vec3<T> matrix_12( select(m_border1,f1_p.x,f1_i.x), select(m_border1,f1_p.y,f1_i.y), select(m_border1,f1_p.z,f1_i.z) );
      const Vec3<T> matrix_22( select(m_border2,f2_p.x,f2_i.x), select(m_border2,f2_p.y,f2_i.y), select(m_border2,f2_p.z,f2_i.z) );
      const Vec3<T> matrix_21( select(m_border3,f3_p.x,f3_i.x), select(m_border3,f3_p.y,f3_i.y), select(m_border3,f3_p.z,f3_i.z) );
#else
      const Vec3<T> matrix_11( select(m_border,f0_p.x,f0_i.x), select(m_border,f0_p.y,f0_i.y), select(m_border,f0_p.z,f0_i.z) );
      const Vec3<T> matrix_12( select(m_border,f1_p.x,f1_i.x), select(m_border,f1_p.y,f1_i.y), select(m_border,f1_p.z,f1_i.z) );
      const Vec3<T> matrix_22( select(m_border,f2_p.x,f2_i.x), select(m_border,f2_p.y,f2_i.y), select(m_border,f2_p.z,f2_i.z) );
      const Vec3<T> matrix_21( select(m_border,f3_p.x,f3_i.x), select(m_border,f3_p.y,f3_i.y), select(m_border,f3_p.z,f3_i.z) );
#endif
      
      const Vec3<T> matrix_00 = Vec3<T>(matrix[0][0].x,matrix[0][0].y,matrix[0][0].z);
      const Vec3<T> matrix_10 = Vec3<T>(matrix[1][0].x,matrix[1][0].y,matrix[1][0].z);
      const Vec3<T> matrix_20 = Vec3<T>(matrix[2][0].x,matrix[2][0].y,matrix[2][0].z);
      const Vec3<T> matrix_30 = Vec3<T>(matrix[3][0].x,matrix[3][0].y,matrix[3][0].z);
      
      const Vec3<T> matrix_01 = Vec3<T>(matrix[0][1].x,matrix[0][1].y,matrix[0][1].z);
      const Vec3<T> matrix_02 = Vec3<T>(matrix[0][2].x,matrix[0][2].y,matrix[0][2].z);
      const Vec3<T> matrix_03 = Vec3<T>(matrix[0][3].x,matrix[0][3].y,matrix[0][3].z);
      
      const Vec3<T> matrix_31 = Vec3<T>(matrix[3][1].x,matrix[3][1].y,matrix[3][1].z);
      const Vec3<T> matrix_32 = Vec3<T>(matrix[3][2].x,matrix[3][2].y,matrix[3][2].z);
      const Vec3<T> matrix_33 = Vec3<T>(matrix[3][3].x,matrix[3][3].y,matrix[3][3].z);
      
      const Vec3<T> matrix_13 = Vec3<T>(matrix[1][3].x,matrix[1][3].y,matrix[1][3].z);
      const Vec3<T> matrix_23 = Vec3<T>(matrix[2][3].x,matrix[2][3].y,matrix[2][3].z);
      
      /* tangentU */
      const Vec3<T> col0 = deCasteljau_t(vv, matrix_00, matrix_10, matrix_20, matrix_30);
      const Vec3<T> col1 = deCasteljau_t(vv, matrix_01, matrix_11, matrix_21, matrix_31);
      const Vec3<T> col2 = deCasteljau_t(vv, matrix_02, matrix_12, matrix_22, matrix_32);
      const Vec3<T> col3 = deCasteljau_t(vv, matrix_03, matrix_13, matrix_23, matrix_33);
      
      const Vec3<T> tangentU = deCasteljau_tangent_t(uu, col0, col1, col2, col3);
      
      /* tangentV */
      const Vec3<T> row0 = deCasteljau_t(uu, matrix_00, matrix_01, matrix_02, matrix_03);
      const Vec3<T> row1 = deCasteljau_t(uu, matrix_10, matrix_11, matrix_12, matrix_13);
      const Vec3<T> row2 = deCasteljau_t(uu, matrix_20, matrix_21, matrix_22, matrix_23);
      const Vec3<T> row3 = deCasteljau_t(uu, matrix_30, matrix_31, matrix_32, matrix_33);
      
      const Vec3<T> tangentV = deCasteljau_tangent_t(vv, row0, row1, row2, row3);
      
      /* normal = tangentU x tangentV */
      const Vec3<T> n = cross(tangentV,tangentU);
      return n;
    }
    
    
#if !defined(__MIC__)
    
#if defined(__AVX__)    
    
    static __forceinline avx3f eval8  (const Vec3fa matrix[4][4], const avxf &uu, const avxf &vv) { return eval_t<avxb,avxf>(matrix,uu,vv); }
    static __forceinline avx3f normal8(const Vec3fa matrix[4][4], const avxf &uu, const avxf &vv) { return normal_t<avxb,avxf>(matrix,uu,vv); }
    
#endif
    
    static __forceinline sse3f eval4  (const Vec3fa matrix[4][4], const ssef &uu, const ssef &vv) { return eval_t<sseb,ssef>(matrix,uu,vv); }
    static __forceinline sse3f normal4(const Vec3fa matrix[4][4], const ssef &uu, const ssef &vv) { return normal_t<sseb,ssef>(matrix,uu,vv); }
    
#else    
    
    
    static __forceinline mic_f extract_f_m_mic_f(const Vec3fa matrix[4][4],
						 const size_t n)
    {
      const mic_f row = load16f(&matrix[n][0]);
      __aligned(64) float xyzw[16];
      compactustore16f_low(0x8888,xyzw,row);
      return broadcast4to16f(xyzw);
    }
    
    static __forceinline mic3f extract_f_m_mic3f(const Vec3fa matrix[4][4],
						 const size_t n)
    {
      return mic3f( extract_f_m(matrix,n,0), extract_f_m(matrix,n,1), extract_f_m(matrix,n,2) );
    }
    
    static __forceinline mic_f eval4(const Vec3fa matrix[4][4],
				     const mic_f uu,
				     const mic_f vv) 
    {
      const mic_m m_border = (uu == 0.0f) | (uu == 1.0f) | (vv == 0.0f) | (vv == 1.0f);
      
      const mic_f f0_p = (Vec3fa_t)matrix[1][1];
      const mic_f f1_p = (Vec3fa_t)matrix[1][2];
      const mic_f f2_p = (Vec3fa_t)matrix[2][2];
      const mic_f f3_p = (Vec3fa_t)matrix[2][1];
      
      const mic_f f0_m = extract_f_m_mic_f(matrix,0);
      const mic_f f1_m = extract_f_m_mic_f(matrix,1);
      const mic_f f2_m = extract_f_m_mic_f(matrix,2);
      const mic_f f3_m = extract_f_m_mic_f(matrix,3);
      
      const mic_f one_minus_uu = mic_f(1.0f) - uu;
      const mic_f one_minus_vv = mic_f(1.0f) - vv;      
      
#if 1
      const mic_f inv0 = rcp(uu+vv);
      const mic_f inv1 = rcp(one_minus_uu+vv);
      const mic_f inv2 = rcp(one_minus_uu+one_minus_vv);
      const mic_f inv3 = rcp(uu+one_minus_vv);
#else
      const mic_f inv0 = 1.0f/(uu+vv);
      const mic_f inv1 = 1.0f/(one_minus_uu+vv);
      const mic_f inv2 = 1.0f/(one_minus_uu+one_minus_vv);
      const mic_f inv3 = 1.0f/(uu+one_minus_vv);
#endif
      
      const mic_f F0 = select(m_border,f0_p, (          uu * f0_p +           vv * f0_m) * inv0);
      const mic_f F1 = select(m_border,f1_p, (one_minus_uu * f1_m +           vv * f1_p) * inv1);
      const mic_f F2 = select(m_border,f2_p, (one_minus_uu * f2_p + one_minus_vv * f2_m) * inv2);
      const mic_f F3 = select(m_border,f3_p, (          uu * f3_m + one_minus_vv * f3_p) * inv3);
      
      const mic_f B0_u = one_minus_uu * one_minus_uu * one_minus_uu;
      const mic_f B0_v = one_minus_vv * one_minus_vv * one_minus_vv;
      const mic_f B1_u = 3.0f * one_minus_uu * one_minus_uu * uu;
      const mic_f B1_v = 3.0f * one_minus_vv * one_minus_vv * vv;
      const mic_f B2_u = 3.0f * one_minus_uu * uu * uu;
      const mic_f B2_v = 3.0f * one_minus_vv * vv * vv;
      const mic_f B3_u = uu * uu * uu;
      const mic_f B3_v = vv * vv * vv;
      
      const mic_f res = 
	(B0_u * (Vec3fa_t)matrix[0][0] + B1_u * (Vec3fa_t)matrix[0][1] + B2_u * (Vec3fa_t)matrix[0][2] + B3_u * (Vec3fa_t)matrix[0][3]) * B0_v + 
	(B0_u * (Vec3fa_t)matrix[1][0] + B1_u *                     F0 + B2_u *                     F1 + B3_u * (Vec3fa_t)matrix[1][3]) * B1_v + 
	(B0_u * (Vec3fa_t)matrix[2][0] + B1_u *                     F3 + B2_u *                     F2 + B3_u * (Vec3fa_t)matrix[2][3]) * B2_v + 
	(B0_u * (Vec3fa_t)matrix[3][0] + B1_u * (Vec3fa_t)matrix[3][1] + B2_u * (Vec3fa_t)matrix[3][2] + B3_u * (Vec3fa_t)matrix[3][3]) * B3_v; 
      return res;
    }
    
    
    static __forceinline mic3f eval16(const Vec3fa matrix[4][4],
				      const mic_f &uu,
				      const mic_f &vv) 
    {
      return eval_t<mic_m,mic_f>(matrix,uu,vv);
    }
    
    
    
    static __forceinline mic3f normal16(const Vec3fa matrix[4][4],
					const mic_f &uu,
					const mic_f &vv) 
    {
      return normal_t<mic_m,mic_f>(matrix,uu,vv);
    }
    
#endif
    
    
    static __forceinline Vec3fa normal(const Vec3fa matrix[4][4],
				       const float uu,
				       const float vv) 
    {
#if defined(__MIC__)
      const mic_f row0 = load16f(&matrix[0][0]);
      const mic_f row1 = load16f(&matrix[1][0]);
      const mic_f row2 = load16f(&matrix[2][0]);
      const mic_f row3 = load16f(&matrix[3][0]);
      
      __aligned(64) Vec3fa f_m[2][2];
      compactustore16f_low(0x8888,(float*)&f_m[0][0],row0);
      compactustore16f_low(0x8888,(float*)&f_m[0][1],row1);
      compactustore16f_low(0x8888,(float*)&f_m[1][1],row2);
      compactustore16f_low(0x8888,(float*)&f_m[1][0],row3);
#else
      __aligned(64) Vec3fa f_m[2][2];
      f_m[0][0] = extract_f_m_Vec3fa(matrix,0);
      f_m[0][1] = extract_f_m_Vec3fa(matrix,1);
      f_m[1][1] = extract_f_m_Vec3fa(matrix,2);
      f_m[1][0] = extract_f_m_Vec3fa(matrix,3);      
#endif      
      return normal(matrix,f_m,uu,vv);
    }
    
    __forceinline Vec3fa eval(const float uu, const float vv) const
    {
      Vec3fa_t v_11, v_12, v_22, v_21;
      computeInnerVertices(v,f,uu,vv,v_11, v_12, v_22, v_21);
      
      const float one_minus_uu = 1.0f - uu;
      const float one_minus_vv = 1.0f - vv;      
      
      const float B0_u = one_minus_uu * one_minus_uu * one_minus_uu;
      const float B0_v = one_minus_vv * one_minus_vv * one_minus_vv;
      const float B1_u = 3.0f * one_minus_uu * one_minus_uu * uu;
      const float B1_v = 3.0f * one_minus_vv * one_minus_vv * vv;
      const float B2_u = 3.0f * one_minus_uu * uu * uu;
      const float B2_v = 3.0f * one_minus_vv * vv * vv;
      const float B3_u = uu * uu * uu;
      const float B3_v = vv * vv * vv;
      
      const Vec3fa_t res = 
	(B0_u * v[0][0] + B1_u * v[0][1] + B2_u * v[0][2] + B3_u * v[0][3]) * B0_v + 
	(B0_u * v[1][0] + B1_u * v_11    + B2_u * v_12    + B3_u * v[1][3]) * B1_v + 
	(B0_u * v[2][0] + B1_u * v_21    + B2_u * v_22    + B3_u * v[2][3]) * B2_v + 
	(B0_u * v[3][0] + B1_u * v[3][1] + B2_u * v[3][2] + B3_u * v[3][3]) * B3_v; 
      
      return res;
      
    }
    
    __forceinline BBox3fa bounds() const
    {
      const Vec3fa *const cv = &v[0][0];
      BBox3fa bounds ( cv[0] );
      for (size_t i = 1; i<16 ; i++)
	bounds.extend( cv[i] );
      bounds.extend(f[0][0]);
      bounds.extend(f[1][0]);
      bounds.extend(f[1][1]);
      bounds.extend(f[1][1]);
      return bounds;
    }
    
    friend std::ostream &operator<<(std::ostream &o, const GregoryPatch &p)
    {
      for (size_t y=0;y<4;y++)
	for (size_t x=0;x<4;x++)
	  o << "v[" << y << "][" << x << "] " << p.v[y][x] << std::endl;
      
      for (size_t y=0;y<2;y++)
	for (size_t x=0;x<2;x++)
	  o << "f[" << y << "][" << x << "] " << p.f[y][x] << std::endl;
      return o;
    } 
  };
}
