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

#include "light.isph"
#include "samplers/shapesampler.isph"
#include "shapes/trianglemesh.isph"

struct TriangleLight
{
  AreaLight base;
  vec3f v0;          //!< First vertex of the triangle
  vec3f v1;          //!< Second vertex of the triangle
  vec3f v2;          //!< Third vertex of the triangle
  vec3f e1;          //!< First edge of the tringle (v0-v1)
  vec3f e2;          //!< Second edge of the triangle (v2-v0)
  vec3f Ng;          //!< Normal of triangle light
  vec3f L;           //!< Radiance (W/(m^2*sr))
  Shape* shape;      //!< Triangle shape
};

/*! Intersects the triangle with a ray of origin O and direction D. */
inline void intersect(const uniform TriangleLight* uniform this, const vec3f O, const vec3f D, float& t, float& u, float& v, float& w)
{
  const vec3f C = sub(this->v0,O);
  const vec3f R = cross(D,C);
  const float det = dot(this->Ng,D);
  const float rcp_det = rcp(det);
  t = dot(C,this->Ng) * rcp_det;
  u = dot(R,this->e2) * rcp_det;
  v = dot(R,this->e1) * rcp_det;
  w = 1.0f-u-v;
}

TriangleLight* uniform TriangleLight__New(const uniform vec3f& v0,
                                          const uniform vec3f& v1,
                                          const uniform vec3f& v2,
                                          const uniform vec3f& L);

uniform Light* uniform TriangleLight__transform(const uniform Light *uniform _this, 
                                               const uniform AffineSpace3f& xfm) 
{
  const uniform TriangleLight* uniform this = (const uniform TriangleLight* uniform)_this;
  return (uniform Light* uniform) TriangleLight__New(xfmPoint(xfm,this->v0),
                                                     xfmPoint(xfm,this->v1),
                                                     xfmPoint(xfm,this->v2),
                                                     this->L);
}

Shape* uniform TriangleLight__shape(const uniform Light *uniform _this) 
{
  const uniform TriangleLight* uniform this = (const uniform TriangleLight* uniform)_this;
  return this->shape;
}

varying vec3f TriangleLight__Le(const uniform AreaLight* uniform _this, 
                                   varying const DifferentialGeometry &dg, 
                                   varying const vec3f& wo) 
{
  const uniform TriangleLight* uniform this = (const uniform TriangleLight* uniform)_this;
  return this->L;
}

varying vec3f TriangleLight__eval(const uniform Light* uniform _this, 
                                     varying const DifferentialGeometry& dg, 
                                     varying const vec3f& wi) 
{
  const uniform TriangleLight* uniform this = (const uniform TriangleLight* uniform)_this;
  float t,u,v,w; intersect(this,dg.P,wi,t,u,v,w);
  if (0.0f <= t & min(min(u,v),w) >= 0) return this->L;
  return make_vec3f(0.0f);
}

varying vec3f TriangleLight__sample(const uniform Light* uniform _this,
                                    varying const DifferentialGeometry& dg, 
                                    varying Sample3f& wi,
                                    varying float& tMax,
                                    varying const vec2f& s) 
{
  const uniform TriangleLight* uniform this = (const uniform TriangleLight* uniform)_this;
  const vec3f d = sub(uniformSampleTriangle(s.x,s.y,this->v0,this->v1,this->v2),dg.P);
  tMax = length(d);
  float dDotNg = dot(d,this->Ng);
  if (dDotNg >= 0) return make_vec3f(0.0f);
  wi = make_Sample3f(mul(d,rcp(tMax)),2.0f*tMax*tMax*tMax*rcp(abs(dDotNg)));
  return this->L;
}

void TriangleLight__Destructor(uniform RefCount* uniform _this)
{ 
  uniform TriangleLight* uniform this = (uniform TriangleLight* uniform) _this;
  RefCount__DecRef(&this->shape->base);
  Light__Destructor(_this);
}

void TriangleLight__Constructor(uniform TriangleLight* uniform this, 
                                const uniform vec3f v0,
                                const uniform vec3f v1,
                                const uniform vec3f v2,
                                const uniform vec3f L) 
{
  AreaLight__Constructor(&this->base,TriangleLight__Destructor,AREA_LIGHT,
                         TriangleLight__transform,TriangleLight__shape,TriangleLight__eval,TriangleLight__sample,
                         TriangleLight__Le);
  this->v0 = v0;
  this->v1 = v1;
  this->v2 = v2;
  this->e1 = sub(v0,v1);
  this->e2 = sub(v2,v0);
  this->Ng = cross(this->e1,this->e2);
  this->L = L;

  /* create triangle shape */
  uniform vec3fa* uniform positions = uniform new uniform vec3fa[3];
  positions[0] = make_vec3fa(v0.x,v0.y,v0.z);
  positions[1] = make_vec3fa(v1.x,v1.y,v1.z);
  positions[2] = make_vec3fa(v2.x,v2.y,v2.z);
  uniform vec4i* uniform triangles = uniform new uniform vec4i[1];
  triangles[0] = make_vec4i(0,1,2,0);
  this->shape = (uniform Shape* uniform) TriangleMesh__new(positions,NULL,NULL,NULL,3,triangles,1);
  RefCount__IncRef(&this->shape->base);
}

TriangleLight* uniform TriangleLight__New(const uniform vec3f& v0,
                                          const uniform vec3f& v1,
                                          const uniform vec3f& v2,
                                          const uniform vec3f& L)
{
  uniform TriangleLight* uniform light = uniform new uniform TriangleLight;
  TriangleLight__Constructor(light,v0,v1,v2,L);
  return light;
}

export void* uniform TriangleLight__new(const uniform vec3f& v0,
                                        const uniform vec3f& v1,
                                        const uniform vec3f& v2,
                                        const uniform vec3f& L)
{
  return TriangleLight__New(v0,v1,v2,L);
}
