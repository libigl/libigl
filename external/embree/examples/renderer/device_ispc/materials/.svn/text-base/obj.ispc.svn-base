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

#include "material.isph"
#include "textures/texture.isph"
#include "brdfs/transmission.isph"
#include "brdfs/lambertian.isph"
#include "brdfs/specular.isph"

struct Obj
{
  Material base;
  Texture* map_d;  float d; 
  Texture* map_Kd; vec3f Kd;
  Texture* map_Ks; vec3f Ks;
  Texture* map_Ns; float Ns;
  Texture* map_Bump;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

void Obj__shade(const uniform Material *uniform _this,
                const Ray&                  ray,           
                const Medium&               currentMedium, 
                const DifferentialGeometry& dg,            
                uniform CompositedBRDF&     brdfs)         
{
  uniform const Obj* uniform this = (uniform const Obj* uniform)_this;


  /*! transmission */
  float d = this->d;  
  if (this->map_d) { const vec3f c = this->map_d->get(this->map_d,dg.st); d *= c.x; }
  if (d < 1.0f) COMPOSITED_BRDF_ADD(brdfs,Transmission,make_vec3f(1.0f-d));

  /*! diffuse component */
  vec3f Kd = mul(d, this->Kd);  
  if (this->map_Kd) Kd = mul(Kd, this->map_Kd->get(this->map_Kd,dg.st));  
  if (ne(Kd,make_vec3f(0.0f))) COMPOSITED_BRDF_ADD(brdfs,Lambertian,Kd);

  /*! specular exponent */
  float Ns = this->Ns;  
  if (this->map_Ns) { const vec3f c = this->map_Ns->get(this->map_Ns,dg.st); Ns *= c.x; }
  
  /*! specular component */
  vec3f Ks = mul(d, this->Ks);  
  if (this->map_Ks) Ks = mul(Ks, this->map_Ks->get(this->map_Ks,dg.st));  
  if (ne(Ks,make_vec3f(0.0f))) COMPOSITED_BRDF_ADD(brdfs,Specular,Ks,Ns);
}

void Obj__Destructor(uniform RefCount* uniform _this)
{
  uniform Obj* uniform this = (uniform Obj* uniform) _this;
  if (this->map_d ) RefCount__DecRef(&this->map_d->base);
  if (this->map_Kd) RefCount__DecRef(&this->map_Kd->base);
  if (this->map_Ks) RefCount__DecRef(&this->map_Ks->base);
  if (this->map_Ns) RefCount__DecRef(&this->map_Ns->base);
  if (this->map_Bump) RefCount__DecRef(&this->map_Bump->base);
  Material__Destructor(_this);
}

void Obj__Constructor(uniform Obj* uniform this,
                      uniform Texture* uniform map_d,  const uniform float d,
                      uniform Texture* uniform map_Kd, const uniform vec3f Kd,
                      uniform Texture* uniform map_Ks, const uniform vec3f Ks,
                      uniform Texture* uniform map_Ns, const uniform float Ns,
                      uniform Texture* uniform map_Bump)
{
  Material__Constructor(&this->base,Obj__Destructor,Obj__shade,NULL);
  this->map_d  = map_d;  this->d  = d;  if (this->map_d ) RefCount__IncRef(&map_d->base);
  this->map_Kd = map_Kd; this->Kd = Kd; if (this->map_Kd) RefCount__IncRef(&map_Kd->base);
  this->map_Ks = map_Ks; this->Ks = Ks; if (this->map_Ks) RefCount__IncRef(&map_Ks->base);
  this->map_Ns = map_Ns; this->Ns = Ns; if (this->map_Ns) RefCount__IncRef(&map_Ns->base);
  this->map_Bump = map_Bump; if (this->map_Bump) RefCount__IncRef(&map_Bump->base);
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform Obj__new(void* uniform map_d,  const uniform float& d,
                              void* uniform map_Kd, const uniform vec3f& Kd,
                              void* uniform map_Ks, const uniform vec3f& Ks,
                              void* uniform map_Ns, const uniform float& Ns,
                              void* uniform map_Bump)
{
  uniform Obj* uniform this = uniform new uniform Obj;
  Obj__Constructor(this,
                   (uniform Texture* uniform)map_d,d,
                   (uniform Texture* uniform)map_Kd,Kd,
                   (uniform Texture* uniform)map_Ks,Ks,
                   (uniform Texture* uniform)map_Ns,Ns,
                   (uniform Texture* uniform)map_Bump);
  return this;
}
