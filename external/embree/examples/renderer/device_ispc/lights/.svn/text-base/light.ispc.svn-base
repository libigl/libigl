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

void Light__Destructor(uniform RefCount* uniform this) {
  LOG(print("Light__Destructor\n"));
  RefCount__Destructor(this);
}

void Light__Constructor(uniform Light* uniform this,
                        uniform DestructorFunc destructor,
                        uniform LightType type,
                        uniform LightTransformFunc transform,
                        uniform ShapeFunc shape,
                        uniform EvalFunc eval,
                        uniform SampleFunc sample)
{
  LOG(print("Light__Constructor\n"));
  RefCount__Constructor(&this->base,destructor);
  this->type = type;
  this->transform = transform;
  this->shape = shape;
  this->eval = eval;
  this->sample = sample;
}

void AreaLight__Constructor(uniform AreaLight* uniform this,
                            uniform DestructorFunc destructor,
                            uniform LightType type,
                            uniform LightTransformFunc transform,
                            uniform ShapeFunc shape,
                            uniform EvalFunc eval,
                            uniform SampleFunc sample_,
                            uniform AreaLeFunc Le)
{
  Light__Constructor(&this->base,destructor,type,transform,shape,eval,sample_);
  this->Le = Le;
}

void EnvironmentLight__Constructor(uniform EnvironmentLight* uniform this,
                                   uniform DestructorFunc destructor,
                                   uniform LightType type,
                                   uniform LightTransformFunc transform,
                                   uniform ShapeFunc shape,
                                   uniform EvalFunc eval,
                                   uniform SampleFunc sample_,
                                   uniform EnvironmentLeFunc Le)
{
  Light__Constructor(&this->base,destructor,type,transform,shape,eval,sample_);
  this->Le = Le;
}

export uniform LightType Light__getType(void* uniform _this)
{
  uniform Light *uniform this = (uniform Light *uniform) _this;
  return this->type;
}

export void* uniform Light__transform(void* uniform _this,
                                      uniform const vec3f& vx, 
                                      uniform const vec3f& vy, 
                                      uniform const vec3f& vz, 
                                      uniform const vec3f& p)
{
  const uniform Light *uniform this = (const uniform Light *uniform) _this;
  const uniform AffineSpace3f xfm = make_AffineSpace3f(vx,vy,vz,p);
  return this->transform(this,xfm);
}

export void* uniform Light__shape(void *uniform _this)
{
  const uniform Light *uniform this = (const uniform Light *uniform) _this;
  if (!this->shape) return NULL;
  return this->shape(this);
}
