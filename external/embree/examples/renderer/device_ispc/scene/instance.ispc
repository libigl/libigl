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

#include "instance.isph"

void Instance__Destructor(uniform RefCount* uniform _this) 
{
  LOG(print("Instance__Destructor\n"));
  uniform Instance* uniform this = (uniform Instance* uniform) _this;
  RefCount__DecRef(&this->shape   ->base);
  RefCount__DecRef(&this->material->base);
  RefCount__DecRef(&this->light->base.base);
  RefCount__Destructor(_this);
}

void Instance__Constructor(uniform Instance* uniform this,
                           uniform Shape* uniform shape,
                           uniform Material* uniform material,
                           uniform AreaLight* uniform light)
{
  LOG(print("Instance__Constructor\n"));
  RefCount__Constructor(&this->base,Instance__Destructor);
  RefCount__IncRef(&shape   ->base); this->shape    = shape;
  RefCount__IncRef(&material->base); this->material = material;
  RefCount__IncRef(&light->base.base); this->light  = light;
}

export void* uniform Instance__new(void* uniform shape,
                                   void* uniform material,
                                   void* uniform light)
{
  uniform Instance *uniform this = uniform new uniform Instance;
  Instance__Constructor(this,
                        (uniform Shape*     uniform) shape, 
                        (uniform Material*  uniform) material, 
                        (uniform AreaLight* uniform) light);
  return this;
}
