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

#include "scene.isph"
#include "../../../../embree_ispc/include/intersector.isph"

void Scene__Destructor(uniform RefCount* uniform _this)
{
  LOG(print("Scene__Destructor\n"));
  uniform Scene* uniform this = (uniform Scene* uniform) _this;

  for (uniform int i=0; i<this->num_allLights; i++)
    RefCount__DecRef(&this->allLights[i]->base);

  for (uniform int i=0; i<this->num_envLights; i++)
    RefCount__DecRef(&this->envLights[i]->base.base);

  for (uniform int i=0; i<this->num_geometries; i++)
    RefCount__DecRef(&this->geometry[i]->base);

  rtcDeleteIntersector(this->intersector);
  rtcDeleteGeometry(this->geom);

  RefCount__Destructor(_this);
}

void Scene__Constructor(uniform Scene* uniform this,
                        void* uniform geom,
                        void* uniform traverserTy,
                        const uniform int num_allLights,
                        void** uniform allLights,
                        const uniform int num_envLights,
                        void** uniform envLights,
                        const uniform int num_geometries,
                        void** uniform geometry)
{
  LOG(print("Scene__Constructor\n"));
  RefCount__Constructor(&this->base,Scene__Destructor);

  /*! add list of lights */
  this->num_allLights = num_allLights;
  this->allLights = (uniform Light**) allLights;

  /*! add list of environment lights */
  this->num_envLights = num_envLights;
  this->envLights = (uniform EnvironmentLight**) envLights;

  /*! add geometries */
  this->num_geometries = num_geometries;
  this->geometry = (uniform Instance**) geometry;

  /*! add acceleration structure */
  this->geom = (uniform RTCGeometry* uniform) geom;
  this->intersector = rtcQueryIntersector(this->geom,traverserTy);
}

export void* uniform Scene__new(void* uniform accel,
                                void* uniform traverserTy,
                                const uniform int num_allLights,
                                void** uniform allLights,
                                const uniform int num_envLights,
                                void** uniform envLights,
                                const uniform int num_geometries,
                                void** uniform geometry)
{
  uniform Scene *uniform this = uniform new uniform Scene;
  Scene__Constructor(this,accel,traverserTy,num_allLights,allLights,num_envLights,envLights,num_geometries,geometry);
  return this;
}

