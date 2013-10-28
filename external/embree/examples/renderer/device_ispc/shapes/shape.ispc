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

#include "shape.isph"

void Shape__Destructor(uniform RefCount* uniform this)
{
  LOG(print("Shape__Destructor\n"));
  RefCount__Destructor(this);
}

void Shape__Constructor(uniform Shape* uniform this,
                        uniform DestructorFunc destructor,
                        uniform GetNumVerticesFunc getNumVertices,
                        uniform GetNumTrianglesFunc getNumTriangles,
                        uniform PostIntersectFunc postIntersect,
                        uniform ShapeTransformFunc transform,
                        uniform ExtractFunc extract_)
{
  LOG(print("Shape__Constructor\n"));
  RefCount__Constructor(&this->base,destructor);
  this->getNumVertices = getNumVertices;
  this->getNumTriangles = getNumTriangles;
  this->postIntersect = postIntersect;
  this->transform = transform;
  this->extract = extract_;
}

export uniform int Shape__getNumVertices(void *uniform _this) {
  const uniform Shape* uniform this = (const uniform Shape* uniform) _this;
  return this->getNumVertices(this);
}

export uniform int Shape__getNumTriangles(void *uniform _this) {
  const uniform Shape* uniform this = (const uniform Shape* uniform) _this;
  return this->getNumTriangles(this);
}

export void* uniform Shape__transform(void* uniform _this,
                                      const uniform vec3f& vx,
                                      const uniform vec3f& vy,
                                      const uniform vec3f& vz,
                                      const uniform vec3f& p) 
{
  const uniform Shape* uniform this = (const uniform Shape* uniform) _this;
  uniform AffineSpace3f xfm = make_AffineSpace3f(vx,vy,vz,p);
  return this->transform(this,xfm);
}

export void Shape__extract(void* uniform _this,
                           const uniform int& id0,
                           uniform RTCTriangle* uniform triangles,
                           uniform int& numTriangles,
                           uniform RTCVertex* uniform vertices,
                           uniform int& numVertices,
                           uniform vec3f& lower_o,
                           uniform vec3f& upper_o)
{
  const uniform Shape* uniform this = (const uniform Shape* uniform) _this;
  this->extract(this,id0,triangles,numTriangles,vertices,numVertices,lower_o,upper_o);
}
