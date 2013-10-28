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

#include "trianglemesh.isph"

struct TriangleMesh
{
  uniform Shape base;
  uniform vec3fa* uniform position;   //!< Position array.
  uniform vec3fa* uniform motion;     //!< Motion array.
  uniform vec3fa* uniform normal;     //!< Normal array (can be empty).
  uniform vec2f* uniform texcoord;   //!< Texture coordinates array (can be empty).
  uniform vec4i*  uniform triangle;   //!< Triangle indices array.

  uniform int numPositions;
  uniform int numTriangles;
};

///////////////////////////////////////////////////////////////////////////////
// Implementation

uniform int TriangleMesh__getNumVertices(uniform const Shape *uniform _this) {
  const uniform TriangleMesh* uniform this = (const uniform TriangleMesh* uniform) _this;
  return this->numPositions;
}

uniform int TriangleMesh__getNumTriangles(uniform const Shape *uniform _this) {
  const uniform TriangleMesh* uniform this = (const uniform TriangleMesh* uniform) _this;
  return this->numTriangles;
}

void TriangleMesh__postIntersect(uniform const Shape *uniform _this,
                                 const varying Ray &ray, 
                                 varying DifferentialGeometry &dg) 
{
  uniform const TriangleMesh* uniform this = (uniform const TriangleMesh *uniform)_this;
  const vec4i tri = this->triangle[ray.id1];

  const float u = ray.u;
  const float v = ray.v;
  const float w = 1.0f-u-v;
  const float t = ray.tfar;

  dg.P  = add(ray.org, mul(t,ray.dir));
  dg.Ng = normalize(ray.Ng);
  dg.error = max(abs(ray.tfar),reduce_max(abs(dg.P)));

  if (this->texcoord) {
    const vec2f st0 = this->texcoord[tri.x];
    const vec2f st1 = this->texcoord[tri.y];
    const vec2f st2 = this->texcoord[tri.z];
    dg.st = add(add(mul(w,st0),mul(u,st1)),mul(v,st2));
  } else {
    dg.st = make_vec2f(u,v);
  }

  
  if (this->normal) {
    const vec3f n0 = make_vec3f(this->normal[tri.x]);
    const vec3f n1 = make_vec3f(this->normal[tri.y]);
    const vec3f n2 = make_vec3f(this->normal[tri.z]);
    vec3f Ns = add(add(mul(w,n0),mul(u,n1)),mul(v,n2));
    float len2 = dot(Ns,Ns);
    cif (len2 > 0.f) {
      Ns = mul(rsqrt(len2), Ns);
    } else {
      Ns = dg.Ng;
    }
    cif (dot(Ns,dg.Ng) < 0.f) Ns = neg(Ns);
    dg.Ns = Ns;
  } else {
    dg.Ns = dg.Ng; 
  }

  /*! face forward normals */
  cif (dot(dg.Ng,ray.dir) >= 0.0f) {
    dg.Ng = neg(dg.Ng); dg.Ns = neg(dg.Ns);
  }
}

void TriangleMesh__extract(const uniform Shape* uniform _this,
                           const uniform int& id0,
                           uniform RTCTriangle* uniform triangles,
                           uniform int& numTriangles,
                           uniform RTCVertex* uniform vertices,
                           uniform int& numVertices,
                           uniform vec3f& lower_o,
                           uniform vec3f& upper_o)
{
  const uniform TriangleMesh* uniform this = (const uniform TriangleMesh* uniform) _this;

  vec3f lower = make_vec3f(+inf);
  vec3f upper = make_vec3f(-inf);
  foreach(j=0 ... this->numTriangles) {
    const vec4i tri = this->triangle[j];
    triangles[numTriangles+j].v0 = numVertices+tri.x;
    triangles[numTriangles+j].v1 = numVertices+tri.y;
    triangles[numTriangles+j].v2 = numVertices+tri.z;
    triangles[numTriangles+j].id0 = id0;
    triangles[numTriangles+j].id1 = j;
  }
  foreach(j=0 ... this->numPositions) {
    const vec3f p = make_vec3f(this->position[j]);
    vertices[numVertices+j].x = p.x;
    vertices[numVertices+j].y = p.y;
    vertices[numVertices+j].z = p.z;
    lower = min(lower,p);
    upper = max(upper,p);
  }
  numTriangles+=this->numTriangles;
  numVertices +=this->numPositions;
  lower_o = __reduce_min(lower);
  upper_o = __reduce_max(upper);
}

void TriangleMesh__Destructor(uniform RefCount* uniform _this)
{
  uniform TriangleMesh* uniform this = (uniform TriangleMesh* uniform) _this;
  delete[] this->position;
  delete[] this->motion;
  delete[] this->normal;
  delete[] this->texcoord;
  delete[] this->triangle;
  RefCount__Destructor(_this);
}

void* uniform TriangleMesh__transform(const uniform Shape* uniform _this,  const uniform AffineSpace3f& xfm)
{
  const uniform TriangleMesh* uniform this = (const uniform TriangleMesh* uniform) _this;
  uniform TriangleMesh *uniform mesh = uniform new uniform TriangleMesh;

  Shape__Constructor(&mesh->base,
                     TriangleMesh__Destructor,
                     TriangleMesh__getNumVertices,
                     TriangleMesh__getNumTriangles,
                     TriangleMesh__postIntersect,
                     TriangleMesh__transform,
                     TriangleMesh__extract);
  
  uniform int numVertices = this->numPositions;
  uniform int numTriangles = this->numTriangles;
  mesh->numPositions = numVertices;
  mesh->numTriangles = numTriangles;
  mesh->position = NULL;
  mesh->motion = NULL;
  mesh->normal = NULL;
  mesh->texcoord = NULL;
  mesh->triangle = NULL;

  if (this->position) {
    mesh->position = uniform new uniform vec3fa[numVertices];
    foreach(i=0 ... numVertices)
      mesh->position[i] = make_vec3fa(xfmPoint(xfm,make_vec3f(this->position[i])));
  }

  if (this->motion) {
    mesh->motion = uniform new uniform vec3fa[numVertices];
    foreach(i=0 ... numVertices)
      mesh->motion[i] = make_vec3fa(xfmVector(xfm,make_vec3f(this->motion[i])));
  }

  if (this->normal) {
    uniform LinearSpace3f normal2world = transposed(rcp(xfm.l));
    mesh->normal = uniform new uniform vec3fa[numVertices];
    foreach(i=0 ... numVertices)
      mesh->normal[i] = make_vec3fa(xfmVector(normal2world,make_vec3f(this->normal[i])));
  }
  
  if (this->texcoord) {
    mesh->texcoord = uniform new uniform vec2f[numVertices];
    foreach(i=0 ... numVertices)
      mesh->texcoord[i] = this->texcoord[i];
  }

  if (this->triangle) {
    mesh->triangle = uniform new uniform vec4i[numTriangles];
    foreach(i=0 ... numTriangles)
      mesh->triangle[i] = this->triangle[i];
  }    
  return mesh;
}

void TriangleMesh__Constructor(uniform TriangleMesh* uniform this,
                               uniform vec3fa *uniform position,
                               uniform vec3fa *uniform motion,
                               uniform vec3fa *uniform normal,
                               uniform vec2f  *uniform texcoord,
                               uniform int     numVertices,
                               uniform vec4i  *uniform triangle,
                               uniform int     numTriangles)
{
  Shape__Constructor(&this->base,
                     TriangleMesh__Destructor,
                     TriangleMesh__getNumVertices,
                     TriangleMesh__getNumTriangles,
                     TriangleMesh__postIntersect,
                     TriangleMesh__transform,
                     TriangleMesh__extract);

  this->numPositions = numVertices;
  this->numTriangles = numTriangles;
  this->position = NULL;
  this->motion = NULL;
  this->normal = NULL;
  this->texcoord = NULL;
  this->triangle = NULL;

  if (position) {
    this->position = uniform new uniform vec3fa[numVertices];
    foreach(i=0 ... numVertices)
      this->position[i] = position[i];
  }

  if (motion) {
    this->motion = uniform new uniform vec3fa[numVertices];
    foreach(i=0 ... numVertices)
      this->motion[i] = motion[i];
  }

  if (normal) {
    this->normal = uniform new uniform vec3fa[numVertices];
    foreach(i=0 ... numVertices)
      this->normal[i] = normal[i];
  }
  
  if (texcoord) {
    this->texcoord = uniform new uniform vec2f[numVertices];
    foreach(i=0 ... numVertices)
      this->texcoord[i] = texcoord[i];
  }

  if (triangle) {
    this->triangle = uniform new uniform vec4i[numTriangles];
    foreach(i=0 ... numTriangles)
      this->triangle[i] = triangle[i];
  }    
}

///////////////////////////////////////////////////////////////////////////////
// External API

export void* uniform TriangleMesh__new( uniform vec3fa *uniform position,
                                        uniform vec3fa *uniform motion,
                                        uniform vec3fa *uniform normal,
                                        uniform vec2f  *uniform texcoord,
                                        uniform int     numPositions,
                                        uniform vec4i  *uniform triangle,
                                        uniform int     numTriangles)
{
  uniform TriangleMesh* uniform this = uniform new uniform TriangleMesh;
  TriangleMesh__Constructor(this,
                            (uniform vec3fa *uniform) position,
                            (uniform vec3fa *uniform) motion,
                            (uniform vec3fa *uniform) normal,
                            texcoord,
                            numPositions,
                            triangle,
                            numTriangles);
  return this;
}
