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

#include "shapes/trianglemesh_normals.h"

namespace embree
{
  TriangleMeshWithNormals::TriangleMeshWithNormals(const Parms& parms)
    : Shape(parms)
  {
    if (Variant v = parms.getData("positions")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong position format");
      vertices.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) vertices[i].p = v.data->getVector3f(i);
    }
    if (Variant v = parms.getData("normals")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong normal format");
      vertices.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) vertices[i].n = v.data->getVector3f(i);
    }
    if (Variant v = parms.getData("indices")) {
      if (!v.data || v.type != Variant::INT3) throw std::runtime_error("wrong triangle format");
      triangles.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) triangles[i] = v.data->getVector3i(i);
    }
  }

  Ref<Shape> TriangleMeshWithNormals::transform(const AffineSpace3f& xfm) const
  {
    /*! do nothing for identity matrix */
    if (xfm == AffineSpace3f(one))
      return (Shape*)this;

    /*! create transformed mesh */
    TriangleMeshWithNormals* mesh = new TriangleMeshWithNormals(ty);
    mesh->vertices.resize(vertices.size());
    for (size_t i=0; i<vertices.size(); i++) mesh->vertices[i].p = xfmPoint (xfm,*(Vector3f*)(void*)&vertices[i].p);
    for (size_t i=0; i<vertices.size(); i++) mesh->vertices[i].n = xfmNormal(xfm,*(Vector3f*)(void*)&vertices[i].n);
    mesh->triangles = triangles;
    return mesh;
  }

  size_t TriangleMeshWithNormals::numTriangles() const {
    return triangles.size();
  }

  size_t TriangleMeshWithNormals::numVertices() const {
    return vertices.size();
  }

  BBox3f TriangleMeshWithNormals::extract(size_t id, RTCTriangle* triangles_o, size_t& numTriangles, Vec3fa* positions_o, size_t& numVertices) const
  {
    BBox3f bounds = empty;
    for (size_t j=0; j<triangles.size(); j++) {
      const TriangleMeshWithNormals::Triangle& tri = triangles[j];
      triangles_o[numTriangles++] = RTCTriangle((int)numVertices+tri.v0,(int)numVertices+tri.v1,(int)numVertices+tri.v2,(int)id,(int)j);
    }
    for (size_t j=0; j<vertices.size(); j++) {
      const Vector3f p = vertices[j].p;
      positions_o[numVertices++] = Vector3f(p.x,p.y,p.z);
      bounds.grow(p);
    }
    return bounds;
  }

  void TriangleMeshWithNormals::postIntersect(const Ray& ray, DifferentialGeometry& dg) const
  {
    const Triangle& tri = triangles[ray.id1];
    const Vertex& v0 = vertices[tri.v0];
    const Vertex& v1 = vertices[tri.v1];
    const Vertex& v2 = vertices[tri.v2];
    float u = ray.u, v = ray.v, w = 1.0f-u-v, t = ray.tfar;
    Vector3f dPdu = v1.p - v0.p, dPdv = v2.p - v0.p;
    dg.P = ray.org+t*ray.dir;
    dg.Ng = normalize(ray.Ng);
    dg.st = Vec2f(u,v);
    Vector3f Ns = w*v0.n + u*v1.n + v*v2.n;
    float len2 = dot(Ns,Ns);
    Ns = len2 > 0 ? Ns*rsqrt(len2) : Vector3f(dg.Ng);
    if (dot(Ns,dg.Ng) < 0) Ns = -Ns;
    dg.Ns = Ns;
    dg.Tx = dPdu;
    dg.Ty = dPdv;
    dg.error = max(abs(ray.tfar),reduce_max(abs(dg.P)));
  }
}
