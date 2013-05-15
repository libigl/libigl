// ======================================================================== //
// Copyright 2009-2012 Intel Corporation                                    //
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

#include "shapes/trianglemesh.h"
#include "shapes/trianglemesh_normals.h"

namespace embree
{
  Ref<Shape> TriangleMesh::Handle::create(const Parms& parms) 
  {
    if (parms.size() == 3 && parms.getData("positions") && parms.getData("normals") && parms.getData("indices"))
      return new TriangleMeshWithNormals(parms);
    else 
      return new TriangleMesh(parms);
  }

  TriangleMesh::TriangleMesh (const Parms& parms)
  {
    if (Variant v = parms.getData("positions")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong position format");
      position.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) position[i] = v.data->getVec3f(i);
    }
    if (Variant v = parms.getData("motions")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong motion vector format");
      motion.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) motion[i] = v.data->getVec3f(i);
    }
    if (Variant v = parms.getData("normals")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong normal format");
      normal.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) normal[i] = v.data->getVec3f(i);
    }
    if (Variant v = parms.getData("texcoords")) {
      if (!v.data || v.type != Variant::FLOAT2) throw std::runtime_error("wrong texcoords0 format");
      texcoord.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) texcoord[i] = v.data->getVec2f(i);
    }
    if (Variant v = parms.getData("texcoords0")) {
      if (!v.data || v.type != Variant::FLOAT2) throw std::runtime_error("wrong texcoords0 format");
      texcoord.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) texcoord[i] = v.data->getVec2f(i);
    }
    if (Variant v = parms.getData("indices")) {
      if (!v.data || v.type != Variant::INT3) throw std::runtime_error("wrong triangle format");
      triangles.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) triangles[i] = v.data->getVec3i(i);
    }
  }

  Ref<Shape> TriangleMesh::transform(const AffineSpace3f& xfm) const
  {
    /*! do nothing for identity matrix */
    if (xfm == AffineSpace3f(one))
      return (Shape*)this;

    /*! create transformed */
    TriangleMesh* mesh = new TriangleMesh;
    mesh->position.resize(position.size());
    for (size_t i=0; i<position.size(); i++) mesh->position[i] = xfmPoint(xfm,position[i]);
    mesh->motion.resize(motion.size());
    for (size_t i=0; i<motion.size(); i++) mesh->motion[i] = xfmVector(xfm,motion[i]);
    mesh->normal.resize(normal.size()  );
    for (size_t i=0; i<normal.size(); i++) mesh->normal[i] = xfmNormal(xfm,normal[i]);
    mesh->texcoord  = texcoord;
    mesh->triangles = triangles;
    return (Shape*)mesh;
  }

  size_t TriangleMesh::numTriangles() const {
    return triangles.size();
  }

  size_t TriangleMesh::numVertices() const {
    return position.size() + motion.size();
  }

  BBox3f TriangleMesh::extract(size_t id, BuildTriangle* triangles_o, size_t& numTriangles, BuildVertex* vertices_o, size_t& numVertices) const
  {
    BBox3f bounds = empty;
    if (motion.size()) 
    {
      for (size_t j=0; j<triangles.size(); j++) {
        const TriangleMesh::Triangle& tri = triangles[j];
        triangles_o[numTriangles++] = BuildTriangle((int)numVertices+2*tri.v0,(int)numVertices+2*tri.v1,(int)numVertices+2*tri.v2,(int)(0x80000000 | id),(int)j);
      }
      for (size_t j=0; j<position.size(); j++) {
        const Vec3f p = position[j];
        const Vec3f dpdt = motion[j];
        vertices_o[numVertices++] = BuildVertex(p.x,p.y,p.z);
        vertices_o[numVertices++] = BuildVertex(dpdt.x,dpdt.y,dpdt.z);
        bounds.grow(p);
        bounds.grow(p+dpdt);
      }
    }
    else
    {
      for (size_t j=0; j<triangles.size(); j++) {
        const TriangleMesh::Triangle& tri = triangles[j];
        triangles_o[numTriangles++] = BuildTriangle((int)numVertices+tri.v0,(int)numVertices+tri.v1,(int)numVertices+tri.v2,(int)id,(int)j);
      }
      for (size_t j=0; j<position.size(); j++) {
        const Vec3f p = position[j];
        vertices_o[numVertices++] = BuildVertex(p.x,p.y,p.z);
        bounds.grow(p);
      }
    }
    return bounds;
  }

  void TriangleMesh::postIntersect(const Ray& ray, DifferentialGeometry& dg) const
  {
    const Triangle& tri = triangles[dg.id1];
    Vec3f p0 = position[tri.v0], p1 = position[tri.v1], p2 = position[tri.v2];
    if (unlikely(motion.size())) {
      p0 += ray.time * motion[tri.v0];
      p1 += ray.time * motion[tri.v1];
      p2 += ray.time * motion[tri.v2];
    }
    const float u = dg.u, v = dg.v, w = 1.0f-u-v, t = dg.t;

    const Vec3f dPdu = p1-p0, dPdv = p2-p0;
    dg.P  = ray.org+t*ray.dir;
    dg.Ng = normalize(cross(dPdv,dPdu));

    /* interpolate texture coordinates */
    if (texcoord.size()) {
      const Vec2f st0 = texcoord[tri.v0];
      const Vec2f st1 = texcoord[tri.v1];
      const Vec2f st2 = texcoord[tri.v2];
      dg.st = st0*w + st1*u + st2*v;
    }
    else {
      dg.st = Vec2f(u,v);
    }

    /* interpolate shading normal */
    if (normal.size())
    {
      const Vec3f n0 = normal[tri.v0], n1 = normal[tri.v1], n2 = normal[tri.v2];
      Vec3f Ns = w*n0 + u*n1 + v*n2;
      float len2 = dot(Ns,Ns);
      Ns = len2 > 0 ? Ns*rsqrt(len2) : dg.Ng;
      if (dot(Ns,dg.Ng) < 0) Ns = -Ns;
      dg.Ns = Ns;
    }
    else
      dg.Ns = dg.Ng;

    dg.error = max(abs(dg.t),reduce_max(abs(dg.P)));
  }
}
