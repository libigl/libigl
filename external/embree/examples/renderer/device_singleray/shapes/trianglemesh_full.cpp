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

#include "shapes/trianglemesh.h"

namespace embree
{
  TriangleMeshFull::TriangleMeshFull (const Parms& parms)
    : Shape(parms)
  {
    if (Variant v = parms.getData("positions")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong position format");
      position.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) position[i] = v.data->getVector3f(i);
    }
    if (Variant v = parms.getData("motions")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong motion vector format");
      motion.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) motion[i] = v.data->getVector3f(i);
    }
    if (Variant v = parms.getData("normals")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong normal format");
      normal.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) normal[i] = v.data->getVector3f(i);
    }
    if (Variant v = parms.getData("tangent_x")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong tangent format");
      tangent_x.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) tangent_x[i] = v.data->getVector3f(i);
    }
    if (Variant v = parms.getData("tangent_y")) {
      if (!v.data || v.type != Variant::FLOAT3) throw std::runtime_error("wrong tangent format");
      tangent_y.resize(v.data->size());
      for (size_t i=0; i<v.data->size(); i++) tangent_y[i] = v.data->getVector3f(i);
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
      for (size_t i=0; i<v.data->size(); i++) triangles[i] = v.data->getVector3i(i);
    }
  }

  Ref<Shape> TriangleMeshFull::transform(const AffineSpace3f& xfm) const
  {
    /*! do nothing for identity matrix */
    if (xfm == AffineSpace3f(one))
      return (Shape*)this;

    /*! create transformed */
    TriangleMeshFull* mesh = new TriangleMeshFull(ty);
    mesh->position.resize(position.size());
    for (size_t i=0; i<position.size(); i++) mesh->position[i] = xfmPoint(xfm,position[i]);
    mesh->motion.resize(motion.size());
    for (size_t i=0; i<motion.size(); i++) mesh->motion[i] = xfmVector(xfm,motion[i]);
    mesh->normal.resize(normal.size()  );
    for (size_t i=0; i<normal.size(); i++) mesh->normal[i] = xfmNormal(xfm,normal[i]);
    mesh->tangent_x.resize(tangent_x.size()  );
    for (size_t i=0; i<tangent_x.size(); i++) mesh->tangent_x[i] = xfmVector(xfm,tangent_x[i]);
    mesh->tangent_y.resize(tangent_y.size()  );
    for (size_t i=0; i<tangent_y.size(); i++) mesh->tangent_y[i] = xfmVector(xfm,tangent_y[i]);
    mesh->texcoord  = texcoord;
    mesh->triangles = triangles;
    return mesh;
  }

  size_t TriangleMeshFull::numTriangles() const {
    return triangles.size();
  }

  size_t TriangleMeshFull::numVertices() const {
    return position.size() + motion.size();
  }

  BBox3f TriangleMeshFull::extract(size_t id, RTCTriangle* triangles_o, size_t& numTriangles, Vec3fa* positions_o, size_t& numVertices) const
  {
    BBox3f bounds = empty;
    if (motion.size()) 
    {
      for (size_t j=0; j<triangles.size(); j++) {
        const TriangleMeshFull::Triangle& tri = triangles[j];
        triangles_o[numTriangles++] = RTCTriangle((int)numVertices+2*tri.v0,(int)numVertices+2*tri.v1,(int)numVertices+2*tri.v2,(int)(0x80000000 | id),(int)j);
      }
      for (size_t j=0; j<position.size(); j++) {
        const Vector3f p = position[j];
        const Vector3f dpdt = motion[j];
        positions_o[numVertices++] = Vector3f(p.x,p.y,p.z);
        positions_o[numVertices++] = Vector3f(dpdt.x,dpdt.y,dpdt.z);
        bounds.grow(p);
        bounds.grow(p+dpdt);
      }
    }
    else
    {
      for (size_t j=0; j<triangles.size(); j++) {
        const TriangleMeshFull::Triangle& tri = triangles[j];
        triangles_o[numTriangles++] = RTCTriangle((int)numVertices+tri.v0,(int)numVertices+tri.v1,(int)numVertices+tri.v2,(int)id,(int)j);
      }
      for (size_t j=0; j<position.size(); j++) {
        const Vector3f p = position[j];
        positions_o[numVertices++] = Vector3f(p.x,p.y,p.z);
        bounds.grow(p);
      }
    }
    return bounds;
  }

  void TriangleMeshFull::postIntersect(const Ray& ray, DifferentialGeometry& dg) const
  {
    const Triangle& tri = triangles[ray.id1];
    Vector3f p0 = position[tri.v0], p1 = position[tri.v1], p2 = position[tri.v2];
    if (unlikely(motion.size())) {
      p0 += ray.time * motion[tri.v0];
      p1 += ray.time * motion[tri.v1];
      p2 += ray.time * motion[tri.v2];
    }
    const float u = ray.u, v = ray.v, w = 1.0f-u-v, t = ray.tfar;

    const Vector3f dPdu = p1-p0, dPdv = p2-p0;
    dg.P  = ray.org+t*ray.dir;
    dg.Ng = normalize(ray.Ng);

    /* interpolate texture coordinates */
    float dsdu, dtdu;
    float dsdv, dtdv;
    if (texcoord.size()) {
      const Vec2f st0 = texcoord[tri.v0];
      const Vec2f st1 = texcoord[tri.v1];
      const Vec2f st2 = texcoord[tri.v2];
      dg.st = st0*w + st1*u + st2*v;
      dsdu = st1.x-st0.x; dtdu = st1.y-st0.y;
      dsdv = st2.x-st0.x; dtdv = st2.y-st0.y;
    }
    else {
      dg.st = Vec2f(u,v);
      dsdu = 1; dtdu = 0;
      dsdv = 0; dtdv = 1;
    }

    /* interpolate shading normal */
    if (normal.size())
    {
      const Vector3f n0 = normal[tri.v0], n1 = normal[tri.v1], n2 = normal[tri.v2];
      Vector3f Ns = w*n0 + u*n1 + v*n2;
      float len2 = dot(Ns,Ns);
      Ns = len2 > 0 ? Ns*rsqrt(len2) : Vector3f(dg.Ng);
      if (dot(Ns,dg.Ng) < 0) Ns = -Ns;
      dg.Ns = Ns;
    }
    else
      dg.Ns = dg.Ng;

    /* interpolate x tangent direction */
    if (tangent_x.size()) { 
      const Vector3f t0 = tangent_x[tri.v0], t1 = tangent_x[tri.v1], t2 = tangent_x[tri.v2];
      dg.Tx = w*t0 + u*t1 + v*t2;
    }
    else {
      const Vector3f dPds = normalize(dPdu*dtdv - dPdv*dtdu);
      dg.Tx = normalize(dPds-dot(dPds,dg.Ns)*dg.Ns);
    }

    /* interpolate y tangent direction */
    if (tangent_y.size()) {
      const Vector3f t0 = tangent_y[tri.v0], t1 = tangent_y[tri.v1], t2 = tangent_y[tri.v2];
      dg.Ty = w*t0 + u*t1 + v*t2;
    } else {
      const Vector3f dPdt = normalize(dPdv*dsdu - dPdu*dsdv);
      dg.Ty = normalize(dPdt-dot(dPdt,dg.Ns)*dg.Ns);
    }

    dg.error = max(abs(ray.tfar),reduce_max(abs(dg.P)));
  }
}
