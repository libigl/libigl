// +-------------------------------------------------------------------------
// | mesh.tpp
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#pragma once

#include <algorithm>
#include "unsafeRayTriIsct.h"
#include <cfloat>
#include <cmath>
#include <sstream>

#include "unionFind.h"

// constructors
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh() {}
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh(Mesh &&cp)
    : tris(cp.tris), verts(cp.verts)
{}
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh(const RawMesh<VertData,TriData> &raw) :
    tris(raw.triangles.size()), verts(raw.vertices)
{
    // fill out the triangles
    for(uint i=0; i<raw.triangles.size(); i++) {
        tris[i].data = raw.triangles[i];
        tris[i].a = raw.triangles[i].a;
        tris[i].b = raw.triangles[i].b;
        tris[i].c = raw.triangles[i].c;
    }
}
template<class VertData, class TriData>
Mesh<VertData,TriData>::~Mesh() {}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::operator=(Mesh &&src)
{
    tris = src.tris;
    verts = src.verts;
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::valid() const
{
    for(uint i=0; i<verts.size(); i++) {
        if(!std::isfinite(verts[i].pos.x) ||
           !std::isfinite(verts[i].pos.y) ||
           !std::isfinite(verts[i].pos.z)) {
            std::ostringstream message;
            message << "vertex #" << i << " has non-finite coordinates: "
                    << verts[i].pos;
            CORK_ERROR(message.str());
            return false;
        }
    }
    
    for(uint i=0; i<tris.size(); i++) {
        if(tris[i].a >= verts.size() ||
           tris[i].b >= verts.size() ||
           tris[i].c >= verts.size()) {
            std::ostringstream message;
            message << "triangle #" << i << " should have indices in "
                    << "the range 0 to " << (verts.size()-1)
                    << ", but it has invalid indices: "
                    << tris[i].a << ", " << tris[i].b << ", " << tris[i].c;
            CORK_ERROR(message.str());
            return false;
        }
    }
    
    return true;
}

template<class VertData, class TriData>
RawMesh<VertData,TriData> Mesh<VertData,TriData>::raw() const
{
    RawMesh<VertData,TriData> result;
    result.vertices = verts;
    result.triangles.resize(tris.size());
    for(uint i=0; i<tris.size(); i++) {
        result.triangles[i]   = tris[i].data;
        result.triangles[i].a = tris[i].a;
        result.triangles[i].b = tris[i].b;
        result.triangles[i].c = tris[i].c;
    }
    return result;
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::disjointUnion(const Mesh &cp)
{
    uint oldVsize = verts.size();
    uint oldTsize = tris.size();
    uint cpVsize  = cp.verts.size();
    uint cpTsize  = cp.tris.size();
    uint newVsize = oldVsize + cpVsize;
    uint newTsize = oldTsize + cpTsize;
    
    std::vector<int> v_remap(cpVsize); // oh this is obvious...
    verts.resize(newVsize);
    tris.resize(newTsize);
    
    for(uint i=0; i<cpVsize; i++)
        verts[oldVsize + i] = cp.verts[i];
    
    for(uint i=0; i<cpTsize; i++) {
        auto &tri = tris[oldTsize + i];
        tri = cp.tris[i];
        tri.a += oldVsize;
        tri.b += oldVsize;
        tri.c += oldVsize;
    }
}





// Picking.
// Dumb Implementation just passes over all triangles w/o any precomputed
// acceleration structure
template<class VertData, class TriData>
typename Mesh<VertData,TriData>::Isct
    Mesh<VertData,TriData>::pick(Ray3d ray)
{
    Isct result;
    result.ray = ray;
    result.exists = false;
    
    double mint = DBL_MAX;
    
    // pass all triangles over ray
    for(uint i=0; i<tris.size(); i++) {
        const Tri  &tri = tris[i];
        
        uint   a = tri.a;
        uint   b = tri.b;
        uint   c = tri.c;
        Vec3d va = verts[a].pos;
        Vec3d vb = verts[b].pos;
        Vec3d vc = verts[c].pos;
        // normalize vertex order (to prevent leaks)
        if(a > b) { std::swap(a, b); std::swap(va, vb); }
        if(b > c) { std::swap(b, c); std::swap(vb, vc); }
        if(a > b) { std::swap(a, b); std::swap(va, vb); }
        
        double t;
        Vec3d  bary;
        if(isct_ray_triangle(ray, va, vb, vc, &t, &bary)) {
            if(t > 0 && t < mint) {
                result.exists = true;
                mint = t;
                result.tri_id = i;
                result.isct = ray.p + t * ray.r;
                result.bary = bary;
            }
        }
    }
    
    return result;
}




template<class VertData, class TriData>
bool Mesh<VertData,TriData>::isClosed()
{
    EGraphCache<int> chains = createEGraphCache<int>();
    chains.for_each([&](uint i, uint j, EGraphEntry<int> &entry) {
        entry.data = 0;
    });
    // count up how many times each edge is encountered in one
    // orientation vs. the other
    for(Tri &tri : tris) {
        chains(tri.a, tri.b).data ++;
        chains(tri.b, tri.a).data --;
        
        chains(tri.b, tri.c).data ++;
        chains(tri.c, tri.b).data --;
        
        chains(tri.c, tri.a).data ++;
        chains(tri.a, tri.c).data --;
    }
    // now go through and see if any of these are non-zero
    bool closed = true;
    chains.for_each([&](uint i, uint j, EGraphEntry<int> &entry) {
        if(entry.data != 0)
            closed = false;
    });
    return closed;
}




static inline
bool contains(const ShortVec<uint, 8> &list, uint item)
{
    for(uint k : list)
        if(k == item)
            return true;
    return false;
}

template<class VertData, class TriData>
typename Mesh<VertData,TriData>::NeighborCache
    Mesh<VertData,TriData>::createNeighborCache()
{
    NeighborCache result;
    result.skeleton.resize(verts.size());
    
    for(uint tid = 0; tid < tris.size(); tid++) {
        const Tri &tri = tris[tid];
        
        result(tri.a, tri.b).tids.push_back(tid);
        result(tri.b, tri.a).tids.push_back(tid);
        
        result(tri.a, tri.c).tids.push_back(tid);
        result(tri.c, tri.a).tids.push_back(tid);
        
        result(tri.b, tri.c).tids.push_back(tid);
        result(tri.c, tri.b).tids.push_back(tid);
    }
    
    return result;
}

// This function signature is an amazing disaster...
#ifdef _WIN32
template<class VertData, class TriData>
template<class Edata>
typename Mesh<VertData,TriData>::EGraphCache<Edata>
#else
template<class VertData, class TriData>
template<class Edata>
typename Mesh<VertData,TriData>::template EGraphCache<Edata>
#endif
    Mesh<VertData,TriData>::createEGraphCache()
{
    EGraphCache<Edata> result;
    result.skeleton.resize(verts.size());
    
    for(uint tid = 0; tid < tris.size(); tid++) {
        const Tri &tri = tris[tid];
        
        result(tri.a, tri.b).tids.push_back(tid);
        result(tri.b, tri.a).tids.push_back(tid);
        
        result(tri.a, tri.c).tids.push_back(tid);
        result(tri.c, tri.a).tids.push_back(tid);
        
        result(tri.b, tri.c).tids.push_back(tid);
        result(tri.c, tri.b).tids.push_back(tid);
    }
    
    return result;
}


template<class VertData, class TriData>
std::vector<uint> Mesh<VertData,TriData>::getComponentIds()
{
    UnionFind uf(verts.size());
    for(const Tri &tri : tris) {
        uf.unionIds(tri.a, tri.b);
        uf.unionIds(tri.a, tri.c);
    }
    
    return uf.dump();
}









