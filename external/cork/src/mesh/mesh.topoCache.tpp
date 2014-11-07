// +-------------------------------------------------------------------------
// | mesh.topoCache.tpp
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

#include "iterPool.h"

/*
 *  Allows for topological algorithms to manipulate
 *  a more familiar pointer data structure based on a simplicial complex.
 *  This structure can be regenerated from the more basic
 *  vertex/triangle arrays using
 *      createTopoCache()
 *  Once manipulations have been satisfactorily performed,
 *  the underlying vertex/triangle arrays can be cleaned up for
 *  further use by topologically insensitive algorithms by
 *      commitTopoCache()
 */

#define INVALID_ID uint(-1)

struct TopoVert {
    uint                    ref;        // index to actual data
    void*                   data;       // algorithm specific handle
                        
    ShortVec<Tptr, 8>       tris;       // triangles this vertex is incident on
    ShortVec<Eptr, 8>       edges;      // edges this vertex is incident on
};

struct TopoEdge {
    void*                   data;       // algorithm specific handle
    
    Vptr                    verts[2];   // endpoint vertices
    ShortVec<Tptr, 2>       tris;       // incident triangles
};

struct TopoTri {
    uint                    ref;        // index to actual data
    void*                   data;       // algorithm specific handle
    
    Vptr                    verts[3];   // vertices of this triangle
    Eptr                    edges[3];   // edges of this triangle
                                        // opposite to the given vertex
};


template<class VertData, class TriData>
struct Mesh<VertData, TriData>::TopoCache {
    IterPool<TopoVert>    verts;
    IterPool<TopoEdge>    edges;
    IterPool<TopoTri>     tris;
    
    Mesh *mesh;
    TopoCache(Mesh *owner);
    virtual ~TopoCache() {}
    
    // until commit() is called, the Mesh::verts and Mesh::tris
    // arrays will still contain garbage entries
    void commit();
    
    bool isValid();
    void print();
    
    // helpers to create bits and pieces
    inline Vptr newVert();
    inline Eptr newEdge();
    inline Tptr newTri();
    
    // helpers to release bits and pieces
    inline void freeVert(Vptr);
    inline void freeEdge(Eptr);
    inline void freeTri(Tptr);
    
    // helper to delete geometry in a structured way
    inline void deleteTri(Tptr);
    
    // helper to flip triangle orientation
    inline void flipTri(Tptr);
    
private:
    void init();
};


template<class VertData, class TriData> inline
Vptr Mesh<VertData, TriData>::TopoCache::newVert()
{
    uint        ref         = mesh->verts.size();
                mesh->verts.push_back(VertData());
    Vptr        v           = verts.alloc(); // cache.verts
                v->ref      = ref;
                return v;
}
template<class VertData, class TriData> inline
Eptr Mesh<VertData, TriData>::TopoCache::newEdge()
{
    Eptr        e           = edges.alloc(); // cache.edges
                return e;
}
template<class VertData, class TriData> inline
Tptr Mesh<VertData, TriData>::TopoCache::newTri()
{
    uint        ref         = mesh->tris.size();
                mesh->tris.push_back(Tri());
    Tptr        t           = tris.alloc(); // cache.tris
                t->ref      = ref;
                return t;
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeVert(Vptr v)
{
    verts.free(v);
}
template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeEdge(Eptr e)
{
    edges.free(e);
}
template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeTri(Tptr t)
{
    tris.free(t);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::deleteTri(Tptr tri)
{
    // first, unhook the triangle from its faces
    for(uint k=0; k<3; k++) {
        Vptr            v                   = tri->verts[k];
                        v->tris.erase(tri);
        Eptr            e                   = tri->edges[k];
                        e->tris.erase(tri);
    }
    // now, let's check for any edges which no longer border triangles
    for(uint k=0; k<3; k++) {
        Eptr            e                   = tri->edges[k];
        if(e->tris.size() == 0) {
            // delete edge
            // unhook from vertices
            Vptr        v0                  = e->verts[0];
                        v0->edges.erase(e);
            Vptr        v1                  = e->verts[1];
                        v1->edges.erase(e);
            freeEdge(e);
        }
    }
    // now, let's check for any vertices which no longer border triangles
    for(uint k=0; k<3; k++) {
        Vptr            v                   = tri->verts[k];
        if(v->tris.size() == 0) {
            freeVert(v);
        }
    }
    // finally, release the triangle
    freeTri(tri);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::flipTri(Tptr t)
{
    std::swap(t->verts[0], t->verts[1]);
    std::swap(t->edges[0], t->edges[1]);
    std::swap(mesh->tris[t->ref].v[0], mesh->tris[t->ref].v[1]);
}



template<class VertData, class TriData>
Mesh<VertData, TriData>::TopoCache::TopoCache(
    Mesh *owner
) : mesh(owner) {
    init();
}



// support structure for cache construction
struct TopoEdgePrototype {
    uint vid;
    ShortVec<Tptr, 2> tris;
    TopoEdgePrototype() {}
    TopoEdgePrototype(uint v) : vid(v) {}
};
inline TopoEdgePrototype& getTopoEdgePrototype(
    uint a, uint b,
    std::vector< ShortVec<TopoEdgePrototype, 8> > &prototypes
) {
    uint N = prototypes[a].size();
    for(uint i=0; i<N; i++) {
        if(prototypes[a][i].vid == b)
            return prototypes[a][i];
    }
    prototypes[a].push_back(TopoEdgePrototype(b));
    return prototypes[a][N];
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::init()
{
    // first lay out vertices
    std::vector<Vptr> temp_verts(mesh->verts.size()); // need temp. reference
    for(uint i=0; i<mesh->verts.size(); i++) {
        Vptr vert = verts.alloc(); // cache.verts.alloc()
        vert->ref = i;
        temp_verts[i] = vert;
    }
    
    // We need to still do the following
    //  * Generate TopoTris
    //  * Generate TopoEdges
    // ---- Hook up references between
    //  * Triangles and Vertices
    //  * Triangles and Edges
    //  * Vertices and Edges
    
    // We handle two of these items in a pass over the triangles,
    //  * Generate TopoTris
    //  * Hook up Triangles and Vertices
    // building a structure to handle the edges as we go:
    std::vector< ShortVec<TopoEdgePrototype, 8> > edgeacc(mesh->verts.size());
    for(uint i=0; i<mesh->tris.size(); i++) {
        Tptr tri = tris.alloc(); // cache.tris.alloc()
        tri->ref = i;
        const Tri &ref_tri = mesh->tris[i];
        
        // triangles <--> verts
        uint vids[3];
        for(uint k=0; k<3; k++) {
            uint vid = vids[k] = ref_tri.v[k];
            tri->verts[k] = temp_verts[vid];
            temp_verts[vid]->tris.push_back(tri);
        }
        // then, put these in arbitrary but globally consistent order
        if(vids[0] > vids[1])   std::swap(vids[0], vids[1]);
        if(vids[1] > vids[2])   std::swap(vids[1], vids[2]);
        if(vids[0] > vids[1])   std::swap(vids[0], vids[1]);
        // and accrue in structure
        getTopoEdgePrototype(vids[0], vids[1], edgeacc).tris.push_back(tri);
        getTopoEdgePrototype(vids[0], vids[2], edgeacc).tris.push_back(tri);
        getTopoEdgePrototype(vids[1], vids[2], edgeacc).tris.push_back(tri);
    }
    
    // Now, we can unpack the edge accumulation to
    //  * Generate TopoEdges
    //  * Hook up Triangles and Edges
    //  * Hook up Vertices and Edges
    for(uint vid0=0; vid0 < edgeacc.size(); vid0++) {
      for(TopoEdgePrototype &proto : edgeacc[vid0]) {
        uint vid1 = proto.vid;
        Vptr v0 = temp_verts[vid0];
        Vptr v1 = temp_verts[vid1];
        
        Eptr edge = edges.alloc(); // cache.edges.alloc()
        // edges <--> verts
        edge->verts[0] = v0;
        v0->edges.push_back(edge);
        edge->verts[1] = v1;
        v1->edges.push_back(edge);
        // edges <--> tris
        for(Tptr tri : proto.tris) {
            edge->tris.push_back(tri);
            for(uint k=0; k<3; k++) {
                if(v0 != tri->verts[k] && v1 != tri->verts[k]) {
                    tri->edges[k] = edge;
                    break;
                }
            }
        }
    }}
    
    //ENSURE(isValid());
    //print();
}




template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::commit()
{
    //ENSURE(isValid());
    
    // record which vertices are live
    std::vector<bool> live_verts(mesh->verts.size(), false);
    verts.for_each([&](Vptr vert) { // cache.verts
        live_verts[vert->ref] = true;
    });
    
    // record which triangles are live, and record connectivity
    std::vector<bool> live_tris(mesh->tris.size(), false);
    tris.for_each([&](Tptr tri) { // cache.tris
        live_tris[tri->ref] = true;
        for(uint k=0; k<3; k++)
            mesh->tris[tri->ref].v[k] = tri->verts[k]->ref;
    });
    
    // compact the vertices and build a remapping function
    std::vector<uint> vmap(mesh->verts.size());
    uint write = 0;
    for(uint read = 0; read < mesh->verts.size(); read++) {
        if(live_verts[read]) {
            vmap[read] = write;
            mesh->verts[write] = mesh->verts[read];
            write++;
        } else {
            vmap[read] = INVALID_ID;
        }
    }
    mesh->verts.resize(write);
    
    // rewrite the vertex reference ids
    verts.for_each([&](Vptr vert) { // cache.verts
        vert->ref = vmap[vert->ref];
    });
    
    std::vector<uint> tmap(mesh->tris.size());
    write = 0;
    for(uint read = 0; read < mesh->tris.size(); read++) {
        if(live_tris[read]) {
            tmap[read] = write;
            mesh->tris[write] = mesh->tris[read];
            for(uint k=0; k<3; k++)
                mesh->tris[write].v[k] = vmap[mesh->tris[write].v[k]];
            write++;
        } else {
            tmap[read] = INVALID_ID;
        }
    }
    mesh->tris.resize(write);
    
    // rewrite the triangle reference ids
    tris.for_each([&](Tptr tri) { // cache.tris
        tri->ref = tmap[tri->ref];
    });
}



// support functions for validity check
template<class T, class Container> inline
bool count(const Container &contain, const T &val) {
    uint c=0;
    for(const T &t : contain)
        if(t == val)    c++;
    return c;
}
template<class T> inline
bool count2(const T arr[], const T &val) {
    return ((arr[0] == val)? 1 : 0) + ((arr[1] == val)? 1 : 0);
}
template<class T> inline
bool count3(const T arr[], const T &val) {
    return ((arr[0] == val)? 1 : 0) + ((arr[1] == val)? 1 : 0)
                                    + ((arr[2] == val)? 1 : 0);
}

template<class VertData, class TriData>
bool Mesh<VertData, TriData>::TopoCache::isValid()
{
    //print();
    std::set<Vptr> vaddr;
    std::set<Eptr> eaddr;
    std::set<Tptr> taddr;
    verts.for_each([&vaddr](Vptr v) { vaddr.insert(v); });
    edges.for_each([&eaddr](Eptr e) { eaddr.insert(e); });
    tris.for_each( [&taddr](Tptr t) { taddr.insert(t); });
    
    // check verts
    verts.for_each([&](Vptr v) {
        ENSURE(v->ref < mesh->verts.size());
        // make sure each edge pointer goes somewhere and that
        // the pointed-to site also points back correctly
        for(Eptr e : v->edges) {
            ENSURE(eaddr.count(e) > 0); // pointer is good
            ENSURE(count2(e->verts, v) == 1); // back-pointer is good
        }
        for(Tptr t : v->tris) {
            ENSURE(taddr.count(t) > 0);
            ENSURE(count3(t->verts, v) == 1);
        }
    });
    
    // check edges
    edges.for_each([&](Eptr e) {
        // check for non-degeneracy
        ENSURE(e->verts[0] != e->verts[1]);
        for(uint k=0; k<2; k++) {
            Vptr v = e->verts[k];
            ENSURE(vaddr.count(v) > 0);
            ENSURE(count(v->edges, e) == 1);
        }
        for(Tptr t : e->tris) {
            ENSURE(taddr.count(t) > 0);
            ENSURE(count3(t->edges, e) == 1);
        }
    });
    
    // check triangles
    tris.for_each([&](Tptr t) {
        // check for non-degeneracy
        ENSURE(t->verts[0] != t->verts[1] && t->verts[1] != t->verts[2]
                                          && t->verts[0] != t->verts[2]);
        for(uint k=0; k<3; k++) {
            Vptr v = t->verts[k];
            ENSURE(vaddr.count(v) > 0);
            ENSURE(count(v->tris, t) == 1);
            
            Eptr e = t->edges[k];
            ENSURE(eaddr.count(e) == 1);
            ENSURE(count(e->tris, t) == 1);
            
            // also need to ensure that the edges are opposite the
            // vertices as expected
            Vptr v0 = e->verts[0];
            Vptr v1 = e->verts[1];
            ENSURE((v0 == t->verts[(k+1)%3] && v1 == t->verts[(k+2)%3])
                || (v0 == t->verts[(k+2)%3] && v1 == t->verts[(k+1)%3]));
        }
    });
    
    return true;
}







std::ostream& operator<<(std::ostream &out, const TopoVert& vert)
{
    out << "ref(" << vert.ref << ") "
        << "e(" << vert.edges.size() << "):";
    for(Eptr e : vert.edges)
        out << e << ";";
    out << " "
        << "t(" << vert.tris.size() << "):";
    for(Tptr t : vert.tris)
        out << t << ";";
    return out;
}

std::ostream& operator<<(std::ostream &out, const TopoEdge& edge)
{
    out << "v(2):" << edge.verts[0] << "(" << edge.verts[0]->ref << ");"
                   << edge.verts[1] << "(" << edge.verts[1]->ref << ");";
    out << " "
        << "t(" << edge.tris.size() << "):";
    for(Tptr t : edge.tris)
        out << t << ";";
    return out;
}

std::ostream& operator<<(std::ostream &out, const TopoTri& tri)
{
    out << "ref(" << tri.ref << ") ";
    out << "v(3):" << tri.verts[0] << "(" << tri.verts[0]->ref << ");"
                   << tri.verts[1] << "(" << tri.verts[1]->ref << ");"
                   << tri.verts[2] << "(" << tri.verts[2]->ref << ");";
    out << " ";
    out << "e(3):" << tri.edges[0] << ";"
                   << tri.edges[1] << ";"
                   << tri.edges[2] << ";";
    return out;
}


template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::print()
{
    using std::cout;
    using std::endl;
    
    cout << "dumping remeshing cache for debug..." << endl;
    cout << "TRIS" << endl;
    int tri_count = 0;
    tris.for_each([&](Tptr t) {
        cout << " " << t << ": " << *t << endl;
        tri_count++;
    });
    cout << "There were " << tri_count << " TRIS" << endl;
    cout << "EDGES" << endl;
    int edge_count = 0;
    edges.for_each([&](Eptr e) {
        cout << " " << e << ": " << endl;
        cout << "  v " << e->verts[0] << "; "
                       << e->verts[1] << endl;
        cout << "  t (" << e->tris.size() << ")" << endl;
        for(Tptr t : e->tris)
        cout << "    " << t << endl;
        edge_count++;
    });
    cout << "There were " << edge_count << " EDGES" << endl;
    cout << "VERTS" << endl;
    int vert_count = 0;
    verts.for_each([&](Vptr v) {
        cout << " " << v << ": ref(" << v->ref << ")" << endl;
        cout << "  e (" << v->edges.size() << ")" << endl;
        for(Eptr e : v->edges)
        cout << "    " << e << endl;
        cout << "  t (" << v->tris.size() << ")" << endl;
        for(Tptr t : v->tris)
        cout << "    " << t << endl;
        vert_count++;
    });
    cout << "There were " << vert_count << " VERTS" << endl;
}




