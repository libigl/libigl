// +-------------------------------------------------------------------------
// | mesh.decl.h
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

#include "rawMesh.h"

#include "prelude.h"
#include <vector>
#include <set>

#include "vec.h"
#include "ray.h"
#include "shortVec.h"

#include "iterPool.h"


struct BoolVertexData {
};
struct BoolTriangleData {
    byte bool_alg_data; // internal use by algorithm
        // please copy value when the triangle is subdivided
};


template<class VertData, class TriData>
struct IsctVertEdgeTriInput
{
    VertData*   e[2];
    VertData*   t[3];
};

template<class VertData, class TriData>
struct IsctVertTriTriTriInput
{
    VertData*   t[3][3];
};

template<class VertData, class TriData>
struct SubdivideTriInput
{
    TriData*    pt;
    VertData*   pv[3];
    VertData*   v[3];
};

// in order to perform intersections, VertData and TriData must support
struct IsctVertexData {
    // specify how to compute new data for vertices formed by intersections
    /*
    // vertices on edge and triangle forming an intersection...
    void isct(IsctVertEdgeTriInput input);
    void isct(IsctVertTriTriTriInput input);
    void isctInterpolate(const VertData &v0, const VertData &v1);
    */
};
struct IsctTriangleData {
    // specify how to compute new data for a triangle in the event
    // that it is merged with another triangle (merge)
    // split into two triangles (split)
    // or that the triangle is moved (move)
    /*
    void subdivide(SubdivideTriInput input);
    */
};



// in order to perform remeshing, VertData and TriData must support
struct RemeshVertexData {
    bool manifold; // whether this point is manifold.
                   // useful for modifying interpolation behavior
    // specify how to compute new data for a vertex in the event of
    // either an edge collapse (via merge) or edge split (via interpolate)
    /*
    void merge(const VertData &v0, const VertData &v1);
    void interpolate(const VertData &v0, const VertData &v1);
    */
};
struct RemeshTriangleData {
    // specify how to compute new data for a triangle in the event
    // that it is merged with another triangle (merge)
    // split into two triangles (split)
    // or that the triangle is moved (move)
    /*
    void merge(const TriData &t0, const TriData &t1);
    static void split(TriData &t0, TriData &t1, const TriData &t_orig);
    void move(const TriData &t_old);
    */
};

struct RemeshOptions
{
    double maxEdgeLength;
    double minEdgeLength;
    double minAngle;
    double maxAngle;
    RemeshOptions() :
        maxEdgeLength(1.0),
        minEdgeLength(0.3),
        minAngle(5.0),
        maxAngle(170.0)
    {}
};


// only for internal use, please do not use as client
    struct TopoVert;
    struct TopoEdge;
    struct TopoTri;
    
    typedef TopoVert* Vptr;
    typedef TopoEdge* Eptr;
    typedef TopoTri*  Tptr;
    //using Vptr = TopoVert*;
    //using Eptr = TopoEdge*;
    //using Tptr = TopoTri*;
// end internal items


template<class VertData, class TriData>
class Mesh
{
public:
    Mesh();
    Mesh(Mesh &&src);
    Mesh(const RawMesh<VertData,TriData> &raw);
    virtual ~Mesh();
    
    void operator=(Mesh &&src);
    
    // validity check:
    //  - all numbers are well-defined and finite
    //  - all triangle vertex indices are in the right range
    bool valid() const;
    
    RawMesh<VertData,TriData> raw() const;
    
    inline int numVerts() const { return verts.size(); }
    inline int numTris() const { return tris.size(); }
    
    inline void for_verts(std::function<void(VertData &)> func);
    inline void for_tris(std::function<void(TriData &,
                          VertData &, VertData &, VertData &)> func);
    inline void for_edges(
        std::function<void(VertData &, VertData &)> start,
        std::function<void(TriData &t,
                           VertData &, VertData &, VertData &)> each_tri
    );
    
    // form the disjoint union of two meshes
    void disjointUnion(const Mesh &cp);
    
    struct Isct {
        Ray3d   ray;
        bool    exists;
        
        uint    tri_id;
        Vec3d   isct;
        Vec3d   bary;
    };
    Isct pick(Ray3d ray);
    inline void accessIsct(const Isct &isct,
        std::function<void(TriData &,
            VertData &, VertData &, VertData &)> func);
    
    // checks if the mesh is closed
    bool isClosed();
    
public: // REMESHING module
    // REQUIRES:
    //  - MinimalData
    //  - RemeshData
    void remesh();
    RemeshOptions remesh_options;
    
public: // ISCT (intersections) module
    void resolveIntersections(); // makes all intersections explicit
    bool isSelfIntersecting(); // is the mesh self-intersecting?
    // TESTING
    void testingComputeStaticIsctPoints(std::vector<Vec3d> *points);
    void testingComputeStaticIsct(std::vector<Vec3d> *points,
               std::vector< std::pair<Vec3d,Vec3d> > *edges);
    
public: // BOOLean operation module
    // all of the form
    //      this = this OP rhs
    void boolUnion(Mesh &rhs);
    void boolDiff(Mesh &rhs);
    void boolIsct(Mesh &rhs);
    void boolXor(Mesh &rhs);
    
private:    // Internal Formats
    struct Tri {
        TriData data;
        union {
            struct {
                uint a, b, c; // vertex ids
            };
            uint v[3];
        };
        
        inline Tri() {}
    };
    
    inline void merge_tris(uint tid_result, uint tid0, uint tid1);
    inline void split_tris(uint t0ref, uint t1ref, uint t_orig_ref);
    inline void move_tri(Tri &t_new, Tri &t_old);
    inline void subdivide_tri(uint t_piece_ref, uint t_parent_ref);
    
private:    // DATA
    std::vector<Tri>        tris;
    std::vector<VertData>   verts;
    
private:    // caches
    struct NeighborEntry {
        uint vid;
        ShortVec<uint, 2> tids;
        inline NeighborEntry() {}
        inline NeighborEntry(uint vid_) : vid(vid_) {}
    };
    struct NeighborCache {
        std::vector< ShortVec<NeighborEntry, 8> > skeleton;
        inline NeighborEntry& operator()(uint i, uint j) {
            uint N = skeleton[i].size();
            for(uint k = 0; k < N; k++) {
                if(skeleton[i][k].vid == j)
                    return skeleton[i][k];
            }
            skeleton[i].push_back(NeighborEntry(j));
            return skeleton[i][N];
        }
    };
    NeighborCache createNeighborCache();
    
    // parallel to vertex array
    std::vector<uint> getComponentIds();
    
    // like the neighbor cache, but more customizable
    template<class Edata>
    struct EGraphEntry {
        uint                vid;
        ShortVec<uint, 2>   tids;
        Edata               data;
        inline EGraphEntry() {}
        inline EGraphEntry(uint vid_) : vid(vid_) {}
    };
    template<class Edata>
    struct EGraphCache {
        std::vector< ShortVec<EGraphEntry<Edata>, 8> > skeleton;
        inline EGraphEntry<Edata> & operator()(uint i, uint j) {
            uint N = skeleton[i].size();
            for(uint k = 0; k < N; k++) {
                if(skeleton[i][k].vid == j)
                    return skeleton[i][k];
            }
            skeleton[i].push_back(EGraphEntry<Edata>(j));
            return skeleton[i][N];
        }
        inline void for_each(std::function<void(
                uint i, uint j, EGraphEntry<Edata> &entry
            )> action
        ) {
            for(uint i=0; i<skeleton.size(); i++) {
                for(auto &entry : skeleton[i]) {
                    action(i, entry.vid, entry);
                }
            }
        }
    };
    template<class Edata>
    EGraphCache<Edata> createEGraphCache();
    
    
private:    // TopoCache Support
    struct TopoCache;
private:    // Isct Support
    class  IsctProblem; // implements intersection functionality
        class TriangleProblem; // support type for IsctProblem
        typedef TriangleProblem* Tprob;
        //using Tprob = TriangleProblem*;
private:    // Bool Support
    class BoolProblem;
    
private:    // Remeshing Support
    struct RemeshScratchpad;
    
    Eptr allocateRemeshEdge(RemeshScratchpad &);
    void deallocateRemeshEdge(RemeshScratchpad &, Eptr);
    
    void edgeSplit(RemeshScratchpad &,
                   Eptr e_split);
    void edgeCollapse(RemeshScratchpad &,
                      Eptr e_collapse,
                      bool collapsing_tetrahedra_disappear);
    
    // Need edge scoring routines...
    void scoreAndEnqueue(
        std::set< std::pair<double, Eptr> > &queue, Eptr edge);
    void dequeue(
        std::set< std::pair<double, Eptr> > &queue, Eptr edge);
    double computeEdgeScore(Eptr edge);
    
    // support functions
    void populateTriFromTopoTri(Tptr t);
    // calls the first function once, then the second once for each triangle
    inline void edgeNeighborhood(
        Eptr edge,
        std::function<void(VertData &v0, VertData &v1)> once,
        std::function<void(VertData &v0, VertData &v1,
                           VertData &vopp, TriData &t)> each_tri
    );
};


// inline functions

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::for_verts(
    std::function<void(VertData &v)> func
) {
    for(auto &v : verts)
        func(v);
}

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::for_tris(
    std::function<void(TriData &, VertData &, VertData &, VertData &)> func
) {
    for(auto &tri : tris) {
        auto &a = verts[tri.a];
        auto &b = verts[tri.b];
        auto &c = verts[tri.c];
        func(tri.data, a, b, c);
    }
}

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::for_edges(
    std::function<void(VertData &, VertData &)> start,
    std::function<void(TriData &t,
                       VertData &, VertData &, VertData &)> each_tri
) {
    NeighborCache cache = createNeighborCache();
    for(uint i=0; i<cache.skeleton.size(); i++) {
        for(auto &entry : cache.skeleton[i]) {
            uint j = entry.vid;
            start(verts[i], verts[j]);
            for(uint tid : entry.tids) {
                Tri &tri = tris[tid];
                each_tri(tri.data, verts[tri.a], verts[tri.b], verts[tri.c]);
            }
        }
    }
}

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::accessIsct(
    const Isct &isct,
    std::function<void(TriData &, VertData &, VertData &, VertData &)> func
) {
    Tri &tri = tris[isct.tri_id];
    auto &a = verts[tri.a];
    auto &b = verts[tri.b];
    auto &c = verts[tri.c];
    func(tri.data, a, b, c);
}



