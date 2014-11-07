// +-------------------------------------------------------------------------
// | mesh.remesh.tpp
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

#include <map>
#include <set>

#include "memPool.h"

//#include "files.h"

#include "mesh.topoCache.tpp"



enum EdgeRemeshOperation {
    EDGE_SPLIT,
    EDGE_COLLAPSE,
    EDGE_NOTHING
};

struct RemeshEdgeAuxiliary {
    double                  score;      // priority for queue
    EdgeRemeshOperation     op;         // operation to perform
};

template<class VertData, class TriData>
struct Mesh<VertData, TriData>::RemeshScratchpad {
    TopoCache                               cache;
    MemPool<RemeshEdgeAuxiliary>            edge_data;
    std::set< std::pair<double, Eptr> >     queue;
    
    RemeshScratchpad(Mesh<VertData, TriData> *mesh) : cache(mesh) {}
};


template<class VertData, class TriData>
void Mesh<VertData, TriData>::remesh()
{
    //std::cout << "remesh initial tri count: " << tris.size() << std::endl;
    if(verts.size() == 0)   return; // pathology guard
    
    /*{ // dump out mesh before performing re-mesh
        Files::FileMesh filemesh;
        filemesh = transduce<Files::FileVertex, Files::FileTriangle,
                             VertData, TriData>
        (raw(),
        [](Files::FileVertex &out, const VertData &in) {
            out.pos = in.pos;
        },
        [](Files::FileTriangle &out, const TriData &in) {
            out.a = in.a; out.b = in.b; out.c = in.c;
        });
        if(Files::writeTriMesh("remesh_in.off", &filemesh) > 0) {
            CORK_ERROR("Unable to write out remesh input to remesh_in.off");
        }
    }*/
    
    // create a scratchpad and set it up
    RemeshScratchpad scratchpad(this);
    
    // compute which vertices are boundary or not
    for(VertData &v : verts) {
        v.manifold = true;
    }
    scratchpad.cache.edges.for_each([this, &scratchpad](Eptr edge) {
        if(edge->tris.size() != 2) {
            verts[edge->verts[0]->ref].manifold = false;
            verts[edge->verts[1]->ref].manifold = false;
        }
    });
    
    // Then, we set up the priority queue
    // (allocating auxiliary data as we go)
    scratchpad.cache.edges.for_each([this, &scratchpad](Eptr edge) {
        edge->data = scratchpad.edge_data.alloc();
        scoreAndEnqueue(scratchpad.queue, edge);
    });
    
    int cutoff = 10000;
    // now let's go into a loop pulling work off of the queue
    while(scratchpad.queue.size() > 0) {
        // pop the end of the queue
        auto        it_top          = scratchpad.queue.end();
                    it_top--;
        //double      top_score       = it_top->first;
        Eptr        top_edge        = it_top->second;
                    scratchpad.queue.erase(it_top);
        
        // Which operation should be performed to this edge?
        EdgeRemeshOperation         op =
                reinterpret_cast<RemeshEdgeAuxiliary*>(top_edge->data)->op;
        if(op == EDGE_SPLIT) {
                    //std::cout << "edge split" << std::endl;
                    edgeSplit(scratchpad, top_edge);
        } else if (op == EDGE_COLLAPSE) {
                    //std::cout << "edge collapse" << std::endl;
                    edgeCollapse(scratchpad, top_edge, true);
        } // else do nothing (should have score <= 0.0 then...)
        cutoff--;
        if(cutoff < 0)
            break;
    }
    
    //std::cout << " cache.tris is "
    //          << scratchpad.cache.tris.size() << std::endl;
    //std::cout << "prefinalize tri count: " << tris.size() << std::endl;
    
    // finally, we commit all of the results of this remeshing operation
    // causing the data storage to defrag and clean out dead items
    scratchpad.cache.commit();
    
    //std::cout << "remesh final tri count: " << tris.size() << std::endl;
}



template<class VertData, class TriData>
void Mesh<VertData, TriData>::scoreAndEnqueue(
    std::set< std::pair<double, Eptr> >     &queue,
    Eptr                                     edge
) {
    double score = computeEdgeScore(edge);
    if(score > 0.0) // only enqueue edges with actual work to do
        queue.insert(std::make_pair(score, edge));
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::dequeue(
    std::set< std::pair<double, Eptr> >     &queue,
    Eptr                                     edge
) {
    double score =
            reinterpret_cast<RemeshEdgeAuxiliary*>(edge->data)->score;
    auto it = queue.find(std::make_pair(score, edge));
    if(it != queue.end())
        queue.erase(it);
}



template<class VertData, class TriData>
inline void Mesh<VertData, TriData>::edgeNeighborhood(
    Eptr edge,
    std::function<void(VertData &v0, VertData &v1)> once,
    std::function<void(VertData &v0, VertData &v1,
                       VertData &vopp, TriData &t)> each_tri
) {
    Vptr                    v0          = edge->verts[0];
    Vptr                    v1          = edge->verts[1];
    VertData                &data0      = verts[v0->ref];
    VertData                &data1      = verts[v1->ref];
    
    once(data0, data1);
    
    for(Tptr tri : edge->tris) {
        TriData             &tdata      = tris[tri->ref].data;
        for(uint i=0; i<3; i++) {
            if(tri->edges[i] == edge) {
                Vptr        vopp        = tri->verts[i];
                VertData    &dataopp    = verts[vopp->ref];
                each_tri(data0, data1, dataopp, tdata);
            }
        }
    }
}

template<class VertData, class TriData>
double Mesh<VertData, TriData>::computeEdgeScore(Eptr edge) {
    double edge_length;
    double min_angle = 360.0; // clearly overkill
    double max_angle = 0.0;
    edgeNeighborhood(edge,
    [&edge_length]
    (VertData &v0, VertData &v1) { // once
        edge_length = len(v1.pos - v0.pos);
    },
    [&min_angle, &max_angle]
    (VertData &v0, VertData &v1, VertData &vopp, TriData &) { // per tri
        Vec3d e0 = v0.pos - vopp.pos;
        Vec3d e1 = v1.pos - vopp.pos;
        double cos_angle = dot(e0, e1) / (len(e0)*len(e1));
        double angle = rad2deg(acos(cos_angle));
        if(angle < min_angle) min_angle = angle;
        if(angle > max_angle) max_angle = angle;
        // CURRENTLY OMITTING DIHEDRAL ANGLE CRITERION
    });
    
    RemeshEdgeAuxiliary *edge_aux =
                reinterpret_cast<RemeshEdgeAuxiliary*>(edge->data);
    
    // extract violation quantities:
    //      if any of these are positive,
    //      then an operation should be performed to this edge
    double excess_length        = edge_length - remesh_options.maxEdgeLength;
    double deficient_length     = remesh_options.minEdgeLength - edge_length;
    double excess_angle         = max_angle - remesh_options.maxAngle;
    double deficient_angle      = remesh_options.minAngle - max_angle;
    if(excess_length <= 0.0 && deficient_length <= 0.0 &&
       excess_angle  <= 0.0 && deficient_angle  <= 0.0) {
        // no op on edge
        edge_aux->score = -1.0;
        edge_aux->op = EDGE_NOTHING;
        return edge_aux->score;
    }
    
    // rescale violation quantities:
    //      make quantities more commensurate for use as priorities 
    excess_length       /= remesh_options.maxEdgeLength;
    deficient_length    /= remesh_options.minEdgeLength;
    excess_angle        /= (180.0 - remesh_options.maxAngle);
    deficient_angle     /= remesh_options.minAngle;
    edge_aux->score = -1.0;
    if(excess_length > edge_aux->score) {
        edge_aux->score = excess_length;
        edge_aux->op = EDGE_SPLIT;
    }
    if(deficient_length > edge_aux->score) {
        edge_aux->score = deficient_length;
        edge_aux->op = EDGE_COLLAPSE;
    }
    if(excess_angle > edge_aux->score) {
        edge_aux->score = excess_angle;
        edge_aux->op = EDGE_SPLIT;
    }
    if(deficient_angle > edge_aux->score) {
        edge_aux->score = deficient_angle;
        edge_aux->op = EDGE_COLLAPSE;
    }
    
    return edge_aux->score;
}



// erase the first occurrence of a value from a short vector,
// !!! KNOWING that it MUST BE PRESENT
template<class T, uint LEN> inline
void remove(ShortVec<T,LEN> &vec, T val)
{
    uint last_i = vec.size()-1;
    for(uint i=0; i<last_i; i++) {
        if(vec[i] == val) {
            std::swap(vec[i], vec[last_i]);
            break;
        }
    }
    vec.resize(last_i);
}

// erase the first occurrence of a value,
// IF THE VALUE OCCURS
template<class T, uint LEN> inline
void maybe_erase(ShortVec<T,LEN> &vec, T val)
{
    uint N = vec.size();
    for(uint i=0; i<N; i++) {
        if(vec[i] == val) {
            std::swap(vec[i], vec[N-1]);
            vec.resize(N-1);
            break;
        }
    }
}



template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::populateTriFromTopoTri(Tptr tri) {
    Tri &tri_ref = tris[tri->ref];
    for(uint k=0; k<3; k++)
        tri_ref.v[k] = tri->verts[k]->ref;
}

template<class VertData, class TriData>
Eptr Mesh<VertData, TriData>::allocateRemeshEdge(
    RemeshScratchpad &scratchpad
) {
    Eptr        e           = scratchpad.cache.newEdge();
                e->data     = scratchpad.edge_data.alloc();
    return      e;
}

template<class VertData, class TriData>
void Mesh<VertData, TriData>::deallocateRemeshEdge(
    RemeshScratchpad &scratchpad,
    Eptr e
) {
                dequeue(scratchpad.queue, e);
                scratchpad.cache.freeEdge(e);
}


template<class VertData, class TriData>
inline void Mesh<VertData, TriData>::merge_tris(
    uint tid_result, uint tid0, uint tid1
) {
    tris[tid_result].data.merge(tris[tid0].data, tris[tid1].data);
}
template<class VertData, class TriData>
inline void Mesh<VertData, TriData>::split_tris(
    uint t0ref, uint t1ref, uint t_orig_ref
) {
    TriData::split(tris[t0ref].data, tris[t1ref].data, tris[t_orig_ref].data);
}
template<class VertData, class TriData>
inline void Mesh<VertData, TriData>::move_tri(Tri &t_new, Tri &t_old)
{
    t_new.data.move(t_old.data);
}






/*
    Strategy:
     1. Construct new geometry in parallel to old geometry w/o
            tampering with old geometry in any way.
     2. Fill out data for new geometry
     3. Determine whether or not to commit or abort the operation
     4. Connect new geometry to old geometry simultaneously to old geometry.
     5. Delete the old geometry appropriately
*/

struct EdgeWedge {
    Eptr e0, e1;
    
    EdgeWedge() : e0(nullptr), e1(nullptr) {}
    bool full() const { return e0 && e1; }
    Eptr one()  const { return (e0)? e0 : e1; }
};

struct TriWedge {
    Tptr t0, t1;
    
    TriWedge() : t0(nullptr), t1(nullptr) {}
    bool full() const { return t0 && t1; }
    Tptr one()  const { return (t0)? t0 : t1; }
};

struct VptrRemap {
    Vptr v0, v1, v_merged;
    VptrRemap(Vptr v0_, Vptr v1_, Vptr v_merged_) :
        v0(v0_), v1(v1_), v_merged(v_merged_)
    {}
    Vptr operator[](Vptr key) {
        if(key == v0 || key == v1)          return v_merged;
        else                                return key;
    }
};

template<class T>
struct PtrRemap {
    std::map<T*, T*> remap;
    
    void set(T* key, T* val) {
        remap[key] = val;
    }
    
    T* operator[](T* key) {
        auto it = remap.find(key);
        if(it == remap.end())       return key;
        else                        return it->second;
    }
};

template<class VertData, class TriData>
void Mesh<VertData, TriData>::edgeCollapse(
    RemeshScratchpad &scratchpad,
    Eptr e_collapse,
    bool collapsing_tetrahedra_disappear
) {
    /*
     *
     *  Two cases: both with and without a triangle filling the
     *      triangular arrangement (i.e. wedge) of edges
     *
     *                     *   
     *                   /   \
     *                 / # # # \
     *               /  # # # #  \  
     *             / # # # # # # # \
     *           /  # # # # # # # #  \
     *         / # # # # # # # # # # # \
     *  vid0 *---------------------------* vid1
     *                eid_collapse
     *
     *                     *   
     *                   /   \
     *                 /       \
     *               /           \  
     *             /               \
     *           /                   \
     *         /                       \
     *  vid0 *---------------------------* vid1
     *                eid_collapse
     *
     *
     *  Additionally, we need to worry about how tetrahedral structures
     *  can collapse.  If the collapse will identify the two
     *  non-collapsing faces, then we call this structure a triangle wedge
     *
     *
     *                     * 
     *                    /|\
     *                   / | \
     *                  /  |  \
     *                 / # | # \
     *                / #  |  # \
     *               /   # | # # \
     *              / # #  *  # # \
     *             /  #  /   \  #  \
     *            / #  /       \  # \
     *           /   /           \   \ 
     *          /  /               \  \
     *         / /                   \ \
     *        //                       \\
     *  vid0 *---------------------------* vid1
     *                 eid_collapse
     *
     *
     *
        In summary, there are three features we should be interested in:
        --  hollow edge wedges that will be collapsed
        --  filled edge wedges that will be collapsed
        --  triangle wedges that will be collapsed
        
        We can identify all edge wedges (hollow or filled) via:
            JOIN FILTER(EDGES(A,X), X!=B) WITH
                 FILTER(EDGES(B,Y), Y!=A)
                 WHERE X=Y
        We can identify all triangle wedges via:
            JOIN FILTER(TRIANGLES(A,X,Y), X!=B && Y!=B) WITH
                 FILTER(TRIANGLES(B,Z,W), Z!=A && W!=A)
                 WHERE (X,Y)=(Z,W) or (X,Y)=(W,Z)
     */
    
    /* Identify all the following components
     * and separate them according to type:
     *  -   e_collapse:     the collapsing edge
     *  -   tris_collapse:  the collapsing triangles
     *  -   v0, v1:         the merging vertices
     *  -   edge_wedges:    the merging edges
     *  -   tri_wedges:     the merging (or disappearing) triangles
     *  -   edges_moving:   the persisting, but moving edges
     *  -   tris_moving:    the persisting, but moving triangles
     */
    ShortVec<Tptr, 2>       tris_collapse   = e_collapse->tris;
    
    Vptr                    v0              = e_collapse->verts[0];
    Vptr                    v1              = e_collapse->verts[1];
    
    ShortVec<EdgeWedge, 2>  edge_wedges;
    ShortVec<Eptr, 10>      edges_moving;
    
    ShortVec<TriWedge, 2>   tri_wedges;
    ShortVec<Tptr, 14>      tris_moving;
    
    // BUILD all the edge wedges
    std::map<Vptr, EdgeWedge>   edge_w_map;
    for(Eptr e : v0->edges) {
        if(e == e_collapse) continue;
        Vptr                ev0                 = e->verts[0];
        Vptr                ev1                 = e->verts[1];
        Vptr                key                 = (ev0 != v0)? ev0 : ev1;
                            edge_w_map[key].e0  = e;
    }
    for(Eptr e : v1->edges) {
        if(e == e_collapse) continue;
        Vptr                ev0                 = e->verts[0];
        Vptr                ev1                 = e->verts[1];
        Vptr                key                 = (ev0 != v1)? ev0 : ev1;
                            edge_w_map[key].e1  = e;
    }
    for(const auto &pair : edge_w_map) {
        if(pair.second.full())
                            edge_wedges.push_back(pair.second);
        else
                            edges_moving.push_back(pair.second.one());
    }
    
    // BUILD all the triangle wedges
    std::map<Eptr, TriWedge>    tri_w_map;
    for(Tptr t : v0->tris) {
        if(t->edges[0] == e_collapse ||
           t->edges[1] == e_collapse ||
           t->edges[2] == e_collapse)   continue;
        
        Eptr                key                 = nullptr;
        for(uint k=0; k<3; k++)
            if(t->verts[k] == v0)
                            key                 = t->edges[k];
        ENSURE(key != NULL);
                            tri_w_map[key].t0   = t;
    }
    for(Tptr t : v1->tris) {
        if(t->edges[0] == e_collapse ||
           t->edges[1] == e_collapse ||
           t->edges[2] == e_collapse)   continue;
        
        Eptr                key                 = nullptr;
        for(uint k=0; k<3; k++)
            if(t->verts[k] == v1)
                            key                 = t->edges[k];
        ENSURE(key != NULL);
                            tri_w_map[key].t1   = t;
    }
    for(const auto &pair : tri_w_map) {
        if(pair.second.full())
                            tri_wedges.push_back(pair.second);
        else
                            tris_moving.push_back(pair.second.one());
    }
    
    /* Next, we need to create the new geometry which will supplant
     * the just enumerated parts:
     *  -   v_merged:       the merged vertex
     *  -   edges_merged:   the merged edges (parallel to edge_wedges)
     *  -   tris_merged:    the merged triangles (parallel to tri_wedges)
     *  -   edges_moved:    the moved edges (parallel to edges_moving)
     *  -   tris_moved:     the moved triangles (parallel to tris_moving)
     */
    /* As part of this, we will build a record of the remapping to
     * aid us when we start to connect all of this geometry
     */
    
    Vptr                    v_merged            = scratchpad.cache.newVert();
    ShortVec<Eptr, 2>       edges_merged(edge_wedges.size());
    ShortVec<Eptr, 10>      edges_moved(edges_moving.size());
    ShortVec<Tptr, 2>       tris_merged(tri_wedges.size());
    ShortVec<Tptr, 14>      tris_moved(tris_moving.size());
    
    VptrRemap               vptr_remap(v0, v1, v_merged);
    PtrRemap<TopoEdge>    eptr_remap;
    PtrRemap<TopoTri>     tptr_remap;
    
    // Create new edges and enter re-mappings
    eptr_remap.set(e_collapse, nullptr); // mark this edge as dying
    for(uint i=0; i<edge_wedges.size(); i++) {
                            edges_merged[i] = allocateRemeshEdge(scratchpad);
                            eptr_remap.set(edge_wedges[i].e0, edges_merged[i]);
                            eptr_remap.set(edge_wedges[i].e1, edges_merged[i]);
    }
    for(uint i=0; i<edges_moving.size(); i++) {
                            edges_moved[i] = allocateRemeshEdge(scratchpad);
                            eptr_remap.set(edges_moving[i], edges_moved[i]);
    }
    
    // Create new triangles and enter re-mappings
    for(Tptr t : tris_collapse)
                            // mark this triangle as dying
                            tptr_remap.set(t, nullptr); 
    for(uint i=0; i<tri_wedges.size(); i++) {
        if(collapsing_tetrahedra_disappear) {
                            tris_merged[i]      = nullptr;
        } else {
                            tris_merged[i]      = scratchpad.cache.newTri();
        }
                            tptr_remap.set(tri_wedges[i].t0, tris_merged[i]);
                            tptr_remap.set(tri_wedges[i].t1, tris_merged[i]);
    }
    for(uint i=0; i<tris_moving.size(); i++) {
                            tris_moved[i]       = scratchpad.cache.newTri();
                            tptr_remap.set(tris_moving[i], tris_moved[i]);
    }
    
    /* Now we want to connect up all of this new geometry to each other
     * and to the existing geometry,
     *  BUT!!!
     * we also want to be careful not to point any existing geometry
     * at the new pieces.  Before doing that, we want to be able to
     * compute data for the new geometry and confirm or deny that we
     * want to commit this operation. (e.g. check for unwanted collisions)
     *
     * We adopt the following strategy:
     *  -   first point the new triangles at edges and vertices
     *  -   next, point the new edges at vertices; and 
     *          point the new edges at new triangles.
     *  -   finally, point the new vertex at
     *  -       the new edges and new triangles
     *
     * If an edge cannot find any valid triangles it is incident to,
     *  then it must be marked for deletion, etc.
     * If the merged vertex cannot find any valid tris/edges it is incident to,
     *  then it must be marked for deletion too!
     */
    
    // First, the triangles
    if(!collapsing_tetrahedra_disappear) {
        for(uint i=0; i<tri_wedges.size(); i++) {
            Tptr            t0                  = tri_wedges[i].t0;
            //Tptr            t1                  = tri_wedges[i].t1;
            Tptr            t_new               = tris_merged[i];
            
            for(uint k=0; k<3; k++) {
                            t_new->verts[k]     = vptr_remap[t0->verts[k]];
                            t_new->edges[k]     = eptr_remap[t0->edges[k]];
            }
                            populateTriFromTopoTri(t_new);
        }
    }
    for(uint i=0; i<tris_moving.size(); i++) {
        Tptr                t_old               = tris_moving[i];
        Tptr                t_new               = tris_moved[i];
        
        for(uint k=0; k<3; k++) {
                            t_new->verts[k]     = vptr_remap[t_old->verts[k]];
                            t_new->edges[k]     = eptr_remap[t_old->edges[k]];
        }
                            populateTriFromTopoTri(t_new);
    }
    
    // Next, the edges
    for(uint i=0; i<edge_wedges.size(); i++) {
        Eptr                e0                  = edge_wedges[i].e0;
        Eptr                e1                  = edge_wedges[i].e1;
        Eptr                e_new               = edges_merged[i];
        
        // plug in all the valid triangles...
        for(Tptr t : e0->tris) {
                            t                   = tptr_remap[t];
            if(t)           e_new->tris.push_back(t);
        }
        for(Tptr t : e1->tris) {
                            t                   = tptr_remap[t];
            if(t)           e_new->tris.push_back(t);
        }
        
        if(e_new->tris.size() == 0) { // if there are no parent triangles left
            // then we need to kill this edge
                            eptr_remap.set(e0, nullptr);
                            eptr_remap.set(e1, nullptr);
                            deallocateRemeshEdge(scratchpad, e_new);
                            edges_merged[i]     = nullptr;
        }
        else { // otherwise, let's go ahead and finish hooking up this edge
            for(uint k=0; k<2; k++)
                            e_new->verts[k]     = vptr_remap[e0->verts[k]];
        }
    }
    for(uint i=0; i<edges_moving.size(); i++) {
        Eptr                e_old               = edges_moving[i];
        Eptr                e_new               = edges_moved[i];
        
        // note: should never have any dead/null triangles
        for(Tptr t : e_old->tris) {
                            t                   = tptr_remap[t];
            ENSURE(t);
                            e_new->tris.push_back(t);
        }
        for(uint k=0; k<2; k++)
                            e_new->verts[k]     = vptr_remap[e_old->verts[k]];
    }
    
    // Finally, the vertex
    {
        // Should do this directly, not via the remap translation.
        // Working via re-maps will lead to duplicates of merged geometry.
        // However, we can exploit the fact that this vertex is unique,
        //    and that we already have lists of all the incident geometry
        
        for(Tptr t : tris_merged)
            if(t)           v_merged->tris.push_back(t);
        for(Tptr t : tris_moved) // cannot be dead
                            v_merged->tris.push_back(t);
        if(v_merged->tris.size() == 0) {
                            scratchpad.cache.freeVert(v_merged);
                            v_merged            = nullptr;
        } else {
            for(Eptr e : edges_merged)
                if(e)           v_merged->edges.push_back(e);
            for(Eptr e : edges_moved) // cannot be dead
                                v_merged->edges.push_back(e);
            // it's impossible to have triangles incident w/o edges too.
            ENSURE(v_merged->edges.size() > 0);
        }
    }
    
    
    // OK, here we get to finally compute data for all the new geometry
    // Once we've done that, we can also check to see whether we actually
    // want to commit this operation or not.
    
    // merge vertices' data
    if(v_merged) { // REMEMBER: vertex could be deleted by now
        VertData            &data_new           = verts[v_merged->ref];
        const VertData      &data0              = verts[v0->ref];
        const VertData      &data1              = verts[v1->ref];
                            data_new.merge(data0, data1);
    }
    // merge triangles' data
    for(uint i=0; i<tri_wedges.size(); i++) {
        Tptr                t_new               = tris_merged[i];
        // REMEMBER: these triangles could be deleted
        if(!t_new)          continue;
        
        Tptr                t0                  = tri_wedges[i].t0;
        Tptr                t1                  = tri_wedges[i].t1;
        
                            merge_tris(t_new->ref, t0->ref, t1->ref);
    }
    // update moved triangles' data
    for(uint i=0; i<tris_moving.size(); i++) {
        // NOTE: moved triangles cannot be deleted
        Tptr                t_new               = tris_moved[i];
        Tri                 &tri_new            = tris[t_new->ref];
        
        Tptr                t_old               = tris_moving[i];
        Tri                 &tri_old            = tris[t_old->ref];
        
                            move_tri(tri_new, tri_old);
    }
    
    // TODO: COMMIT OPTION will be left as a stub for now!
    if(false) {
        // DESTROY THE GEOMETRY WE JUST CREATED
        // AND RETURN
    }
    
    // Find and Store all of the existing, unchanged edges
    // that border this operation.  We will need these references
    // when we account for changes to edge operation priorities later
    ShortVec<Eptr, 16>      borderEdges;
    for(const TriWedge &e_wedge : tri_wedges) {
        Tptr                t0                  = e_wedge.t0; 
        for(uint k=0; k<3; k++) {
            if(t0->verts[k] == v0) {
                            borderEdges.push_back(t0->edges[k]);
                            break;
            }
        }
    }
    for(Tptr t : tris_moved) {
        for(uint k=0; k<3; k++) {
            if(t->verts[k] == v_merged) {
                            borderEdges.push_back(t->edges[k]);
                            break;
            }
        }
    }
    
    /* Now that we've got the go ahead, let's hook in the new geometry to
     *  the existing geometry!
     * We can do this in the following order:
     *  -   take all the new edges, and add them to their existing endpoint's
     *      edge list.  (all new edges must have exactly one such endpoint)
     *  -   take all the new triangles, and add them to their two existing
     *      endpoints' and one existing edge's triangle lists.
     */
    // first, consolidate arrays
    ShortVec<Eptr, 12>      new_edges;
    ShortVec<Tptr, 16>      new_tris;
    {
        for(Eptr e : edges_merged)
            if(e)           new_edges.push_back(e);
        for(Eptr e : edges_moved) // cannot be deleted
                            new_edges.push_back(e);
        for(Tptr t : tris_merged)
            if(t)           new_tris.push_back(t);
        for(Tptr t : tris_moved) // cannot be deleted
                            new_tris.push_back(t);
    }
    // hook up edges to existing geometry
    for(Eptr edge : new_edges) {
        Vptr                ev0                 = edge->verts[0];
        Vptr                ev1                 = edge->verts[1];
        Vptr                v_old               = (ev0 != v_merged)? ev0 : ev1;
                            v_old->edges.push_back(edge);
    }
    // hook up triangles to existing geometry
    for(Tptr tri : new_tris) {
        for(uint k=0; k<3; k++) {
            if(tri->verts[k] != v_merged)       continue;
            Eptr            e                   = tri->edges[k];
            Vptr            tv0                 = e->verts[0];
            Vptr            tv1                 = e->verts[1];
                            e->tris.push_back(tri);
                            tv0->tris.push_back(tri);
                            tv1->tris.push_back(tri);
                            break;
        }
    }
    
    /* Now, let's kill all the old geometry.  We can use the checklist
     * we built at the begining:
     *  -   e_collapse:     the collapsing edge
     *  -   tris_collapse:  the collapsing triangles
     *  -   v0, v1:         the merging vertices
     *  -   edge_wedges:    the merging edges
     *  -   tri_wedges:     the merging (or disappearing) triangles
     *  -   edges_moving:   the persisting, but moving edges
     *  -   tris_moving:    the persisting, but moving triangles
     *
     * We need to be careful to free this geometry in top down order;
     * starting with the triangles and moving towards the vertices.
     * If we furthermore guarantee that any singular edges or vertices
     * created by a triangle deletion are also deleted, then we can
     * focus all of our attention on just deleting triangles 
     */
    ShortVec<Tptr, 16>      dead_tris;
    ShortVec<Eptr, 16>      dead_edges;
    ShortVec<Vptr, 2>       dead_verts;
    
    // assemble the list of triangles to kill
    for(Tptr t : tris_collapse)
                            dead_tris.push_back(t);
    for(const auto &t_wedge : tri_wedges) {
                            dead_tris.push_back(t_wedge.t0);
                            dead_tris.push_back(t_wedge.t1);
    }
    for(Tptr t : tris_moving)
                            dead_tris.push_back(t);
    
    // process the list of triangles
    for(Tptr tri : dead_tris) {
        // Let's unhook this triangle from its faces first
        for(uint k=0; k<3; k++) {
            Vptr            v                   = tri->verts[k];
                            remove(v->tris, tri);
            if(v->tris.size() == 0)
                            dead_verts.push_back(v);
            
            Eptr            e                   = tri->edges[k];
                            remove(e->tris, tri);
            if(e->tris.size() == 0)
                            dead_edges.push_back(e);
        }
        // now that we're disconnected, go ahead and jettison the triangle
                            scratchpad.cache.freeTri(tri);
    }
    
    // now, we can process the list of edges
    for(Eptr edge : dead_edges) {
        // Let's unhook this edge from its vertices
        for(uint k=0; k<2; k++) {
            Vptr            v                   = edge->verts[k];
                            remove(v->edges, edge);
            // the triangle removal was enough to
            // determine which vertices should die.
            // re-adding them here would lead to duplicates
        }
        // and then jetisson the edge
                            deallocateRemeshEdge(scratchpad, edge);
        // If this edge is in the border edge list,
        // then we need to remove it right away!
                            maybe_erase(borderEdges, edge);
    }
    
    // Finally, polish off by getting rid of any vertices that talked too much
    for(Vptr vert : dead_verts) {
                            scratchpad.cache.freeVert(vert);
                            if(vert == v_merged)    v_merged = nullptr;
    }
    
    // We pause a moment here to update the manifoldness of any
    // vertices for which it might have changed
    // ONLY do if the merged vertex is still alive...
    if(v_merged) {
        verts[v_merged->ref].manifold = true;
        for(Eptr e : v_merged->edges) {
            if(e->tris.size() != 2)
                verts[v_merged->ref].manifold = false;
            
            // process neighboring point
            Vptr v = e->verts[0];
            if(v == v_merged)   v = e->verts[1];
            verts[v->ref].manifold = true;
            for(Eptr ee : v->edges) {
                if(ee->tris.size() != 2) {
                    verts[v->ref].manifold = false;
                    break;
                }
            }
        }
    }
    
    // Before we're completely done, we will go through and
    // adjust priorities for edges which might have been effected by this op.
    // Only explicitly dequeue pre-existing edges we did not delete!
    for(Eptr e : edges_merged) { if(e) { // might be deleted
                            scoreAndEnqueue(scratchpad.queue, e);
    }}
    for(Eptr e : edges_moved) { // def. not deleted
                            scoreAndEnqueue(scratchpad.queue, e);
    }
    for(Eptr e : borderEdges) {
            // border edges could have been deleted...
                            dequeue(scratchpad.queue, e);
                            scoreAndEnqueue(scratchpad.queue, e);
    }
    
    // that should more or less complete an edge collapse
}


/*
    Strategy:
     1. Construct new geometry in parallel to old geometry w/o
            tampering with old geometry in any way.
     2. Fill out data for new geometry
     3. Determine whether or not to commit or abort the operation
     4. Connect new geometry to old geometry simultaneously to old geometry.
     5. Delete the old geometry appropriately
*/


template<class VertData, class TriData>
void Mesh<VertData,TriData>::edgeSplit(
    RemeshScratchpad &scratchpad,
    Eptr e_split
) {
    /*
        Clean-up picture please...
            
     *                     *  vs_opp
     *                   /   \
     *                 /       \
     *       e0s     /           \     e1s  
     *             /               \             
     *           /      ts_orig      \           
     *         /                       \         
     *    v0 *---------------------------* v1  
     *                  e_split                    
     *                     _
     *                     |
     *                     v
     *
     *                     *   vs_opp
     *                   / | \
     *                 /   |   \
     *        e0s    /     | <----------- es_mid
     *             /       |       \
     *           /         |         \  
     *         /  t0s_new  |  t1s_new  \
     *    v0 *-------------*-------------* v1
     *           e0_new  v_new  e1_new
     *
        
        BEGIN
            new vertex
                -- invoke interpolation callback
        FOR EACH TRIANGLE SUB-PROBLEM:
            new next triangle (from split)
            new prev triangle (from split)
                -- invoke triangle split callback
            delete triangle
            FIXUP all of the Tri-Edges
    */
    
    /* Identify all the following components
     * and separate them according to type:
     *  -   e_split:        the edge to be split
     *  -   ts_orig:        the triangles to be split
     *  -   v0, v1:         the two endpoint vertices of the edge
     *  (the following is helper data; do not delete these vertices)
     *  -   vs_opp:         the vertices opposite eid_split for each triangle
     */
    
    ShortVec<Tptr, 2>       ts_orig             = e_split->tris;
    
    Vptr                    v0                  = e_split->verts[0];
    Vptr                    v1                  = e_split->verts[1];
    
    ShortVec<Vptr, 2>       vs_opp(ts_orig.size());
    for(uint i=0; i<ts_orig.size(); i++) {
        Tptr                t_orig              = ts_orig[i];
        for(uint k=0; k<3; k++) {
            if(t_orig->edges[k] == e_split) {
                            vs_opp[i]           = t_orig->verts[k];
                            break;
            }
        }
    }
    
    /* Next, we need to create the new geometry which will supplant
     * the just enumerated parts:
     *  -   v_new:          the vertex introduced by the split
     *  -   e0_new,
     *      e1_new:         the two pieces of e_split
     *  -   t0s_new,
     *      t1s_new:        the two pieces of each triangle in tids_orig
     *  -   es_mid:         the new edges splitting each triangle
     */
    Vptr                    v_new       = scratchpad.cache.newVert();
    Eptr                    e0_new      = allocateRemeshEdge(scratchpad);
    Eptr                    e1_new      = allocateRemeshEdge(scratchpad);
    ShortVec<Tptr, 2>       t0s_new(ts_orig.size());
    ShortVec<Tptr, 2>       t1s_new(ts_orig.size());
    ShortVec<Eptr, 2>       es_mid(ts_orig.size());
    
    for(uint i=0; i<ts_orig.size(); i++) {
                            t0s_new[i]  = scratchpad.cache.newTri();
                            t1s_new[i]  = scratchpad.cache.newTri();
                            es_mid[i]   = allocateRemeshEdge(scratchpad);
    }
    
    /* Now we want to connect up all of this new geometry to each other
     * and to the existing geometry,
     *  BUT!!!
     * we also want to be careful not to point any existing geometry
     * at the new pieces.  Before doing that, we want to be able to
     * compute data for the new geometry and confirm or deny that we
     * want to commit this operation. (e.g. check for unwanted collisions)
     *
     * We adopt the following strategy:
     *  -   first point the new triangles at edges and vertices
     *  -   next, point the new edges at vertices; and 
     *          point the new edges at new triangles.
     *  -   finally, point the new vertex at
     *  -       the new edges and new triangles
     */
    
    // hook up t0s_new and t1s_new
    // also go ahead and hook up es_mid
    for(uint i=0; i<ts_orig.size(); i++) {
        Tptr                t_orig              = ts_orig[i];
        Tptr                t0                  = t0s_new[i];
        Tptr                t1                  = t1s_new[i];
        Eptr                e_mid               = es_mid[i];
        
        // replace every edge and vertex appropriately for the two variants
        for(uint k=0; k<3; k++) {
            Vptr            v_orig              = t_orig->verts[k];
            Eptr            e_orig              = t_orig->edges[k];
            if(v_orig == v0) {
                            t0->verts[k]        = v_orig;
                            t0->edges[k]        = e_mid;
                            t1->verts[k]        = v_new;
                            t1->edges[k]        = e_orig;
            } else if(v_orig == v1) {
                            t0->verts[k]        = v_new;
                            t0->edges[k]        = e_orig;
                            t1->verts[k]        = v_orig;
                            t1->edges[k]        = e_mid;
            } else {
                            t0->verts[k]        = v_orig;
                            t0->edges[k]        = e0_new;
                            t1->verts[k]        = v_orig;
                            t1->edges[k]        = e1_new;
            }
        }
                            populateTriFromTopoTri(t0);
                            populateTriFromTopoTri(t1);
        // set up the mid edge from the split
                            e_mid->verts[0]     = v_new;
                            e_mid->verts[1]     = vs_opp[i];
                            e_mid->tris.resize(2);
                            e_mid->tris[0]      = t0;
                            e_mid->tris[1]      = t1;
    }
    // hook up e0_new and e1_new
                            e0_new->verts[0]    = v0;
                            e0_new->verts[1]    = v_new;
                            e1_new->verts[0]    = v_new;
                            e1_new->verts[1]    = v1;
    for(uint i=0; i<ts_orig.size(); i++) {
                            e0_new->tris.push_back(t0s_new[i]);
                            e1_new->tris.push_back(t1s_new[i]);
    }
    // hook up v_new
                            v_new->edges.push_back(e0_new);
                            v_new->edges.push_back(e1_new);
    for(uint i=0; i<ts_orig.size(); i++) {
                            v_new->edges.push_back(es_mid[i]);
                            v_new->tris.push_back(t0s_new[i]);
                            v_new->tris.push_back(t1s_new[i]);
    }
    
    
    // OK, here we get to finally compute data for all the new geometry
    // Once we've done that, we can also check to see whether we actually
    // want to commit this operation or not.
    
    // interpolate data onto the new vertex
    {
        VertData            &data_new           = verts[v_new->ref];
        const VertData      &data0              = verts[v0->ref];
        const VertData      &data1              = verts[v1->ref];
                            data_new.interpolate(data0, data1);
                            data_new.manifold   = (ts_orig.size() == 2);
    }
    
    // split triangles' data
    for(uint i=0; i<ts_orig.size(); i++) {
        Tptr                t_orig              = ts_orig[i];
        
        Tptr                t0                  = t0s_new[i];
        Tptr                t1                  = t1s_new[i];
                            
                            split_tris(t0->ref, t1->ref, t_orig->ref);
    }
    
    // TODO: COMMIT OPTION will be left as a stub for now!
    if(false) {
        // DESTROY THE GEOMETRY WE JUST CREATED AND RETURN
    }
    
    // record border edges for later priority updates
    ShortVec<Eptr, 4>       borderEdges;
    for(Tptr t : ts_orig) { // TODO: evacuate all of this to the end...
        for(uint k=0; k<3; k++) {
            if(t->edges[k] == e_split)  continue;
                            borderEdges.push_back(t->edges[k]);
        }
    }
    
    /* Now that we've got the go ahead, let's hook in the new geometry to
     *  the existing geometry!
     * We can do this in the following order:
     *  -   take all the new edges, and add them to their existing endpoint's
     *      edge list.  (all new edges must have exactly one such endpoint)
     *  -   take all the new triangles, and add them to their two existing
     *      endpoints' and one existing edge's triangle lists.
     */
    
    // add new edges to v0 and v1
                            v0->edges.push_back(e0_new);
                            v1->edges.push_back(e1_new);
    
    // now, let's tackle the other edges and triangles in tandem...
    for(uint i=0; i<ts_orig.size(); i++) {
        Tptr                t_orig              = ts_orig[i];
        Vptr                v_opp               = vs_opp[i];
        
        // add mid edge and two tris to v_opp
                            v_opp->edges.push_back(es_mid[i]);
                            v_opp->tris.push_back(t0s_new[i]);
                            v_opp->tris.push_back(t1s_new[i]);
        // add resp. tris to v0 and v1
                            v0->tris.push_back(t0s_new[i]);
                            v1->tris.push_back(t1s_new[i]);
        // find the two non-split edges and add resp. tris
        for(uint k=0; k<3; k++) {
            if(t_orig->verts[k] == v0) {
                Eptr        e1                  = t_orig->edges[k];
                            e1->tris.push_back(t1s_new[i]);
            } else if(t_orig->verts[k] == v1) {
                Eptr        e0                  = t_orig->edges[k];
                            e0->tris.push_back(t0s_new[i]);
            }
        }
    }
    
    /* Now, let's kill all the old geometry.  This consists of:
     *  -   e_split:        the edge to be split
     *  -   ts_orig:        the triangles to be split
     *
     * Luckily, in triangle splitting we know exactly which things
     * must be deleted.  A split cannot make any geometry newly singular.
     */
    
    // kill triangles
    for(Tptr t : ts_orig) {
        // First, unhook this triangle from its faces
        for(uint k=0; k<3; k++) {
            Vptr            v               = t->verts[k];
                            remove(v->tris, t);
            
            Eptr            e               = t->edges[k];
                            remove(e->tris, t);
        }
        // now that we're disconnected, jettison the triangle
                            scratchpad.cache.freeTri(t);
    }
    
    // now, kill the edge that we split
                            remove(v0->edges, e_split);
                            remove(v1->edges, e_split);
                            deallocateRemeshEdge(scratchpad, e_split);
    
    
    // recompute edge scores for all edges whose scores might be effected
    // Don't need to dequeue newly created edges...
                            scoreAndEnqueue(scratchpad.queue, e0_new);
                            scoreAndEnqueue(scratchpad.queue, e1_new);
    for(Eptr e : es_mid) {
                            scoreAndEnqueue(scratchpad.queue, e);
    }
    for(Eptr e : borderEdges) {
                            dequeue(scratchpad.queue, e);
                            scoreAndEnqueue(scratchpad.queue, e);
    }
}

