// +-------------------------------------------------------------------------
// | mesh.bool.tpp
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

#include <queue>

template<class VertData, class TriData>
class Mesh<VertData,TriData>::BoolProblem
{
public:
    BoolProblem(Mesh *owner) : mesh(owner)
    {}
    virtual ~BoolProblem() {}
    
    // do things
    void doSetup(Mesh &rhs);
    
    // choose what to remove
    enum TriCode { KEEP_TRI, DELETE_TRI, FLIP_TRI };
    void doDeleteAndFlip(
        std::function<TriCode(byte bool_alg_data)> classify
    );

private: // methods
    struct BoolEdata {
        bool is_isct;
    };
    
    inline byte& boolData(uint tri_id) {
        return mesh->tris[tri_id].data.bool_alg_data;
    }
    
    void populateECache()
    {
        ecache = mesh->createEGraphCache<BoolEdata>();
        
        // label some of the edges as intersection edges and others as not
        ecache.for_each([&](uint i, uint j, EGraphEntry<BoolEdata> &entry) {
            entry.data.is_isct = false;
            byte operand = boolData(entry.tids[0]);
            for(uint k=1; k<entry.tids.size(); k++) {
                if(boolData(entry.tids[k]) != operand) {
                    entry.data.is_isct = true;
                    break;
                }
            }
        });
    }
    
    inline void for_ecache(
        std::function<void(uint i, uint j,
                           bool isisct,
                           const ShortVec<uint, 2> &tids)> action
    ) {
        ecache.for_each([&](uint i, uint j, EGraphEntry<BoolEdata> &entry) {
            if(entry.data.is_isct) {
                ShortVec<uint, 2> tid0s;
                ShortVec<uint, 2> tid1s;
                for(uint tid : entry.tids) {
                    if(boolData(tid) & 1)
                        tid1s.push_back(tid);
                    else
                        tid0s.push_back(tid);
                }
                action(i,j, true, tid1s);
                action(i,j, true, tid0s);
            } else {
                action(i,j, false, entry.tids);
            }
        });
    }
    
    void prepInsideOutsideTests()
    {
    }
    
    bool isInside(uint tid, byte operand) {
        // find the point to trace outward from...
        Vec3d p(0,0,0);
        p += mesh->verts[mesh->tris[tid].a].pos;
        p += mesh->verts[mesh->tris[tid].b].pos;
        p += mesh->verts[mesh->tris[tid].c].pos;
        p /= 3.0;
        // ok, we've got the point, now let's pick a direction
        Ray3d r;
        r.p = p;
        r.r = Vec3d(drand(0.5,1.5), drand(0.5,1.5), drand(0.5, 1.5));
        
        
        int winding = 0;
        // pass all triangles over ray
        for(Tri &tri : mesh->tris) {
            // ignore triangles from the same operand surface
            if((tri.data.bool_alg_data & 1) == operand)   continue;
            
            double flip = 1.0;
            uint   a = tri.a;
            uint   b = tri.b;
            uint   c = tri.c;
            Vec3d va = mesh->verts[a].pos;
            Vec3d vb = mesh->verts[b].pos;
            Vec3d vc = mesh->verts[c].pos;
            // normalize vertex order (to prevent leaks)
            if(a > b) { std::swap(a, b); std::swap(va, vb); flip = -flip; }
            if(b > c) { std::swap(b, c); std::swap(vb, vc); flip = -flip; }
            if(a > b) { std::swap(a, b); std::swap(va, vb); flip = -flip; }
            
            double t;
            Vec3d bary;
            if(isct_ray_triangle(r, va, vb, vc, &t, &bary)) {
                Vec3d normal = flip * cross(vb - va, vc - va);
                if(dot(normal, r.r) > 0.0) { // UNSAFE
                    winding++;
                } else {
                    winding--;
                }
            }
        }
        
        // now, we've got a winding number to work with...
        return winding > 0;
    }
    
private: // data
    Mesh                        *mesh;
    EGraphCache<BoolEdata>      ecache;
};


static inline double triArea(
    Vec3d a, Vec3d b, Vec3d c
) {
    return len(cross(b-a, c-a));
}


template<class VertData, class TriData>
void Mesh<VertData,TriData>::BoolProblem::doSetup(
    Mesh &rhs
) {
    // Label surfaces...
    mesh->for_tris([](TriData &tri, VertData&, VertData&, VertData&) {
        tri.bool_alg_data = 0;
    });
    rhs.for_tris([](TriData &tri, VertData&, VertData&, VertData&) {
        tri.bool_alg_data = 1;
    });
    
    mesh->disjointUnion(rhs);
    mesh->resolveIntersections();
    
    populateECache();
    
    // form connected components;
    // we get one component for each connected component in one
    // of the two input meshes.
    // These components are not necessarily uniformly inside or outside
    // of the other operand mesh.
    UnionFind uf(mesh->tris.size());
    for_ecache([&](uint, uint, bool, const ShortVec<uint, 2> &tids) {
        uint tid0 = tids[0];
        for(uint k=1; k<tids.size(); k++)
            uf.unionIds(tid0, tids[k]);
    });
    
    // we re-organize the results of the union find as follows:
    std::vector<uint> uq_ids(mesh->tris.size(), uint(-1));
    std::vector< std::vector<uint> > components;
    for(uint i=0; i<mesh->tris.size(); i++) {
        uint ufid = uf.find(i);
        if(uq_ids[ufid] == uint(-1)) { // unassigned
            uint N = components.size();
            components.push_back(std::vector<uint>());
            
            uq_ids[ufid] = uq_ids[i] = N;
            components[N].push_back(i);
        } else { // assigned already
            uq_ids[i] = uq_ids[ufid]; // propagate assignment
            components[uq_ids[i]].push_back(i);
        }
    }
    
    std::vector<bool> visited(mesh->tris.size(), false);
    
    // find the "best" triangle in each component,
    // and ray cast to determine inside-ness vs. outside-ness
    for(auto &comp : components) {
        // find max according to score
        uint best_tid = comp[0];
        double best_area = 0.0;
        // SEARCH
        for(uint tid : comp) {
            Vec3d va = mesh->verts[mesh->tris[tid].a].pos;
            Vec3d vb = mesh->verts[mesh->tris[tid].b].pos;
            Vec3d vc = mesh->verts[mesh->tris[tid].c].pos;
            
            double area = triArea(va, vb, vc);
            if(area > best_area) {
                best_area = area;
                best_tid = tid;
            }
        }
        
        byte operand = boolData(best_tid);
        bool inside = isInside(best_tid, operand);
        
        // NOW PROPAGATE classification throughout the component.
        // do a breadth first propagation
        std::queue<uint> work;
        
        // begin by tagging the first triangle
        boolData(best_tid) |= (inside)? 2 : 0;
        visited[best_tid] = true;
        work.push(best_tid);
        
        while(!work.empty()) {
            uint curr_tid = work.front();
            work.pop();
            
            for(uint k=0; k<3; k++) {
                uint a = mesh->tris[curr_tid].v[k];
                uint b = mesh->tris[curr_tid].v[(k+1)%3];
                auto &entry = ecache(a,b);
                byte inside_sig = boolData(curr_tid) & 2;
                if(entry.data.is_isct)  inside_sig ^= 2;
                for(uint tid : entry.tids) {
                    if(visited[tid])                    continue;
                    if((boolData(tid)&1) != operand)    continue;
                    
                    boolData(tid) |= inside_sig;
                    visited[tid] = true;
                    work.push(tid);
                }
            }
        }
    }
}


template<class VertData, class TriData>
void Mesh<VertData,TriData>::BoolProblem::doDeleteAndFlip(
    std::function<TriCode(byte bool_alg_data)> classify
) {
    TopoCache topocache(mesh);
    
    std::vector<Tptr> toDelete;
    topocache.tris.for_each([&](Tptr tptr) {
        TriCode code = classify(boolData(tptr->ref));
        switch(code) {
        case DELETE_TRI:
            toDelete.push_back(tptr);
            break;
        case FLIP_TRI:
            topocache.flipTri(tptr);
            break;
        case KEEP_TRI:
        default:
            break;
        }
    });
    
    for(Tptr tptr : toDelete) {
        topocache.deleteTri(tptr);
    }
    
    topocache.commit();
}





template<class VertData, class TriData>
void Mesh<VertData,TriData>::boolUnion(Mesh &rhs)
{
    BoolProblem bprob(this);
    
    bprob.doSetup(rhs);
    
    bprob.doDeleteAndFlip([](byte data) -> typename BoolProblem::TriCode {
        if((data & 2) == 2)     // part of op 0/1 INSIDE op 1/0
            return BoolProblem::DELETE_TRI;
        else                    // part of op 0/1 OUTSIDE op 1/0
            return BoolProblem::KEEP_TRI;
    });
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::boolDiff(Mesh &rhs)
{
    BoolProblem bprob(this);
    
    bprob.doSetup(rhs);
    
    bprob.doDeleteAndFlip([](byte data) -> typename BoolProblem::TriCode {
        if(data == 2 ||         // part of op 0 INSIDE op 1
           data == 1)           // part of op 1 OUTSIDE op 0
            return BoolProblem::DELETE_TRI;
        else if(data == 3)      // part of op 1 INSIDE op 1
            return BoolProblem::FLIP_TRI;
        else                    // part of op 0 OUTSIDE op 1
            return BoolProblem::KEEP_TRI;
    });
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::boolIsct(Mesh &rhs)
{
    BoolProblem bprob(this);
    
    bprob.doSetup(rhs);
    
    bprob.doDeleteAndFlip([](byte data) -> typename BoolProblem::TriCode {
        if((data & 2) == 0)     // part of op 0/1 OUTSIDE op 1/0
            return BoolProblem::DELETE_TRI;
        else                    // part of op 0/1 INSIDE op 1/0
            return BoolProblem::KEEP_TRI;
    });
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::boolXor(Mesh &rhs)
{
    BoolProblem bprob(this);
    
    bprob.doSetup(rhs);
    
    bprob.doDeleteAndFlip([](byte data) -> typename BoolProblem::TriCode {
        if((data & 2) == 0)     // part of op 0/1 OUTSIDE op 1/0
            return BoolProblem::KEEP_TRI;
        else                    // part of op 0/1 INSIDE op 1/0
            return BoolProblem::FLIP_TRI;
    });
}















