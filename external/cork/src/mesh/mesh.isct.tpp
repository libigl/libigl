// +-------------------------------------------------------------------------
// | mesh.isct.tpp
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

#include "mesh.topoCache.tpp"

#include "bbox.h"
#include "quantization.h"
#include "empty3d.h"

#include "aabvh.h"

#define REAL double
extern "C" {
#include "triangle.h"
}

// Alec:
#include <iostream>
#include <stdexcept>

struct GenericVertType;
    struct IsctVertType;
    struct OrigVertType;
struct GenericEdgeType;
    struct IsctEdgeType;
    struct OrigEdgeType;
    struct SplitEdgeType;
struct GenericTriType;

struct GluePointMarker;

//using GVptr             = GenericVertType*;
//    using IVptr         = IsctVertType*;
//    using OVptr         = OrigVertType*;
//using GEptr             = GenericEdgeType*;
//    using IEptr         = IsctEdgeType*;
//    using OEptr         = OrigEdgeType*;
//    using SEptr         = SplitEdgeType*;
//using GTptr             = GenericTriType*;

//using GluePt = GluePointMarker*;

typedef GenericVertType*    GVptr;
typedef IsctVertType*           IVptr;
typedef OrigVertType*	        OVptr;
typedef GenericEdgeType*    GEptr;
typedef IsctEdgeType*           IEptr;
typedef OrigEdgeType*	        OEptr;
typedef SplitEdgeType*	        SEptr;
typedef GenericTriType*     GTptr;

typedef GluePointMarker*    GluePt;

struct GenericVertType
{
    virtual ~GenericVertType() {}
    Vptr                    concrete;
    Vec3d                   coord;
    
    bool                    boundary;
    uint                    idx; // temporary for triangulation marshalling
    
    ShortVec<GEptr,2>       edges;
};
struct IsctVertType : public GenericVertType
{
    GluePt                  glue_marker;
};
struct OrigVertType : public GenericVertType {};

struct GenericEdgeType
{
    virtual ~GenericEdgeType() {}
    Eptr                    concrete;
    
    bool                    boundary;
    uint                    idx; // temporary for triangulation marshalling
    
    GVptr                   ends[2];
    ShortVec<IVptr, 1>      interior;       
};
struct IsctEdgeType : public GenericEdgeType
{
public:
    // use to detect duplicate instances within a triangle
    Tptr                    other_tri_key;
};
struct OrigEdgeType : public GenericEdgeType {};
struct SplitEdgeType : public GenericEdgeType {};

struct GenericTriType
{
    Tptr                    concrete;
    
    GVptr                   verts[3];
};

struct GluePointMarker
{
    // list of all the vertices to be glued...
    ShortVec<IVptr, 3>      copies;
    bool                    split_type; // splits are introduced
                                        // manually, not via intersection
                                        // and therefore use only e pointer
    bool                    edge_tri_type; // true if edge-tri intersection
                                           // false if tri-tri-tri
    Eptr                    e;
    Tptr                    t[3];
};



template<uint LEN> inline
IEptr find_edge(ShortVec<IEptr,LEN> &vec, Tptr key)
{
    for(IEptr ie : vec) {
        if(ie->other_tri_key == key)
            return ie;
    }
    return nullptr;
}

inline Vptr commonVert(Tptr t0, Tptr t1)
{
    for(uint i=0; i<3; i++) {
      for(uint j=0; j<3; j++) {
        if(t0->verts[i] == t1->verts[j])
            return t0->verts[i];
      }
    }
    return nullptr;
}

inline bool hasCommonVert(Tptr t0, Tptr t1)
{
    return (t0->verts[0] == t1->verts[0] ||
            t0->verts[0] == t1->verts[1] ||
            t0->verts[0] == t1->verts[2] ||
            t0->verts[1] == t1->verts[0] ||
            t0->verts[1] == t1->verts[1] ||
            t0->verts[1] == t1->verts[2] ||
            t0->verts[2] == t1->verts[0] ||
            t0->verts[2] == t1->verts[1] ||
            t0->verts[2] == t1->verts[2]);
}

inline bool hasCommonVert(Eptr e, Tptr t)
{
    return (e->verts[0] == t->verts[0] ||
            e->verts[0] == t->verts[1] ||
            e->verts[0] == t->verts[2] ||
            e->verts[1] == t->verts[0] ||
            e->verts[1] == t->verts[1] ||
            e->verts[1] == t->verts[2]);
}

inline void disconnectGE(GEptr ge)
{
    ge->ends[0]->edges.erase(ge);
    ge->ends[1]->edges.erase(ge);
    for(IVptr iv : ge->interior)
        iv->edges.erase(ge);
}

// should deal with via pointers
template<class VertData, class TriData>
class Mesh<VertData, TriData>::TriangleProblem
{
public:
    TriangleProblem() {}
    ~TriangleProblem() {}
    
    inline void init(IsctProblem *iprob, Tptr t) {
        the_tri             = t;
        // extract original edges/verts
        for(uint k=0; k<3; k++)
            overts[k]       = iprob->newOrigVert(the_tri->verts[k]);
        for(uint k=0; k<3; k++) {
            oedges[k]       = iprob->newOrigEdge(the_tri->edges[k],
                                                 overts[(k+1)%3],
                                                 overts[(k+2)%3]);
        }
    }
    
private: // may actually not add edge, but instead just hook up endpoint
    inline void addEdge(
        IsctProblem *iprob, IVptr iv, Tptr tri_key
    ) {
        IEptr       ie              = find_edge(iedges, tri_key);
        if(ie) { // if the edge is already present
                    ie->ends[1]     = iv;
                    iv->edges.push_back(ie);
        } else { // if the edge is being added
                    ie              = iprob->newIsctEdge(iv, tri_key);
                    iedges.push_back(ie);
        }
    }
    void addBoundaryHelper(
        Eptr edge, IVptr iv
    ) {
                    iv->boundary    = true;
                    iverts.push_back(iv);
        // hook up point to boundary edge interior!
        for(uint k=0; k<3; k++) {
            OEptr   oe              = oedges[k];
            if(oe->concrete == edge) {
                    oe->interior.push_back(iv);
                    iv->edges.push_back(oe);
                    break;
            }
        }
    }
public:
    // specify reference glue point and edge piercing this triangle.
    IVptr addInteriorEndpoint(
        IsctProblem *iprob, Eptr edge, GluePt glue
    ) {
        IVptr       iv              = iprob->newIsctVert(edge, the_tri, glue);
                    iv->boundary    = false;
                    iverts.push_back(iv);
        for(Tptr tri_key : edge->tris) {
                    addEdge(iprob, iv, tri_key);
        }
        return iv;
    }
    // specify the other triangle cutting this one, the edge cut,
    // and the resulting point of intersection
    void addBoundaryEndpoint(
        IsctProblem *iprob, Tptr tri_key, Eptr edge, IVptr iv
    ) {
                    iv              = iprob->copyIsctVert(iv);
                    addBoundaryHelper(edge, iv);
        // handle edge extending into interior
                    addEdge(iprob, iv, tri_key);
    }
    IVptr addBoundaryEndpoint(
        IsctProblem *iprob, Tptr tri_key, Eptr edge, Vec3d coord, GluePt glue
    ) {
        IVptr       iv              = iprob->newSplitIsctVert(coord, glue);
                    addBoundaryHelper(edge, iv);
        // handle edge extending into interior
                    addEdge(iprob, iv, tri_key);
        return iv;
    }
    // Should only happen for manually inserted split points on
    // edges, not for points computed via intersection...
    IVptr addBoundaryPointAlone(
        IsctProblem *iprob, Eptr edge, Vec3d coord, GluePt glue
    ) {
        IVptr       iv              = iprob->newSplitIsctVert(coord, glue);
                    addBoundaryHelper(edge, iv);
        return iv;
    }
    void addInteriorPoint(
        IsctProblem *iprob, Tptr t0, Tptr t1, GluePt glue
    ) {
        // note this generates wasted re-computation of coordinates 3X
        IVptr       iv              = iprob->newIsctVert(the_tri, t0, t1, glue);
                    iv->boundary    = false;
                    iverts.push_back(iv);
        // find the 2 interior edges
        for(IEptr ie : iedges) {
            if(ie->other_tri_key == t0 ||
               ie->other_tri_key == t1) {
                    ie->interior.push_back(iv);
                    iv->edges.push_back(ie);
            }
        }
    }
    
    // run after we've accumulated all the elements
    void consolidate(IsctProblem *iprob) {
        // identify all intersection edges missing endpoints
        // and check to see if we can assign an original vertex
        // as the appropriate endpoint.
        for(IEptr ie : iedges) {
            if(ie->ends[1] == nullptr) {
                // try to figure out which vertex must be the endpoint...
                Vptr vert = commonVert(the_tri, ie->other_tri_key);
                if(!vert) {
                    std::cout << "the  edge is "
                              << ie->ends[0] << ",  "
                              << ie->ends[1] << std::endl;
                    IVptr iv = dynamic_cast<IVptr>(ie->ends[0]);
                    std::cout << "   "
                              << iv->glue_marker->edge_tri_type
                              << std::endl;
                    std::cout << "the   tri is " << the_tri << ": "
                              << *the_tri << std::endl;
                    std::cout << "other tri is " << ie->other_tri_key << ": "
                              << *(ie->other_tri_key) << std::endl;
                    std::cout << "coordinates for triangles" << std::endl;
                    std::cout << "the tri" << std::endl;
                    for(uint k=0; k<3; k++)
                        std::cout << iprob->vPos(the_tri->verts[k])
                                  << std::endl;
                    for(uint k=0; k<3; k++)
                        std::cout << iprob->vPos(ie->other_tri_key->verts[k])
                                  << std::endl;
                    std::cout << "degen count:"
                              << Empty3d::degeneracy_count << std::endl;
                    std::cout << "exact count: "
                              << Empty3d::exact_count << std::endl;
                }
                ENSURE(vert); // bad if we can't find a common vertex
                // then, find the corresponding OVptr, and connect
                for(uint k=0; k<3; k++) {
                    if(overts[k]->concrete == vert) {
                        ie->ends[1] = overts[k];
                        overts[k]->edges.push_back(ie);
                        break;
                    }
                }
            }
        }
        
        ENSURE(isValid());
    }
    
    bool isValid() const {
        ENSURE(the_tri);
        
        return true;
    }
    
    void subdivide(IsctProblem *iprob) {
        // collect all the points, and create more points as necessary
        ShortVec<GVptr, 7> points;
        for(uint k=0; k<3; k++) {
            points.push_back(overts[k]);
            //std::cout << k << ": id " << overts[k]->concrete->ref << std::endl;
        }
        for(IVptr iv : iverts) {
            //iprob->buildConcreteVert(iv);
            points.push_back(iv);
            /*std::cout << "  " << points.size() - 1
                          << " (" << iv->glue_marker->edge_tri_type
                          << ") ";
            if(iv->glue_marker->edge_tri_type) {
                Eptr e = iv->glue_marker->e;
                Vec3d p0 = iprob->vPos(e->verts[0]);
                Vec3d p1 = iprob->vPos(e->verts[1]);
                std::cout << " "
                          << e->verts[0]->ref << p0
                          << "  " << e->verts[1]->ref << p1;
                
                Tptr t = iv->glue_marker->t[0];
                p0 = iprob->vPos(t->verts[0]);
                p1 = iprob->vPos(t->verts[1]);
                Vec3d p2 = iprob->vPos(t->verts[2]);
                std::cout << "        "
                          << t->verts[0]->ref << p0
                          << "  " << t->verts[1]->ref << p1
                          << "  " << t->verts[2]->ref << p2;
            }
            std::cout
                          << std::endl;*/
        }
        for(uint i=0; i<points.size(); i++)
            points[i]->idx = i;
        
        // split edges and marshall data
        // for safety, we zero out references to pre-subdivided edges,
        // which may have been destroyed
        ShortVec<GEptr, 8> edges;
        for(uint k=0; k<3; k++) {
            //std::cout << "oedge:  "
            //          << oedges[k]->ends[0]->idx << "; "
            //          << oedges[k]->ends[1]->idx << std::endl;
            //for(IVptr iv : oedges[k]->interior)
            //    std::cout << "  " << iv->idx << std::endl;
            subdivideEdge(iprob, oedges[k], edges);
            oedges[k]       = nullptr;
        }
        //std::cout << "THE TRI: " << the_tri->verts[0]->ref
        //          << "; " << the_tri->verts[1]->ref
        //          << "; " << the_tri->verts[2]->ref
        //          << std::endl;
        for(IEptr &ie : iedges) {
            //std::cout << "iedge:  "
            //          << ie->ends[0]->idx << "; "
            //          << ie->ends[1]->idx << std::endl;
            //std::cout << "other tri: " << ie->other_tri_key->verts[0]->ref
            //          << "; " << ie->other_tri_key->verts[1]->ref
            //          << "; " << ie->other_tri_key->verts[2]->ref
            //          << std::endl;
            //for(IVptr iv : ie->interior)
            //    std::cout << "  " << iv->idx 
            //              << " (" << iv->glue_marker->edge_tri_type
            //              << ") " << std::endl;
            subdivideEdge(iprob, ie, edges);
            ie              = nullptr;
        }
        for(uint i=0; i<edges.size(); i++)
            edges[i]->idx = i;
        
        // find 2 dimensions to project onto
        // get normal
        Vec3d normal = cross( overts[1]->coord - overts[0]->coord,
                              overts[2]->coord - overts[0]->coord );
        uint normdim = maxDim(abs(normal));
        uint dim0 = (normdim+1)%3;
        uint dim1 = (normdim+2)%3;
        double sign_flip = (normal.v[normdim] < 0.0)? -1.0 : 1.0;
        
        struct triangulateio in, out;
        
        /* Define input points. */
        in.numberofpoints           = points.size();
        in.numberofpointattributes  = 0;
        in.pointlist                = new REAL[in.numberofpoints * 2];
        in.pointattributelist       = nullptr;
        in.pointmarkerlist          = new int[in.numberofpoints];
        for(int k=0; k<in.numberofpoints; k++) {
            in.pointlist[k*2 + 0] = points[k]->coord.v[dim0];
            in.pointlist[k*2 + 1] = points[k]->coord.v[dim1] * sign_flip;
            in.pointmarkerlist[k] = (points[k]->boundary)? 1 : 0;
        }
        
        /* Define the input segments */
        in.numberofsegments = edges.size();
        in.numberofholes = 0;// yes, zero
        in.numberofregions = 0;// not using regions
        in.segmentlist = new int[in.numberofsegments * 2];
        in.segmentmarkerlist = new int[in.numberofsegments];
        for(int k=0; k<in.numberofsegments; k++) {
            in.segmentlist[k*2 + 0] = edges[k]->ends[0]->idx;
            in.segmentlist[k*2 + 1] = edges[k]->ends[1]->idx;
            in.segmentmarkerlist[k] = (edges[k]->boundary)? 1 : 0;
        }
        
        // to be safe... declare 0 triangle attributes on input
        in.numberoftriangles = 0;
        in.numberoftriangleattributes = 0;
        
        /* set for flags.... */
        out.pointlist = nullptr;
        out.pointattributelist = nullptr; // not necessary if using -N or 0 attr
        out.pointmarkerlist = nullptr;
        out.trianglelist = nullptr; // not necessary if using -E
        //out.triangleattributelist = null; // not necessary if using -E or 0 attr
        //out.trianglearealist = // only needed with -r and -a
        //out.neighborlist = null; // only neccesary if -n is used
        out.segmentlist = nullptr; // NEED THIS; output segments go here
        out.segmentmarkerlist = nullptr; // NEED THIS for OUTPUT SEGMENTS
        //out.edgelist = null; // only necessary if -e is used
        //out.edgemarkerlist = null; // only necessary if -e is used
                
        // solve the triangulation problem
        char *params = (char*)("pzQYY");
        //char *debug_params = (char*)("pzYYVC");
        triangulate(params, &in, &out, nullptr);
        
        if(out.numberofpoints != in.numberofpoints) {
            std::cout << "out.numberofpoints: "
                      << out.numberofpoints << std::endl;
            std::cout << "points.size(): " << points.size() << std::endl;
            std::cout << "dumping out the points' coordinates" << std::endl;
            for(uint k=0; k<points.size(); k++) {
                GVptr gv = points[k];
                std::cout << "  " << gv->coord
                          << "  " << gv->idx << std::endl;
            }
            
            std::cout << "dumping out the segments" << std::endl;
            for(int k=0; k<in.numberofsegments; k++)
                std::cout << "  " << in.segmentlist[k*2 + 0]
                          << "; " << in.segmentlist[k*2 + 1]
                          << " (" << in.segmentmarkerlist[k]
                          << ") " << std::endl;
            
            std::cout << "dumping out the solved for triangles now..."
                      << std::endl;
            for(int k=0; k<out.numberoftriangles; k++) {
                std::cout << "  "
                          << out.trianglelist[(k*3)+0] << "; "
                          << out.trianglelist[(k*3)+1] << "; "
                          << out.trianglelist[(k*3)+2] << std::endl;
            }
        }
        ENSURE(out.numberofpoints == in.numberofpoints);
        
        //std::cout << "number of points in: " << in.numberofpoints
        //          << std::endl;
        //std::cout << "number of edges in: " << in.numberofsegments
        //          << std::endl;
        //std::cout << "number of triangles out: " << out.numberoftriangles
        //          << std::endl;
        
        gtris.resize(out.numberoftriangles);
        for(int k=0; k<out.numberoftriangles; k++) {
            GVptr       gv0         = points[out.trianglelist[(k*3)+0]];
            GVptr       gv1         = points[out.trianglelist[(k*3)+1]];
            GVptr       gv2         = points[out.trianglelist[(k*3)+2]];
                        gtris[k]    = iprob->newGenericTri(gv0, gv1, gv2);
        }
        
        
        // clean up after triangulate...
            // in free
        free(in.pointlist);
        free(in.pointmarkerlist);
        free(in.segmentlist);
        free(in.segmentmarkerlist);
            // out free
        free(out.pointlist);
        //free(out.pointattributelist);
        free(out.pointmarkerlist);
        free(out.trianglelist);
        //free(out.triangleattributelist);
        //free(out.trianglearealist);
        //free(out.neighborlist);
        free(out.segmentlist);
        free(out.segmentmarkerlist);
        //free(out.edgelist);
        //free(out.edgemarkerlist);
    }

private:
    void subdivideEdge(IsctProblem *iprob, GEptr ge, ShortVec<GEptr, 8> &edges)
    {
        if(ge->interior.size() == 0) {
            //if(typeid(ge) == typeid(IEptr)) { // generate new edge
            //    iprob->buildConcreteEdge(ge);
            //}
            edges.push_back(ge);
        } else if(ge->interior.size() == 1) { // common case
            SEptr       se0     = iprob->newSplitEdge(ge->ends[0],
                                                      ge->interior[0],
                                                      ge->boundary);
            SEptr       se1     = iprob->newSplitEdge(ge->interior[0],
                                                      ge->ends[1],
                                                      ge->boundary);
                        //iprob->buildConcreteEdge(se0);
                        //iprob->buildConcreteEdge(se1);
                        edges.push_back(se0);
                        edges.push_back(se1);
            
            // get rid of old edge
                        iprob->releaseEdge(ge);
        } else { // sorting is the uncommon case
            // determine the primary dimension and direction of the edge
            Vec3d       dir     = ge->ends[1]->coord - ge->ends[0]->coord;
            uint        dim     = (fabs(dir.x) > fabs(dir.y))?
                                    ((fabs(dir.x) > fabs(dir.z))? 0 : 2) :
                                    ((fabs(dir.y) > fabs(dir.z))? 1 : 2);
            double      sign    = (dir.v[dim] > 0.0)? 1.0 : -1.0;
            
            // pack the interior vertices into a vector for sorting
            std::vector< std::pair<double,IVptr> > verts;
            for(IVptr iv : ge->interior) {
                        verts.push_back(std::make_pair(
                            // if the sort is ascending, then we're good...
                            sign * iv->coord.v[dim],
                            iv
                        ));
            }
            // ... and sort the vector
                        std::sort(verts.begin(), verts.end());
            // then, write the verts into a new container with the endpoints
            std::vector<GVptr>  allv(verts.size()+2);
                        allv[0]             = ge->ends[0];
                        allv[allv.size()-1] = ge->ends[1];
            for(uint k=0; k<verts.size(); k++)
                        allv[k+1]           = verts[k].second;
            
            // now create and accumulate new split edges
            for(uint i=1; i<allv.size(); i++) {
                SEptr   se      = iprob->newSplitEdge(allv[i-1],
                                                      allv[i],
                                                      ge->boundary);
                        edges.push_back(se);
            }
            // get rid of old edge
                        iprob->releaseEdge(ge);
        }
    }
    
public: // data
    ShortVec<IVptr, 4>      iverts;
    ShortVec<IEptr, 2>      iedges;
    // original triangle elements
    OVptr                   overts[3];
    OEptr                   oedges[3];
    
    ShortVec<GTptr, 8>      gtris;
    
    Tptr                    the_tri;
};


template<class VertData, class TriData>
class Mesh<VertData,TriData>::IsctProblem : public TopoCache
{
public:
    IsctProblem(Mesh *owner) : TopoCache(owner)
    {
        // initialize all the triangles to NOT have an associated tprob
        TopoCache::tris.for_each([](Tptr t) {
            t->data = nullptr;
        });
        
        // Callibrate the quantization unit...
        double maxMag = 0.0;
        for(VertData &v : TopoCache::mesh->verts) {
            maxMag = std::max(maxMag, max(abs(v.pos)));
        }
        Quantization::callibrate(maxMag);
        
        // and use vertex auxiliary data to store quantized vertex coordinates
        uint N = TopoCache::mesh->verts.size();
        quantized_coords.resize(N);
        uint write = 0;
        TopoCache::verts.for_each([&](Vptr v) {
#ifdef _WIN32
            Vec3d raw = mesh->verts[v->ref].pos;
#else
            Vec3d raw = TopoCache::mesh->verts[v->ref].pos;
#endif
            quantized_coords[write].x = Quantization::quantize(raw.x);
            quantized_coords[write].y = Quantization::quantize(raw.y);
            quantized_coords[write].z = Quantization::quantize(raw.z);
            v->data = &(quantized_coords[write]);
            write++;
        });
    }
    
    virtual ~IsctProblem() {}
    
    // access auxiliary quantized coordinates
    inline Vec3d vPos(Vptr v) const {
        return *(reinterpret_cast<Vec3d*>(v->data));
    }
    
    Tprob getTprob(Tptr t) {
        Tprob prob = reinterpret_cast<Tprob>(t->data);
        if(!prob) {
            t->data = prob = tprobs.alloc();
            prob->init(this, t);
        }
        return prob;
    }
    GluePt newGluePt() {
        GluePt glue = glue_pts.alloc();
        glue->split_type = false;
        return glue;
    }
    
    inline IVptr newIsctVert(Eptr e, Tptr t, GluePt glue) {
        IVptr       iv                  = ivpool.alloc();
                    iv->concrete        = nullptr;
                    iv->coord           = computeCoords(e, t);
                    iv->glue_marker     = glue;
                    glue->copies.push_back(iv);
        return      iv;
    }
    inline IVptr newIsctVert(Tptr t0, Tptr t1, Tptr t2, GluePt glue) {
        IVptr       iv                  = ivpool.alloc();
                    iv->concrete        = nullptr;
                    iv->coord           = computeCoords(t0, t1, t2);
                    iv->glue_marker     = glue;
                    glue->copies.push_back(iv);
        return      iv;
    }
    inline IVptr newSplitIsctVert(Vec3d coords, GluePt glue) {
        IVptr       iv                  = ivpool.alloc();
                    iv->concrete        = nullptr;
                    iv->coord           = coords;
                    iv->glue_marker     = glue;
                    glue->copies.push_back(iv);
        return      iv;
    }
    inline IVptr copyIsctVert(IVptr orig) {
        IVptr       iv                  = ivpool.alloc();
                    iv->concrete        = nullptr;
                    iv->coord           = orig->coord;
                    iv->glue_marker     = orig->glue_marker;
                    orig->glue_marker->copies.push_back(iv);
        return      iv;
    }
    inline IEptr newIsctEdge(IVptr endpoint, Tptr tri_key) {
        IEptr       ie                  = iepool.alloc();
                    ie->concrete        = nullptr;
                    ie->boundary        = false;
                    ie->ends[0]         = endpoint;
                    endpoint->edges.push_back(ie);
                    ie->ends[1]         = nullptr; // other end null
                    ie->other_tri_key   = tri_key;
        return      ie;
    }
    
    inline OVptr newOrigVert(Vptr v) {
        OVptr       ov                  = ovpool.alloc();
                    ov->concrete        = v;
                    ov->coord           = vPos(v);
                    ov->boundary        = true;
        return      ov;
    }
    inline OEptr newOrigEdge(Eptr e, OVptr v0, OVptr v1) {
        OEptr       oe                  = oepool.alloc();
                    oe->concrete        = e;
                    oe->boundary        = true;
                    oe->ends[0]         = v0;
                    oe->ends[1]         = v1;
                    v0->edges.push_back(oe);
                    v1->edges.push_back(oe);
        return      oe;
    }
    inline SEptr newSplitEdge(GVptr v0, GVptr v1, bool boundary) {
        SEptr       se                  = sepool.alloc();
                    se->concrete        = nullptr;
                    se->boundary        = boundary;
                    se->ends[0]         = v0;
                    se->ends[1]         = v1;
                    v0->edges.push_back(se);
                    v1->edges.push_back(se);
        return      se;
    }
    
    inline GTptr newGenericTri(GVptr v0, GVptr v1, GVptr v2) {
        GTptr       gt                  = gtpool.alloc();
                    gt->verts[0]        = v0;
                    gt->verts[1]        = v1;
                    gt->verts[2]        = v2;
                    gt->concrete        = nullptr;
        return      gt;
    }
    
    inline void releaseEdge(GEptr ge) {
        disconnectGE(ge);
        IEptr       ie      = dynamic_cast<IEptr>(ge);
        if(ie) {
                    iepool.free(ie);
        } else {
            OEptr   oe      = dynamic_cast<OEptr>(ge);
                    ENSURE(oe);
                    oepool.free(oe);
        }
    }
    
    inline void killIsctVert(IVptr iv) {
        iv->glue_marker->copies.erase(iv);
        if(iv->glue_marker->copies.size() == 0)
            glue_pts.free(iv->glue_marker);
        
        for(GEptr ge : iv->edges) {
            // disconnect
            ge->interior.erase(iv);
            if(ge->ends[0] == iv)   ge->ends[0] = nullptr;
            if(ge->ends[1] == iv)   ge->ends[1] = nullptr;
        }
        
        ivpool.free(iv);
    }
    
    inline void killIsctEdge(IEptr ie) {
        // an endpoint may be an original vertex
        if(ie->ends[1])
            ie->ends[1]->edges.erase(ie);
        iepool.free(ie);
    }
    
    inline void killOrigVert(OVptr ov) {
        ovpool.free(ov);
    }
    inline void killOrigEdge(OEptr oe) {
        oepool.free(oe);
    }
    
    /*
    inline void buildConcreteVert(GVptr gv) {
                    gv->concrete        = newVert();
                    // NEED DATA SOMEHOW??
    }
    // make sure the endpoints have concrete versions first!
    inline void buildConcreteEdge(GEptr ge) {
        Eptr        e   = ge->concrete  = newEdge();
        Vptr        v0  = e->verts[0]   = ge->ends[0]->concrete;
        Vptr        v1  = e->verts[1]   = ge->ends[1]->concrete;
                    v0->edges.push_back(e);
                    v1->edges.push_back(e);
    }
    inline void buildConcreteTri(GVptr gv0, GVptr gv1, GVptr gv2) {
        // create edges as necessary...
    }*/
    
    bool hasIntersections(); // test for iscts, exit if one is found
    
    void findIntersections();
    void resolveAllIntersections();
private:
    // if we encounter ambiguous degeneracies, then this
    // routine returns false, indicating that the computation aborted.
    bool tryToFindIntersections();
    // In that case, we can perturb the positions of points
    void perturbPositions();
    // in order to give things another try, discard partial work
    void reset();
public:
    
    void dumpIsctPoints(std::vector<Vec3d> *points);
    void dumpIsctEdges(std::vector< std::pair<Vec3d,Vec3d> > *edges);
    
protected: // DATA
    IterPool<GluePointMarker>   glue_pts;
    IterPool<TriangleProblem>   tprobs;
    
    IterPool<IsctVertType>      ivpool;
    IterPool<OrigVertType>      ovpool;
    IterPool<IsctEdgeType>      iepool;
    IterPool<OrigEdgeType>      oepool;
    IterPool<SplitEdgeType>     sepool;
    IterPool<GenericTriType>    gtpool;
private:
    std::vector<Vec3d>          quantized_coords;
private:
    inline void for_edge_tri(std::function<bool(Eptr e, Tptr t)>);
    inline void bvh_edge_tri(std::function<bool(Eptr e, Tptr t)>);

    inline GeomBlob<Eptr> edge_blob(Eptr e);
    inline BBox3d bboxFromTptr(Tptr t);
    
    inline BBox3d buildBox(Eptr e) const;
    inline BBox3d buildBox(Tptr t) const;
    
    inline void marshallArithmeticInput(Empty3d::TriIn &input, Tptr t) const;
    inline void marshallArithmeticInput(Empty3d::EdgeIn &input, Eptr e) const;
    inline void marshallArithmeticInput(
        Empty3d::TriEdgeIn &input, Eptr e, Tptr t) const;
    inline void marshallArithmeticInput(
        Empty3d::TriTriTriIn &input, Tptr t0, Tptr t1, Tptr t2) const;
    
    bool checkIsct(Eptr e, Tptr t) const;
    bool checkIsct(Tptr t0, Tptr t1, Tptr t2) const;
    
    Vec3d computeCoords(Eptr e, Tptr t) const;
    Vec3d computeCoords(Tptr t0, Tptr t1, Tptr t2) const;
    
    void fillOutVertData(GluePt glue, VertData &data);
    void fillOutTriData(Tptr tri, Tptr parent);
private:
    class EdgeCache;
    
private: // functions here to get around a GCC bug...
    void createRealPtFromGluePt(GluePt glue);
    void createRealTriangles(Tprob tprob, EdgeCache &ecache);
};

template<class T, uint LEN> inline
void for_pairs(
    ShortVec<T,LEN> &vec,
    std::function<void(T&,T&)> func
) {
    for(uint i=0; i<vec.size(); i++)
        for(uint j=i+1; j<vec.size(); j++)
            func(vec[i], vec[j]);
}

struct TriTripleTemp
{
    Tptr t0, t1, t2;
    TriTripleTemp(Tptr tp0, Tptr tp1, Tptr tp2) :
        t0(tp0), t1(tp1), t2(tp2)
    {}
};

template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::for_edge_tri(
    std::function<bool(Eptr e, Tptr t)> func
) {
    bool aborted = false;
    TopoCache::tris.for_each([&](Tptr t) {
        TopoCache::edges.for_each([&](Eptr e) {
            if(!aborted) {
                if(!func(e, t))
                    aborted = true;
            }
        });
    });
}

template<class VertData, class TriData> inline
GeomBlob<Eptr> Mesh<VertData,TriData>::IsctProblem::edge_blob(
    Eptr e
) {
    GeomBlob<Eptr>  blob;
    blob.bbox = buildBox(e);
    blob.point = (blob.bbox.minp + blob.bbox.maxp) / 2.0;
    blob.id = e;
    return blob;
}

template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::bvh_edge_tri(
    std::function<bool(Eptr e, Tptr t)> func
) {
    std::vector< GeomBlob<Eptr> > edge_geoms;
    TopoCache::edges.for_each([&](Eptr e) {
        edge_geoms.push_back(edge_blob(e));
    });
    AABVH<Eptr> edgeBVH(edge_geoms);
    
    // use the acceleration structure
    bool aborted = false;
    TopoCache::tris.for_each([&](Tptr t) {
        // compute BBox
        BBox3d bbox = buildBox(t);
        if(!aborted) {
            edgeBVH.for_each_in_box(bbox, [&](Eptr e) {
                if(!func(e,t))
                    aborted = true;
            });
        }
    });
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::IsctProblem::tryToFindIntersections()
{
    Empty3d::degeneracy_count = 0;
    // Find all edge-triangle intersection points
    //for_edge_tri([&](Eptr eisct, Tptr tisct)->bool{
    bvh_edge_tri([&](Eptr eisct, Tptr tisct)->bool{
      if(checkIsct(eisct,tisct)) {
        GluePt      glue                    = newGluePt();
                    glue->edge_tri_type     = true;
                    glue->e                 = eisct;
                    glue->t[0]              = tisct;
        // first add point and edges to the pierced triangle
        IVptr iv = getTprob(tisct)->addInteriorEndpoint(this, eisct, glue);
        for(Tptr tri : eisct->tris) {
            getTprob(tri)->addBoundaryEndpoint(this, tisct, eisct, iv);
        }
      }
      if(Empty3d::degeneracy_count > 0)
        return false; // break
      else
        return true; // continue
    });
    if(Empty3d::degeneracy_count > 0) {
        return false;   // restart / abort
    }
    
    // we're going to peek into the triangle problems in order to
    // identify potential candidates for Tri-Tri-Tri intersections
    std::vector<TriTripleTemp> triples;
    tprobs.for_each([&](Tprob tprob) {
        Tptr t0 = tprob->the_tri;
        // Scan pairs of existing edges to create candidate triples
        for_pairs<IEptr,2>(tprob->iedges, [&](IEptr &ie1, IEptr &ie2){
            Tptr t1 = ie1->other_tri_key;
            Tptr t2 = ie2->other_tri_key;
            // This triple might be considered three times,
            // one for each triangle it contains.
            // To prevent duplication, only proceed if this is
            // the least triangle according to an arbitrary ordering
            if(t0 < t1 && t0 < t2) {
                // now look for the third edge.  We're not
                // sure if it exists...
                Tprob prob1 = reinterpret_cast<Tprob>(t1->data);
                for(IEptr ie : prob1->iedges) {
                    if(ie->other_tri_key == t2) {
                        // ADD THE TRIPLE
                        triples.push_back(TriTripleTemp(t0, t1, t2));
                    }
                }
            }
        });
    });
    // Now, we've collected a list of Tri-Tri-Tri intersection candidates.
    // Check to see if the intersections actually exist.
    for(TriTripleTemp t : triples) {
        if(!checkIsct(t.t0, t.t1, t.t2))    continue;
        
        // Abort if we encounter a degeneracy
        if(Empty3d::degeneracy_count > 0)   break;
        
        GluePt      glue                    = newGluePt();
                    glue->edge_tri_type     = false;
                    glue->t[0]              = t.t0;
                    glue->t[1]              = t.t1;
                    glue->t[2]              = t.t2;
        getTprob(t.t0)->addInteriorPoint(this, t.t1, t.t2, glue);
        getTprob(t.t1)->addInteriorPoint(this, t.t0, t.t2, glue);
        getTprob(t.t2)->addInteriorPoint(this, t.t0, t.t1, glue);
    }
    if(Empty3d::degeneracy_count > 0) {
        return false;   // restart / abort
    }
    
    return true;
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::perturbPositions()
{
  // Alec
    //const double EPSILON = 1.0e-5; // perturbation epsilon
    const double EPSILON = 2.5e-9; // perturbation epsilon
    for(Vec3d &coord : quantized_coords) {
        Vec3d perturbation(Quantization::quantize(drand(-EPSILON, EPSILON)),
                           Quantization::quantize(drand(-EPSILON, EPSILON)),
                           Quantization::quantize(drand(-EPSILON, EPSILON)));
        coord += perturbation;
    }
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::reset()
{
    // the data pointer in the triangles points to tproblems
    // that we're about to destroy,
    // so zero out all those pointers first!
    tprobs.for_each([](Tprob tprob) {
        Tptr t = tprob->the_tri;
        t->data = nullptr;
    });
    
    glue_pts.clear();
    tprobs.clear();
    
    ivpool.clear();
    ovpool.clear();
    iepool.clear();
    oepool.clear();
    sepool.clear();
    gtpool.clear();
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::findIntersections()
{
  // Alec:
    int nTrys = 300;
    //perturbPositions(); // always perturb for safety...
    while(nTrys > 0) {
        if(!tryToFindIntersections()) {
            reset();
            perturbPositions();
            nTrys--;
        } else {
            break;
        }
    }
    if(nTrys <= 0) {
        CORK_ERROR("Ran out of tries to perturb the mesh");
        //exit(1);
        throw std::runtime_error("Ran out of tries to perturb the mesh");
    }
    
    // ok all points put together,
    // all triangle problems assembled.
    // Some intersection edges may have original vertices as endpoints
    // we consolidate the problems to check for cases like these.
    tprobs.for_each([&](Tprob tprob) {
        tprob->consolidate(this);
    });
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::IsctProblem::hasIntersections()
{
    bool foundIsct = false;
    Empty3d::degeneracy_count = 0;
    // Find some edge-triangle intersection point...
    bvh_edge_tri([&](Eptr eisct, Tptr tisct)->bool{
      if(checkIsct(eisct,tisct)) {
        foundIsct = true;
        return false; // break;
      }
      if(Empty3d::degeneracy_count > 0) {
        return false; // break;
      }
      return true; // continue
    });
    
    if(Empty3d::degeneracy_count > 0 || foundIsct) {
        std::cout << "This self-intersection might be spurious. "
                     "Degeneracies were detected." << std::endl;
        return true;
    } else {
        return false;
    }
}


template<class VertData, class TriData> inline
BBox3d Mesh<VertData,TriData>::IsctProblem::buildBox(Eptr e) const
{
    Vec3d p0 = vPos(e->verts[0]);
    Vec3d p1 = vPos(e->verts[1]);
    return BBox3d(min(p0, p1), max(p0, p1));
}
template<class VertData, class TriData> inline
BBox3d Mesh<VertData,TriData>::IsctProblem::buildBox(Tptr t) const
{
    Vec3d p0 = vPos(t->verts[0]);
    Vec3d p1 = vPos(t->verts[1]);
    Vec3d p2 = vPos(t->verts[2]);
    return BBox3d(min(p0, min(p1, p2)), max(p0, max(p1, p2)));
}

template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::marshallArithmeticInput(
    Empty3d::EdgeIn &input, Eptr e
) const {
    input.p[0] = vPos(e->verts[0]);
    input.p[1] = vPos(e->verts[1]);
}
template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::marshallArithmeticInput(
    Empty3d::TriIn &input, Tptr t
) const {
    input.p[0] = vPos(t->verts[0]);
    input.p[1] = vPos(t->verts[1]);
    input.p[2] = vPos(t->verts[2]);
}
template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::marshallArithmeticInput(
    Empty3d::TriEdgeIn &input,
    Eptr e, Tptr t
) const {
    marshallArithmeticInput(input.edge, e);
    marshallArithmeticInput(input.tri, t);
}
template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::marshallArithmeticInput(
    Empty3d::TriTriTriIn &input,
    Tptr t0, Tptr t1, Tptr t2
) const {
    marshallArithmeticInput(input.tri[0], t0);
    marshallArithmeticInput(input.tri[1], t1);
    marshallArithmeticInput(input.tri[2], t2);
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::IsctProblem::checkIsct(Eptr e, Tptr t) const
{
    // simple bounding box cull; for acceleration, not correctness
    BBox3d      ebox        = buildBox(e);
    BBox3d      tbox        = buildBox(t);
    if(!hasIsct(ebox, tbox))
                return      false;
    
    // must check whether the edge and triangle share a vertex
    // if so, then trivially we know they intersect in exactly that vertex
    // so we discard this case from consideration.
    if(hasCommonVert(e, t))
                return      false;
    
    Empty3d::TriEdgeIn input;
    marshallArithmeticInput(input, e, t);
    //bool empty = Empty3d::isEmpty(input);
    bool empty = Empty3d::emptyExact(input);
    return !empty;
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::IsctProblem::checkIsct(
    Tptr t0, Tptr t1, Tptr t2
) const {
    // This function should only be called if we've already
    // identified that the intersection edges
    //      (t0,t1), (t0,t2), (t1,t2)
    // exist.
    // From this, we can conclude that each pair of triangles
    // shares no more than a single vertex in common.
    //  If each of these shared vertices is different from each other,
    // then we could legitimately have a triple intersection point,
    // but if all three pairs share the same vertex in common, then
    // the intersection of the three triangles must be that vertex.
    // So, we must check for such a single vertex in common amongst
    // the three triangles
    Vptr common = commonVert(t0, t1);
    if(common) {
        for(uint i=0; i<3; i++)
            if(common == t2->verts[i])
                return      false;
    }
    
    Empty3d::TriTriTriIn input;
    marshallArithmeticInput(input, t0, t1, t2);
    //bool empty = Empty3d::isEmpty(input);
    bool empty = Empty3d::emptyExact(input);
    return !empty;
}

template<class VertData, class TriData>
Vec3d Mesh<VertData,TriData>::IsctProblem::computeCoords(Eptr e, Tptr t) const
{
    Empty3d::TriEdgeIn input;
    marshallArithmeticInput(input, e, t);
    Vec3d coords = Empty3d::coordsExact(input);
    return coords;
}

template<class VertData, class TriData>
Vec3d Mesh<VertData,TriData>::IsctProblem::computeCoords(
    Tptr t0, Tptr t1, Tptr t2
) const {
    Empty3d::TriTriTriIn input;
    marshallArithmeticInput(input, t0, t1, t2);
    Vec3d coords = Empty3d::coordsExact(input);
    return coords;
}


template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::fillOutVertData(
    GluePt glue, VertData &data
) {
    if(glue->split_type) { // manually inserted split point
        uint v0i = glue->e->verts[0]->ref;
        uint v1i = glue->e->verts[1]->ref;
        data.isctInterpolate(TopoCache::mesh->verts[v0i],
                             TopoCache::mesh->verts[v1i]);
    } else
    if(glue->edge_tri_type) { // edge-tri type
        IsctVertEdgeTriInput<VertData,TriData>      input;
        for(uint k=0; k<2; k++) {
            uint    vid                 = glue->e->verts[k]->ref;
                    input.e[k]          = &(TopoCache::mesh->verts[vid]);
        }
        for(uint k=0; k<3; k++) {
            uint    vid                 = glue->t[0]->verts[k]->ref;
                    input.t[k]          = &(TopoCache::mesh->verts[vid]);
        }
        data.isct(input);
    } else { // tri-tri-tri type
        IsctVertTriTriTriInput<VertData,TriData>    input;
        for(uint i=0; i<3; i++) {
          for(uint j=0; j<3; j++) {
            uint    vid                 = glue->t[i]->verts[j]->ref;
                    input.t[i][j]       = &(TopoCache::mesh->verts[vid]);
        }}
        data.isct(input);
    }
}

template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::subdivide_tri(
    uint t_piece_ref, uint t_parent_ref
) {
    SubdivideTriInput<VertData,TriData>     input;
                input.pt        = &(tris[t_parent_ref].data);
    for(uint k=0; k<3; k++) {
                input.pv[k]     = &(verts[tris[t_parent_ref].v[k]]);
                input.v[k]      = &(verts[tris[t_piece_ref].v[k]]);
    }
    tris[t_piece_ref].data.subdivide(input);
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::fillOutTriData(
    Tptr piece, Tptr parent
) {
    TopoCache::mesh->subdivide_tri(piece->ref, parent->ref);
}


template<class VertData, class TriData>
class Mesh<VertData,TriData>::IsctProblem::EdgeCache
{
public:
    EdgeCache(IsctProblem *ip) : iprob(ip), edges(ip->mesh->verts.size()) {}
    
    Eptr operator()(Vptr v0, Vptr v1) {
        uint i = v0->ref;
        uint j = v1->ref;
        if(i > j) std::swap(i,j);
        
        uint N = edges[i].size();
        for(uint k=0; k<N; k++)
            if(edges[i][k].vid == j)
                return edges[i][k].e;
        // if not existing, create it
        edges[i].push_back(EdgeEntry(j));
        Eptr e = edges[i][N].e = iprob->newEdge();
        e->verts[0] = v0;
        e->verts[1] = v1;
        v0->edges.push_back(e);
        v1->edges.push_back(e);
        
        return e;
    }
    
    // k = 0, 1, or 2
    Eptr getTriangleEdge(GTptr gt, uint k, Tptr big_tri)
    {
        GVptr   gv0             = gt->verts[(k+1)%3];
        GVptr   gv1             = gt->verts[(k+2)%3];
        Vptr    v0              = gv0->concrete;
        Vptr    v1              = gv1->concrete;
            // if neither of these are intersection points,
            // then this is a pre-existing edge...
        Eptr    e               = nullptr;
        if(typeid(gv0) == typeid(OVptr) &&
           typeid(gv1) == typeid(OVptr)
        ) {
            // search through edges of original triangle...
            for(uint c=0; c<3; c++) {
                Vptr corner0 = big_tri->verts[(c+1)%3];
                Vptr corner1 = big_tri->verts[(c+2)%3];
                if((corner0 == v0 && corner1 == v1) ||
                   (corner0 == v1 && corner1 == v0)) {
                    e   = big_tri->edges[c];
                }
            }
            ENSURE(e); // Yell if we didn't find an edge
        }
            // otherwise, we need to check the cache to find this edge
        else
        {
            e = operator()(v0, v1);
        }
        return e;
    }
    
    Eptr maybeEdge(GEptr ge)
    {
        uint i = ge->ends[0]->concrete->ref;
        uint j = ge->ends[1]->concrete->ref;
        if(i > j) std::swap(i,j);
        
        uint N = edges[i].size();
        for(uint k=0; k<N; k++)
            if(edges[i][k].vid == j)
                return edges[i][k].e;
        // if we can't find it
        return nullptr;
    }
    
private:
    struct EdgeEntry {
        EdgeEntry(uint id) : vid(id) {}
        EdgeEntry() {}
        uint vid;
        // things
        Eptr e;
    };
    
    IsctProblem *iprob;
    std::vector< ShortVec<EdgeEntry, 8> >   edges;
};

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::createRealPtFromGluePt(GluePt glue) {
    ENSURE(glue->copies.size() > 0);
    Vptr        v               = TopoCache::newVert();
    VertData    &data           = TopoCache::mesh->verts[v->ref];
                data.pos        = glue->copies[0]->coord;
                fillOutVertData(glue, data);
    for(IVptr iv : glue->copies)
                iv->concrete    = v;
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::createRealTriangles(
    Tprob tprob, EdgeCache &ecache
) {
    for(GTptr gt : tprob->gtris) {
        Tptr        t               = TopoCache::newTri();
                    gt->concrete    = t;
        Tri         &tri            = TopoCache::mesh->tris[t->ref];
        for(uint k=0; k<3; k++) {
            Vptr    v               = gt->verts[k]->concrete;
                    t->verts[k]     = v;
                    v->tris.push_back(t);
                    tri.v[k]        = v->ref;
            
            Eptr    e = ecache.getTriangleEdge(gt, k, tprob->the_tri);
                    e->tris.push_back(t);
                    t->edges[k] = e;
        }
                    fillOutTriData(t, tprob->the_tri);
    }
    // Once all the pieces are hooked up, let's kill the old triangle!
    TopoCache::deleteTri(tprob->the_tri);
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::resolveAllIntersections()
{
    // solve a subdivision problem in each triangle
    tprobs.for_each([&](Tprob tprob) {
        tprob->subdivide(this);
    });
    
    // now we have diced up triangles inside each triangle problem
    
    // Let's go through the glue points and create a new concrete
    // vertex object for each of these.
    glue_pts.for_each([&](GluePt glue) {
        createRealPtFromGluePt(glue);
    });
    
    EdgeCache ecache(this);
    
    // Now that we have concrete vertices plugged in, we can
    // go through the diced triangle pieces and create concrete triangles
    // for each of those.
    // Along the way, let's go ahead and hook up edges as appropriate
    tprobs.for_each([&](Tprob tprob) {
        createRealTriangles(tprob, ecache);
    });
    
    // mark all edges as normal by zero-ing out the data pointer
    TopoCache::edges.for_each([](Eptr e) {
        e->data = 0;
    });
    // then iterate over the edges formed by intersections
    // (i.e. those edges without the boundary flag set in each triangle)
    // and mark those by setting the data pointer
    iepool.for_each([&](IEptr ie) {
        // every ie must be non-boundary
        Eptr e = ecache.maybeEdge(ie);
        ENSURE(e);
        e->data = (void*)1;
    });
    sepool.for_each([&](SEptr se) {
        //if(se->boundary)    return; // continue
        Eptr e = ecache.maybeEdge(se);
        ENSURE(e);
        e->data = (void*)1;
    });
    
    // This basically takes care of everything EXCEPT one detail
    // *) The base mesh data structures still need to be compacted
    
    // This detail should be handled by the calling code...
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::dumpIsctPoints(
    std::vector<Vec3d> *points
) {
    points->resize(glue_pts.size());
    uint write = 0;
    glue_pts.for_each([&](GluePt glue) {
        ENSURE(glue->copies.size() > 0);
        IVptr       iv                  = glue->copies[0];
                    (*points)[write]    = iv->coord;
                    write++;
    });
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::dumpIsctEdges(
    std::vector< std::pair<Vec3d,Vec3d> > *edges
) {
    edges->clear();
    tprobs.for_each([&](Tprob tprob) {
        for(IEptr ie : tprob->iedges) {
            GVptr gv0 = ie->ends[0];
            GVptr gv1 = ie->ends[1];
            edges->push_back(std::make_pair(gv0->coord, gv1->coord));
        }
    });
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::testingComputeStaticIsctPoints(
    std::vector<Vec3d> *points
) {
    IsctProblem iproblem(this);
    
    iproblem.findIntersections();
    
    iproblem.dumpIsctPoints(points);
}


template<class VertData, class TriData>
void Mesh<VertData,TriData>::testingComputeStaticIsct(
    std::vector<Vec3d> *points,
    std::vector< std::pair<Vec3d,Vec3d> > *edges
) {
    IsctProblem iproblem(this);
    
    iproblem.findIntersections();
    
    iproblem.dumpIsctPoints(points);
    iproblem.dumpIsctEdges(edges);
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::resolveIntersections()
{
    IsctProblem iproblem(this);
    
    iproblem.findIntersections();
    
    iproblem.resolveAllIntersections();
    
    iproblem.commit();
    
    //iproblem.print();
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::isSelfIntersecting()
{
    IsctProblem iproblem(this);
    
    return iproblem.hasIntersections();
}



