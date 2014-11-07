// +-------------------------------------------------------------------------
// | aabvh.h
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

#include "bbox.h"

#include <stack>

// maximum leaf size
static const uint LEAF_SIZE = 8;

template<class GeomIdx>
struct GeomBlob
{
    BBox3d  bbox;
    Vec3d   point; // representative point, usually the box midpoint
    GeomIdx id;
};

template<class GeomIdx>
struct AABVHNode
{
    BBox3d                          bbox;
    AABVHNode                       *left;
    AABVHNode                       *right;
    ShortVec<uint, LEAF_SIZE>       blobids;
    inline bool isLeaf() const { return left == nullptr; }
};

template<class GeomIdx>
class AABVH
{
public:
    AABVH(const std::vector< GeomBlob<GeomIdx> > &geoms) :
        root(nullptr), blobs(geoms), tmpids(geoms.size())
    {
        ENSURE(blobs.size() > 0);
        
        for(uint k=0; k<tmpids.size(); k++)
            tmpids[k] = k;
        
        root = constructTree(0, tmpids.size(), 2);
    }
    ~AABVH() {}
    
    // query?
    inline void for_each_in_box(
        const BBox3d                        &bbox,
        std::function<void(GeomIdx idx)>    action
    ) {
        // do a recursive search and invoke the action at each
        // piece of geometry
        std::stack< AABVHNode<GeomIdx>* >  nodes;
        nodes.push(root);
        
        while(!nodes.empty()) {
            AABVHNode<GeomIdx> *node    = nodes.top();
                                        nodes.pop();
            
            // check bounding box isct
            if(!hasIsct(node->bbox, bbox))  continue;
            
            // otherwise...
            if(node->isLeaf()) {
                for(uint bid : node->blobids) {
                    if(hasIsct(bbox, blobs[bid].bbox))
                        action(blobs[bid].id);
                }
            } else {
                nodes.push(node->left);
                nodes.push(node->right);
            }
        }
    }
    
private:
    // process range of tmpids including begin, excluding end
    // last_dim provides a hint by saying which dimension a
    // split was last made along
    AABVHNode<GeomIdx>* constructTree(uint begin, uint end, uint last_dim)
    {
        ENSURE(end - begin > 0); // don't tell me to build a tree from nothing
        // base case
        if(end-begin <= LEAF_SIZE) {
            AABVHNode<GeomIdx> *node = node_pool.alloc();
            node->left = nullptr;
            node->blobids.resize(end-begin);
            for(uint k=0; k<end-begin; k++) {
                uint blobid = node->blobids[k] = tmpids[begin + k];
                node->bbox = convex(node->bbox, blobs[blobid].bbox);
            }
            return node;
        }
        // otherwise, let's try to split this geometry up
        
        uint dim = (last_dim+1)%3;
        uint mid = (begin + end) / 2;
        quickSelect(mid, begin, end, dim);
        
        // now recurse
        AABVHNode<GeomIdx> *node = node_pool.alloc();
        node->left = constructTree(begin, mid, dim);
        node->right = constructTree(mid, end, dim);
        node->bbox = convex(node->left->bbox, node->right->bbox);
        return node;
    }
    
    // precondition: begin <= select < end
    void quickSelect(uint select, uint begin, uint end, uint dim)
    {
        // NOTE: values equal to the pivot
        //       may appear on either side of the split
        if(end-1 == select)     return;
        
        // p(ivot)i(ndex) and p(ivot)v(alue)
        uint pi = randMod(end-begin) + begin;
        double pv = blobs[tmpids[pi]].point[dim];
        
        int front = begin;
        int back  = end-1;
        while(front < back) {
            if(blobs[tmpids[front]].point[dim] < pv) {
                front++;
            }
            else if(blobs[tmpids[back]].point[dim] > pv) {
                back--;
            }
            else {
                std::swap(tmpids[front], tmpids[back]);
                front++;
                back--;
            }
        }
        if(front == back && blobs[tmpids[front]].point[dim] <= pv) {
            front++;
        }
        
        if(select < uint(front)) {
            quickSelect(select, begin, front, dim);
        } else {
            quickSelect(select, front, end, dim);
        }
    }
private:
    AABVHNode<GeomIdx>                  *root;
    
    IterPool< AABVHNode<GeomIdx> >      node_pool;
    std::vector< GeomBlob<GeomIdx> >    blobs;
    std::vector<uint>                   tmpids; // used during construction
};























