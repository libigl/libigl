// +-------------------------------------------------------------------------
// | unionFind.h
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

#include "prelude.h"

class UnionFind {
public:
    UnionFind(uint N) : ids(N), rank(N, 0) {
        for(uint i=0; i<N; i++)
            ids[i] = i;
    }
    
    uint find(uint i) {
        uint id = i;
        while(ids[id] != id)
            id = ids[id];
        ids[i] = id; // path compression optimization
        return id;
    }
    uint unionIds(uint i, uint j) {
        uint iid = find(i);
        uint jid = find(j);
        if(iid == jid) // no need to merge the same tree with itself
            return iid;
        // simple implementation (not used)
        // return ids[iid] = jid;
        // instead we attempt to rank balance
        if(rank[iid] > rank[jid]) {
            return ids[jid] = iid;
        } else if(rank[iid] < rank[jid]) {
            return ids[iid] = jid;
        } else { // equal ranks
            rank[jid]++;
            return ids[jid] = iid;
        }
    }
    std::vector<uint> dump() {
        std::vector<uint> result(ids.size());
        for(uint i=0; i<ids.size(); i++)
            result[i] = find(i);
        return result;
    }
    
private:
    std::vector<uint> ids;
    std::vector<uint> rank;
};
