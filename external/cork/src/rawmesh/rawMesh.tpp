// +-------------------------------------------------------------------------
// | rawMesh.tpp
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


template<class VertDataOut, class TriDataOut,
         class VertDataIn,  class TriDataIn>
inline RawMesh<VertDataOut,TriDataOut> transduce(
    const RawMesh<VertDataIn,TriDataIn> &input,
    std::function<void(VertDataOut &, const VertDataIn &)> vertTransduce,
    std::function<void(TriDataOut  &, const TriDataIn  &)> triTransduce
) {
    RawMesh<VertDataOut, TriDataOut> output;
    
    uint nVert = input.vertices.size();
    uint nTri  = input.triangles.size();
    output.vertices.resize(nVert);
    output.triangles.resize(nTri);
    
    for(uint i=0; i<nVert; i++)
        vertTransduce(output.vertices[i], input.vertices[i]);
    for(uint i=0; i<nTri; i++)
        triTransduce(output.triangles[i], input.triangles[i]);
    
    return output;
}

