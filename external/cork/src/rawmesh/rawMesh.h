// +-------------------------------------------------------------------------
// | rawMesh.h
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

#include <functional>

#include "vec.h"

#include <vector>

// vertex and triangle data must satisfy the following minimal specifications!
// These exact field names must be used!
struct MinimalVertexData
{
    Vec3d pos; // required for both RawMesh and Mesh
};

struct MinimalTriangleData
{
    // Vertex Ids: (only used in raw mesh form)
    int a, b, c;
};

// Raw mesh presents an exposed interface to a mesh.
// This allows for easy input/output of data to and from the more powerful
// Mesh data structure, as well as supporting more basic mesh applications
template<class VertData, class TriData>
struct RawMesh
{
    std::vector<VertData> vertices;
    std::vector<TriData> triangles;
};

// allow for data to move between different kinds of raw meshes
template<class VertDataOut, class TriDataOut,
         class VertDataIn,  class TriDataIn>
inline RawMesh<VertDataOut,TriDataOut> transduce(
    const RawMesh<VertDataIn,TriDataIn> &input,
    std::function<void(VertDataOut &, const VertDataIn &)> vertTransduce,
    std::function<void(TriDataOut  &, const TriDataIn  &)> triTransduce
);


#include "rawMesh.tpp"
