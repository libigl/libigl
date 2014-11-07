// +-------------------------------------------------------------------------
// | off.cpp
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
#include "files.h"

#include <iostream>
#include <fstream>
using std::ifstream;
using std::ofstream;
using std::endl;


namespace Files {

using std::string;
using std::vector;

int readOFF(string filename, FileMesh *data)
{
    if(!data) return 1;
    
    ifstream in; // will close on exit from this read function
    in.open(filename.c_str());
    if(!in) return 1;
    
    // "OFF"
    string filetype;
    in >> filetype;
    if(filetype != "OFF") return 1;
    
    // counts of things
    int numvertices, numfaces, numedges;
    in >> numvertices >> numfaces >> numedges;
    if(!in) return 1;
    data->vertices.resize(numvertices);
    data->triangles.resize(numfaces);
    
    // vertex data
    for(auto &v : data->vertices) {
        Vec3d &p = v.pos;
        in >> p.x >> p.y >> p.z;
    }
    if(!in) return 1;
    
    // face data
    for(auto &tri : data->triangles) {
        int polysize;
        in >> polysize;
        if(polysize != 3)   return 1;
        
        in >> tri.a >> tri.b >> tri.c;
    }
    if(!in) return 1;
    
    return 0;
}

int writeOFF(string filename, FileMesh *data)
{
    if(!data) return 1;
    
    ofstream out;
    out.open(filename.c_str());
    if(!out) return 1;
    
    // "OFF"
    out << "OFF" << endl;
    
    // numvertices, numfaces, numedges=0
    int numvertices = data->vertices.size();
    int numfaces = data->triangles.size();
    out << numvertices << ' ' << numfaces << ' ' << 0 << endl;
    
    // vertex data
    for(const auto &v : data->vertices) {
        const Vec3d &p = v.pos;
        out << p.x << ' ' << p.y << ' ' << p.z << endl;
    }
    
    // face data
    for(const auto &tri : data->triangles) {
        out << "3 " << tri.a << ' ' << tri.b << ' ' << tri.c << endl;
    }
    if(!out) return 1;
    
    return 0;
}



} // end namespace Files
