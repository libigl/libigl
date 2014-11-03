// +-------------------------------------------------------------------------
// | files.cpp
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
using std::cout;
using std::cerr;
using std::endl;

namespace Files {

using std::string;

// private module data
namespace {
}

int readTriMesh(string filename, FileMesh *mesh)
{
    int lastdot = filename.find_last_of('.');
    if(lastdot < 0) return 1;
    string suffix = filename.substr(lastdot, filename.length()-lastdot);
    if(suffix == ".ifs")
        return readIFS(filename, mesh);
    else if (suffix == ".off")
        return readOFF(filename, mesh);
    else
        return 1;
}

int writeTriMesh(string filename, FileMesh *mesh)
{
    int lastdot = filename.find_last_of('.');
    if(lastdot < 0) return 1;
    string suffix = filename.substr(lastdot, filename.length()-lastdot);
    if(suffix == ".ifs")
        return writeIFS(filename, mesh);
    else if (suffix == ".off")
        return writeOFF(filename, mesh);
    else
        return 1;
}

} // end namespace Files
