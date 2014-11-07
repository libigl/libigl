// +-------------------------------------------------------------------------
// | ifs.cpp
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
using std::cout;
using std::endl;
#include <ios>


namespace Files {

using std::string;

// private module data
namespace {
    typedef float float32;
    typedef unsigned int uint32;
    
    bool swapBytes;
    bool writeFlipping = false;
}

// All functions return true
// if there are no errors

inline void flip(char *bytes, int num)
{
    char temp;
    for(int i=0; i<num/2; i++) {
        temp = bytes[i];
        bytes[i] = bytes[num-i-1];
        bytes[num-i-1] = temp;
    }
}

inline bool read_uint32(ifstream &in, uint32 &data)
{
    union {
        uint32 val;
        char raw[4];
    };
    
    in.read(raw, 4);
    if(!in)
        return false;
    if(swapBytes)
        flip(raw,4);
    data = val;
    return true;
}

inline bool write_uint32(ofstream &out, uint32 data)
{
    union {
        uint32 val;
        char raw[4];
    };
    
    val = data;
    if(swapBytes)
        flip(raw,4);
    out.write(raw,4);
    if(!out)
        return false;
    return true;
}

inline bool read_float32(ifstream &in, float32 &data)
{
    union {
        float32 val;
        char raw[4];
    };
    
    in.read(raw, 4);
    if(!in)
        return false;
    if(swapBytes)
        flip(raw,4);
    data = val;
    return true;
}

inline bool write_float32(ofstream &out, float32 data)
{
    union {
        float32 val;
        char raw[4];
    };
    
    val = data;
    if(swapBytes)
        flip(raw,4);
    out.write(raw,4);
    if(!out)
        return false;
    return true;
}

inline bool read_vertex(ifstream &in, Vec3d &data)
{
    Vec3f temp;
    bool result =   read_float32(in, temp.x) &&
                    read_float32(in, temp.y) &&
                    read_float32(in, temp.z);
    if(result)
        data = temp; // coercion from float to double here
    return result;
}

inline bool write_vertex(ofstream &out, const Vec3d &data)
{
    // data is coerced from double to float automatically here
    bool result =   write_float32(out, data.x) &&
                    write_float32(out, data.y) &&
                    write_float32(out, data.z);
    return result;
}

inline bool read_triangle(ifstream &in, FileTriangle &data)
{
    uint32 a, b, c;
    if(!read_uint32(in, a))     return false;
    if(!read_uint32(in, b))     return false;
    if(!read_uint32(in, c))     return false;
    data.a = a;
    data.b = b;
    data.c = c;
                                return true;
}

inline bool write_triangle(ofstream &out, const FileTriangle &data)
{
    return  write_uint32(out, data.a) &&
            write_uint32(out, data.b) &&
            write_uint32(out, data.c);
}

/*inline bool read_texturecoord(ifstream &in, TextureCoord &data)
{
    return  read_float32(in, data.u) &&
            read_float32(in, data.v);
}

inline bool write_texturecoord(ofstream &out, const TextureCoord &data)
{
    return  write_float32(out, data.u) &&
            write_float32(out, data.v);;
}*/

// strings are formatted in a minorly funky way in IFS files
// first comes a uint32 with the length of the string
// folowed by the actual string (including null termination)

// return true if the provided string was successfully read
// do not advance the file pointer
bool checkString(ifstream &in, string toCheck)
{
    uint32 length = toCheck.length()+1; // +1 for null-terminator
    uint32 read_length;
    if(!read_uint32(in, read_length)) return false;
    if(read_length != length) return false;
    
    char *buf = new char[length];
    in.read(buf, length);
    if(!in) {
        delete[] buf;
        return false;
    }
    
    // rollback on having read the length
    in.seekg(in.tellg()-(std::streamoff)(4+length));
    bool result = (toCheck.compare(buf) == 0);
    delete[] buf;
    return result;
}

bool readString(ifstream &in, string &data)
{
    uint32 length;
    if(!read_uint32(in, length)) return false;
    
    char *buf = new char[length];
    in.read(buf, length);
    if(!in) {
        delete[] buf;
        return false;
    }
    
    data = string(buf);
    delete[] buf;
    return true;
}

bool writeString(ofstream &out, const string &data)
{
    uint32 length = data.length()+1;
    if(!write_uint32(out, length)) return false;
    
    out.write(data.c_str(), length);
    if(!out) {
        return false;
    }
    
    return true;
}

int readIFS(string filename, FileMesh *data)
{
    string junk_str;
    swapBytes = false;
    if(!data) return 1;
    
    ifstream in;
    in.open(filename.c_str(), std::ios::in | std::ios::binary);
    if(!in) return 1;
    
    // "IFS"
    if(!checkString(in,"IFS")) {
        swapBytes = true;
        if(!checkString(in,"IFS")) {
            swapBytes = false;
            return 1;
        }
    }
    if(!readString(in,junk_str)) return 1;
    
    // version #
    float32 version;
    bool    hasTextures;
    if(!read_float32(in,version)) return 1;
    if(version == 1.0) {
        hasTextures = false;
    } else if(version == 1.1) {
        hasTextures = true;
    } else {
        return 1;
    }
    
    // model name
    if(!readString(in,junk_str)) return 1;
    //string name = junk_str;
    
    // VERTICES
    if(!checkString(in,"VERTICES")) return 1;
    if(!readString(in,junk_str)) return 1;
    uint32 num_vertices;
    if(!read_uint32(in, num_vertices)) return 1;
    data->vertices.resize(num_vertices);
    for(auto &v : data->vertices) {
        if(!read_vertex(in, v.pos)) return 1;
    }
    
    // TRIANGLES
    if(!checkString(in,"TRIANGLES")) return 1;
    if(!readString(in,junk_str)) return 1;
    uint32 num_tris;
    if(!read_uint32(in, num_tris)) return 1;
    data->triangles.resize(num_tris);
    for(auto &tri : data->triangles) {
        if(!read_triangle(in, tri)) return 1;
    }
    
    // TEXTURECOORD
    if(hasTextures) {
        // noop
    }
    
    return 0;
}

int writeIFS(string filename, FileMesh *data)
{
    swapBytes = writeFlipping;
    string junk_str;
    if(!data) return 1;
    
    ofstream out;
    out.open(filename.c_str(), std::ios::out | std::ios::binary);
    if(!out) return 1;
    
    // "IFS"
    if(!writeString(out,"IFS")) return 1;
    // version #
    float32 version = 1.0; // i.e. no textures
    if(!write_float32(out,version)) return 1;
    // model name
    if(!writeString(out,"no_name")) return 1;
    
    // VERTICES
    if(!writeString(out,"VERTICES")) return 1;
    uint32 num_vertices = data->vertices.size();
    if(!write_uint32(out, num_vertices)) return 1;
    for(const auto &v : data->vertices) {
        if(!write_vertex(out, v.pos)) return 1;
    }
    
    // TRIANGLES
    if(!writeString(out,"TRIANGLES")) return 1;
    uint32 num_tris = data->triangles.size();
    if(!write_uint32(out, num_tris)) return 1;
    for(const auto &tri : data->triangles) {
        if(!write_triangle(out, tri)) return 1;
    }
    
    // TEXTURECOORD
    // not used
    
    return 0;
}

} // end namespace Files
