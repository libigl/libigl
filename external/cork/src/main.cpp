// +-------------------------------------------------------------------------
// | main.cpp
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

// This file contains a command line program that can be used
// to exercise Cork's functionality without having to write
// any code.

#include "files.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <sstream>
using std::stringstream;
using std::string;

using std::ostream;

#include "cork.h"


void file2corktrimesh(
    const Files::FileMesh &in, CorkTriMesh *out
) {
    out->n_vertices = in.vertices.size();
    out->n_triangles = in.triangles.size();
    
    out->triangles = new uint[(out->n_triangles) * 3];
    out->vertices  = new double[(out->n_vertices) * 3];
    
    for(uint i=0; i<out->n_triangles; i++) {
        (out->triangles)[3*i+0] = in.triangles[i].a;
        (out->triangles)[3*i+1] = in.triangles[i].b;
        (out->triangles)[3*i+2] = in.triangles[i].c;
    }
    
    for(uint i=0; i<out->n_vertices; i++) {
        (out->vertices)[3*i+0] = in.vertices[i].pos.x;
        (out->vertices)[3*i+1] = in.vertices[i].pos.y;
        (out->vertices)[3*i+2] = in.vertices[i].pos.z;
    }
}

void corktrimesh2file(
    CorkTriMesh in, Files::FileMesh &out
) {
    out.vertices.resize(in.n_vertices);
    out.triangles.resize(in.n_triangles);
    
    for(uint i=0; i<in.n_triangles; i++) {
        out.triangles[i].a = in.triangles[3*i+0];
        out.triangles[i].b = in.triangles[3*i+1];
        out.triangles[i].c = in.triangles[3*i+2];
    }
    
    for(uint i=0; i<in.n_vertices; i++) {
        out.vertices[i].pos.x = in.vertices[3*i+0];
        out.vertices[i].pos.y = in.vertices[3*i+1];
        out.vertices[i].pos.z = in.vertices[3*i+2];
    }
}

void loadMesh(string filename, CorkTriMesh *out)
{
    Files::FileMesh filemesh;
    
    if(Files::readTriMesh(filename, &filemesh) > 0) {
        cerr << "Unable to load in " << filename << endl;
        exit(1);
    }
    
    file2corktrimesh(filemesh, out);
}
void saveMesh(string filename, CorkTriMesh in)
{
    Files::FileMesh filemesh;
    
    corktrimesh2file(in, filemesh);
    
    if(Files::writeTriMesh(filename, &filemesh) > 0) {
        cerr << "Unable to write to " << filename << endl;
        exit(1);
    }
}



class CmdList {
public:
    CmdList();
    ~CmdList() {}
    
    void regCmd(
        string name,
        string helptxt,
        std::function< void(std::vector<string>::iterator &,
                            const std::vector<string>::iterator &) > body
    );
    
    void printHelp(ostream &out);
    void runCommands(std::vector<string>::iterator &arg_it,
                     const std::vector<string>::iterator &end_it);
    
private:
    struct Command {
        string name;    // e.g. "show" will be invoked with option "-show"
        string helptxt; // lines to be displayed
        std::function< void(std::vector<string>::iterator &,
                       const std::vector<string>::iterator &) > body;
    };
    std::vector<Command> commands;
};

CmdList::CmdList()
{
    regCmd("help",
    "-help                  show this help message",
    [this](std::vector<string>::iterator &,
           const std::vector<string>::iterator &) {
        printHelp(cout);
        exit(0);
    });
}

void CmdList::regCmd(
    string name,
    string helptxt,
    std::function< void(std::vector<string>::iterator &,
                        const std::vector<string>::iterator &) > body
) {
    Command cmd = {
        name,
        helptxt,
        body
    };
    commands.push_back(cmd);
}

void CmdList::printHelp(ostream &out)
{
    out <<
    "Welcome to Cork.  Usage:" << endl <<
    "  > cork [-command arg0 arg1 ... argn]*" << endl <<
    "for example," << endl <<
    "  > cork -union box0.off box1.off result.off" << endl <<
    "Options:" << endl;
    for(auto &cmd : commands)
        out << cmd.helptxt << endl;
    out << endl;
}

void CmdList::runCommands(std::vector<string>::iterator &arg_it,
                          const std::vector<string>::iterator &end_it)
{
    while(arg_it != end_it) {
        string arg_cmd = *arg_it;
        if(arg_cmd[0] != '-') {
            cerr << arg_cmd << endl;
            cerr << "All commands must begin with '-'" << endl;
            exit(1);
        }
        arg_cmd = arg_cmd.substr(1);
        arg_it++;
        
        bool found = true;
        for(auto &cmd : commands) {
            if(arg_cmd == cmd.name) {
                cmd.body(arg_it, end_it);
                found = true;
                break;
            }
        }
        if(!found) {
            cerr << "Command -" + arg_cmd + " is not recognized" << endl;
            exit(1);
        }
    }
}


std::function< void(
    std::vector<string>::iterator &,
    const std::vector<string>::iterator &
) >
genericBinaryOp(
    std::function< void(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out) >
        binop
) {
    return [binop]
    (std::vector<string>::iterator &args,
     const std::vector<string>::iterator &end) {
        // data...
        CorkTriMesh in0;
        CorkTriMesh in1;
        CorkTriMesh out;
        
        if(args == end) { cerr << "too few args" << endl; exit(1); }
        loadMesh(*args, &in0);
        args++;
        
        if(args == end) { cerr << "too few args" << endl; exit(1); }
        loadMesh(*args, &in1);
        args++;
        
        binop(in0, in1, &out);
        
        if(args == end) { cerr << "too few args" << endl; exit(1); }
        saveMesh(*args, out);
        args++;
        
        freeCorkTriMesh(&out);
        
        delete[] in0.vertices;
        delete[] in0.triangles;
        delete[] in1.vertices;
        delete[] in1.triangles;
    };
}


int main(int argc, char *argv[])
{
    initRand(); // that's useful
    
    if(argc < 2) {
        cout << "Please type 'cork -help' for instructions" << endl;
        exit(0);
    }
    
    // store arguments in a standard container
    std::vector<string> args(argc);
    for(int k=0; k<argc; k++) {
        args[k] = argv[k];
    }
    
    auto arg_it = args.begin();
    // *arg_it is the program name to begin with, so advance!
    arg_it++;
    
    CmdList cmds;
    
    // add cmds
    cmds.regCmd("solid",
    "-solid in              Determine whether the input mesh represents\n"
    "                       a solid object.  (aka. watertight) (technically\n"
    "                         solid == closed and non-self-intersecting)",
    [](std::vector<string>::iterator &args,
       const std::vector<string>::iterator &end) {
        CorkTriMesh in;
        if(args == end) { cerr << "too few args" << endl; exit(1); }
        string filename = *args;
        loadMesh(*args, &in);
        args++;
        
        bool solid = isSolid(in);
        cout << "The mesh " << filename << " is: " << endl;
        cout << "    " << ((solid)? "SOLID" : "NOT SOLID") << endl;
        
        delete[] in.vertices;
        delete[] in.triangles;
    });
    cmds.regCmd("union",
    "-union in0 in1 out     Compute the Boolean union of in0 and in1,\n"
    "                       and output the result",
    genericBinaryOp(computeUnion));
    cmds.regCmd("diff",
    "-diff in0 in1 out      Compute the Boolean difference of in0 and in1,\n"
    "                       and output the result",
    genericBinaryOp(computeDifference));
    cmds.regCmd("isct",
    "-isct in0 in1 out      Compute the Boolean intersection of in0 and in1,\n"
    "                       and output the result",
    genericBinaryOp(computeIntersection));
    cmds.regCmd("xor",
    "-xor in0 in1 out       Compute the Boolean XOR of in0 and in1,\n"
    "                       and output the result\n"
    "                       (aka. the symmetric difference)",
    genericBinaryOp(computeSymmetricDifference));
    cmds.regCmd("resolve",
    "-resolve in0 in1 out   Intersect the two meshes in0 and in1,\n"
    "                       and output the connected mesh with those\n"
    "                       intersections made explicit and connected",
    genericBinaryOp(resolveIntersections));
    
    
    cmds.runCommands(arg_it, args.end());
    
    return 0;
}









