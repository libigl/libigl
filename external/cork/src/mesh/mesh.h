// +-------------------------------------------------------------------------
// | mesh.h
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

// SIMPLE USAGE:
//  In order to get standard template inclusion behavior/usage
//  just include this file.  Then the entire template code
//  will be included in the compilation unit.
// ADVANCED USAGE:
//  Only include "mesh.decl.h" where-ever you would normally include a
//  header file.  This will avoid including the implementation code in
//  the current compilation unit.
//  Then, create a seperate cpp file which includes "mesh.h" and
//  explicitly instantiates the template with the desired template
//  parameters.
//  By following this scheme, you can prevent re-compiling the entire
//  template implementation in every usage compilation unit and every
//  time those compilation units are recompiled during development.

// Declaration
#include "mesh.decl.h"

// Implementation
#include "mesh.tpp"
#include "mesh.remesh.tpp"
#include "mesh.isct.tpp"
#include "mesh.bool.tpp"



