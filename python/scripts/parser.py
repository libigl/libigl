# This file is part of libigl, a simple c++ geometry processing library.
#
# Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public License
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.
import sys
import os
from threading import Thread

import clang.cindex
import ccsyspath
import itertools
from mako.template import Template







def get_annotations(node):
    return [c.displayname for c in node.get_children()
            if c.kind == clang.cindex.CursorKind.ANNOTATE_ATTR]


class Function(object):
    def __init__(self, cursor):
        self.name = cursor.spelling
        self.annotations = get_annotations(cursor)
        self.access = cursor.access_specifier
#        template_pars = [c.extent for c in cursor.get_children() if c.kind == clang.cindex.CursorKind.TEMPLATE_TYPE_PARAMETER]
        parameter_dec = [c for c in cursor.get_children() if c.kind == clang.cindex.CursorKind.PARM_DECL]
        
        parameters = []
        for p in parameter_dec:
            children = []
            for c in p.get_children():
#                print(c.spelling)
                children.append(c.spelling)
            parameters.append((p.spelling, p.type.spelling, children))

        self.parameters = parameters
        self.documentation = cursor.raw_comment


class Enum(object):
    def __init__(self, cursor):
        self.name = cursor.spelling
        self.constants = [c.spelling for c in cursor.get_children() if c.kind ==
                       clang.cindex.CursorKind.ENUM_CONSTANT_DECL]
        self.documentation = cursor.raw_comment

class Class(object):
    def __init__(self, cursor):
        self.name = cursor.spelling
#        self.functions = []
        self.annotations = get_annotations(cursor)

#        for c in cursor.get_children():
#            if (c.kind == clang.cindex.CursorKind.CXX_METHOD and
#                c.access_specifier == clang.cindex.AccessSpecifier.PUBLIC):
#                f = Function(c)
#                self.functions.append(f)



def traverse(c, path, objects):
    if c.location.file and not c.location.file.name.endswith(path):
        return

    if c.spelling == "PARULA_COLOR_MAP": # Fix to prevent python stack overflow from infinite recursion
        return

#    print(c.kind, c.spelling)


    if c.kind == clang.cindex.CursorKind.TRANSLATION_UNIT or c.kind == clang.cindex.CursorKind.UNEXPOSED_DECL:
        # Ignore  other cursor kinds
        pass

    elif c.kind == clang.cindex.CursorKind.NAMESPACE:
        objects["namespaces"].append(c.spelling)
        # print("Namespace", c.spelling, c.get_children())
        pass

    elif c.kind == clang.cindex.CursorKind.FUNCTION_TEMPLATE:
#        print("Function Template", c.spelling, c.raw_comment)
        objects["functions"].append(Function(c))
        return

    elif c.kind == clang.cindex.CursorKind.FUNCTION_DECL:
        # print("FUNCTION_DECL", c.spelling, c.raw_comment)
        objects["functions"].append(Function(c))
        return

    elif c.kind == clang.cindex.CursorKind.ENUM_DECL:
        # print("ENUM_DECL", c.spelling, c.raw_comment)
        objects["enums"].append(Enum(c))
        return

    elif c.kind == clang.cindex.CursorKind.CLASS_DECL:
        objects["classes"].append(Class(c))
        return

    elif c.kind == clang.cindex.CursorKind.CLASS_TEMPLATE:
        objects["classes"].append(Class(c))
        return

    elif c.kind == clang.cindex.CursorKind.STRUCT_DECL:
        objects["structs"].append(Class(c))
        return

    else:
        # print("Unknown", c.kind, c.spelling)
        pass

    for child_node in c.get_children():
        traverse(child_node, path, objects)


def parse(path):
    index = clang.cindex.Index.create()
    # Clang can't parse files with missing definitions, add static library definition or not?
    args = ['-x', 'c++', '-std=c++11', '-fparse-all-comments', '-DIGL_STATIC_LIBRARY']
    args.append('-I/usr/include/eigen3/') # TODO Properly add all needed includes
    syspath = ccsyspath.system_include_paths('clang++') # Add the system libraries
    incargs = [(b'-I' + inc).decode("utf-8") for inc in syspath]
    args.extend(incargs)

    tu = index.parse(path, args)
    objects = {"functions": [], "enums": [], "namespaces": [], "classes": [], "structs": []}
    traverse(tu.cursor, path, objects)
    return objects

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 parser.py <headerfile_path>")
        exit(-1)
    parse(sys.argv[1])
