import sys
import os
from threading import Thread

import clang.cindex
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
        #        parameter_dec = [c for c in cursor.get_children() if c.kind == clang.cindex.CursorKind.PARM_DECL]
        #        print(parameter_dec, template_pars)
        #        print(cursor.get_num_template_arguments(), cursor.get_template_argument_type(0), cursor.get_template_argument_value(0), template_pars, parameter_dec)
        self.parameters = []
        self.parnames = []
        self.documentation = cursor.raw_comment


class Enum(object):
    def __init__(self, cursor):
        self.name = cursor.spelling
        self.annotations = get_annotations(cursor)
        self.access = cursor.access_specifier
        #        template_pars = [c.extent for c in cursor.get_children() if c.kind == clang.cindex.CursorKind.TEMPLATE_TYPE_PARAMETER]
        #        parameter_dec = [c for c in cursor.get_children() if c.kind == clang.cindex.CursorKind.PARM_DECL]
        #        print(parameter_dec, template_pars)
        #        print(cursor.get_num_template_arguments(), cursor.get_template_argument_type(0), cursor.get_template_argument_value(0), template_pars, parameter_dec)
        self.parameters = []
        self.parnames = []
        self.documentation = cursor.raw_comment

# class Class(object):
#    def __init__(self, cursor):
#        self.name = cursor.spelling
#        self.functions = []
#        self.annotations = get_annotations(cursor)

#        for c in cursor.get_children():
#            if (c.kind == clang.cindex.CursorKind.CXX_METHOD and
#                c.access_specifier == clang.cindex.AccessSpecifier.PUBLIC):
#                f = Function(c)
#                self.functions.append(f)

def find_namespace_node(c):
    if (c.kind == clang.cindex.CursorKind.NAMESPACE and c.spelling == "igl"):
        return c
    else:
        for child_node in c.get_children():
            return find_namespace_node(child_node)


def traverse(c, path, objects):
    if c.location.file and not c.location.file.name.endswith(path):
        return

    # print(c.kind, c.spelling)


    if c.kind == clang.cindex.CursorKind.TRANSLATION_UNIT or c.kind == clang.cindex.CursorKind.UNEXPOSED_DECL:
        # Ignore  other cursor kinds
        pass

    elif c.kind == clang.cindex.CursorKind.NAMESPACE:
        objects["namespaces"].append(c.spelling)
        # print("Namespace", c.spelling, c.get_children())
        pass

    elif c.kind == clang.cindex.CursorKind.FUNCTION_TEMPLATE:
        # print("Function Template", c.spelling, c.raw_comment)
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

    else:
        # print("Unknown", c.kind, c.spelling)
        pass

    for child_node in c.get_children():
        traverse(child_node, path, objects)


def parse(path):
    index = clang.cindex.Index.create()
    tu = index.parse(path, ['-x', 'c++', '-std=c++11', '-fparse-all-comments', '-DIGL_STATIC_LIBRARY'])
    # Clang can't parse files with missing definitions, add static library definition
    objects = {"functions": [], "enums": [], "namespaces": [], "classes": []}
    traverse(tu.cursor, path, objects)

    #    tpl = Template(filename='bind.mako')
    #    rendered = tpl.render(functions=functions)

    #    OUTPUT_DIR = 'generated'

    #    if not os.path.isdir(OUTPUT_DIR): os.mkdir(OUTPUT_DIR)

    #    with open("generated/{}.bind.cc".format(sys.argv[1]), "w") as f:
    #        f.write(rendered)
    return objects

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 parser.py <headerfile_path>")
        exit(-1)
    parse(sys.argv[1])
