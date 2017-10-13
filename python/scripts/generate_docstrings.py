# This file is part of libigl, a simple c++ geometry processing library.
#
# Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public License
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.
#!/usr/bin/env python
#
#  Syntax: generate_docstrings.py <path to libigl C++ header_files> <path to python binding C++ files>
#
#  Extract documentation from C++ header files to use it in libiglPython bindings
#

import os, sys, glob
from joblib import Parallel, delayed
from multiprocessing import cpu_count
from mako.template import Template
from parser import parse


# http://stackoverflow.com/questions/3207219/how-to-list-all-files-of-a-directory-in-python
def get_filepaths(directory):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.
    root_file_paths = []

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

            if root.endswith(directory): # Add only the files in the root directory
                root_file_paths.append(filepath)

    return file_paths, root_file_paths  # file_paths contains all file paths, core_file_paths only the ones in <directory>


def get_name_from_path(path, basepath, prefix, postfix):
    f_clean = os.path.relpath(path, basepath)
    f_clean = f_clean.replace(postfix, "")
    f_clean = f_clean.replace(prefix, "")
    f_clean = f_clean.replace("/", "_")
    f_clean = f_clean.replace("\\", "_")
    f_clean = f_clean.replace(" ", "_")
    f_clean = f_clean.replace(".", "_")
    return f_clean


if __name__ == '__main__':

    if len(sys.argv) != 3:
        print('Syntax: %s generate_docstrings.py <path to libigl C++ header_files> <path to python binding C++ files>' % sys.argv[0])
        exit(-1)

    # List all files in the given folder and subfolders
    cpp_base_path = sys.argv[1]
    py_base_path = sys.argv[2]
    cpp_file_paths, cpp_root_file_paths = get_filepaths(cpp_base_path)
    py_file_paths, py_root_file_paths = get_filepaths(py_base_path)

    # Add all the .h filepaths to a dict
    mapping = {}
    for f in cpp_file_paths:
        if f.endswith(".h"):
            name = get_name_from_path(f, cpp_base_path, "", ".h")
            mapping[name] = f

    # Add all python binding files to a list
    implemented_names = []
    core_implemented_names = []
    for f in py_file_paths:
        if f.endswith(".cpp"):
            name = get_name_from_path(f, py_base_path, "py_", ".cpp")
            implemented_names.append(name)
            if f in py_root_file_paths:
                core_implemented_names.append(name)

    implemented_names.sort()
    core_implemented_names.sort()

    # Create a list of cpp header files for which a python binding file exists
    files_to_parse = []
    for n in implemented_names:
        if n not in mapping:
            print("No cpp header file for python function %s found." % n)
            continue
        files_to_parse.append(mapping[n])
        # print(mapping[n])

    # Parse c++ header files
    job_count = cpu_count()
    dicts = Parallel(n_jobs=job_count)(delayed(parse)(path) for path in files_to_parse)

    hpplines = []
    cpplines = []

    for idx, n in enumerate(implemented_names):
        d = dicts[idx]
        contained_elements = sum(map(lambda x: len(x), d.values()))
        # Check for files that don't contain functions/enums/classes
        if contained_elements == 0:
            print("Function %s contains no parseable content in cpp header. Something might be wrong." % n)
            continue
        else:
            names = []
            namespaces = "_".join(d["namespaces"])  # Assumption that all entities lie in deepest namespace
            for f in d["functions"]:
                h_string = "extern const char *__doc_" + namespaces + "_" + f.name + ";\n"
                docu_string = "See " + f.name + " for the documentation."
                if f.documentation:
                    docu_string = f.documentation
                cpp_string = "const char *__doc_" + namespaces + "_" + f.name + " = R\"igl_Qu8mg5v7(" + docu_string + ")igl_Qu8mg5v7\";\n"

                if f.name not in names:  # Prevent multiple additions of declarations, TODO: Possible fix is to merge comments and add them to all functions
                    hpplines.append(h_string)
                    cpplines.append(cpp_string)
                names.append(f.name)

    # Change directory to become independent of execution directory
    path = os.path.dirname(__file__)
    if path != "":
        os.chdir(path)

    # Update the two files py_doc.h and py_doc.cpp
    with open('../py_doc.h', 'w') as fh:
        fh.writelines(hpplines)
    with open('../py_doc.cpp', 'w') as fc:
        fc.writelines(cpplines)

    # Write python_shared_cpp file
    tpl = Template(filename='python_shared.mako')
    rendered = tpl.render(functions=implemented_names)
    with open("../python_shared.cpp", 'w') as fs:
        fs.write(rendered)

    # Write py_igl_cpp file with all core library files
    tpl = Template(filename='py_igl.mako')
    rendered = tpl.render(functions=core_implemented_names)
    with open("../py_igl.cpp", 'w') as fs:
        fs.write(rendered)
