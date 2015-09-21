#!/usr/bin/env python3
#
#  Syntax: mkdoc.py <path_of_header_files>
#
#  Extract documentation from C++ header files to use it in libiglPython bindings
#

import os, sys, glob

#http://stackoverflow.com/questions/3207219/how-to-list-all-files-of-a-directory-in-python
def get_filepaths(directory):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Syntax: %s <path_of_header_files>' % sys.argv[0])
        exit(-1)

    # Open two files, py_doc.h and py_doc.cpp
    fh = open('py_doc.h', 'w')
    fc = open('py_doc.cpp', 'w')

    # List all files in the given folder and subfolders
    base_path = sys.argv[1]
    full_file_paths = get_filepaths(base_path)

    # Add all the .h files
    for f in full_file_paths:
      if f.endswith(".h"):
        f_clean = f[len(base_path):]
        f_clean = f_clean.replace(base_path, "")
        f_clean = f_clean.replace(".h", "")
        f_clean = f_clean.replace("/", "_")
        f_clean = f_clean.replace("\\", "_")
        f_clean = f_clean.replace(" ", "_")
        f_clean = f_clean.replace(".", "_")

        #tmp = open(f, 'r', encoding="utf8")
        tmp_string = f.replace("../include/", "libigl/") # " " # tmp.read()
        tmp_string = "See " + tmp_string + " for the documentation."
        #tmp.close()

        h_string = "extern const char *__doc_" + f_clean + ";\n"
        cpp_string = "const char *__doc_" + f_clean + " = R\"igl_Qu8mg5v7(" + tmp_string + ")igl_Qu8mg5v7\";\n"

        fh.write(h_string)
        fc.write(cpp_string)

    # Close files
    fh.close()
    fc.close()
