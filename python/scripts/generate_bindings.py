# This file is part of libigl, a simple c++ geometry processing library.
#
# Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public License
# v. 2.0. If a copy of the MPL was not distributed with this file, You can
# obtain one at http://mozilla.org/MPL/2.0/.
#!/usr/bin/env python3
#
#  Syntax: generate_docstrings.py <path_to_c++_header_files> <path_to_python_files>
#
#  Extract documentation from C++ header files to use it in libiglPython bindings
#

import os, sys, glob
import pickle

import shutil
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

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.


def get_name_from_path(path, basepath, prefix, postfix):
    f_clean = path[len(basepath):]
    f_clean = f_clean.replace(basepath, "")
    f_clean = f_clean.replace(postfix, "")
    f_clean = f_clean.replace(prefix, "")
    f_clean = f_clean.replace("/", "_")
    f_clean = f_clean.replace("\\", "_")
    f_clean = f_clean.replace(" ", "_")
    f_clean = f_clean.replace(".", "_")
    return f_clean


def map_parameter_types(name, cpp_type, parsed_types, errors, enum_types):
    # TODO Replace with proper regex matching and derive types from templates, comment parsing, names in cpp files
    # CAUTION: This is work in progress mapping code to get a grip of the problem
    # Types to map
    #    const int dim -> const int& dim ?
    result = []

    if cpp_type.startswith("const"):
        result.append("const ")
        cpp_type = cpp_type[6:]  # Strip const part

    # Handle special types
    skip_parsing = False
    if cpp_type.startswith("MatY"):
        result.append("Eigen::SparseMatrix<double>&")
        skip_parsing = True
    if cpp_type.startswith("Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>"):
        result.append("Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>")
        skip_parsing = True
    if cpp_type == "std::vector<std::vector<Scalar> > &":
        result.append("std::vector<std::vector<double> > &")
        skip_parsing = True
    if cpp_type == "std::vector<std::vector<Index> > &":
        result.append("std::vector<std::vector<int> > &")
        skip_parsing = True
    for constant in enum_types:
        if cpp_type.endswith(constant):
            result.append(cpp_type)
            skip_parsing = True

    if len(parsed_types) == 0:
        errors.append("Empty typechain: %s" % cpp_type)
        if cpp_type == "int" or cpp_type == "bool" or cpp_type == "unsigned int":
            return cpp_type, True
        else:
            return cpp_type, False

    # print(parsed_types, cpp_type)
    if not skip_parsing:
        for i, t in enumerate(parsed_types):

            if t == "Eigen":
                result.append("Eigen::")
                continue
            if t == "std":
                result.append("std::")
                continue

            if t == "PlainObjectBase" or t == "MatrixBase":
                if name == "F":
                    result.append("MatrixXi&")
                elif name == "V":
                    result.append("MatrixXd&")
                else:
                    result.append("MatrixXd&")
                break
            if t == "MatrixXi" or t == "VectorXi":
                result.append("MatrixXi&")
                break
            if t == "MatrixXd" or t == "VectorXd":
                result.append("MatrixXd&")
                break
            if t == "SparseMatrix" and len(parsed_types) >= i + 2 and (
                    parsed_types[i + 1] == "Scalar" or parsed_types[i + 1] == "T"):
                result.append("SparseMatrix<double>&")
                break
            if t == "SparseVector" and len(parsed_types) >= i + 2 and (parsed_types[i + 1] == "Scalar" or parsed_types[
                    i + 1] == "T"):
                result.append("SparseMatrix<double>&")
                break

            if t == "bool" or t == "int" or t == "double" or t == "unsigned" or t == "string":
                if cpp_type.endswith("&"):
                    result.append(t + " &")
                else:
                    result.append(t)
                break

            else:
                errors.append("Unknown typechain: %s" % cpp_type)
                return cpp_type, False


    return "".join(result), True


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Syntax: %s <path_to_c++_files>' % sys.argv[0])
        exit(-1)

    errors = {"missing": [], "empty": [], "others": [], "incorrect": [], "render": [], "various": []}
    files = {"complete": [], "partial": [], "errors": [], "others": [], "empty": []}

    # List all files in the given folder and subfolders
    cpp_base_path = sys.argv[1]
    cpp_file_paths = get_filepaths(cpp_base_path)

    # Add all the .h filepaths to a dict
    print("Collecting cpp files for parsing...")
    mapping = {}
    cppmapping = {}
    for f in cpp_file_paths:
        if f.endswith(".h"):
            name = get_name_from_path(f, cpp_base_path, "", ".h")
            mapping[name] = f

        if f.endswith(".cpp"):
            name = get_name_from_path(f, cpp_base_path, "", ".cpp")
            cppmapping[name] = f

    # Add all python binding files to a list
    implemented_names = list(mapping.keys())  # ["point_mesh_squared_distance"]
    implemented_names.sort()
    single_postfix = ""
    single_prefix = ""

    # Create a list of all cpp header files
    files_to_parse = []
    cppfiles_to_parse = []
    for n in implemented_names:
        files_to_parse.append(mapping[n])

        if n not in cppmapping:
            errors["missing"].append("No cpp source file for function %s found." % n)
        else:
            cppfiles_to_parse.append(cppmapping[n])

    # Parse c++ header files
    print("Parsing header files...")
    load_headers = False
    if load_headers:
        with open("headers.dat", 'rb') as fs:
            dicts = pickle.load(fs)
    else:
        job_count = cpu_count()
        dicts = Parallel(n_jobs=job_count)(delayed(parse)(path) for path in files_to_parse)

    if not load_headers:
        print("Saving parsed header files...")
        with open("headers.dat", 'wb') as fs:
            pickle.dump(dicts, fs)

    # Not yet needed, as explicit template parsing does not seem to be supported in clang
    # Parse c++ source files
    # cppdicts = Parallel(n_jobs=job_count)(delayed(parse)(path) for path in cppfiles_to_parse)

    # Change directory to become independent of execution directory
    print("Generating directory tree for binding files...")
    path = os.path.dirname(__file__)
    if path != "":
        os.chdir(path)
    try:
        shutil.rmtree("generated")
    except:
        pass # Ignore missing generated directory
    os.makedirs("generated/complete")
    os.mkdir("generated/partial")

    print("Generating and writing binding files...")
    for idx, n in enumerate(implemented_names):
        d = dicts[idx]
        contained_elements = sum(map(lambda x: len(x), d.values()))

        # Skip files that don't contain functions/enums/classes
        if contained_elements == 0:
            errors["empty"].append("Function %s contains no parseable content in cpp header. Something might be wrong." % n)
            files["empty"].append(n)
            continue

        # Add functions with classes to others
        if len(d["classes"]) != 0 or len(d["structs"]) != 0:
            errors["others"].append("Function %s contains classes/structs in cpp header. Skipping" % n)
            files["others"].append(n)
            continue

        # Work on files that contain only functions/enums and namespaces
        if len(d["functions"]) + len(d["namespaces"]) + len(d["enums"]) == contained_elements:
            correct_functions = []
            incorrect_functions = []

            # Collect enums to generate binding files
            enums = []
            enum_types = []
            for e in d["enums"]:
                enums.append({"name": e.name, "namespaces": d["namespaces"], "constants": e.constants})
                enum_types.append(e.name)

            # Collect functions to generate binding files
            for f in d["functions"]:
                parameters = []
                correct_function = True
                f_errors = []
                for p in f.parameters:
                    typ, correct = map_parameter_types(p[0], p[1], p[2], f_errors, enum_types)
                    correct_function &= correct
                    parameters.append({"name": p[0], "type": typ})

                if correct_function and len(parameters) > 0: #TODO add constants like EPS
                    correct_functions.append({"parameters": parameters, "namespaces": d["namespaces"], "name": f.name})
                elif len(parameters) > 0:
                    incorrect_functions.append({"parameters": parameters, "namespaces": d["namespaces"], "name": f.name})
                    errors["incorrect"].append("Incorrect function in %s: %s, %s\n" % (n, f.name, ",".join(f_errors)))
                else:
                    errors["various"].append("Function without pars in %s: %s, %s\n" % (n, f.name, ","
                                                                                                     "".join(f_errors)))

            # Write binding files
            try:
                tpl = Template(filename='basic_function.mako')
                rendered = tpl.render(functions=correct_functions, enums=enums)
                tpl1 = Template(filename='basic_function.mako')
                rendered1 = tpl.render(functions=incorrect_functions, enums=enums)
                path = "generated/"
                if len(incorrect_functions) == 0 and (len(correct_functions) != 0 or len(enums) != 0):
                    path += "complete/"
                    with open(path + single_prefix + "py_" + n + ".cpp", 'w') as fs:
                        fs.write(rendered)
                    files["complete"].append(n)
                else:
                    path += "partial/"
                    with open(path + single_prefix + "py_" + n + ".cpp", 'w') as fs:
                        fs.write("// COMPLETE BINDINGS ========================\n")
                        fs.write(rendered)
                        fs.write("\n\n\n\n// INCOMPLETE BINDINGS ========================\n")
                        fs.write(rendered1)

                    if len(correct_functions) != 0:
                        files["partial"].append(n)
                    else:
                        files["errors"].append(n)

            except Exception as e:
                files["errors"].append(n)
                errors["render"].append("Template rendering failed:" + n + " " + str(correct_functions) + ", incorrect "
                                                                                                "functions are " + str(
                    incorrect_functions) + str(e) + "\n")

    print("Writing error and overview files...")
    with open("errors.txt" + single_postfix, 'w') as fs:
        l = list(errors.keys())
        l.sort()
        for k in l:
            fs.write("%s: %i \n" %(k, len(errors[k])))
            fs.writelines("\n".join(errors[k]))
            fs.write("\n\n\n")

    with open("files.txt" + single_postfix, 'w') as fs:
        l = list(files.keys())
        l.sort()
        for k in l:
            fs.write("%s: %i \n" %(k, len(files[k])))
            fs.writelines("\n".join(files[k]))
            fs.write("\n\n\n")
