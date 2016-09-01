import pyigl as igl
import sys


TUTORIAL_SHARED_PATH = "../../tutorial/shared/"

def check_dependencies(deps):
    available = [hasattr(igl, m) for m in deps]
    all_available = True
    for i, d in enumerate(available):
        if not d:
            all_available = False
            print("The libigl python bindings were compiled without %s support. Please recompile with the CMAKE flag LIBIGL_WITH_%s." %(deps[i], deps[i].upper()))

    if not all_available:
        sys.exit(-1)


def print_usage(key_dict):
    print("Usage:")
    for k in key_dict.keys():
        print("%s : %s" %(k, key_dict[k]))
