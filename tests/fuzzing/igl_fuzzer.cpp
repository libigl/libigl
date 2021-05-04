// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2021 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.


#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <igl/MshLoader.h>
#include <iostream>
#include <fstream>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size){
        char *new_str = (char *)malloc(size+1);
        if (new_str == NULL){
                return 0;
        }
        memcpy(new_str, data, size);
        new_str[size] = '\0';
		
        std::ofstream mshfile;
        mshfile.open ("file.msh");
        mshfile << new_str;
        mshfile.close();

        igl::MshLoader msh_loader("example.msh");

        free(new_str);
        return 0;
}
