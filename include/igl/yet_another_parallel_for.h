// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2019 [Hsien-Yu Meng] [spur398@gmail.com]
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/

#ifndef YET_ANOTHER_IGL_PARALLEL_FOR_H
#define YET_ANOTHER_IGL_PARALLEL_FOR_H
#include "igl_inline.h"
#include <omp.h>

namespace yc_igl{
    template<typename FunctionType, class ...Ts>
    void parallel_for(
                    size_t loop_size,
                    FunctionType&& func,
                    Ts&& ...ts)
    {
        size_t i = 0;
#pragma omp parallel for private(i) default(shared)
        for(i = 0; i < loop_size; i++)
        {
            std::forward<FunctionType>(func)(i, std::forward<Ts>(ts)...);
        }
    }


    template<typename PreFunctionType,
        typename FunctionType,
        typename AccumFunctionType>
    void parallel_for(
                      size_t loop_size,
                      size_t min_parallel,
                      PreFunctionType&& prep_func,
                      FunctionType&& func,
                      AccumFunctionType&& accum_func)
    {


        size_t sthc = std::thread::hardware_concurrency();
        size_t nthreads =
#ifdef IGL_PARALLEL_FOR_FORCE_SERIAL
            0;
#else
        (loop_size < min_parallel) ? 0 : (sthc == 0 ? 8 : sthc);
#endif
        //size_t slice = std::max(std::round((static_cast<double>(loop_size)+1)/static_cast<double>(nthreads)),1);
        size_t slice = std::ceil(static_cast<double>(loop_size) / static_cast<double>(nthreads));
        if(nthreads == 0)
        {
            std::forward<PreFunctionType>(prep_func)(1);
            for(size_t i = 0; i < loop_size; i++)
            {
                std::forward<FunctionType>(func)(i, 0);
                std::forward<AccumFunctionType>(accum_func)(0);
            }
        }
        else{
            std::forward<PreFunctionType>(prep_func)(nthreads);

            size_t thread_id = 0;
#pragma omp parallel for private(thread_id) default(shared)
            for(thread_id = 0; thread_id < nthreads; ++thread_id)
            {
                for(size_t iter_id = thread_id * slice;
                    iter_id < (thread_id + 1) * slice && iter_id < loop_size;
                    ++iter_id)
                {
                    std::forward<FunctionType>(func)(iter_id, thread_id);
                }
            }
#pragma omp barrier
            for(thread_id = 0; thread_id < nthreads;++thread_id)
            {
                std::forward<AccumFunctionType>(accum_func)(thread_id);
            }
        }
    }
}
#endif
