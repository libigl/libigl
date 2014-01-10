/*
 *  Copyright (c) 2004-2010, Bruno Levy
 *  All rights reserved.
 *
 *  CNC: Concurrent Number Cruncher, original code by Luc Buatois
 *  Copyright (C) 2008-2010 GOCAD/ASGA, INRIA/ALICE
 *
 *  Sparse matrix-vector multiplication (SpMV) CUDA kernels based on code
 *  by Nathan Bell and Michael Garland at NVIDIA.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef CNC_TIMER_H
#define CNC_TIMER_H

#ifdef OS_WIN
#include <windows.h>
typedef __int64 i64 ;



class CNCTimer {

public:
	CNCTimer () ;
	double GetTime () ;
	double GetElapsedTime ( double old_time ) ;

private:
	i64 _freq ;
	i64 _clocks ;
};


CNCTimer::CNCTimer () : _clocks(0) {
	QueryPerformanceFrequency((LARGE_INTEGER *)&_freq);
}


double CNCTimer::GetTime () {
    QueryPerformanceCounter((LARGE_INTEGER *)&_clocks);
	return (double)_clocks / (double)_freq;
}


double CNCTimer::GetElapsedTime ( double old_time ) {
    QueryPerformanceCounter((LARGE_INTEGER *)&_clocks);
	return ((double)_clocks / (double)_freq - old_time) ;
}

#endif


#ifdef OS_LINUX
typedef long long i64;
typedef long long LARGE_INTEGER;


class CNCTimer {

public:
    CNCTimer () ;
    double GetTime () ;
    double GetElapsedTime ( double old_time ) ;

private:
	i64 _freq ;
	i64 _clocks ;

	cudaEvent_t start;
	cudaEvent_t stop;
};


CNCTimer::CNCTimer () : _clocks(0) {
    //QueryPerformanceFrequency((LARGE_INTEGER *)&_freq);
    //    lpFrequency
    //    [out] Pointer to a variable that receives the current 
    //    performance-counter frequency, in counts per second. 
    //    If the installed hardware does not support a high-resolution 
    //    performance counter, this parameter can be zero. 
    _freq = 1000000; 
    cudaEventCreate(&start);
    cudaEventRecord(start, 0);
}


double CNCTimer::GetTime () {
    float elaps;   
    cudaEventCreate(&stop);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    /* return milliseconds from start to stop event */
    cudaEventElapsedTime(&elaps, start, stop);
    cudaEventDestroy(stop);
    _clocks = (long long)(elaps * 1000.0); // store microseconds
    return (double)_clocks / (double)_freq; // returns seconds
}


double CNCTimer::GetElapsedTime ( double old_time ) {
    CNCTimer::GetTime();
    return ((double)_clocks / (double)_freq - old_time) ; // returns seconds
}


#endif



#endif // CNC_TIMER_H


