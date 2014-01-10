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
#ifndef CNC_BASIC_ARRAYS_H
#define CNC_BASIC_ARRAYS_H

#include <NL/nl.h>
#include "cnc_kernels.h"


template <class T> class CNCArray1d {
public:
    typedef CNCArray1d<T> thisclass ;

    CNCArray1d(unsigned int size = 0, unsigned int alignment = 1) {
        data_ = NULL ;
        base_mem_ = NULL ;
        size_ = 0 ;
        allocate(size, alignment) ;
    } 

    inline ~CNCArray1d() { deallocate() ; }


    /** does not preserve previous values stored in the array */
    void allocate(unsigned int size, unsigned int alignment = 1) {
        deallocate() ;
        if(size != 0) {
            base_mem_ = (char*)malloc(size * sizeof(T) + alignment -1) ;
	    char* p = base_mem_ ;
	    // GMY 20090825 original: while(unsigned __int64(p) % alignment){p++;}
	    while (((unsigned long long) p) % alignment) { ++p; }
            data_ = (T*)p ;
            for(unsigned int i=0; i<size; i++) {
                // Direct call to the constructor, see dlist.h for more explanations.
                new(&data_[i])T() ;                    
            }
        }
        size_ = size ;
    }

    void set_all(const T& value) {
        for(unsigned int i=0; i<size_; i++) {
            data_[i] = value ;
        }
    }

    T& operator()(unsigned int i) {
        return data_[i] ;
    }

    const T& operator()(unsigned int i) const {
        return data_[i] ;
    }

    T& operator[](unsigned int index) {
        return data_[index] ;
    }

    const T& operator[](unsigned int index) const {
        return data_[index] ;
    }

    T& from_linear_index(unsigned int index) {
        return data_[index] ;
    }

    const T& from_linear_index(unsigned int index) const {
        return data_[index] ;
    }

    unsigned int size() const { return size_ ; }

	unsigned int alignment() const { return alignment_ ; }

    void clear() { allocate(0) ; }

    /** low-level access, for experts only. */
    const T* data() const { return data_ ; }

    /** low-level access, for experts only. */
    T* data() { return data_ ; }

    unsigned int mem_usage() const {
        return size_ * sizeof(T) + sizeof(thisclass) ;
    }

	void print () const {
		for ( unsigned int index = 0; index<size_; index++ ) {
			printf ( "array[%d]=%lg\n", index, data_[index] ) ;
		}
	}

protected:
    T* data_ ;
    unsigned int size_ ;
    char* base_mem_ ;
    unsigned int alignment_ ;

    void deallocate() {
        if(size_ != 0) {
            for(unsigned int i=0; i<size_; i++) {
                // direct call to the destructor
                data_[i].~T() ;
            }
            free(base_mem_) ;
            data_ = NULL ;
            base_mem_ = NULL ;
            size_ = 0 ;
        }
    }

private:
    CNCArray1d(const thisclass& rhs) ;
    thisclass& operator=(const thisclass& rhs) ;

} ;

typedef CNCArray1d<float>  floatArray  ;
typedef CNCArray1d<double> doubleArray ;


#endif

