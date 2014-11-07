// +-------------------------------------------------------------------------
// | fixint.h
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#pragma once

/*
 *
 *  fixint.H
 *
 *      Support for known size big integers
 *
 *      See the BitInt wrapper class towards the end of the file
 *  for the type you want to use.  Specifically you want to declare
 *
 *      BitInt<128>::Rep myInt;
 *
 *  to get an integer guaranteed to hold at least 128 bits.
 *  Provided the size of all integers used in this fashion
 *  is declared appropriately, the template magic will
 *  statically decide which specialized function to call.
 *
 */



#ifdef _WIN32
#include <mpir.h>
#else
#include <gmp.h>
#endif

#include <iostream>
#include <iomanip>
#include <string>

namespace FixInt {

// static assertion using templates.  Usage:
//      ASSERT_STATIC< expr_to_assert >::test();
template<bool test_val>
struct ASSERT_STATIC {};
template<>
struct ASSERT_STATIC<true> { static void test() {}; };


const static unsigned int LIMB_BIT_SIZE  = GMP_NUMB_BITS;
const static unsigned int SIGN_BIT_OFFSET = LIMB_BIT_SIZE-1;
// mask defining the position of the sign bit in a limb
#define LIMB_SIGN_MASK (mp_limb_t(1) << SIGN_BIT_OFFSET)
// function which tests sign and returns 1 if it's set; 0 otherwise
#define SIGN_BOOL(limbs,n) \
        (((limbs)[(n)-1] & LIMB_SIGN_MASK) >> SIGN_BIT_OFFSET)
// function which tests sign and returns -1 if it's set; 1 otherwise
#define SIGN_INT(limbs,n) (((limbs)[(n)-1] & LIMB_SIGN_MASK)? -1 : 1)
// pattern of all 1s filling up a limb
#define ONES_PATTERN (mp_limb_t(-1))
#define ZERO_PATTERN (mp_limb_t(0))
// function which tests sign and returns 
#define SIGN_LIMB(limbs,n) (((limbs)[(n)-1] & LIMB_SIGN_MASK)? \
                                ONES_PATTERN : ZERO_PATTERN)
#define BITS_TO_LIMBS(n) ((((n)-1)/ LIMB_BIT_SIZE)+1)

template<int Nlimbs>
class LimbInt {
public:
    mp_limb_t limbs[Nlimbs];
    
    inline LimbInt(int init) {
        ASSERT_STATIC<(sizeof(int) <= sizeof(mp_limb_t))>::test();
        limbs[0] = init;
        if(Nlimbs > 1) {
            mp_limb_t fill = SIGN_LIMB(limbs,1);
            for(int i=1; i<Nlimbs; i++)
                limbs[i] = fill;
        }
    }
    
    inline LimbInt() {}
};

template<int N>
std::ostream& operator<<(std::ostream &out, const LimbInt<N> &num)
{
    out << "[" << std::hex << num.limbs[0];
    for(int k=0; k<N; k++)
        out << ";" << std::hex << num.limbs[k];
    out << "]";
    return out;
}


template<int Nout, int Nin>
inline void promote(LimbInt<Nout> &out, const LimbInt<Nin> &in)
{
    ASSERT_STATIC<(Nout >= Nin)>::test();
    
    mpn_copyi(out.limbs, in.limbs, Nin);
    
    if(Nout > Nin) { // fill out the higher order bits...
        mp_limb_t fill = SIGN_LIMB(out.limbs, Nin);
        for(int i=Nin; i<Nout; i++)
            out.limbs[i] = fill;
    }
}


template<int Nout, int Nlhs, int Nrhs>
inline void add(LimbInt<Nout> &out,
                const LimbInt<Nlhs> &lhs,
                const LimbInt<Nrhs> &rhs)
{
    ASSERT_STATIC<(Nout >= Nlhs && Nout >= Nrhs)>::test();
    int Nmax = Nlhs;
    if(Nrhs > Nlhs) Nmax = Nrhs;
    mp_limb_t carry;
    
    if(Nlhs == Nrhs) {
        carry = mpn_add_n(out.limbs, lhs.limbs, rhs.limbs, Nlhs);
    } else if(Nlhs > Nrhs) {
        mp_limb_t rhs_is_neg = SIGN_BOOL(rhs.limbs, Nrhs);
        carry = mpn_add(out.limbs, lhs.limbs, Nlhs, rhs.limbs, Nrhs);
        mp_limb_t borrow = mpn_sub_1(out.limbs+Nrhs,
                                     out.limbs+Nrhs, Nlhs-Nrhs,
                                     rhs_is_neg);
        if(Nout > Nmax) carry = carry | (mp_limb_t(1) - borrow);
    } else { // Nrhs > Nlhs
        mp_limb_t lhs_is_neg = SIGN_BOOL(lhs.limbs, Nlhs);
        carry = mpn_add(out.limbs, rhs.limbs, Nrhs, lhs.limbs, Nlhs);
        mp_limb_t borrow = mpn_sub_1(out.limbs+Nlhs,
                                     out.limbs+Nlhs, Nrhs-Nlhs,
                                     lhs_is_neg);
        if(Nout > Nmax) carry = carry | (mp_limb_t(1) - borrow);
    }
    
    // fill out any new higher order bits
    if(Nout > Nmax) {
        mp_limb_t rhs_is_neg = SIGN_BOOL(rhs.limbs, Nrhs);
        mp_limb_t lhs_is_neg = SIGN_BOOL(lhs.limbs, Nlhs);
        mp_limb_t fill_bit = rhs_is_neg ^ lhs_is_neg ^ carry;
        mp_limb_t fill = (fill_bit)? ONES_PATTERN : ZERO_PATTERN;
        for(int i=Nmax; i<Nout; i++)
            out.limbs[i] = fill;
    }
}

template<int Nout, int Nlhs, int Nrhs>
inline void sub(LimbInt<Nout> &out,
                const LimbInt<Nlhs> &lhs,
                const LimbInt<Nrhs> &rhs)
{
    // for testing...
    LimbInt<Nrhs> tempright;
    mpn_neg(tempright.limbs, rhs.limbs, Nrhs);
    add(out, lhs, tempright);
    /*ASSERT_STATIC<(Nout >= Nlhs && Nout >= Nrhs)>::test();
    const mp_limb_t *left   = lhs.limbs;
    const mp_limb_t *right  = rhs.limbs;
    int Nleft               = Nlhs;
    int Nright              = Nrhs;
    if(Nrhs > Nlhs) {
        const mp_limb_t *temp = left;   left = right;       right = temp;
        int tempi = Nleft;              Nleft = Nright;     Nright = tempi;
    }
    
    mpn_copyi(out.limbs, left, Nleft);
    if(Nout != Nleft) { // make sure higher order limbs are filled
        mp_limb_t fill = SIGN_LIMB(left, Nleft);
        for(int i=Nleft; i<Nout; i++)
            out.limbs[i] = fill;
    }
        
    if(Nout == Nright)
        mpn_sub_n(out.limbs, out.limbs, right, Nout);
    else
        mpn_sub(out.limbs, out.limbs, Nout, right, Nright);
        
        
    
    ASSERT_STATIC<(Nout >= Nlhs && Nout >= Nrhs)>::test();
    int Nmax = Nlhs;
    if(Nrhs > Nlhs) Nmax = Nrhs;
    mp_limb_t carry;
    
    if(Nlhs == Nrhs) {
        carry = mpn_add_n(out.limbs, lhs.limbs, rhs.limbs, Nlhs);
    } else if(Nlhs > Nrhs) {
        mp_limb_t rhs_is_neg = SIGN_BOOL(rhs.limbs, Nrhs);
        carry = mpn_add(out.limbs, lhs.limbs, Nlhs, rhs.limbs, Nrhs);
        mp_limb_t borrow = mpn_sub_1(out.limbs+Nrhs,
                                     out.limbs+Nrhs, Nlhs-Nrhs,
                                     rhs_is_neg);
        if(Nout > Nmax) carry = carry | (mp_limb_t(1) - borrow);
    } else { // Nrhs > Nlhs
        mp_limb_t lhs_is_neg = SIGN_BOOL(lhs.limbs, Nlhs);
        carry = mpn_add(out.limbs, rhs.limbs, Nrhs, lhs.limbs, Nlhs);
        mp_limb_t borrow = mpn_sub_1(out.limbs+Nlhs,
                                     out.limbs+Nlhs, Nrhs-Nlhs,
                                     lhs_is_neg);
        if(Nout > Nmax) carry = carry | (mp_limb_t(1) - borrow);
    }
    
    // fill out any new higher order bits  (DOES THIS HAVE TO USE THE CARRY?)
    if(Nout > Nmax) {
        mp_limb_t rhs_is_neg = SIGN_BOOL(rhs.limbs, Nrhs);
        mp_limb_t lhs_is_neg = SIGN_BOOL(lhs.limbs, Nlhs);
        mp_limb_t fill_bit = rhs_is_neg ^ lhs_is_neg ^ carry;
        mp_limb_t fill = -fill_bit;
        for(int i=Nmax; i<Nout; i++)
            out.limbs[i] = fill;
    }*/
}

// be careful of negating a value which becomes simply itself again...
// You should leave yourself one bit of safety
// WARNING: I DON'T HANDLE THE 1000000... bit pattern here!
template<int Nout, int Nin>
inline void neg(LimbInt<Nout> &out,
                const LimbInt<Nin> &in)
{
    ASSERT_STATIC<(Nout >= Nin)>::test();
    
    mpn_neg(out.limbs, in.limbs, Nin);
    
    if(Nout > Nin) {
        mp_limb_t fill = SIGN_LIMB(out.limbs, Nin);
        for(int i=Nin; i<Nout; i++)
            out.limbs[i] = fill;
    }
}

template<int Nout, int Nlhs, int Nrhs>
inline void mul(LimbInt<Nout> &out,
                const LimbInt<Nlhs> &lhs,
                const LimbInt<Nrhs> &rhs)
{
    ASSERT_STATIC<(Nout >= Nlhs && Nout >= Nrhs)>::test();
    // handle the possibility that we need to use more space
    // than we have...
    mp_limb_t *res = out.limbs;
    LimbInt<Nlhs + Nrhs> tempresult;
    if(Nout < Nlhs + Nrhs)
        res = tempresult.limbs;
    
    // multiply
    if(Nlhs == Nrhs)
        mpn_mul_n(res, lhs.limbs, rhs.limbs, Nlhs);
    else if(Nlhs > Nrhs)
        mpn_mul(res, lhs.limbs, Nlhs, rhs.limbs, Nrhs);
    else // need to flip in order to satisfy calling condition...
        mpn_mul(res, rhs.limbs, Nrhs, lhs.limbs, Nlhs);
    
    mp_limb_t lhs_sign = SIGN_BOOL(lhs.limbs,Nlhs);
    mp_limb_t rhs_sign = SIGN_BOOL(rhs.limbs,Nrhs);
    
    mpn_submul_1((res+Nlhs), rhs.limbs, Nrhs, lhs_sign);
    mpn_submul_1((res+Nrhs), lhs.limbs, Nlhs, rhs_sign);
    
    // transfer large result if we had one...
    if(Nout < Nlhs + Nrhs)
        mpn_copyi(out.limbs, res, Nout);
    
    // if we have more limbs than needed for the multiply,
    // fill out the extra higher order limbs...
    if(Nout > Nlhs + Nrhs) {
        mp_limb_t fill = SIGN_LIMB(out.limbs, Nlhs+Nrhs);
        for(int i=Nlhs+Nrhs; i<Nout; i++)
            out.limbs[i] = fill;
    }
}

template<int N>
inline
int sign(const LimbInt<N> &in)
{
    bool nonzero = false;
    for(int i=0; i<N; i++) {
        bool limbnonzero = (in.limbs[i] != 0);
        nonzero = nonzero || limbnonzero;
    }
    return SIGN_INT(in.limbs,N) * int(nonzero);
}

//template<int N>
//inline
//double approximate(const LimbInt<N> &in)
//{
//    std::cout << in << std::endl;
//    // not the most efficient implementation, but it should work
//    LimbInt<N>  tmp;
//    double      sign;
//    if(SIGN_BOOL(in.limbs, N-1)) {
//        neg(tmp, in);
//        sign = -1.0;
//    } else {
//        tmp = in;
//        sign = 1.0;
//    }
//    std::cout << tmp << std::endl;
//    std::cout << typeid(mp_limb_t).name() << std::endl;
//    std::cout << double(in.limbs[2]) << std::endl;
//    double result = 0.0;
//    for(int i=0; i<N; i++) {
//        result += ldexp(double(in.limbs[i]), i * LIMB_BIT_SIZE);
//    }
//    return sign * result;
//}

template<int N>
inline
std::string toString(const LimbInt<N> &num)
{
    char cbuf[(N*LIMB_BIT_SIZE*3)/10 + 3];
    LimbInt<N> garbage = num;
    bool neg = SIGN_BOOL(num.limbs,N);
    if(neg)
        mpn_neg(garbage.limbs, garbage.limbs, N);
    
    int count = mpn_get_str(reinterpret_cast<unsigned char*>(cbuf),
                            10, garbage.limbs, N);
    
    std::string result = "";
    
    if(neg) result += '-';
    int i=0;
    for(;i<count-1; i++)
        if(cbuf[i] != char(0)) break;
    for(; i<count; i++)
        result += (cbuf[i] + '0');
    
    return result;
}


// In order to declare an integer in terms of # of bits,
// use the idiom
//      BitInt<128>::Rep my128bitNumberVariable;
template<int Nbits>
class BitInt {
public:
    typedef LimbInt<BITS_TO_LIMBS(Nbits)> Rep;
private:
    BitInt();
};

} // end namespace FIXINT
