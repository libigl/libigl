// +-------------------------------------------------------------------------
// | empty3d.cpp
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
#include "empty3d.h"

#include "ext4.h"
#include "absext4.h"
#include "fixext4.h"
#include "gmpext4.h"

#include "quantization.h"

#include <cfloat>

namespace Empty3d {

// externalized counters...
int degeneracy_count = 0;
int exact_count = 0;
int callcount = 0;

using namespace Ext4;
using namespace AbsExt4;
using namespace FixExt4;
using namespace GMPExt4;

void toExt(Ext4_1 &out, const Vec3d &in)
{
    out.e0 = in.x;
    out.e1 = in.y;
    out.e2 = in.z;
    out.e3 = 1.0f;
}

void toAbsExt(AbsExt4_1 &out, const Vec3d &in)
{
    Vec3d vec = abs(in);
    out.e0 = vec.x;
    out.e1 = vec.y;
    out.e2 = vec.z;
    out.e3 = 1.0f;
}

void toVec3d(Vec3d &out, const Ext4_1 &in)
{
    // Warning: beware of division by zero!
    out.x = in.e0 / in.e3;
    out.y = in.e1 / in.e3;
    out.z = in.e2 / in.e3;
}

const static int IN_BITS = Quantization::BITS + 1; // +1 for sign bit

void toFixExt(FixExt4_1<IN_BITS> &out, const Vec3d &in)
{
    out.e0 = BitInt<IN_BITS>::Rep(Quantization::quantize2int(in[0]));
    out.e1 = BitInt<IN_BITS>::Rep(Quantization::quantize2int(in[1]));
    out.e2 = BitInt<IN_BITS>::Rep(Quantization::quantize2int(in[2]));
    out.e3 = BitInt<IN_BITS>::Rep(1);
}

void toGmpExt(GmpExt4_1 &out, const Vec3d &in)
{
    out.e0 = Quantization::quantize2int(in.x);
    out.e1 = Quantization::quantize2int(in.y);
    out.e2 = Quantization::quantize2int(in.z);
    out.e3 = 1;
}

void toVec3d(Vec3d &out, const GmpExt4_1 &in)
{
    Vec4d tmp;
    tmp.x = in.e0.get_d();
    tmp.y = in.e1.get_d();
    tmp.z = in.e2.get_d();
    tmp.w = in.e3.get_d();
    tmp /= tmp.w;
    for(uint k=0; k<3; k++)
        out.v[k] = Quantization::RESHRINK * tmp.v[k];
}

//template<int BITS>
//void appxFixExt(Vec3d &out, FixExt4_1<BITS> &in)
//{
//    Vec4d temp;
//    temp[0] = approximate(in.e0);
//    temp[1] = approximate(in.e1);
//    temp[2] = approximate(in.e2);
//    temp[3] = approximate(in.e3);
//    std::cout << "appx " << temp << std::endl;
//    for(uint k=0; k<3; k++) {
//        out[k] = temp[k];
//        out[k] /= temp[3];
//        out[k] *= Quantization::RESHRINK;
//    }
//}

const static double EPS                 = DBL_EPSILON;
const static double EPS2                = EPS * EPS;

inline bool filterCheck(double val, double absval, double coeff) {
    return fabs(val) > absval*coeff;
}


bool isEmpty(const TriEdgeIn &input)
{
    callcount++;
    
    Ext4_2 temp_e2;
    
    // construct the triangle
    Ext4_1 tp0, tp1, tp2;
    toExt(tp0, input.tri.p[0]);
    toExt(tp1, input.tri.p[1]);
    toExt(tp2, input.tri.p[2]);
    Ext4_3 t_ext3;
    join(temp_e2, tp0, tp1);
    join(t_ext3, temp_e2, tp2);
    
    // construct the edge
    Ext4_1 ep0, ep1;
    toExt(ep0, input.edge.p[0]);
    toExt(ep1, input.edge.p[1]);
    Ext4_2 e_ext2;
    join(e_ext2, ep0, ep1);
    
    // compute the point of intersection
    Ext4_1 p_isct;
    meet(p_isct, e_ext2, t_ext3);
    // need to adjust for negative w-coordinate
    if(p_isct.e3 < 0.0)
        neg(p_isct, p_isct);
    
    // a_t0 is the variation of t_ext3 with tp0 replaced with p_isct
    // a_t1 is the variation of t_ext3 with tp1 replaced with p_isct
    // and so on...
    Ext4_3 a_t0, a_t1, a_t2;
    join(temp_e2, p_isct, tp1);     join(a_t0, temp_e2, tp2);
    join(temp_e2, tp0,    p_isct);  join(a_t1, temp_e2, tp2);
    join(temp_e2, tp0,    tp1);     join(a_t2, temp_e2, p_isct);
    Ext4_2 a_e0, a_e1;
    join(a_e0, p_isct, ep1);
    join(a_e1, ep0,    p_isct);
    
    return (inner(t_ext3, a_t0) < 0.0) ||
           (inner(t_ext3, a_t1) < 0.0) ||
           (inner(t_ext3, a_t2) < 0.0) ||
           (inner(e_ext2, a_e0) < 0.0) ||
           (inner(e_ext2, a_e1) < 0.0);
}

Vec3d coords(const TriEdgeIn &input)
{
    Ext4_2 temp_e2;
    
    // construct the triangle
    Ext4_1 tp0, tp1, tp2;
    toExt(tp0, input.tri.p[0]);
    toExt(tp1, input.tri.p[1]);
    toExt(tp2, input.tri.p[2]);
    Ext4_3 t_ext3;
    join(temp_e2, tp0, tp1);
    join(t_ext3, temp_e2, tp2);
    
    // construct the edge
    Ext4_1 ep0, ep1;
    toExt(ep0, input.edge.p[0]);
    toExt(ep1, input.edge.p[1]);
    Ext4_2 e_ext2;
    join(e_ext2, ep0, ep1);
    
    // compute the point of intersection
    Ext4_1 p_isct;
    meet(p_isct, e_ext2, t_ext3);
    // no need to adjust for negative w-coordinate.
    // will drop out in divide.
    
    Vec3d result;
    toVec3d(result, p_isct);
    return result;
}

const static double COEFF_IT12_PISCT    = 10.0*EPS + 64.0*EPS2;
const static double COEFF_IT12_S1       = 20.0*EPS + 256.0*EPS2;
const static double COEFF_IT12_S2       = 24.0*EPS + 512.0*EPS2;

int emptyFilter(const TriEdgeIn &input)
{
    Ext4_2 temp2;                           AbsExt4_2 ktemp2;
    Ext4_1 ep[2];                           AbsExt4_1 kep[2];
    Ext4_1 tp[3];                           AbsExt4_1 ktp[3];
    Ext4_2 e_ext2;                          AbsExt4_2 ke_ext2;
    Ext4_3 t_ext3;                          AbsExt4_3 kt_ext3;
    
    // load the points
    for(int i=0; i<2; i++) {
        toExt(ep[i], input.edge.p[i]);      abs(kep[i], ep[i]);
    }
    for(int i=0; i<3; i++) {
        toExt(tp[i], input.tri.p[i]);       abs(ktp[i], tp[i]);
    }
    // form the edge and triangle
    join(e_ext2, ep[0], ep[1]);             join(ke_ext2, kep[0], kep[1]);
    join(temp2,  tp[0], tp[1]);             join(ktemp2,  ktp[0], ktp[1]);
    join(t_ext3, temp2, tp[2]);             join(kt_ext3, ktemp2, ktp[2]);
    
    // compute the point of intersection
    Ext4_1 pisct;                           AbsExt4_1 kpisct;
    meet(pisct, e_ext2, t_ext3);            meet(kpisct, ke_ext2, kt_ext3);
    // We perform one of the filter exit tests here...
    if(!filterCheck(pisct.e3, kpisct.e3, COEFF_IT12_PISCT))
        return 0; // i.e. uncertain
    // need to adjust for negative w-coordinate
    if(pisct.e3 < 0.0)
        neg(pisct, pisct);
    
    bool uncertain = false;
    // process edge
    for(int i=0; i<2; i++) {
        Ext4_2 a;                           AbsExt4_2 ka;
        join(a, (i==0)? pisct : ep[0],
                (i==1)? pisct : ep[1]);
                                            join(ka, (i==0)? kpisct : kep[0],
                                                     (i==1)? kpisct : kep[1]);
        double dot = inner(e_ext2, a);      double kdot = inner(ke_ext2, ka);
        // now figure out what to do...
            bool outside = dot < 0.0;
            bool reliable = filterCheck(dot, kdot, COEFF_IT12_S1);
            if(reliable && outside)
                return 1; // i.e. true
            if(!reliable)   uncertain = true;
    }
    // process triangle
    for(int i=0; i<3; i++) {
        Ext4_3 a;                           AbsExt4_3 ka;
        join(temp2, (i==0)? pisct : tp[0],
                    (i==1)? pisct : tp[1]);
        join(a,     temp2,
                    (i==2)? pisct : tp[2]);
                                    join(ktemp2, (i==0)? kpisct : ktp[0],
                                                 (i==1)? kpisct : ktp[1]);
                                    join(ka,     ktemp2,
                                                 (i==2)? kpisct : ktp[2]);
        double dot = inner(t_ext3, a);      double kdot = inner(kt_ext3, ka);
        // now figure out what to do...
            bool outside = dot < 0.0;
            bool reliable = filterCheck(dot, kdot, COEFF_IT12_S2);
            if(reliable && outside)
                return 1; // i.e. true
            if(!reliable)   uncertain = true;
    }
    if(uncertain)
        return 0;
    else
        return -1; // i.e. false (the intersection is not empty)
}

bool exactFallback(const TriEdgeIn &input)
{
    // How many bits do we need for various intermediary values?
    // Here we label the amount with the relevant type (i.e. EXT2)
    // and the relevant role
    const static int LINE_BITS       = 2*IN_BITS + 1;
    const static int TRI_BITS        = LINE_BITS + IN_BITS + 2;
    const static int ISCT_BITS       = TRI_BITS + LINE_BITS + 2;
    const static int LINE_A_BITS     = ISCT_BITS + IN_BITS + 1;
    const static int TRI_A_BITS      = LINE_A_BITS + IN_BITS + 2;
    const static int INNER_LINE_BITS = LINE_A_BITS + LINE_BITS + 3;
    const static int INNER_TRI_BITS  = TRI_A_BITS + TRI_BITS + 2;
    
    // pull in points
    FixExt4_1<IN_BITS>                  ep[2];
    FixExt4_1<IN_BITS>                  tp[3];
    for(uint i=0; i<2; i++)
        toFixExt(ep[i], input.edge.p[i]);
    for(uint i=0; i<3; i++)
        toFixExt(tp[i], input.tri.p[i]);
    
    // construct geometry
    FixExt4_2<LINE_BITS>                e;
    join(e, ep[0], ep[1]);
    FixExt4_2<LINE_BITS>                temp_up;
    FixExt4_3<TRI_BITS>                 t;
    join(temp_up, tp[0], tp[1]);
    join(t,     temp_up, tp[2]);
    
    // compute the point of intersection
    FixExt4_1<ISCT_BITS>                pisct;
    meet(pisct, e, t);
    // need to adjust for negative w-coordinate
    int e3sign = sign(pisct.e3);
    if(e3sign < 0) {
        neg(pisct, pisct);
    } else if(e3sign == 0) {
        degeneracy_count++;
        return true;
    }
    
    // process edge
    FixExt4_2<LINE_A_BITS>              ae0, ae1;
    BitInt<INNER_LINE_BITS>::Rep        test_e0, test_e1;
    int                                 sign_e0, sign_e1;
    join(ae0, pisct, ep[1]);
    join(ae1, ep[0], pisct);
    inner(test_e0, e, ae0);
    inner(test_e1, e, ae1);
    sign_e0 = sign(test_e0);
    sign_e1 = sign(test_e1);
    
    // process triangle
    FixExt4_3<TRI_A_BITS>               at0, at1, at2;
    FixExt4_2<LINE_A_BITS>              temp0, temp1;
    FixExt4_2<LINE_BITS>                temp2;
    BitInt<INNER_TRI_BITS>::Rep         test_t0, test_t1, test_t2;
    int                                 sign_t0, sign_t1, sign_t2;
    join(temp0, pisct, tp[1]);      join(at0, temp0, tp[2]);
    join(temp1, tp[0], pisct);      join(at1, temp1, tp[2]);
    join(temp2, tp[0], tp[1]);      join(at2, temp2, pisct);
    inner(test_t0, t, at0);
    inner(test_t1, t, at1);
    inner(test_t2, t, at2);
    sign_t0 = sign(test_t0);
    sign_t1 = sign(test_t1);
    sign_t2 = sign(test_t2);
    
    if(sign_e0 < 0 || sign_e1 < 0 ||
       sign_t0 < 0 || sign_t1 < 0 || sign_t2 < 0)
        return true;
    
    if(sign_e0 == 0 || sign_e1 == 0 ||
       sign_t0 == 0 || sign_t1 == 0 || sign_t2 == 0)
    {
        degeneracy_count++;
    }
    return false;
}

bool emptyExact(const TriEdgeIn &input)
{
    callcount++;
    int filter = emptyFilter(input);
    if(filter == 0) {
        exact_count++;
        return exactFallback(input);
    }
    else
        return filter > 0;
}

Vec3d coordsExact(const TriEdgeIn &input)
{
    // How many bits do we need for various intermediary values?
    // Here we label the amount with the relevant type (i.e. EXT2)
    // and the relevant role
    //const static int LINE_BITS       = 2*IN_BITS + 1;
    //const static int TRI_BITS        = LINE_BITS + IN_BITS + 2;
    //const static int ISCT_BITS       = TRI_BITS + LINE_BITS + 2;
    //const static int LINE_A_BITS     = ISCT_BITS + IN_BITS + 1;
    //const static int TRI_A_BITS      = LINE_A_BITS + IN_BITS + 2;
    //const static int INNER_LINE_BITS = LINE_A_BITS + LINE_BITS + 3;
    //const static int INNER_TRI_BITS  = TRI_A_BITS + TRI_BITS + 2;
    
    // pull in points
    GmpExt4_1                           ep[2];
    GmpExt4_1                           tp[3];
    for(uint i=0; i<2; i++)
        toGmpExt(ep[i], input.edge.p[i]);
    for(uint i=0; i<3; i++)
        toGmpExt(tp[i], input.tri.p[i]);
    
    // construct geometry
    GmpExt4_2                           e;
    join(e, ep[0], ep[1]);
    GmpExt4_2                           temp_up;
    GmpExt4_3                           t;
    join(temp_up, tp[0], tp[1]);
    join(t,     temp_up, tp[2]);
    
    // compute the point of intersection
    GmpExt4_1                           pisct;
    meet(pisct, e, t);
    
    // convert to double
    Vec3d result;
    toVec3d(result, pisct);
    //std::cout << result << std::endl;
    return result;
}








bool isEmpty(const TriTriTriIn &input)
{
    callcount++;
    
    Ext4_2 temp_e2;
    
    // construct the triangles
    Ext4_1 tps[3][3];
    Ext4_3 t_ext3s[3];
    for(uint ti=0; ti<3; ti++) {
        for(uint pi=0; pi<3; pi++) {
            toExt(tps[ti][pi], input.tri[ti].p[pi]);
        }
        Ext4_2 temp_e2; // help out liveness analysis
        join(temp_e2, tps[ti][0], tps[ti][1]);
        join(t_ext3s[ti], temp_e2, tps[ti][2]);
    }
    
    // compute the point of intersection
    Ext4_1 p_isct;
    meet(temp_e2, t_ext3s[0], t_ext3s[1]);
    meet(p_isct, temp_e2, t_ext3s[2]);
    // need to adjust for negative w-coordinate
    if(p_isct.e3 < 0.0)
        neg(p_isct, p_isct);
    
    // Test whether p_isct is inside each triangle
    // For each triangle, this is done by creating three modified
    // versions, replacing each point with p_isct.
    // If none of these modifications flips the orientation of
    // the resulting plane, then the point must lie inside the triangle
    for(uint ti=0; ti<3; ti++) {
        for(uint pi=0; pi<3; pi++) { // three copies...
            Ext4_3 a;
            Ext4_2 temp_e2; // help out liveness analysis
            join(temp_e2,       ((pi==0)? p_isct : tps[ti][0]),
                                ((pi==1)? p_isct : tps[ti][1]));
            join(a, temp_e2,    ((pi==2)? p_isct : tps[ti][2]));
            double test = inner(t_ext3s[ti], a);
            if(test < 0.0) // AHA, p_isct IS outside this triangle
                return true;
        }
    }
    // well, p_isct must be inside all of the triangles.
    return false;
}

Vec3d coords(const TriTriTriIn &input)
{
    Ext4_2 temp_e2;
    
    // construct the triangles
    Ext4_1 tps[3][3];
    Ext4_3 t_ext3s[3];
    for(uint ti=0; ti<3; ti++) {
        for(uint pi=0; pi<3; pi++) {
            toExt(tps[ti][pi], input.tri[ti].p[pi]);
        }
        Ext4_2 temp_e2; // help out liveness analysis
        join(temp_e2, tps[ti][0], tps[ti][1]);
        join(t_ext3s[ti], temp_e2, tps[ti][2]);
    }
    
    // compute the point of intersection
    Ext4_1 p_isct;
    meet(temp_e2, t_ext3s[0], t_ext3s[1]);
    meet(p_isct, temp_e2, t_ext3s[2]);
    // no need to adjust for negative w-coordinate.
    // will come out in the divide.
    
    Vec3d result;
    toVec3d(result, p_isct);
    return result;
}

const static double COEFF_IT222_PISCT   = 20.0*EPS + 256.0*EPS2;
const static double COEFF_IT222_S2      = 34.0*EPS + 1024.0*EPS2;

int emptyFilter(const TriTriTriIn &input)
{
    Ext4_2 temp2;                           AbsExt4_2 ktemp2;
    Ext4_1 p[3][3];                         AbsExt4_1 kp[3][3];
    Ext4_3 t[3];                            AbsExt4_3 kt[3];
    
    // load the points and form triangles
    for(uint i=0; i<3; i++) {
      for(uint j=0; j<3; j++) {
        toExt(p[i][j], input.tri[i].p[j]);  abs(kp[i][j], p[i][j]);
      }
      join(temp2, p[i][0], p[i][1]);        join(ktemp2, kp[i][0], kp[i][1]);
      join(t[i],  temp2,   p[i][2]);        join(kt[i],  ktemp2,   kp[i][2]);
    }
    
    // compute the point of intersection
    Ext4_1 pisct;                           AbsExt4_1 kpisct;
    meet(temp2, t[0],  t[1]);               meet(ktemp2, kt[0],  kt[1]);
    meet(pisct, temp2, t[2]);               meet(kpisct, ktemp2, kt[2]);
    // We perform one of the filter exit tests here...
    if(!filterCheck(pisct.e3, kpisct.e3, COEFF_IT222_PISCT))
        return 0; // i.e. uncertain
    // need to adjust for negative w-coordinate
    if(pisct.e3 < 0.0)
        neg(pisct, pisct);
    
    bool uncertain = false;
    for(int i=0; i<3; i++) {
      for(int j=0; j<3; j++) {
        Ext4_2 b;                           AbsExt4_2 kb;
        Ext4_3 a;                           AbsExt4_3 ka;
        join(b, (j==0)? pisct : p[i][0],
                (j==1)? pisct : p[i][1]);
        join(a, b,
                (j==2)? pisct : p[i][2]);
                                        join(kb, (j==0)? kpisct : kp[i][0],
                                                 (j==1)? kpisct : kp[i][1]);
                                        join(ka, kb,
                                                 (j==2)? kpisct : kp[i][2]);
        double dot = inner(t[i], a);        double kdot = inner(kt[i], ka);
        // then consider the filtering...
            bool outside = dot < 0.0;
            bool reliable = filterCheck(dot, kdot, COEFF_IT222_S2);
            if(reliable && outside)
                return 1; // i.e. true
            if(!reliable)   uncertain = true;
      }
    }
    if(uncertain)
        return 0;
    else
        return -1; // i.e. false (the intersection is not empty)
}

bool exactFallback(const TriTriTriIn &input)
{
    // How many bits do we need for various intermediary values?
    // Here we label the amount with the relevant type (i.e. EXT2)
    // and the relevant role
    const static int EXT2_UP_BITS = 2*IN_BITS + 1;
    const static int EXT3_UP_BITS = EXT2_UP_BITS + IN_BITS + 2;
    const static int EXT2_DN_BITS = 2*EXT3_UP_BITS + 1;
    const static int ISCT_BITS    = EXT2_DN_BITS + EXT3_UP_BITS + 2;
    const static int EXT2_TA_BITS = ISCT_BITS + IN_BITS + 1;
    const static int EXT3_TA_BITS = EXT2_TA_BITS + IN_BITS + 2;
    const static int INNER_BITS   = EXT3_TA_BITS + EXT3_UP_BITS + 2;
    
    FixExt4_1<IN_BITS>                  p[3][3];
    FixExt4_3<EXT3_UP_BITS>             t[3];
    for(uint i=0; i<3; i++) {
        for(uint j=0; j<3; j++) {
            toFixExt(p[i][j], input.tri[i].p[j]);
        }
        FixExt4_2<EXT2_UP_BITS>         temp;
        join(temp, p[i][0], p[i][1]);
        join(t[i], temp,    p[i][2]);
    }
    
    // compute the point of intersection
    FixExt4_1<ISCT_BITS>                pisct;
    {
        FixExt4_2<EXT2_DN_BITS>         temp;
        meet(temp,  t[0], t[1]);
        meet(pisct, temp, t[2]);
    }
    // need to adjust for negative w-coordinate
    int e3sign = sign(pisct.e3);
    if(e3sign < 0) {
        neg(pisct, pisct);
    } else if(e3sign == 0) {
        degeneracy_count++;
        return true;
    }
    
    bool uncertain = false;
    for(uint i=0; i<3; i++) {
        FixExt4_3<EXT3_TA_BITS>         a[3];
        FixExt4_2<EXT2_TA_BITS>         temp;
        FixExt4_2<EXT2_UP_BITS>         tmp2;
        
        join(temp,   pisct, p[i][1]);   join(a[0], temp, p[i][2]);
        join(temp, p[i][0],   pisct);   join(a[1], temp, p[i][2]);
        join(tmp2, p[i][0], p[i][1]);   join(a[2], tmp2, pisct);
        for(uint j=0; j<3; j++) {
            BitInt<INNER_BITS>::Rep     test;
            int                         testsign;
            inner(test, a[j], t[i]);
            testsign = sign(test);
            if(testsign < 0)
                return true;
            if(testsign == 0)
                uncertain = true;
        }
    }
    if(uncertain) {
        degeneracy_count++;
    }
    return false;
}

bool emptyExact(const TriTriTriIn &input)
{
    callcount++;
    int filter = emptyFilter(input);
    if(filter == 0) {
        exact_count++;
        return exactFallback(input);
    }
    else
        return filter > 0;
}

Vec3d coordsExact(const TriTriTriIn &input)
{
    // How many bits do we need for various intermediary values?
    // Here we label the amount with the relevant type (i.e. EXT2)
    // and the relevant role
    //const static int EXT2_UP_BITS = 2*IN_BITS + 1;
    //const static int EXT3_UP_BITS = EXT2_UP_BITS + IN_BITS + 2;
    //const static int EXT2_DN_BITS = 2*EXT3_UP_BITS + 1;
    //const static int ISCT_BITS    = EXT2_DN_BITS + EXT3_UP_BITS + 2;
    //const static int EXT2_TA_BITS = ISCT_BITS + IN_BITS + 1;
    //const static int EXT3_TA_BITS = EXT2_TA_BITS + IN_BITS + 2;
    //const static int INNER_BITS   = EXT3_TA_BITS + EXT3_UP_BITS + 2;
    
    GmpExt4_1                           p[3][3];
    GmpExt4_3                           t[3];
    for(uint i=0; i<3; i++) {
        for(uint j=0; j<3; j++) {
            toGmpExt(p[i][j], input.tri[i].p[j]);
        }
        GmpExt4_2                       temp;
        join(temp, p[i][0], p[i][1]);
        join(t[i], temp,    p[i][2]);
    }
    
    // compute the point of intersection
    GmpExt4_1                           pisct;
    {
        GmpExt4_2                       temp;
        meet(temp,  t[0], t[1]);
        meet(pisct, temp, t[2]);
    }
    
    // convert to double
    Vec3d result;
    toVec3d(result, pisct);
    return result;
}



} // end namespace Empty3d



