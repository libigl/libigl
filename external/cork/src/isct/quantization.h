// +-------------------------------------------------------------------------
// | quantization.h
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

#include<cmath>

namespace Quantization {

// NOTE: none of these values should be modified by the clients
static const int BITS = 30;
// MAGNIFY * RESHRINK == 1
extern double MAGNIFY;
extern double RESHRINK;
inline int quantize2int(double number) {
    return int(number * MAGNIFY);
}
inline double quantizedInt2double(int number) {
    return RESHRINK * double(number);
}
inline double quantize(double number) {
    return RESHRINK * double(int(number * MAGNIFY));
}

// given the specified number of bits,
// and bound on the coordinate values of points,
// fit as fine-grained a grid as possible over the space.
inline void callibrate(double maximumMagnitude)
{
    int max_exponent;
    std::frexp(maximumMagnitude, &max_exponent);
    max_exponent++; // ensure that 2^max_exponent > maximumMagnitude
    
    // set constants
    MAGNIFY = std::pow(2.0, BITS - max_exponent);
    // we are guaranteed that maximumMagnitude * MAGNIFY < 2.0^BITS
    RESHRINK = std::pow(2.0, max_exponent - BITS);
}


} // end namespace Quantization


