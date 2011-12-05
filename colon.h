#ifndef IGL_COLON_H
#define IGL_COLON_H
#include <Eigen/Dense>
namespace igl
{
  // Note:
  // This should be potentially replaced with eigen's LinSpaced() function

  // Colon operator like matlab's colon operator. Enumerats values between low
  // and hi with step step.
  // Templates:
  //   L  should be a eigen matrix primitive type like int or double
  //   S  should be a eigen matrix primitive type like int or double
  //   H  should be a eigen matrix primitive type like int or double
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   low  starting value if step is valid then this is *always* the first
  //     element of I
  //   step  step difference between sequential elements returned in I,
  //     remember this will be cast to template T at compile time. If low<hi
  //     then step must be positive. If low>hi then step must be negative.
  //     Otherwise I will be set to empty.
  //   hi  ending value, if (hi-low)%step is zero then this will be the last
  //     element in I. If step is positive there will be no elements greater
  //     than hi, vice versa if hi<low
  // Output:
  //   I  list of values from low to hi with step size step
  template <typename L,typename S,typename H,typename T>
  inline void colon(
    const L low, 
    const S step, 
    const H hi, 
    Eigen::Matrix<T,Eigen::Dynamic,1> & I);
  // Same as above but step == (T)1
  template <typename L,typename H,typename T>
  inline void colon(
    const L low, 
    const H hi, 
    Eigen::Matrix<T,Eigen::Dynamic,1> & I);
  // Return output rather than set in reference
  template <typename T,typename L,typename H>
  inline Eigen::Matrix<T,Eigen::Dynamic,1> colon(
    const L low, 
    const H hi);
}

// Implementation
#include <cstdio>

template <typename L,typename S,typename H,typename T>
inline void igl::colon(
  const L low, 
  const S step, 
  const H hi, 
  Eigen::Matrix<T,Eigen::Dynamic,1> & I)
{
  if(low < hi)
  {
    if(step < 0)
    {
      I.resize(0);
      fprintf(stderr,"Error: colon() low(%g)<hi(%g) but step(%g)<0\n",
        (double)low,
        (double)hi,
        (double)step);
      return;
    }
  }
  if(low > hi)
  {
    if(step > 0)
    {
      I.resize(0);
      fprintf(stderr,"Error: colon() low(%g)>hi(%g) but step(%g)<0\n",
        (double)low,
        (double)hi,
        (double)step);
      return;
    }
  }
  // resize output
  int n = floor(double((hi-low)/step))+1;
  I.resize(n);
  int i = 0;
  T v = (T)low;
  while((low<hi && (H)v<=hi) || (low>hi && (H)v>=hi))
  {
    I(i) = v;
    v = v + (T)step;
    i++;
  }
  assert(i==n);
}

template <typename L,typename H,typename T>
inline void igl::colon(
  const L low, 
  const H hi, 
  Eigen::Matrix<T,Eigen::Dynamic,1> & I)
{
  return igl::colon(low,(T)1,hi,I);
}

template <typename T,typename L,typename H>
inline Eigen::Matrix<T,Eigen::Dynamic,1> igl::colon(
  const L low, 
  const H hi)
{
  Eigen::Matrix<T,Eigen::Dynamic,1> I;
  igl::colon(low,hi,I);
  return I;
}

#endif

