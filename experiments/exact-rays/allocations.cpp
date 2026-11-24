#include "api.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <new>

namespace
{
bool track_allocations=false;
std::size_t allocation_count=0;
}
void *operator new(std::size_t size)
{
  if(track_allocations) { ++allocation_count; }
  if(void *memory=std::malloc(size==0 ? 1 : size)) { return memory; }
  throw std::bad_alloc();
}
void *operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void *memory) noexcept { std::free(memory); }
void operator delete[](void *memory) noexcept { std::free(memory); }
void operator delete(void *memory,std::size_t) noexcept { std::free(memory); }
void operator delete[](void *memory,std::size_t) noexcept { std::free(memory); }

int main()
{
  const double epsilon=std::ldexp(1.0,-52);
  const Ray ray{Point(0,0,0),Point(0,0,1)};
  const Triangle first{Point(0,0,1),Point(1,1+epsilon,1),Point(1-epsilon,1,1)};
  const Triangle second{Point(0,0,2),Point(1,1+epsilon,2),Point(1-epsilon,1,2)};
  shewchuk::classify(ray,first); // Initialize the backend before measuring.
  track_allocations=true;
  const int original_class=baseline::classify(ray,first);
  const int original_compare=baseline::compare(ray,first,second);
  const std::size_t original_allocations=allocation_count;
  allocation_count=0;
  const int changed_class=shewchuk::classify(ray,first);
  const int changed_compare=shewchuk::compare(ray,first,second);
  const std::size_t changed_allocations=allocation_count;
  track_allocations=false;
  std::cout<<"baseline_allocations,shewchuk_allocations\n"
    <<original_allocations<<","<<changed_allocations<<"\n";
  return original_class!=1 || changed_class!=1 || original_compare!=-1
    || changed_compare!=-1 || original_allocations==0 || changed_allocations!=0;
}
