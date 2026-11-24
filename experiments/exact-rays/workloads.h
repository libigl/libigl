#pragma once
#include "api.h"
#include <cmath>
#include <random>
#include <string>
#include <utility>

namespace exact_ray_experiment
{
struct Case { Ray ray; Triangle first, second; };
inline std::vector<Case> cases(const std::string &kind, int count)
{
  std::mt19937 generator(2549);
  std::uniform_real_distribution<double> uniform(-1,1);
  const auto point = [&]() { return Point(uniform(generator),uniform(generator),uniform(generator)); };
  std::vector<Case> result;
  for(int i=0; i<count; ++i)
  {
    Case sample{{point()*2, point()}, {point(),point(),point()}, {point(),point(),point()}};
    if(kind=="coplanar")
    {
      sample.ray={Point(0.25,0.25,0),Point(1,0,0)};
      sample.first={Point(0,0,0),Point(1,0,0),Point(0,1,0)};
    }
    if(kind=="thin")
    {
      const double epsilon=std::ldexp(1.0,-52);
      sample.ray={Point(0,0,0),Point(0,0,1)};
      sample.first={Point(0,0,1),Point(1,1+epsilon,1),Point(1-epsilon,1,1)};
      sample.second=sample.first;
      sample.second.a.z()=sample.second.b.z()=sample.second.c.z()=2;
    }
    if(kind=="ordered" || kind=="dense_ordered")
    {
      sample.ray={Point(0,0,0),Point(0,0,1)};
      sample.first={Point(-1,-1,2+uniform(generator)),Point(1,-1,2+uniform(generator)),Point(0,1,2+uniform(generator))};
      sample.second={Point(-2,-2,2+uniform(generator)),Point(2,-2,2+uniform(generator)),Point(0,2,2+uniform(generator))};
      if(i%4==0) { sample.second=sample.first; }
      if(i%3==0) { std::swap(sample.second.b,sample.second.c); }
    }
    if(kind=="dense_ordered")
    {
      // Use an invertible transform so hit triangles have no systematic zero
      // determinant entries; exercise the cost of checking nonzero factors.
      const auto transform=[](const Point &point) {
        return Point(point.x()+2*point.y()+3*point.z(),
          3*point.x()+5*point.y()+7*point.z(),
          11*point.x()+13*point.y()+17*point.z());
      };
      sample.ray.origin=transform(sample.ray.origin);
      sample.ray.direction=transform(sample.ray.direction);
      for(Triangle *triangle : {&sample.first,&sample.second})
      {
        triangle->a=transform(triangle->a);
        triangle->b=transform(triangle->b);
        triangle->c=transform(triangle->c);
      }
    }
    const double scale=std::ldexp(1.0,(i%7-3)*20);
    sample.ray.origin*=scale; sample.ray.direction*=scale;
    for(Triangle *triangle : {&sample.first,&sample.second})
    {
      triangle->a*=scale; triangle->b*=scale; triangle->c*=scale;
    }
    result.push_back(sample);
  }
  return result;
}

}
