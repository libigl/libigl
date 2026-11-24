#include <cstdint>
#include "api.h"
#include "workloads.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace
{
using exact_ray_experiment::cases;
volatile std::int64_t checksum = 0;

void require(bool condition, const char *message)
{
  if(!condition) { throw std::runtime_error(message); }
}

void mesh(int faces_count, bool layered, Eigen::MatrixXd &vertices, Eigen::MatrixXi &faces)
{
  vertices.resize(3*faces_count,3); faces.resize(faces_count,3);
  std::mt19937 generator(712);
  std::uniform_real_distribution<double> uniform(-1,1);
  for(int face=0; face<faces_count; ++face)
  {
    if(layered)
    {
      double height=double((face*7919)%faces_count+1)/faces_count;
      vertices.row(3*face)<<-1,-1,height;
      vertices.row(3*face+1)<<1,-1,height;
      vertices.row(3*face+2)<<0,1,height;
    }else
    {
      for(int corner=0; corner<3; ++corner)
        for(int axis=0; axis<3; ++axis) vertices(3*face+corner,axis)=uniform(generator);
    }
    faces.row(face)<<3*face,3*face+1,3*face+2;
  }
}

void accuracy()
{
  std::size_t classifications=0, comparisons=0, mesh_queries=0;
  for(const std::string kind : {"random","coplanar","thin","ordered","dense_ordered"})
  {
    for(const auto &sample : cases(kind,25000))
    {
      const int original=baseline::classify(sample.ray,sample.first);
      require(original==shewchuk::classify(sample.ray,sample.first),"Classification differs from baseline");
      ++classifications;
      if(original==1 && baseline::classify(sample.ray,sample.second)==1)
      {
        require(baseline::compare(sample.ray,sample.first,sample.second)
          ==shewchuk::compare(sample.ray,sample.first,sample.second),"Ordering differs from baseline");
        ++comparisons;
      }
    }
  }
  for(bool layered : {false,true})
  {
    for(int size : {1,16,128})
    {
      Eigen::MatrixXd vertices; Eigen::MatrixXi faces;
      mesh(size,layered,vertices,faces);
      for(const auto &sample : cases("random",100))
      {
        Ray ray=layered ? Ray{Point(0,0,-1),Point(0,0,1)} : sample.ray;
        std::vector<int> original, changed;
        baseline::all(ray,vertices,faces,original);
        shewchuk::all(ray,vertices,faces,changed);
        require(original==changed,"All-hit ordering differs from baseline");
        int first_original=-1, first_changed=-1;
        require(baseline::first(ray,vertices,faces,first_original)
          ==shewchuk::first(ray,vertices,faces,first_changed),"First-hit status differs from baseline");
        require(first_original==first_changed,"First-hit index differs from baseline");
        mesh_queries+=2;
      }
    }
  }
  std::cout<<"classifications,hit_comparisons,mesh_queries,mismatches\n"
    <<classifications<<","<<comparisons<<","<<mesh_queries<<",0\n";
}

template<typename Operation>
double duration(const Operation &operation, const int repetitions)
{
  const auto start=std::chrono::steady_clock::now();
  std::int64_t sum=0;
  for(int i=0; i<repetitions; ++i) { sum+=operation(); }
  checksum=sum;
  return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/repetitions;
}

template<typename First, typename Second>
void benchmark(const std::string &name, const int operations, First original, Second changed)
{
  // Warm both paths, calibrate to >=50ms per sample, then alternate ordering.
  duration(original,1); duration(changed,1);
  int repeats_original=1, repeats_changed=1;
  while(duration(original,repeats_original)*repeats_original<0.05) { repeats_original*=2; }
  while(duration(changed,repeats_changed)*repeats_changed<0.05) { repeats_changed*=2; }
  std::vector<double> original_times,changed_times;
  for(int repetition=0; repetition<7; ++repetition)
  {
    if(repetition%2)
    {
      changed_times.push_back(duration(changed,repeats_changed));
      original_times.push_back(duration(original,repeats_original));
    }else
    {
      original_times.push_back(duration(original,repeats_original));
      changed_times.push_back(duration(changed,repeats_changed));
    }
  }
  std::sort(original_times.begin(),original_times.end());
  std::sort(changed_times.begin(),changed_times.end());
  std::cout<<name<<","<<operations<<","<<original_times[3]*1e9/operations
    <<","<<changed_times[3]*1e9/operations<<","<<original_times[3]/changed_times[3]<<"\n"<<std::flush;
}

void performance()
{
  std::cout<<"case,operations,baseline_ns,shewchuk_ns,speedup\n";
  for(const std::string kind : {"random","ordered","dense_ordered","coplanar","thin"})
  {
    const auto samples=cases(kind,1024);
    benchmark("classify_"+kind,samples.size(),[&](){
      int sum=0; for(const auto &sample:samples) sum+=baseline::classify(sample.ray,sample.first); return sum;
    },[&](){
      int sum=0; for(const auto &sample:samples) sum+=shewchuk::classify(sample.ray,sample.first); return sum;
    });
  }
  for(const std::string kind : {"ordered","dense_ordered","thin"})
  {
    const auto samples=cases(kind,512);
    benchmark("compare_"+kind,samples.size(),[&](){
      int sum=0; for(const auto &sample:samples) sum+=baseline::compare(sample.ray,sample.first,sample.second); return sum;
    },[&](){
      int sum=0; for(const auto &sample:samples) sum+=shewchuk::compare(sample.ray,sample.first,sample.second); return sum;
    });
  }
  for(bool layered : {false,true})
  {
    for(int size : {128,1024})
    {
      Eigen::MatrixXd vertices; Eigen::MatrixXi faces;
      mesh(size,layered,vertices,faces);
      const auto samples=cases("random",8);
      const std::string label=std::string(layered?"layers_":"random_")+std::to_string(size);
      benchmark("mesh_first_"+label,samples.size(),[&](){
        int sum=0; for(const auto &sample:samples) { int face=-1;
          baseline::first(layered?Ray{Point(0,0,-1),Point(0,0,1)}:sample.ray,vertices,faces,face); sum+=face; } return sum;
      },[&](){
        int sum=0; for(const auto &sample:samples) { int face=-1;
          shewchuk::first(layered?Ray{Point(0,0,-1),Point(0,0,1)}:sample.ray,vertices,faces,face); sum+=face; } return sum;
      });
      benchmark("mesh_all_"+label,samples.size(),[&](){
        int sum=0; for(const auto &sample:samples) { std::vector<int> hits;
          baseline::all(layered?Ray{Point(0,0,-1),Point(0,0,1)}:sample.ray,vertices,faces,hits); sum+=hits.size(); } return sum;
      },[&](){
        int sum=0; for(const auto &sample:samples) { std::vector<int> hits;
          shewchuk::all(layered?Ray{Point(0,0,-1),Point(0,0,1)}:sample.ray,vertices,faces,hits); sum+=hits.size(); } return sum;
      });
    }
  }
}

void known_failures()
{
  const Triangle triangle{Point(-1,-1,0),Point(1,-1,0),Point(0,1,0)};
  Eigen::MatrixXd vertices(3,3); vertices<<triangle.a,triangle.b,triangle.c;
  Eigen::MatrixXi faces(1,3); faces<<0,1,2;
  const Point direction(1,49,1);
  std::cout<<"case,baseline,shewchuk,mathematical_expectation\n";
  for(const auto &named : std::vector<std::pair<std::string,Ray>>{
    {"distant_origin",{-std::ldexp(1.0,60)*direction,direction}},
    {"tiny_direction",{Point(0,0,-1),Point(0,0,std::numeric_limits<double>::min())}}})
  {
    int face=-1;
    std::cout<<named.first<<","<<baseline::first(named.second,vertices,faces,face)
      <<","<<shewchuk::first(named.second,vertices,faces,face)<<",1\n";
  }
  for(int exponent : {-400,-220,220,400})
  {
    const double scale=std::ldexp(1.0,exponent);
    const Ray ray{Point(0,0,0),Point(0,0,1)};
    Triangle first{Point(-scale,-scale,scale),Point(scale,-scale,scale),Point(0,scale,scale)};
    Triangle second=first; second.a.z()=second.b.z()=second.c.z()=2*scale;
    std::cout<<"range_compare_"<<exponent<<","<<baseline::compare(ray,first,second)
      <<","<<shewchuk::compare(ray,first,second)<<",-1\n";
    first.a.z()=first.b.z()=first.c.z()=-scale;
    std::cout<<"range_behind_"<<exponent<<","<<baseline::classify(ray,first)
      <<","<<shewchuk::classify(ray,first)<<",0\n";
  }
  const double epsilon=std::ldexp(1.0,-52);
  const Triangle thin{Point(0,0,1),Point(1,1+epsilon,1),Point(1-epsilon,1,1)};
  const Ray ray{Point(0,0,0),Point(0,0,1)};
  std::cout<<"approximate_t,"<<baseline::approximate(ray,thin)<<","<<shewchuk::approximate(ray,thin)<<",1\n";
}
}

int main(int argc,char **argv)
{
  try
  {
    std::cout<<std::setprecision(10);
    if(argc!=2) { throw std::runtime_error("Use --accuracy, --performance or --known-failures"); }
    const std::string mode=argv[1];
    if(mode=="--accuracy") accuracy();
    else if(mode=="--performance") performance();
    else if(mode=="--known-failures") known_failures();
    else throw std::runtime_error("Unknown mode");
  }catch(const std::exception &error)
  {
    std::cerr<<error.what()<<"\n";
    return EXIT_FAILURE;
  }
}
