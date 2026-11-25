// Compare certified fast-path variants with the original exact traversal.
#include "filter_ablation.h"
#include <igl/predicates/exact_point_mesh_closest_triangle.h>
#include <igl/default_num_threads.h>
#include <igl/eytzinger_aabb.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>
#include <algorithm>
#include <chrono>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>

using Point=Eigen::RowVector3d;
struct Method
{
  std::string name;
  std::function<Eigen::Index(const Point &)> evaluate;
  std::function<void()> audit;
  std::vector<Eigen::Index> output;
  std::vector<double> samples;
};
struct Experiment
{
  Eigen::MatrixXd vertices,lower,upper,queries;
  Eigen::MatrixXi faces;
  Eigen::VectorXi leaves;
  std::vector<Method> methods;
  template<bool FastUpper,bool CacheDistance,bool Fused,bool Enclosures>
  void add(const char *name)
  {
    methods.push_back({name,[&](const Point &q)
    {
      filter_profile::Counts unused;
      return filter_profile::search<FastUpper,CacheDistance,Fused,Enclosures>(
        q,vertices,faces,lower,upper,leaves,unused);
    },[&,name]()
    {
      filter_profile::Counts counts;
      for(Eigen::Index i=0;i<queries.rows();++i)
      {
        const Point q=queries.row(i);
        const auto found=filter_profile::search<FastUpper,CacheDistance,Fused,Enclosures,true>(
          q,vertices,faces,lower,upper,leaves,counts);
        if(found!=methods.front().output[i]) throw std::runtime_error("Audit disagrees with baseline");
      }
      std::cout<<"# audit "<<name<<" nodes="<<counts.nodes<<" leaves="<<counts.leaves
        <<" updates="<<counts.updates<<" feature_fallbacks="<<counts.feature_fallbacks
        <<" interval_bounds="<<counts.interval_bounds<<" filter_compares="<<counts.filter_compares
        <<" enclosure_compares="<<counts.enclosure_compares<<" interval_compares="<<counts.interval_compares
        <<" expansion_compares="<<counts.expansion_compares<<'\n';
    },std::vector<Eigen::Index>(queries.rows()),{}});
  }
  void run(Method &method)
  {
    std::vector<std::exception_ptr> errors;
    igl::parallel_for(queries.rows(),[&](std::size_t n) { errors.resize(n); },
      [&](Eigen::Index i,std::size_t thread)
      {
        try { const Point q=queries.row(i); method.output[i]=method.evaluate(q); }
        catch(...) { if(!errors[thread]) errors[thread]=std::current_exception(); }
      },[&](std::size_t thread) { if(errors[thread]) std::rethrow_exception(errors[thread]); },10000);
  }
};
int positive(const char *text)
{
  std::size_t consumed=0;
  const long long value=std::stoll(text,&consumed);
  if(text[consumed]!='\0' || value<1 || value>std::numeric_limits<int>::max())
    throw std::invalid_argument("Counts must be positive integers fitting in int");
  return static_cast<int>(value);
}
int main(int argc,char **argv)
{
  if(argc<2 || argc>5) { std::cerr<<"Usage: profile_point_filters mesh.obj [points=10000] [repetitions=9] [threads=1]\n"; return 1; }
  try
  {
    const int count=argc>2 ? positive(argv[2]) : 10000;
    const int repetitions=argc>3 ? positive(argv[3]) : 9;
    const int threads=argc>4 ? positive(argv[4]) : 1;
    igl::default_num_threads(threads); Eigen::setNbThreads(1);
    igl::predicates::exactinit();
    Experiment e;
    if(!igl::read_triangle_mesh(argv[1],e.vertices,e.faces) || e.vertices.rows()==0 || e.faces.rows()==0 ||
       e.vertices.cols()!=3 || e.faces.cols()!=3 || !e.vertices.allFinite() ||
       e.faces.minCoeff()<0 || e.faces.maxCoeff()>=e.vertices.rows())
      throw std::invalid_argument("Expected a nonempty finite triangle mesh with valid indices");
    Eigen::MatrixXd lo(e.faces.rows(),3),hi(e.faces.rows(),3);
    for(Eigen::Index i=0;i<e.faces.rows();++i)
    {
      lo.row(i)=e.vertices.row(e.faces(i,0)); hi.row(i)=lo.row(i);
      for(int j=1;j<3;++j)
      { lo.row(i)=lo.row(i).cwiseMin(e.vertices.row(e.faces(i,j))); hi.row(i)=hi.row(i).cwiseMax(e.vertices.row(e.faces(i,j))); }
    }
    igl::eytzinger_aabb(lo,hi,e.lower,e.upper,e.leaves);
    const Point minimum=e.vertices.colwise().minCoeff(),extent=e.vertices.colwise().maxCoeff()-minimum;
    std::mt19937_64 generator(896);
    e.queries.resize(count,3);
    for(int i=0;i<count;++i) for(int axis=0;axis<3;++axis)
    {
      const double unit=(generator()>>11)*(1.0/9007199254740992.0);
      e.queries(i,axis)=minimum(axis)+extent(axis)*(-.1+1.2*unit);
    }
    e.add<false,false,false,false>("baseline");
    e.add<true,false,false,false>("filtered_upper");
    e.add<false,true,false,false>("reuse_distance");
    e.add<true,true,false,false>("reuse_filtered_upper");
    e.add<true,true,true,false>("fused_reuse_upper");
    e.add<false,true,false,true>("enclosures");
    e.add<false,true,true,true>("fused_enclosures");
    e.methods.push_back({"production",[&](const Point &q)
    { return igl::predicates::exact_point_mesh_closest_triangle(q,e.vertices,e.faces,e.lower,e.upper,e.leaves); },
      [](){},std::vector<Eigen::Index>(count),{}});
    for(auto &method:e.methods) e.run(method);
    for(int repetition=0;repetition<repetitions;++repetition)
      for(std::size_t offset=0;offset<e.methods.size();++offset)
      {
        auto &method=e.methods[(offset+repetition)%e.methods.size()];
        const auto start=std::chrono::steady_clock::now(); e.run(method);
        method.samples.push_back(std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count());
      }
    std::cout<<std::fixed<<std::setprecision(6)<<"# points="<<count<<" seed=896 padding=0.1 threads="<<threads
      <<" repetitions="<<repetitions<<" min_parallel=10000\nmethod,median_ms,min_ms,max_ms\n";
    for(auto &method:e.methods)
    {
      if(method.output!=e.methods.front().output) throw std::runtime_error("Timing variant disagrees with baseline");
      std::sort(method.samples.begin(),method.samples.end());
      std::cout<<method.name<<','<<method.samples[method.samples.size()/2]<<','<<method.samples.front()<<','<<method.samples.back()<<'\n';
    }
    for(auto &method:e.methods) method.audit();
    std::cout<<"# all_triangle_ids_identical=1\n";
  }
  catch(const std::exception &error) { std::cerr<<error.what()<<'\n'; return 1; }
}
