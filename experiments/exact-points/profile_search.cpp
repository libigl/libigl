// Isolate tree representation, traversal and exact-metric costs on tutorial data.
#include "search_ablation.h"
#include <igl/AABB.h>
#include <igl/default_num_threads.h>
#include <igl/eytzinger_aabb.h>
#include <igl/eytzinger_aabb_sdf.h>
#include <igl/fast_winding_number.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/predicates/approximate_signed_distance_exact_sign.h>
#include <algorithm>
#include <chrono>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

#ifdef PROFILE_EXPANSIONS
// Only the separate audit executable interposes arithmetic. Timed builds have
// no counter branches or linker wrappers. Audit queries run serially.
std::size_t expansion_sums=0;
extern "C" int __real_predicates_expansion_sum(int,const double*,int,const double*,double*);
extern "C" int __wrap_predicates_expansion_sum(int n,const double *a,int m,const double *b,double *out)
{
  ++expansion_sums;
  return __real_predicates_expansion_sum(n,a,m,b,out);
}
#endif

namespace
{
namespace p=igl::predicates;
namespace k=p::exact_point_triangle_kernel;
namespace m=p::exact_point_mesh_kernel;
using Point=Eigen::RowVector3d;
using Tree=igl::AABB<Eigen::MatrixXd,3>;
struct FlatTree
{
  Eigen::MatrixXd lower,upper;
  Eigen::VectorXi leaves;
  void init(const Eigen::MatrixXd &vertices,const Eigen::MatrixXi &faces)
  {
    Eigen::MatrixXd minimum(faces.rows(),3),maximum(faces.rows(),3);
    for(Eigen::Index face=0;face<faces.rows();++face)
    {
      minimum.row(face)=vertices.row(faces(face,0));
      maximum.row(face)=minimum.row(face);
      for(int corner=1;corner<3;++corner)
      {
        minimum.row(face)=minimum.row(face).cwiseMin(vertices.row(faces(face,corner)));
        maximum.row(face)=maximum.row(face).cwiseMax(vertices.row(faces(face,corner)));
      }
    }
    igl::eytzinger_aabb(minimum,maximum,lower,upper,leaves);
  }
};
struct Output
{
  double squared=0,distance=0;
  int triangle=-1;
  Point closest=Point::Zero();
  p::PointMeshClosestPoint exact;
  p::PointMeshLocation location=p::PointMeshLocation::OUTSIDE;
};
struct Method
{
  std::string name;
  std::function<void(const Point &,Output &)> evaluate;
  std::vector<Output> output;
  std::vector<double> samples;
};

template<typename Evaluate> bool parallel_queries(const Eigen::MatrixXd &queries,const Evaluate &evaluate)
{
  std::vector<std::exception_ptr> errors;
  return igl::parallel_for(queries.rows(),
    [&](std::size_t threads) { errors.resize(threads); },
    [&](Eigen::Index i,std::size_t thread)
    {
      try { const Point query=queries.row(i); evaluate(i,query); }
      catch(...) { if(!errors[thread]) errors[thread]=std::current_exception(); }
    },
    [&](std::size_t thread) { if(errors[thread]) std::rethrow_exception(errors[thread]); },10000);
}
void ordinary_flat(const Point &query,const Eigen::MatrixXd &vertices,
  const Eigen::MatrixXi &faces,const FlatTree &tree,Output &output)
{
  double best=std::numeric_limits<double>::infinity();
  igl::eytzinger_aabb_sdf<true>(query,[&](int face)
  {
    double squared; Point closest;
    igl::point_simplex_squared_distance<3>(query,vertices,faces,face,squared,closest);
    if(squared<best)
    {
      best=squared; output.triangle=face; output.closest=closest;
    }
    return squared;
  },tree.lower,tree.upper,tree.leaves,output.squared);
}
std::array<k::Point,3> triangle(const Eigen::MatrixXd &vertices,const Eigen::MatrixXi &faces,int face)
{
  return {{k::point(vertices.row(faces(face,0))),k::point(vertices.row(faces(face,1))),
    k::point(vertices.row(faces(face,2)))}};
}
int positive_argument(const char *text)
{
  std::size_t consumed=0;
  const long long value=std::stoll(text,&consumed);
  if(text[consumed]!='\0' || value<1 || value>std::numeric_limits<int>::max())
    throw std::invalid_argument("Counts must be positive integers fitting in int");
  return static_cast<int>(value);
}
// Benchmark-only local-sign prototype. For a closest vertex, ring[0] is
// already the incumbent edge; begin at ring[1] to avoid comparing it to itself.
int location_skip_self(const Point &query,const Eigen::MatrixXd &vertices,
  const Eigen::MatrixXi &faces,const Eigen::MatrixXi &neighbors,
  const Eigen::MatrixXi &neighbor_edges,const p::PointMeshClosestPoint &closest)
{
  const auto q=k::point(query),vertex=m::corner(vertices,faces,closest.triangle,closest.feature_id);
  const auto corners=triangle(vertices,faces,closest.triangle);
  if(closest.feature!=p::ClosestPointFeature::VERTEX || q==vertex)
    return m::local_location(vertices,faces,neighbors,neighbor_edges,q,corners,
      closest.triangle,static_cast<int>(closest.feature),closest.feature_id);
  m::Corner edge{closest.triangle,closest.feature_id};
  const auto ring=m::one_ring(faces,neighbors,neighbor_edges,edge);
  for(std::size_t i=1;i<ring.size();++i)
  {
    const auto &candidate=ring[i];
    const auto a=m::corner(vertices,faces,candidate.face,(candidate.corner+1)%3);
    const auto b=m::corner(vertices,faces,edge.face,(edge.corner+1)%3);
    if(k::perpendicular_compare(q,vertex,a,b)<0) edge=candidate;
  }
  const int fold_label=m::fold(vertices,faces,neighbors,neighbor_edges,edge);
  if(fold_label!=0) return fold_label;
  return -k::orientation(m::corner(vertices,faces,edge.face,0),
    m::corner(vertices,faces,edge.face,1),m::corner(vertices,faces,edge.face,2),q);
}
void audit_ordinary(const char *name,const Eigen::MatrixXd &queries,
  const Eigen::MatrixXd &vertices,const Eigen::MatrixXi &faces,const FlatTree &tree)
{
  std::size_t leaves=0,updates=0;
  for(Eigen::Index i=0;i<queries.rows();++i)
  {
    const Point q=queries.row(i);
    double best=std::numeric_limits<double>::infinity(),result;
    igl::eytzinger_aabb_sdf<true>(q,[&](int face)
    {
      ++leaves;
      double squared; Point closest;
      igl::point_simplex_squared_distance<3>(q,vertices,faces,face,squared,closest);
      if(squared<best) { ++updates; best=squared; }
      return squared;
    },tree.lower,tree.upper,tree.leaves,result);
  }
  std::cout<<"# audit "<<name<<" leaves="<<leaves<<" updates="<<updates<<'\n';
}
template<bool ExactCompare,bool IntervalUpper,bool KeepTieUpper=false,
  search_profile::FeatureCache Cache=search_profile::FeatureCache::None>
void audit(const char *name,const Eigen::MatrixXd &queries,
  const Eigen::MatrixXd &vertices,const Eigen::MatrixXi &faces,const FlatTree &tree,
  const std::vector<Output> &reference)
{
  search_profile::Counts counts;
  std::size_t different=0,worse=0;
  for(Eigen::Index i=0;i<queries.rows();++i)
  {
    const Point query=queries.row(i);
    const int found=search_profile::search<ExactCompare,IntervalUpper,true,KeepTieUpper,Cache>(
      query,vertices,faces,tree.lower,tree.upper,tree.leaves,counts);
    if(found!=reference[i].triangle)
    {
      ++different;
      if(k::compare(k::point(query),triangle(vertices,faces,found),
          triangle(vertices,faces,reference[i].triangle))!=0) ++worse;
    }
  }
  std::cout<<"# audit "<<name<<" popped_nodes="<<counts.nodes<<" leaves="<<counts.leaves
    <<" updates="<<counts.updates<<" comparisons="<<counts.comparisons
    <<" ties="<<counts.ties<<" shared_vertex_ties="<<counts.shared_vertex_ties
    <<" shared_edge_ties="<<counts.shared_edge_ties<<" tie_updates="<<counts.tie_updates
    <<" compare_interval_fallbacks="<<counts.compare_intervals
    <<" compare_expansion_fallbacks="<<counts.compare_expansions
    <<" interval_upper_bounds="<<counts.upper_intervals
    <<" upper_feature_interval_fallbacks="<<counts.upper_feature_intervals
    <<" upper_feature_expansion_fallbacks="<<counts.upper_feature_expansions
    <<" different_triangle="<<different<<" non_tie="<<worse<<'\n';
  if constexpr(ExactCompare && IntervalUpper)
    if(different) throw std::runtime_error("Audit copy disagrees with production exact traversal");
}

void verify_search_prototypes()
{
  p::exactinit();
  // Include shared features, coordinate-coincident vertices with distinct
  // indices, reversed face orderings, near ties and exact on-surface queries.
  // Check the deterministic winner against a brute-force exact scan.
  const double epsilon=std::ldexp(1.0,-52),far=std::ldexp(1.0,50);
  const std::array<Point,10> queries{{Point(-1,-1,-1),Point(.5,-1,-1),Point(.25,.25,-1),
    Point(0,0,0),Point(.5,0,0),Point(.25,.25,0),Point(.25,.25,epsilon/2),
    Point(far,far,far),Point(.5,epsilon,-1),Point(0,0,epsilon)}};
  std::size_t checks=0;
  for(double separation:{0.0,epsilon,0.5})
  {
    Eigen::MatrixXd vertices(8,3);
    vertices.topRows(4)<<0,0,0, 1,0,0, 0,1,0, 0,0,1;
    vertices.bottomRows(4)=vertices.topRows(4);
    vertices.bottomRows(4).col(2).array()+=separation;
    Eigen::MatrixXi faces(8,3);
    faces.topRows(4)<<0,2,1, 0,1,3, 0,3,2, 1,2,3;
    for(int i=0;i<4;++i) for(int j=0;j<3;++j) faces(i+4,j)=faces(i,2-j)+4;
    FlatTree tree; tree.init(vertices,faces);
    for(const auto &q:queries)
    {
      int expected=0;
      for(int face=1;face<faces.rows();++face)
        if(k::compare(k::point(q),triangle(vertices,faces,face),triangle(vertices,faces,expected))<0)
          expected=face;
      search_profile::Counts unused;
      const std::array<Eigen::Index,6> found{{
        p::exact_point_mesh_closest_triangle(q,vertices,faces,tree.lower,tree.upper,tree.leaves),
        search_profile::search<true,true,false,true>(q,vertices,faces,tree.lower,tree.upper,tree.leaves,unused),
        search_profile::search<true,true,false,false,search_profile::FeatureCache::Coordinates>(q,vertices,faces,tree.lower,tree.upper,tree.leaves,unused),
        search_profile::search<true,true,false,false,search_profile::FeatureCache::Topology>(q,vertices,faces,tree.lower,tree.upper,tree.leaves,unused),
        search_profile::search<true,true,false,true,search_profile::FeatureCache::Coordinates>(q,vertices,faces,tree.lower,tree.upper,tree.leaves,unused),
        search_profile::search<true,true,false,true,search_profile::FeatureCache::Topology>(q,vertices,faces,tree.lower,tree.upper,tree.leaves,unused)}};
      for(auto value:found)
      {
        ++checks;
        if(value!=expected) throw std::runtime_error("Search prototype disagrees with brute-force exact scan");
      }
    }
  }
  std::cout<<"# prototype_brute_force_checks="<<checks<<'\n';
}
}

int main(int argc,char *argv[])
{
  if(argc<2 || argc>5 || std::string(argv[1])=="--help")
  {
    std::cout<<"Usage: "<<argv[0]<<" mesh.obj [points=10000] [repetitions=7] [threads=1]\n";
    return argc<2 || argc>5 ? 1 : 0;
  }
  try
  {
    const int count=argc>2 ? positive_argument(argv[2]) : 10000;
    const int repetitions=argc>3 ? positive_argument(argv[3]) : 7;
    const int threads=argc>4 ? positive_argument(argv[4]) : 1;
#ifdef PROFILE_EXPANSIONS
    if(threads!=1) throw std::invalid_argument("The expansion audit must run with one thread");
#endif
    igl::default_num_threads(threads);
    Eigen::setNbThreads(1);
    verify_search_prototypes();
    Eigen::MatrixXd vertices; Eigen::MatrixXi faces;
    if(!igl::read_triangle_mesh(argv[1],vertices,faces) || vertices.rows()==0 || faces.rows()==0 ||
       vertices.cols()!=3 || faces.cols()!=3 || !vertices.allFinite() ||
       faces.minCoeff()<0 || faces.maxCoeff()>=vertices.rows())
      throw std::invalid_argument("Expected a nonempty finite triangle mesh with valid indices");
    Tree aabb; aabb.init(vertices,faces);
    FlatTree native,serialized;
    native.init(vertices,faces);
    aabb.serialize(serialized.lower,serialized.upper,serialized.leaves);
    Tree restored; restored.init(vertices,faces,native.lower,native.upper,native.leaves);
    igl::FastWindingNumberBVH winding;
    igl::fast_winding_number(vertices.cast<float>().eval(),faces,2,winding);
    Eigen::MatrixXi neighbors,neighbor_edges;
    igl::triangle_triangle_adjacency(faces,neighbors,neighbor_edges);
    p::MeshFeatureExposure exposure;
    p::exact_mesh_feature_exposure(vertices,faces,neighbors,neighbor_edges,exposure);
    const Point minimum=vertices.colwise().minCoeff(),extent=vertices.colwise().maxCoeff()-minimum;
    Eigen::MatrixXd queries(count,3);
    std::mt19937_64 generator(896);
    for(int i=0;i<count;++i) for(int axis=0;axis<3;++axis)
    {
      const double unit=(generator()>>11)*(1.0/9007199254740992.0);
      queries(i,axis)=minimum(axis)+extent(axis)*(-.1+1.2*unit);
    }
    std::vector<Method> methods;
    const auto add=[&](const char *name,auto evaluate)
    { methods.push_back({name,evaluate,std::vector<Output>(count),{}}); };
    add("ordinary_aabb",[&](const Point &q,Output &o)
    { o.squared=aabb.squared_distance(vertices,faces,q,o.triangle,o.closest); });
    add("ordinary_eytzinger_aabb_hierarchy",[&](const Point &q,Output &o)
    { ordinary_flat(q,vertices,faces,serialized,o); });
    add("ordinary_eytzinger_native",[&](const Point &q,Output &o)
    { ordinary_flat(q,vertices,faces,native,o); });
    add("ordinary_aabb_eytzinger_hierarchy",[&](const Point &q,Output &o)
    { o.squared=restored.squared_distance(vertices,faces,q,o.triangle,o.closest); });
    add("ablation_rounded_compare_rounded_upper",[&](const Point &q,Output &o)
    { search_profile::Counts unused; o.triangle=search_profile::search<false,false>(q,vertices,faces,native.lower,native.upper,native.leaves,unused); });
    add("ablation_exact_compare_rounded_upper",[&](const Point &q,Output &o)
    { search_profile::Counts unused; o.triangle=search_profile::search<true,false>(q,vertices,faces,native.lower,native.upper,native.leaves,unused); });
    add("ablation_rounded_compare_interval_upper",[&](const Point &q,Output &o)
    { search_profile::Counts unused; o.triangle=search_profile::search<false,true>(q,vertices,faces,native.lower,native.upper,native.leaves,unused); });
    add("ablation_exact_compare_interval_upper",[&](const Point &q,Output &o)
    { search_profile::Counts unused; o.triangle=search_profile::search<true,true>(q,vertices,faces,native.lower,native.upper,native.leaves,unused); });
    add("prototype_keep_tie_upper",[&](const Point &q,Output &o)
    { search_profile::Counts unused; o.triangle=search_profile::search<true,true,false,true>(q,vertices,faces,native.lower,native.upper,native.leaves,unused); });
    add("prototype_cache_features_coordinates",[&](const Point &q,Output &o)
    { search_profile::Counts unused; o.triangle=search_profile::search<true,true,false,false,search_profile::FeatureCache::Coordinates>(q,vertices,faces,native.lower,native.upper,native.leaves,unused); });
    add("prototype_cache_features_topology",[&](const Point &q,Output &o)
    { search_profile::Counts unused; o.triangle=search_profile::search<true,true,false,false,search_profile::FeatureCache::Topology>(q,vertices,faces,native.lower,native.upper,native.leaves,unused); });
    add("prototype_cache_features_coordinates_keep_tie_upper",[&](const Point &q,Output &o)
    { search_profile::Counts unused; o.triangle=search_profile::search<true,true,false,true,search_profile::FeatureCache::Coordinates>(q,vertices,faces,native.lower,native.upper,native.leaves,unused); });
    add("prototype_cache_features_topology_keep_tie_upper",[&](const Point &q,Output &o)
    { search_profile::Counts unused; o.triangle=search_profile::search<true,true,false,true,search_profile::FeatureCache::Topology>(q,vertices,faces,native.lower,native.upper,native.leaves,unused); });
    const int exact_method=methods.size();
    add("exact_search",[&](const Point &q,Output &o)
    { o.triangle=p::exact_point_mesh_closest_triangle(q,vertices,faces,native.lower,native.upper,native.leaves); });
    add("closest_point",[&](const Point &q,Output &o)
    {
      o.exact=p::point_mesh_closest_point(q,vertices,faces,native.lower,native.upper,native.leaves);
      o.triangle=o.exact.triangle; o.squared=o.exact.squared_distance;
    });
    add("exact_sign_cached",[&](const Point &q,Output &o)
    {
      o.distance=p::approximate_signed_distance_exact_sign(q,vertices,faces,native.lower,native.upper,native.leaves,exposure,o.exact,o.location);
      o.triangle=o.exact.triangle; o.squared=o.exact.squared_distance;
    });
    add("exact_sign_uncached",[&](const Point &q,Output &o)
    {
      o.distance=p::approximate_signed_distance_exact_sign(q,vertices,faces,native.lower,native.upper,native.leaves,neighbors,neighbor_edges,o.exact,o.location);
      o.triangle=o.exact.triangle; o.squared=o.exact.squared_distance;
    });
    add("prototype_uncached_skip_self",[&](const Point &q,Output &o)
    {
      m::dimensions(vertices,faces,neighbors,neighbor_edges);
      o.exact=p::point_mesh_closest_point(q,vertices,faces,native.lower,native.upper,native.leaves);
      o.location=static_cast<p::PointMeshLocation>(location_skip_self(q,vertices,faces,neighbors,neighbor_edges,o.exact));
      o.distance=o.location==p::PointMeshLocation::ON_SURFACE ? 0.0 : static_cast<int>(o.location)*o.exact.distance;
      o.triangle=o.exact.triangle; o.squared=o.exact.squared_distance;
    });
    add("fast_winding",[&](const Point &q,Output &o)
    {
      o.squared=aabb.squared_distance(vertices,faces,q,o.triangle,o.closest);
      o.distance=std::sqrt(o.squared)*(1-2*std::abs(igl::fast_winding_number(winding,2,q.cast<float>().eval())));
    });
    std::cout<<std::fixed<<std::setprecision(6)
      <<"# mesh="<<argv[1]<<" vertices="<<vertices.rows()<<" faces="<<faces.rows()
      <<" points="<<count<<" seed=896 box_padding_per_side=0.1 threads="<<threads
      <<" min_parallel=10000 repetitions="<<repetitions<<'\n';
#ifndef PROFILE_EXPANSIONS
    const auto run=[&](Method &method)
    {
      const bool parallel=parallel_queries(queries,[&](Eigen::Index i,const Point &q)
      { method.evaluate(q,method.output[i]); });
      if(parallel!=(count>=10000 && threads>1))
        throw std::runtime_error("Unexpected parallel_for scheduling");
    };
    for(auto &method:methods) run(method);
    long long checksum=0;
    for(int repetition=0;repetition<repetitions;++repetition)
      for(std::size_t offset=0;offset<methods.size();++offset)
      {
        auto &method=methods[(offset+repetition)%methods.size()];
        const auto start=std::chrono::steady_clock::now(); run(method);
        method.samples.push_back(std::chrono::duration<double,std::milli>(
          std::chrono::steady_clock::now()-start).count());
        for(const auto &output:method.output) checksum+=output.triangle;
      }
    std::cout<<"method,median_ms,min_ms,max_ms,us_per_point\n";
    for(auto &method:methods)
    {
      std::sort(method.samples.begin(),method.samples.end());
      const double median=method.samples[method.samples.size()/2];
      std::cout<<method.name<<','<<median<<','<<method.samples.front()<<','<<method.samples.back()
        <<','<<1000*median/count<<'\n';
    }
    std::cout<<"# triangle_checksum="<<checksum<<'\n';
#else
    for(auto &method:methods)
    {
      std::size_t sums=0,queries_with_expansions=0;
      for(int i=0;i<count;++i)
      {
        expansion_sums=0; const Point q=queries.row(i);
        method.evaluate(q,method.output[i]);
        sums+=expansion_sums; queries_with_expansions+=expansion_sums>0;
      }
      std::cout<<"# expansions "<<method.name<<" sums="<<sums
        <<" queries="<<queries_with_expansions<<'\n';
    }
#endif
    const auto &reference=methods[exact_method].output;
    const auto &signed_reference=methods[exact_method+2].output;
    for(const auto &method:methods)
    {
      std::size_t different=0,non_tie=0; double max_squared_difference=0;
      for(int i=0;i<count;++i)
      {
        const auto &output=method.output[i];
        if(output.triangle<0 || output.triangle>=faces.rows() || !std::isfinite(output.squared) ||
           !std::isfinite(output.distance)) throw std::runtime_error("Invalid query output");
        if(method.name.rfind("exact_sign_",0)==0 || method.name=="prototype_uncached_skip_self")
          if(output.location!=signed_reference[i].location || output.distance!=signed_reference[i].distance)
            throw std::runtime_error("Cached, uncached or skip-self signed outputs disagree");
        if(output.triangle!=reference[i].triangle)
        {
          ++different;
          const Point q=queries.row(i);
          non_tie+=k::compare(k::point(q),triangle(vertices,faces,output.triangle),
            triangle(vertices,faces,reference[i].triangle))!=0;
        }
        if(method.name.rfind("ordinary_",0)==0)
          max_squared_difference=std::max(max_squared_difference,
            std::abs(output.squared-methods[0].output[i].squared));
      }
      std::cout<<"# agreement "<<method.name<<" different_triangle="<<different
        <<" non_tie="<<non_tie<<" max_squared_difference_from_aabb="
        <<std::scientific<<max_squared_difference<<std::fixed<<'\n';
      if((method.name=="ablation_exact_compare_interval_upper" || method.name.rfind("prototype_",0)==0) && different)
        throw std::runtime_error("Timing copy disagrees with production exact traversal");
      if(non_tie) throw std::runtime_error("A sampled approximate result is not an exact nearest triangle");
    }
    audit<false,false>("rounded_rounded",queries,vertices,faces,native,reference);
    audit<true,false>("exact_rounded",queries,vertices,faces,native,reference);
    audit<false,true>("rounded_interval",queries,vertices,faces,native,reference);
    audit<true,true>("exact_interval",queries,vertices,faces,native,reference);
    audit<true,true,true>("keep_tie_upper",queries,vertices,faces,native,reference);
    audit<true,true,false,search_profile::FeatureCache::Coordinates>("cache_coordinates",queries,vertices,faces,native,reference);
    audit<true,true,false,search_profile::FeatureCache::Topology>("cache_topology",queries,vertices,faces,native,reference);
    audit<true,true,true,search_profile::FeatureCache::Coordinates>("cache_coordinates_keep_tie_upper",queries,vertices,faces,native,reference);
    audit<true,true,true,search_profile::FeatureCache::Topology>("cache_topology_keep_tie_upper",queries,vertices,faces,native,reference);
    audit_ordinary("ordinary_eytzinger_native",queries,vertices,faces,native);
    audit_ordinary("ordinary_eytzinger_aabb_hierarchy",queries,vertices,faces,serialized);
    std::cout<<"# cached_uncached_skip_self_signed_outputs_identical=1\n";
  }
  catch(const std::exception &error) { std::cerr<<error.what()<<'\n'; return 1; }
}
