// Profile unsigned, pseudonormal, fast-winding and exact-sign distances.
#include <igl/AABB.h>
#include <igl/centroid.h>
#include <igl/default_num_threads.h>
#include <igl/eytzinger_aabb.h>
#include <igl/box_simplices.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/parallel_for.h>
#include <igl/per_edge_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/read_triangle_mesh.h>
#include <igl/signed_distance.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/predicates/approximate_signed_distance_exact_sign.h>
#include <igl/predicates/find_self_intersections.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
namespace predicates=igl::predicates;
using Point=Eigen::RowVector3d;
using Clock=std::chrono::steady_clock;
constexpr int min_parallel=10000;
constexpr unsigned int seed=896;
constexpr int method_count=5;
const std::array<const char *,method_count> names{{
  "point_mesh_squared_distance", "signed_distance_pseudonormal",
  "signed_distance_fast_winding_number", "exact_sign_uncached", "exact_sign_cached"}};

struct EytzingerTree
{
  Eigen::MatrixXd lower,upper;
  Eigen::VectorXi leaves;
  void init(const Eigen::MatrixXd &vertices,const Eigen::MatrixXi &faces)
  {
    Eigen::MatrixXd minimum,maximum;
    igl::box_simplices(vertices,faces,minimum,maximum);
    igl::eytzinger_aabb(minimum,maximum,lower,upper,leaves);
  }
};
struct Results
{
  Eigen::VectorXd distances;
  Eigen::VectorXi triangles;
  Eigen::MatrixXd closest_points,normals,barycentrics;
};

template<typename Evaluate> bool parallel_queries(const Eigen::MatrixXd &queries,const Evaluate &evaluate)
{
  // All five prepared query paths use this same scheduling, threshold, thread
  // count and exception handling. The scalar callees do not launch workers.
  std::vector<std::exception_ptr> errors;
  return igl::parallel_for(queries.rows(),
    [&](std::size_t threads) { errors.resize(threads); },
    [&](Eigen::Index i,std::size_t thread)
    {
      try { const Point query=queries.row(i); evaluate(i,query); }
      catch(...) { if(!errors[thread]) errors[thread]=std::current_exception(); }
    },
    [&](std::size_t thread) { if(errors[thread]) std::rethrow_exception(errors[thread]); },
    min_parallel);
}
template<typename Function> double elapsed_ms(const Function &function)
{
  const auto start=Clock::now(); function();
  return std::chrono::duration<double,std::milli>(Clock::now()-start).count();
}
void print_timing(const char *phase,const char *method,std::vector<double> samples,int query_count)
{
  std::sort(samples.begin(),samples.end());
  const double median=samples[samples.size()/2];
  std::cout<<phase<<','<<method<<','<<median<<','<<samples.front()<<','<<samples.back()
    <<','<<median*1000/query_count<<'\n';
}
void profile(const char *phase,const std::array<std::function<void()>,method_count> &methods,
  const std::array<Results,method_count> &results,int repetitions,int query_count)
{
  for(const auto &method:methods) method(); // Warm data, code and output allocations.
  std::array<std::vector<double>,method_count> samples;
  double checksum=0;
  for(int repetition=0;repetition<repetitions;++repetition)
    for(int offset=0;offset<method_count;++offset)
    {
      const int method=(offset+repetition)%method_count;
      samples[method].push_back(elapsed_ms(methods[method]));
      if(!results[method].distances.allFinite())
        throw std::runtime_error(std::string("Nonfinite distances from ")+names[method]);
      checksum+=results[method].distances.sum(); // Consume outputs outside the timing.
    }
  for(int method=0;method<method_count;++method)
    print_timing(phase,names[method],samples[method],query_count);
  std::cout<<"# "<<phase<<" checksum="<<std::setprecision(17)<<checksum<<std::setprecision(6)<<'\n';
}
int positive_argument(const char *text)
{
  std::size_t consumed=0;
  const long long value=std::stoll(text,&consumed);
  if(text[consumed]!='\0' || value<1 || value>std::numeric_limits<int>::max())
    throw std::invalid_argument("Counts must be positive integers fitting in int");
  return static_cast<int>(value);
}
}

int main(int argc,char *argv[])
{
  if(argc>5 || (argc>1 && std::string(argv[1])=="--help"))
  {
    std::cout<<"Usage: "<<argv[0]<<" [mesh=armadillo.obj] [points=10000] [repetitions=5] [threads=1]\n";
    return argc>5 ? 1 : 0;
  }
  try
  {
    const std::string filename=argc>1 ? argv[1] : TUTORIAL_SHARED_PATH "/armadillo.obj";
    const int query_count=argc>2 ? positive_argument(argv[2]) : 10000;
    const int repetitions=argc>3 ? positive_argument(argv[3]) : 5;
    const int threads=argc>4 ? positive_argument(argv[4]) : 1;
    // Set this before any libigl call; the value is cached on first use.
    igl::default_num_threads(threads);
    Eigen::setNbThreads(1);
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    if(!igl::read_triangle_mesh(filename,vertices,faces) || vertices.rows()==0 || faces.rows()==0 ||
       vertices.cols()!=3 || faces.cols()!=3 || !vertices.allFinite() ||
       faces.minCoeff()<0 || faces.maxCoeff()>=vertices.rows())
      throw std::invalid_argument("Expected a nonempty finite triangle mesh with valid indices");
    if(!igl::is_edge_manifold(faces) || !igl::is_vertex_manifold(faces))
      throw std::invalid_argument("Expected an edge- and vertex-manifold mesh");
    Point centroid; double volume;
    igl::centroid(vertices,faces,centroid,volume);
    if(!(volume>0)) throw std::invalid_argument("Expected an outward-oriented mesh with positive volume");
    Eigen::MatrixXi intersections; Eigen::Array<bool,Eigen::Dynamic,1> coplanar;
    if(predicates::find_self_intersections(vertices,faces,true,intersections,coplanar))
      throw std::invalid_argument("Expected a mesh without self-intersections");

    const Point minimum=vertices.colwise().minCoeff(),maximum=vertices.colwise().maxCoeff();
    const Point extent=maximum-minimum;
    Eigen::MatrixXd queries(query_count,3);
    std::mt19937_64 generator(seed);
    for(int i=0;i<query_count;++i) for(int axis=0;axis<3;++axis)
    {
      const double unit=(generator()>>11)*(1.0/9007199254740992.0);
      queries(i,axis)=minimum(axis)+extent(axis)*(-.1+1.2*unit);
    }
    std::cout<<std::setprecision(6)<<std::fixed
      <<"# mesh="<<filename<<" vertices="<<vertices.rows()<<" faces="<<faces.rows()<<'\n'
      <<"# points="<<query_count<<" seed="<<seed<<" box_padding_per_side=0.1 repetitions="<<repetitions<<'\n'
      <<"# threads="<<igl::default_num_threads()<<" min_parallel="<<min_parallel
      <<" eigen_threads="<<Eigen::nbThreads()<<'\n';
#ifdef IGL_PARALLEL_FOR_FORCE_SERIAL
    std::cout<<"# IGL_PARALLEL_FOR_FORCE_SERIAL overrides the requested thread count\n";
#endif

    igl::AABB<Eigen::MatrixXd,3> tree;
    const double tree_ms=elapsed_ms([&]{tree.init(vertices,faces);});
    Eigen::MatrixXd face_normals,vertex_normals,edge_normals;
    Eigen::MatrixXi edges; Eigen::VectorXi edge_map;
    const double normals_ms=elapsed_ms([&]
    {
      igl::per_face_normals(vertices,faces,face_normals);
      igl::per_vertex_normals(vertices,faces,igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,face_normals,vertex_normals);
      igl::per_edge_normals(vertices,faces,igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,face_normals,edge_normals,edges,edge_map);
    });
    igl::FastWindingNumberBVH winding_tree;
    const double winding_ms=elapsed_ms([&]{igl::fast_winding_number(vertices.cast<float>().eval(),faces,2,winding_tree);});
    EytzingerTree exact_tree;
    const double exact_tree_ms=elapsed_ms([&]{exact_tree.init(vertices,faces);});
    Eigen::MatrixXi neighbors,neighbor_edges;
    const double adjacency_ms=elapsed_ms([&]{igl::triangle_triangle_adjacency(faces,neighbors,neighbor_edges);});
    if((neighbors.array()<0).any()) throw std::invalid_argument("Expected a watertight mesh without boundary edges");
    predicates::MeshFeatureExposure exposure;
    // Also validates nondegenerate faces and reciprocal, oriented adjacency.
    const double exposure_ms=elapsed_ms([&]{predicates::exact_mesh_feature_exposure(vertices,faces,neighbors,neighbor_edges,exposure);});

    std::cout<<"# validation=closed_manifold_oriented_nondegenerate_no_self_intersections\n"
      <<"# preprocessing rows sum separately timed components; mesh loading, validation and sampling excluded\n"
      <<"phase,method,median_ms,min_ms,max_ms,us_per_point\n";
    const std::array<double,method_count> preprocessing{{tree_ms,tree_ms+normals_ms,tree_ms+winding_ms,
      exact_tree_ms+adjacency_ms,exact_tree_ms+adjacency_ms+exposure_ms}};
    for(int method=0;method<method_count;++method)
      print_timing("preprocess",names[method],{preprocessing[method]},query_count);
    std::cout<<"# components_ms: aabb="<<tree_ms<<" normals="<<normals_ms<<" winding_bvh="<<winding_ms
      <<" eytzinger="<<exact_tree_ms<<" adjacency="<<adjacency_ms<<" exposure="<<exposure_ms<<'\n';

    std::array<Results,method_count> prepared,complete;
    std::array<bool,method_count> used_parallel{};
    for(int method=0;method<method_count;++method)
    {
      auto &result=prepared[method]; result.distances.resize(query_count);
      result.triangles.resize(query_count); result.closest_points.resize(query_count,3);
      if(method==1) result.normals.resize(query_count,3);
      if(method>=3) result.barycentrics.resize(query_count,3);
    }
    std::array<std::function<void()>,method_count> query_methods;
    query_methods[0]=[&]
    {
      auto &result=prepared[0];
      used_parallel[0]=parallel_queries(queries,[&](Eigen::Index i,const Point &query)
      {
        int triangle; Point point;
        result.distances(i)=tree.squared_distance(vertices,faces,query,triangle,point);
        result.triangles(i)=triangle; result.closest_points.row(i)=point;
      });
    };
    query_methods[1]=[&]
    {
      auto &result=prepared[1];
      used_parallel[1]=parallel_queries(queries,[&](Eigen::Index i,const Point &query)
      {
        double sign,squared; int triangle; Point point,normal;
        igl::signed_distance_pseudonormal(tree,vertices,faces,face_normals,vertex_normals,
          edge_normals,edge_map,query,sign,squared,triangle,point,normal);
        result.distances(i)=sign*std::sqrt(squared); result.triangles(i)=triangle;
        result.closest_points.row(i)=point; result.normals.row(i)=normal;
      });
    };
    query_methods[2]=[&]
    {
      auto &result=prepared[2];
      used_parallel[2]=parallel_queries(queries,[&](Eigen::Index i,const Point &query)
      {
        int triangle; Point point;
        const double squared=tree.squared_distance(vertices,faces,query,triangle,point);
        const double winding=igl::fast_winding_number(winding_tree,2,query.cast<float>().eval());
        // Preserve signed_distance's continuous factor rather than thresholding w.
        result.distances(i)=std::sqrt(squared)*(1.-2.*std::abs(winding));
        result.triangles(i)=triangle; result.closest_points.row(i)=point;
      });
    };
    for(int method=3;method<5;++method) query_methods[method]=[&,method]
    {
      auto &result=prepared[method];
      used_parallel[method]=parallel_queries(queries,[&](Eigen::Index i,const Point &query)
      {
        Eigen::Index triangle; Point point,bary;
        result.distances(i)=method==4 ? predicates::approximate_signed_distance_exact_sign(
          query,vertices,faces,exact_tree.lower,exact_tree.upper,exact_tree.leaves,exposure,triangle,point,bary) :
          predicates::approximate_signed_distance_exact_sign(query,vertices,faces,exact_tree.lower,
            exact_tree.upper,exact_tree.leaves,neighbors,neighbor_edges,triangle,point,bary);
        result.triangles(i)=static_cast<int>(triangle); result.closest_points.row(i)=point;
        result.barycentrics.row(i)=bary;
      });
    };
    profile("query",query_methods,prepared,repetitions,query_count);
    if(!std::all_of(used_parallel.begin(),used_parallel.end(),[&](bool value){return value==used_parallel[0];}))
      throw std::logic_error("Inconsistent query parallelism");
    std::cout<<"# all_prepared_methods_used_thread_pool="<<used_parallel[0]<<'\n';

    // Also time the actual convenience APIs including their preprocessing.
    // All five batch APIs use the same 10,000-query parallel_for threshold.
    std::array<std::function<void()>,method_count> full_methods;
    full_methods[0]=[&]
    { auto &r=complete[0]; igl::point_mesh_squared_distance(queries,vertices,faces,r.distances,r.triangles,r.closest_points); };
    for(int method=1;method<3;++method) full_methods[method]=[&,method]
    {
      auto &r=complete[method];
      igl::signed_distance(queries,vertices,faces,method==1 ? igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL :
        igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,r.distances,r.triangles,r.closest_points,r.normals);
    };
    for(int method=3;method<5;++method) full_methods[method]=[&,method]
    {
      auto &r=complete[method];
      EytzingerTree fresh_tree; fresh_tree.init(vertices,faces);
      Eigen::MatrixXi fresh_neighbors,fresh_edges;
      igl::triangle_triangle_adjacency(faces,fresh_neighbors,fresh_edges);
      if(method==4)
      {
        predicates::MeshFeatureExposure fresh_exposure;
        predicates::exact_mesh_feature_exposure(vertices,faces,fresh_neighbors,fresh_edges,fresh_exposure);
        predicates::approximate_signed_distance_exact_sign(queries,vertices,faces,fresh_tree.lower,
          fresh_tree.upper,fresh_tree.leaves,fresh_exposure,r.distances,r.triangles,r.closest_points,r.barycentrics);
      }else predicates::approximate_signed_distance_exact_sign(queries,vertices,faces,fresh_tree.lower,
        fresh_tree.upper,fresh_tree.leaves,fresh_neighbors,fresh_edges,r.distances,r.triangles,r.closest_points,r.barycentrics);
    };
    profile("total",full_methods,complete,repetitions,query_count);

    const auto &reference=prepared[4];
    if(prepared[3].distances!=reference.distances || prepared[3].triangles!=reference.triangles ||
       prepared[3].barycentrics!=reference.barycentrics)
      throw std::logic_error("Cached and uncached exact-sign results disagree");
    for(int method=0;method<method_count;++method)
    {
      const double api_difference=(prepared[method].distances-complete[method].distances).cwiseAbs().maxCoeff();
      const double scale=std::max(1.,prepared[method].distances.cwiseAbs().maxCoeff());
      if(api_difference>1e-12*scale)
        throw std::logic_error(std::string("Prepared and convenience API results disagree: ")+names[method]);
      int sign_mismatches=0; double magnitude_error=0;
      for(int i=0;i<query_count;++i)
      {
        const double value=prepared[method].distances(i);
        if(method>0 && ((value>0)-(value<0))!=((reference.distances(i)>0)-(reference.distances(i)<0))) ++sign_mismatches;
        const double magnitude=method==0 ? std::sqrt(value) : std::abs(value);
        magnitude_error=std::max(magnitude_error,std::abs(magnitude-std::abs(reference.distances(i))));
      }
      std::cout<<"# agreement,"<<names[method]<<",sign_mismatches=";
      if(method==0) std::cout<<"n/a"; else std::cout<<sign_mismatches;
      std::cout<<",max_abs_magnitude_difference="<<std::scientific<<magnitude_error
        <<",max_prepared_vs_api_difference="<<api_difference<<std::fixed<<'\n';
    }
    std::array<int,3> features{{0,0,0}},locations{{0,0,0}};
    for(int i=0;i<query_count;++i)
    {
      ++features[(reference.barycentrics.row(i).array()>0).count()-1];
      ++locations[((reference.distances(i)>0)-(reference.distances(i)<0))+1];
    }
    std::cout<<"# exact_locations: inside="<<locations[0]<<" surface="<<locations[1]<<" outside="<<locations[2]<<'\n'
      <<"# exact_features: vertex="<<features[0]<<" edge="<<features[1]<<" face="<<features[2]<<'\n';
    return 0;
  }catch(const std::exception &error)
  {
    std::cerr<<"Signed-distance profile failed: "<<error.what()<<'\n'; return 1;
  }
}
