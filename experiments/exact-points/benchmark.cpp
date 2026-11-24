#include <igl/icosahedron.h>
#include <igl/upsample.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/eytzinger_aabb.h>
#include <igl/eytzinger_aabb_sdf.h>
#include <igl/point_simplex_squared_distance.h>
#include <igl/predicates/exact_point_mesh_closest_triangle.h>
#include <igl/predicates/exact_point_mesh_classify.h>
#include <igl/predicates/approximate_signed_distance_exact_sign.h>
#include <igl/predicates/exact_point_triangle_closest_feature.h>
#include <igl/predicates/exact_point_triangle_squared_distance_compare.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

using Point=Eigen::RowVector3d;
namespace p=igl::predicates;
volatile long long checksum=0;
template<typename Function> void measure(const char *name,int faces,int queries,const Function &function)
{
  std::vector<double> samples;
  for(int sample=0;sample<7;++sample)
  {
    const auto start=std::chrono::steady_clock::now();
    int iterations=0;
    double seconds;
    do
    {
      checksum+=function(); ++iterations;
      seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }while(seconds<.02);
    samples.push_back(seconds*1e9/(iterations*queries));
  }
  std::sort(samples.begin(),samples.end());
  std::cout<<name<<","<<faces<<","<<samples[3]<<","<<samples.front()<<","<<samples.back()<<"\n";
}
int main()
{
  std::cout<<"operation,faces,median_ns,min_ns,max_ns\n";
  std::mt19937 generator(896);
  std::uniform_real_distribution<double> sample(-1.5,1.5);
  for(int subdivisions:{2,4})
  {
    Eigen::MatrixXd vertices; Eigen::MatrixXi faces;
    igl::icosahedron(vertices,faces);
    for(int i=0;i<subdivisions;++i)
    {
      igl::upsample(vertices,faces);
      vertices.rowwise().normalize();
    }
    Eigen::MatrixXd minimum(faces.rows(),3),maximum(faces.rows(),3),lower,upper;
    Eigen::VectorXi leaves;
    for(int f=0;f<faces.rows();++f)
    {
      minimum.row(f)=vertices.row(faces(f,0)); maximum.row(f)=minimum.row(f);
      for(int i=1;i<3;++i)
      {
        minimum.row(f)=minimum.row(f).cwiseMin(vertices.row(faces(f,i)));
        maximum.row(f)=maximum.row(f).cwiseMax(vertices.row(faces(f,i)));
      }
    }
    igl::eytzinger_aabb(minimum,maximum,lower,upper,leaves);
    Eigen::MatrixXi neighbors,neighbor_edges;
    igl::triangle_triangle_adjacency(faces,neighbors,neighbor_edges);
    p::MeshFeatureExposure exposure;
    p::exact_mesh_feature_exposure(vertices,faces,neighbors,neighbor_edges,exposure);
    std::vector<Point> queries;
    std::vector<Eigen::Index> closest;
    for(int i=0;i<256;++i)
    {
      queries.emplace_back(sample(generator),sample(generator),sample(generator));
      closest.push_back(p::exact_point_mesh_closest_triangle(queries.back(),vertices,faces,lower,upper,leaves));
    }
    const auto approximate=[&](const Point &q,int f)
    {
      double squared; Point point;
      igl::point_simplex_squared_distance<3>(q,vertices,faces,f,squared,point);
      return squared;
    };
    measure("approximate_tree",faces.rows(),queries.size(),[&]()
    {
      long long result=0;
      for(const Point &q:queries)
      {
        double squared;
        igl::eytzinger_aabb_sdf<true>(q,[&](int f){return approximate(q,f);},lower,upper,leaves,squared);
        result+=squared*1e6;
      }
      return result;
    });
    measure("exact_tree",faces.rows(),queries.size(),[&]()
    {
      long long result=0;
      for(const Point &q:queries) result+=p::exact_point_mesh_closest_triangle(q,vertices,faces,lower,upper,leaves);
      return result;
    });
    measure("closest_point_wrapper",faces.rows(),queries.size(),[&]()
    {
      long long result=0;
      for(const Point &q:queries)
      {
        const auto closest_point=p::point_mesh_closest_point(q,vertices,faces,lower,upper,leaves);
        result+=closest_point.triangle+closest_point.distance*1e6+closest_point.closest_point.sum()*1e6;
      }
      return result;
    });
    Eigen::MatrixXd query_matrix(queries.size(),3);
    for(std::size_t i=0;i<queries.size();++i) query_matrix.row(i)=queries[i];
    Eigen::VectorXd signed_distances; Eigen::VectorXi locations;
    std::vector<p::PointMeshClosestPoint> closest_points;
    measure("signed_distance_batch",faces.rows(),queries.size(),[&]()
    {
      p::approximate_signed_distance_exact_sign(query_matrix,vertices,faces,lower,upper,leaves,
        exposure,signed_distances,closest_points,locations);
      long long result=0;
      for(Eigen::Index i=0;i<signed_distances.size();++i)
        result+=signed_distances(i)*1e6+locations(i)+closest_points[i].triangle;
      return result;
    });
    for(bool cached:{false,true})
    {
      measure(cached ? "cached_sign" : "local_sign",faces.rows(),queries.size(),[&]()
      {
        long long result=0;
        for(std::size_t i=0;i<queries.size();++i)
          result+=int(cached ? p::exact_point_mesh_classify(queries[i],vertices,faces,closest[i],exposure) :
            p::exact_point_mesh_classify(queries[i],vertices,faces,closest[i],neighbors,neighbor_edges));
        return result;
      });
    }
    measure("closest_feature",faces.rows(),queries.size(),[&]()
    {
      long long result=0;
      for(std::size_t i=0;i<queries.size();++i)
      {
        const int f=closest[i]; int id;
        const Point a=vertices.row(faces(f,0)),b=vertices.row(faces(f,1)),c=vertices.row(faces(f,2));
        result+=int(p::exact_point_triangle_closest_feature(queries[i],a,b,c,id));
      }
      return result;
    });
    measure("adjacent_distance_compare",faces.rows(),queries.size(),[&]()
    {
      long long result=0;
      for(std::size_t i=0;i<queries.size();++i)
      {
        const int f=closest[i],g=neighbors(f,0);
        const Point a=vertices.row(faces(f,0)),b=vertices.row(faces(f,1)),c=vertices.row(faces(f,2));
        const Point d=vertices.row(faces(g,0)),e=vertices.row(faces(g,1)),h=vertices.row(faces(g,2));
        result+=p::exact_point_triangle_squared_distance_compare(queries[i],a,b,c,d,e,h);
      }
      return result;
    });
    measure("exposure_preprocess",faces.rows(),1,[&]()
    {
      p::MeshFeatureExposure fresh;
      p::exact_mesh_feature_exposure(vertices,faces,neighbors,neighbor_edges,fresh);
      return fresh.vertices.sum();
    });
  }
  std::cerr<<"checksum="<<checksum<<"\n";
}
