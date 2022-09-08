#pragma once

// These are not directly used but would otherwise be included in most files.
// Leaving them included here.
#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>

#include <igl/find.h>

#include <Eigen/Core>
#include <catch2/catch.hpp>

#include <cctype>
#include <string>
#include <functional>
#include <algorithm>
#include <tuple>


// Disable lengthy tests in debug mode
#ifdef NDEBUG
#define IGL_DEBUG_OFF ""
#else
#define IGL_DEBUG_OFF "[!hide]"
#endif


#include <igl/STR.h>
template<>
struct Catch::StringMaker<std::tuple<int,int,int> > 
{
  static std::string convert(std::tuple<int,int,int> const& t)
  {
    return 
      STR("("<<std::get<0>(t)<<","<<std::get<1>(t)<<","<<std::get<2>(t)<<")");
  }
};
template<>
struct Catch::StringMaker<std::tuple<int,int,double> > 
{
  static std::string convert(std::tuple<int,int,double> const& t)
  {
    return 
      STR("("<<std::get<0>(t)<<","<<std::get<1>(t)<<","<<std::get<2>(t)<<")");
  }
};


namespace test_common
{
  template<typename Param, typename Fun>
  void run_test_cases(const std::vector<Param> &params,  Fun test_case)
  {
    for(const auto &p : params)
    {
      // Can't use INFO( p ) because we're not sure how to print p
      test_case(p);
    }
  }

  template<typename Fun>
  void run_test_cases(const std::vector<std::string> &params,  Fun test_case)
  {
    for(const auto &p : params)
    {
      INFO( p );
      test_case(p);
    }
  }

  inline std::vector<std::string> closed_genus_0_meshes()
  {
    return
    {
      "cube.obj",
      "decimated-knight.obj",
      "boolean_minus_test_cube.obj",
      "boolean_minus_test_green.obj",
    };
  };
  inline std::vector<std::string> closed_manifold_meshes()
  {
    std::vector<std::string> meshes = closed_genus_0_meshes();
    meshes.insert(meshes.end(),
    {
      "TinyTorus.obj",
    });
    return meshes;
  };
  inline std::vector<std::string> manifold_meshes()
  {
    std::vector<std::string> meshes = closed_manifold_meshes();
    meshes.insert(meshes.end(),
    {
      "bunny.off",
      "elephant.off",
      "hemisphere.obj",
    });
    return meshes;
  };
  inline std::vector<std::string> tet_meshes()
  {
    return
    {
      "decimated-knight.mesh"
    };
  };
  inline std::vector<std::string> all_meshes()
  {
    std::vector<std::string> meshes = manifold_meshes();
    meshes.insert(meshes.end(),
    {
      "truck.obj",
    });
    return meshes;
  };
  inline std::string data_path(std::string s)
  {
    return std::string(LIBIGL_DATA_DIR) + "/" + s;
  };

  template <typename DerivedA, typename DerivedB>
  void assert_eq(
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedB> & B)
  {
    // Sizes should match
    REQUIRE(A.rows() == B.rows());
    REQUIRE(A.cols() == B.cols());
    for(int i = 0;i<A.rows();i++)
    {
      for(int j = 0;j<A.cols();j++)
      {
        // Create an ijv tuple to trick GoogleTest into printing (i,j) so we
        // know where the disagreement is.
        std::tuple<int,int,typename DerivedA::Scalar> Aijv {i,j,A(i,j)};
        std::tuple<int,int,typename DerivedB::Scalar> Bijv {i,j,B(i,j)};
        REQUIRE(Aijv == Bijv);
      }
    }
  }

  template <typename DerivedA, typename DerivedB>
  void assert_neq(
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedB> & B)
  {
    // Sizes should match
    REQUIRE(A.rows() == B.rows());
    REQUIRE(A.cols() == B.cols());
    bool all_equals = true;
    for(int i = 0;i<A.rows();i++)
    {
      for(int j = 0;j<A.cols();j++)
      {
        if (A(i,j) != B(i,j))
        {
          all_equals = false;
        }
      }
    }
    REQUIRE_FALSE(all_equals);
  }

  template <typename DerivedA, typename DerivedB>
  void assert_eq(
    const Eigen::SparseMatrix<DerivedA> & A,
    const Eigen::SparseMatrix<DerivedB> & B)
  {
    // Sizes should match
    REQUIRE(A.rows() == B.rows());
    REQUIRE(A.cols() == B.cols());
    Eigen::Matrix<long int,Eigen::Dynamic, 1> AI,AJ;
    Eigen::Matrix<long int,Eigen::Dynamic, 1> BI,BJ;
    Eigen::Matrix<DerivedA,Eigen::Dynamic, 1> AV;
    Eigen::Matrix<DerivedB,Eigen::Dynamic, 1> BV;
    // Assumes A and B are in same Major Ordering
    igl::find(A,AI,AJ,AV);
    igl::find(B,BI,BJ,BV);
    // This doesn't generalized to assert_near nicely, and it makes it hard to
    // tell which entries are different:
    assert_eq(AI,BI);
    assert_eq(AJ,BJ);
    assert_eq(AV,BV);
  }

  template <typename DerivedA, typename DerivedB, typename EpsType>
  void assert_near(
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedB> & B,
    const EpsType & eps)
  {
    // Sizes should match
    REQUIRE(A.rows() == B.rows());
    REQUIRE(A.cols() == B.cols());
    for(int i = 0;i<A.rows();i++)
    {
      for(int j = 0;j<A.cols();j++)
      {
        // Create an ijv tuple to trick GoogleTest into printing (i,j) so we
        // know where the disagreement is.
        //
        // Equivalent to ASSERT_NEAR(Aijv,Bijv)

        CAPTURE( i );
        CAPTURE( j );
        {
          // std::tuple<int,int,typename DerivedA::Scalar> Aijv {i,j,A(i,j)};
          // std::tuple<int,int,typename DerivedB::Scalar> Bijv {i,j,B(i,j)+eps};
          REQUIRE(A(i,j) < B(i,j)+eps);
        }
        {
          // std::tuple<int,int,typename DerivedA::Scalar> Aijv {i,j,A(i,j)+eps};
          // std::tuple<int,int,typename DerivedB::Scalar> Bijv {i,j,B(i,j)};
          REQUIRE(A(i,j)+eps > B(i,j));
        }
      }
    }
  }

}
