// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_BOOLEAN_CSG_TREE_H
#define IGL_BOOLEAN_CSG_TREE_H

#include <igl/boolean/string_to_mesh_boolean_type.h>
#include <igl/boolean/MeshBooleanType.h>
#include <igl/boolean/mesh_boolean.h>

namespace igl
{
  namespace boolean
  {
    // Class for defining and computing a constructive solid geometry result
    // out of a tree of boolean operations on "solid" triangle meshes.
    template <typename DerivedF>
    class CSGTree
    {
      private:
        typedef CGAL::Epeck::FT ExactScalar;
        typedef Eigen::Matrix<ExactScalar,Eigen::Dynamic,3> MatrixX3E;
        typedef Eigen::PlainObjectBase<DerivedF> POBF;
        typedef Eigen::Matrix<typename DerivedF::Index,Eigen::Dynamic,1> 
          VectorJ;
        // Resulting mesh
        MatrixX3E m_V;
        POBF m_F;
        VectorJ m_J;
        // Number of birth faces in A + those in B. I.e. sum of original "leaf"
        // faces involved in result.
        size_t m_number_of_birth_faces;
      public:
        CSGTree()
        {
        }
        //typedef Eigen::MatrixXd MatrixX3E;
        //typedef Eigen::MatrixXi POBF;
        // http://stackoverflow.com/a/3279550/148668
        CSGTree(const CSGTree & other)
          // copy things
        {
        }
        // copy-swap idiom
        friend void swap(CSGTree& first, CSGTree& second)
        {
          using std::swap;
        }
        // Pass-by-value (aka copy)
        CSGTree& operator=(CSGTree other)
        {
          swap(*this,other);
          return *this;
        }
        CSGTree(CSGTree&& other):
          // initialize via default constructor
          CSGTree() 
        {
          swap(*this,other);
        }
        CSGTree(
          const CSGTree & A,
          const CSGTree & B,
          const MeshBooleanType & type)
        {
          // conduct boolean operation
          mesh_boolean(A.V(),A.F(),B.V(),B.F(),type,m_V,m_F,m_J);
          // reindex m_J
          std::for_each(m_J.data(),m_J.data()+m_J.size(),
            [&](typename VectorJ::Scalar & j)
            {
              if(j < A.F().rows())
              {
                j = A.J()(j);
              }else
              {
                assert(j<(A.F().rows()+B.F().rows()));
                j = A.number_of_birth_faces()+(B.J()(j-A.F().rows()));
              }
            });
          m_number_of_birth_faces = 
            A.number_of_birth_faces() + B.number_of_birth_faces();
        }
        // Overload using string
        CSGTree(
          const CSGTree & A,
          const CSGTree & B,
          const std::string & s):
          CSGTree(A,B,string_to_mesh_boolean_type(s))
        {
          // do nothing (all done in constructor).
        }
        // "Leaf" node with identity operation
        template <typename DerivedV>
        CSGTree(const Eigen::PlainObjectBase<DerivedV> & V, const POBF & F)//:
        // Possible Eigen bug:
        // https://forum.kde.org/viewtopic.php?f=74&t=128414
          //m_V(V.template cast<ExactScalar>()),m_F(F)
        {
          m_V = V.template cast<ExactScalar>();
          m_F = F;
          // number of faces
          m_number_of_birth_faces = m_F.rows();
          // identity birth index
          m_J = VectorJ::LinSpaced(
            m_number_of_birth_faces,0,m_number_of_birth_faces-1);
        }
        // Returns reference to resulting mesh vertices m_V
        const MatrixX3E & V() const
        {
          return m_V;
        }
        // Returns reference to resulting mesh faces m_F
        template <typename DerivedV>
        Eigen::PlainObjectBase<DerivedV> cast_V() const
        {
          Eigen::PlainObjectBase<DerivedV> dV;
          dV.resize(m_V.rows(),m_V.cols());
          for(int i = 0;i<m_V.size();i++)
          {
            *(dV.data()+i) = CGAL::to_double((*(m_V.data()+i)));
          }
          return dV;
        }
        const POBF & F() const
        {
          return m_F;
        }
        const VectorJ & J() const
        {
          return m_J;
        }
        const size_t & number_of_birth_faces() const
        {
          return m_number_of_birth_faces;
        }
    };
  }
}


#endif
