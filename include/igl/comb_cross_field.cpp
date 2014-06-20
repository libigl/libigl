#include "comb_cross_field.h"

#include <vector>
#include <deque>
#include "per_face_normals.h"
#include "is_border_vertex.h"
#include "tt.h"

namespace igl {
  template <typename DerivedV, typename DerivedF>
  class Comb
  {
  public:
    
    const Eigen::PlainObjectBase<DerivedV> &V;
    const Eigen::PlainObjectBase<DerivedF> &F;
    const Eigen::PlainObjectBase<DerivedV> &PD1;
    const Eigen::PlainObjectBase<DerivedV> &PD2;
    Eigen::PlainObjectBase<DerivedV> N;
    
  private:
    // internal
    Eigen::PlainObjectBase<DerivedF> TT;
    Eigen::PlainObjectBase<DerivedF> TTi;
    
    
  private:
    
    
    static double Sign(double a){return (double)((a>0)?+1:-1);}
    
    ///given 2 vector centered into origin calculate the rotation matrix from first to the second
    static Eigen::Matrix<typename DerivedV::Scalar, 3, 3> RotationMatrix(Eigen::Matrix<typename DerivedV::Scalar, 3, 1> v0,
                                                                         Eigen::Matrix<typename DerivedV::Scalar, 3, 1> v1,
                                                                         bool normalized=true)
    {
      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> rotM;
      const double epsilon=0.00001;
      if (!normalized)
      {
        v0.normalize();
        v1.normalize();
      }
      typename DerivedV::Scalar dot=v0.dot(v1);
      ///control if there is no rotation
      if (dot>((double)1-epsilon))
      {
        rotM = Eigen::Matrix<typename DerivedV::Scalar, 3, 3>::Identity();
        return rotM;
      }
      
      ///find the axis of rotation
      Eigen::Matrix<typename DerivedV::Scalar, 3, 1> axis;
      axis=v0.cross(v1);
      axis.normalize();
      
      ///construct rotation matrix
      typename DerivedV::Scalar u=axis(0);
      typename DerivedV::Scalar v=axis(1);
      typename DerivedV::Scalar w=axis(2);
      typename DerivedV::Scalar phi=acos(dot);
      typename DerivedV::Scalar rcos = cos(phi);
      typename DerivedV::Scalar rsin = sin(phi);
      
      rotM(0,0) =      rcos + u*u*(1-rcos);
      rotM(1,0) =  w * rsin + v*u*(1-rcos);
      rotM(2,0) = -v * rsin + w*u*(1-rcos);
      rotM(0,1) = -w * rsin + u*v*(1-rcos);
      rotM(1,1) =      rcos + v*v*(1-rcos);
      rotM(2,1) =  u * rsin + w*v*(1-rcos);
      rotM(0,2) =  v * rsin + u*w*(1-rcos);
      rotM(1,2) = -u * rsin + v*w*(1-rcos);
      rotM(2,2) =      rcos + w*w*(1-rcos);
      
      return rotM;
    }
    
    
  public:
    ///rotate a given vector from the tangent space
    ///of f0 to the tangent space of f1 by considering the difference of normals
    static Eigen::Matrix<typename DerivedV::Scalar, 3, 1> Rotate(Eigen::Matrix<typename DerivedV::Scalar, 3, 1> N0,
                                                                 Eigen::Matrix<typename DerivedV::Scalar, 3, 1> N1,
                                                                 const Eigen::Matrix<typename DerivedV::Scalar, 3, 1>& dir3D)
    {
      ///find the rotation matrix that maps between normals
      //    vcg::Matrix33<ScalarType> rotation=vcg::RotationMatrix(N0,N1);
      Eigen::Matrix<typename DerivedV::Scalar, 3, 3> rotation = RotationMatrix(N0,N1);
      Eigen::Matrix<typename DerivedV::Scalar, 3, 1> rotated=rotation*dir3D;
      return rotated;
    }
    
  private:
    
    // returns the 90 deg rotation of a (around n) most similar to target b
    /// a and b should be in the same plane orthogonal to N
    static Eigen::Matrix<typename DerivedV::Scalar, 3, 1> K_PI_new(const Eigen::Matrix<typename DerivedV::Scalar, 3, 1>& a,
                                                                   const Eigen::Matrix<typename DerivedV::Scalar, 3, 1>& b,
                                                                   const Eigen::Matrix<typename DerivedV::Scalar, 3, 1>& n)
    {
      Eigen::Matrix<typename DerivedV::Scalar, 3, 1> c = (a.cross(n)).normalized();
      typename DerivedV::Scalar scorea = a.dot(b);
      typename DerivedV::Scalar scorec = c.dot(b);
      if (fabs(scorea)>=fabs(scorec))
        return a*Sign(scorea);
      else
        return c*Sign(scorec);
    }
    
 
    
  public:
    Comb(const Eigen::PlainObjectBase<DerivedV> &_V,
         const Eigen::PlainObjectBase<DerivedF> &_F,
         const Eigen::PlainObjectBase<DerivedV> &_PD1,
         const Eigen::PlainObjectBase<DerivedV> &_PD2
         ):
    V(_V),
    F(_F),
    PD1(_PD1),
    PD2(_PD2)
    {
      igl::per_face_normals(V,F,N);
      igl::tt(V,F,TT,TTi);
    }
    void comb(Eigen::PlainObjectBase<DerivedV> &PD1out,
              Eigen::PlainObjectBase<DerivedV> &PD2out)
    {
//      PD1out = PD1;
//      PD2out = PD2;
      PD1out.setZero(F.rows(),3);PD1out<<PD1;
      PD2out.setZero(F.rows(),3);PD2out<<PD2;
      
      Eigen::VectorXi mark = Eigen::VectorXi::Constant(F.rows(),false);
      
      std::deque<int> d;
      
      d.push_back(0);
      mark(0) = true;
      
      while (!d.empty())
      {
        int f0 = d.at(0);
        d.pop_front();
        for (int k=0; k<3; k++)
        {
          int f1 = TT(f0,k);
          if (f1==-1) continue;
          if (mark(f1)) continue;
          
          Eigen::Matrix<typename DerivedV::Scalar, 3, 1> dir0    = PD1out.row(f0);
          Eigen::Matrix<typename DerivedV::Scalar, 3, 1> dir1    = PD1out.row(f1);
          Eigen::Matrix<typename DerivedV::Scalar, 3, 1> n0    = N.row(f0);
          Eigen::Matrix<typename DerivedV::Scalar, 3, 1> n1    = N.row(f1);

          Eigen::Matrix<typename DerivedV::Scalar, 3, 1> dir0Rot = Rotate(n0,n1,dir0);
          dir0Rot.normalize();
          Eigen::Matrix<typename DerivedV::Scalar, 3, 1> targD   = K_PI_new(dir1,dir0Rot,n1);
          
          PD1out.row(f1)  = targD;
          PD2out.row(f1)  = n1.cross(targD).normalized();
          
          mark(f1) = true;
          d.push_back(f1);

        }
      }
      
      // everything should be marked
      for (int i=0; i<F.rows(); i++)
      {
        assert(mark(i));
      }
    }
    
    
    
  };
}
template <typename DerivedV, typename DerivedF>
IGL_INLINE void igl::comb_cross_field(const Eigen::PlainObjectBase<DerivedV> &V,
                                      const Eigen::PlainObjectBase<DerivedF> &F,
                                      const Eigen::PlainObjectBase<DerivedV> &PD1,
                                      const Eigen::PlainObjectBase<DerivedV> &PD2,
                                      Eigen::PlainObjectBase<DerivedV> &PD1out,
                                      Eigen::PlainObjectBase<DerivedV> &PD2out)
{
  igl::Comb<DerivedV, DerivedF> cmb(V, F, PD1, PD2);
  cmb.comb(PD1out, PD2out);
}