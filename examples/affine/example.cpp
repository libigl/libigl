#include <Eigen/Geometry>
#include <iostream>
using namespace std;

typedef Eigen::Transform<double,3,Eigen::Affine> Tform3;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Quaterniond Quat;

int main(int argc, char * argv[])
{
  Tform3 T;
  T =Tform3::Identity();
  cout<<"T =Tform3::Identity();"<<endl;
  cout<<"T=["<<endl<<T.affine().matrix()<<endl<<"];"<<endl;
  cout<<"T.translate(Vec3(1,2,3));"<<endl;
  T.translate(Vec3(1,2,3));
  cout<<"T=["<<endl<<T.affine().matrix()<<endl<<"];"<<endl;
  cout<<"T.rotate(Quat(0,1,0,0));"<<endl;
  T.rotate(Quat(0,1,0,0));
  cout<<"T=["<<endl<<T.affine().matrix()<<endl<<"];"<<endl;
  cout<<endl;

  T =Tform3::Identity();
  cout<<"T =Tform3::Identity();"<<endl;
  cout<<"T=["<<endl<<T.affine().matrix()<<endl<<"];"<<endl;
  cout<<"T.rotate(Quat(0,1,0,0));"<<endl;
  T.rotate(Quat(0,1,0,0));
  cout<<"T=["<<endl<<T.affine().matrix()<<endl<<"];"<<endl;
  cout<<"T.translate(Vec3(1,2,3));"<<endl;
  T.translate(Vec3(1,2,3));
  cout<<"T=["<<endl<<T.affine().matrix()<<endl<<"];"<<endl;
  cout<<endl;
}
