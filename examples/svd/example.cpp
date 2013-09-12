int main(){}
//#include <igl/svd.h>
//#include <cstdlib>
//#include <Accelerate/Accelerate.h>
//#include <cstdio>
//
///* Auxiliary routines prototypes */
//extern void print_matrix( char* desc, int m, int n, double* a, int lda );
//
///* Parameters */
//
//void print3x3(const char * s, double * a)
//{
//  printf("%s =\n",s);
//  for(int i = 0;i<3;i++)
//  {
//    for(int j = 0;j<3;j++)
//    {
//      printf("%g ",a[j*3+i]);
//    }
//    printf("\n");
//  }
//  printf("\n");
//}
//
//int main(int argc, char * argv[])
//{
//  //// List of rest positions
//  ////        (0,1)
//  ////         / \
//  ////        /   \
//  ////       /     \
//  ////      /       \
//  ////  (-1,0)-----(1,0)
//  ////
//  //double rest[3][3] = {
//  //  {-1,0,0},
//  //  {1,0,0},
//  //  {0,1,0}};
//  //// List of pose positions
//  //// 
//  //// (0,1)
//  ////  |   \
//  ////  |    \
//  ////  |     (1,0)
//  ////  |    /
//  ////  |   /
//  //// (0,-1)
//  //double pose[3][3] = {
//  //  {0,1,0},
//  //  {0,-1,0},
//  //  {1,0,0}};
//  //// Compute covariance matrix C
//  //double C[3*3];
//  //// Initialize to zero
//  //for(int i = 0;i<3*3;i++)
//  //{
//  //  C[i] = 0;
//  //}
//
//  //// Loop over vertices
//  //for(int i = 0;i<3;i++)
//  //{
//  //  // Compute outer product rest[i] * pose[i]
//  //  // Loop over coordinates
//  //  for(int j = 0;j<3;j++)
//  //  {
//  //    // Loop over coordinates
//  //    for(int k = 0;k<3;k++)
//  //    {
//  //      C[k*3+j] = rest[i][j] * pose[i][k];
//  //    }
//  //  }
//  //}
//  //print3x3("C",C);
//
//
//  //
//  //double C[3*3] = {8,3,4,1,5,9,6,7,2};
//  double C[3*3] = {5242.55,3364,-0,-8170.15,-5242.56,0,-0,-0,0};
//  double u[3*3],s[3],vt[3*3];
//  print3x3("C",C);
//  // Compute SVD of C
//  igl::svd3x3(C,u,s,vt);
//  print3x3("u",u);
//  print3x3("vt",vt);
//
//  // Compute R = u*vt
//  double R[3*3];
//  const double _3 = 3;
//  const double _1 = 1;
//  cblas_dgemm(CblasColMajor, CblasNoTrans,CblasNoTrans,3,3,3,1,u,3,vt,3,1,R,3);
//  print3x3("RT (transposed to be row-major)",R);
//
//
//  return 0;
//}
