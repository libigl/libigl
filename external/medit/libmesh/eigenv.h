#ifdef __cplusplus
extern "C" {
#endif

int eigenv(int symmat,double *mat,double lambda[3],double v[3][3]);
int eigen2(double *mm,double *lambda,double vp[2][2]);

#ifdef __cplusplus
}
#endif
