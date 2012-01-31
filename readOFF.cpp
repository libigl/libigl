#include "readOFF.h"

IGL_INLINE bool igl::readOFF (const std::string meshfile, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    int vnum, fnum;
    FILE *fp = fopen (meshfile.c_str(), "r");
    
    if (!fp)
    {
      fprintf (stderr, "readOFF(): could not open file %s", meshfile.c_str());
      return false;
    }
    
    fscanf (fp, "OFF\n%d %d 0\n",  &vnum, &fnum);
    
    V = Eigen::MatrixXd (vnum, 3);
    F = Eigen::MatrixXi (fnum, 3);
    
    for (unsigned i = 0; i < V.rows(); i++)
        fscanf (fp, "%lf %lf %lf\n", &V(i,0), &V(i,1), &V(i,2));
    
    for (unsigned i = 0; i < F.rows(); i++)
        fscanf (fp, "3 %d %d %d\n", &F(i,0), &F(i,1), &F(i,2));
    
    fclose (fp);
    return true;
}
