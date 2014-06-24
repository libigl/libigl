
#define IGL_HEADER_ONLY
#include <igl/principal_curvature.h>
#include <igl/read.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>


using namespace std;
bool cont=true;
bool benchmark=false;

void print_help()
{
    cout << "Usage: Curvature -i meshfile [-o curvfile] [-S scalefile] "
            "[-k N | -e F] [-a] [-l] [-p] [-m N] [-s] [-z] [-b] [-n N] [-M N] [-h] \n"
            "    -k N: compute using the N-ring of vertexes.\n"
            "    -e S: compute using a sphere of radius R = S * mesh_average_edge.\n"
            "    -a: no projection plane, use normals average.\n"
            "    -l: local mode, use vcg normal (when using -a) or incident faces normal (default)\n"
            "    -p: check projection plane calculation\n"
            "    -m N: use montecarlo, extract N vertices maximum.\n"
            "    -b: benchmark mode\n"
            "    -n N: step for benchmarking\n"
            "    -M N: max size for benchmarking\n"
            "    -s: use svd calculation, not pseudoinverse.\n"
            "    -z: check if determinant is almost zero.\n"
            "    -h: print this help and exit.\n";
}

bool isDir(char* file)
{
    std::string s(file);
    return (s.find_last_of('/')==s.size()-1);
}

void app_init(int argc, char* argv[], CurvatureCalculator& c, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    if (argc<2)
    {
        print_help();
        exit(0);
    }
    char* tmp;
    const char* meshName;
    char * scaleFile;
    for (argc--, argv++; argc--; argv++)
    {
        if( (*argv)[0] == '-')
        {
            switch( (*argv)[1] )
            {
            case 'b':
                benchmark=true;
                break;
            case 'n':
                tmp=*++argv;
                if (!strcmp(tmp,"exp"))
                    c.expStep=true;
                else
                    c.step=atoi(tmp);
                --argc;
                break;
            case 'M':
                c.maxSize=atoi(*++argv);
                argc--;
                break;
            case 'i':
                meshName=*++argv;
                igl::read_triangle_mesh(meshName,V,F);
                argc--;
                break;
            case 'k':
                c.st=K_RING_SEARCH;
                c.kRing = atoi (*++argv);
                argc--;
                break;

            case 'e':
                c.st=SPHERE_SEARCH;
                c.sphereRadius = atof (*++argv);
                argc--;
                break;

            case 'a':
                c.nt=AVERAGE;
                break;

            case 'l':
                c.localMode = true;
                break;

            case 'p':
                c.projectionPlaneCheck = true;
                break;

            case 'm':
                c.montecarlo=true;
                c.montecarloN = atoi (*++argv);
                argc--;
                break;

            case 's':
                c.svd = true;
                break;

            case 'z':
                c.zeroDetCheck = true;
                break;

            case 'f':
                /* dc.fitCreaseMode = atoi (*++argv);
                    argc--;*/
                break;

            case 'o':
                tmp=*++argv;
                cerr << "Tmp: " << tmp << endl;
                if (isDir(tmp))
                {
                    std::string tmpS(tmp);
                    std::string tmpM(meshName);
                    size_t pos=tmpM.find_last_of('/');
                    tmpM.assign(tmpM.begin()+pos+1,tmpM.end());
                    c.lastMeshName=tmpS+tmpM;
                }
                else
                    c.lastMeshName=tmp;
                argc--;
                break;

            case 'h':
                print_help();
                exit(0);

            default:
                print_help();
                exit(0);
            }
        }
    }

    return;
}

int main(int argc, char *argv[])
{
    CurvatureCalculator c;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    string filename;

    app_init(argc,argv,c,V,F);

    c.init(V,F);

    c.computeCurvature();

    c.printCurvature(c.lastMeshName);

}
