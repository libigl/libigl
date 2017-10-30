#include <igl/unproject_onto_mesh.h>
#include <igl/viewer/Viewer.h>
#include <hedra/polygonal_read_OFF.h>
#include <hedra/triangulate_mesh.h>
#include <hedra/polygonal_edge_topology.h>
#include <hedra/point_spheres.h>
#include <hedra/shapeup.h>
#include <hedra/planarity.h>
#include <hedra/concyclity.h>
#include <hedra/regularity.h>
#include <hedra/scalar2RGB.h>

enum ViewingMode{PLANARITY, CONCYCLITY, REGULARITY} ViewingMode;


std::vector<int> Handles;
std::vector<Eigen::RowVector3d> HandlePoses;
int CurrentHandle;
Eigen::MatrixXd VOrig, V, C;
Eigen::MatrixXi F, T;
Eigen::VectorXi D, TF;
Eigen::MatrixXi EV, EF, FE, EFi;
Eigen::MatrixXd FEs;
Eigen::VectorXi innerEdges;
Eigen::Vector3d spans;
bool Editing=false;
bool ChoosingHandleMode=false;
double CurrWinZ;

Eigen::VectorXd planarity, concyclity, regularity;

hedra::ShapeupData sudata;
Eigen::MatrixXi SFaceRegular;
Eigen::VectorXi SDFaceRegular;



bool UpdateCurrentView(igl::viewer::Viewer & viewer)
{
    using namespace Eigen;
    using namespace std;
    
    hedra::planarity(V,D,F,planarity);
    hedra::concyclity(V,D,F,concyclity);
    hedra::regularity(V,D,F,regularity);
    
    MatrixXd sphereV;
    MatrixXi sphereT;
    MatrixXd sphereTC;
    Eigen::MatrixXd bc(Handles.size(),V.cols());
    for (int i=0;i<Handles.size();i++)
        bc.row(i)=HandlePoses[i].transpose();
    
    switch(ViewingMode){
        case PLANARITY: hedra::scalar2RGB(planarity, 0.0,1.0, C); break;
        case CONCYCLITY:  hedra::scalar2RGB(concyclity, 0.0,5.0, C); break;
        case REGULARITY:  hedra::scalar2RGB(regularity, 0.0,1.0, C); break;
    }
 
    double sphereRadius=spans.sum()/200.0;
    MatrixXd sphereGreens(Handles.size(),3);
    sphereGreens.col(0).setZero();
    sphereGreens.col(1).setOnes();
    sphereGreens.col(2).setZero();

    hedra::point_spheres(bc, sphereRadius, sphereGreens, 10, false, sphereV, sphereT, sphereTC);
    
    Eigen::MatrixXd bigV(V.rows()+sphereV.rows(),3);
    Eigen::MatrixXi bigT(T.rows()+sphereT.rows(),3);
    Eigen::MatrixXd bigTC(C.rows()+sphereTC.rows(),3);
    for (int i=0;i<T.rows();i++)
        bigTC.row(i)=C.row(TF(i));
    bigTC.block(T.rows(),0,sphereTC.rows(),3)=sphereTC;
    if (sphereV.rows()!=0){
        bigV<<V, sphereV;
        bigT<<T, sphereT+Eigen::MatrixXi::Constant(sphereT.rows(), sphereT.cols(), V.rows());
    } else{
        bigV<<V;
        bigT<<T;
    }
    
    viewer.core.show_lines=false;
    Eigen::MatrixXd OrigEdgeColors(EV.rows(),3);
    OrigEdgeColors.col(0)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    OrigEdgeColors.col(1)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    OrigEdgeColors.col(2)=Eigen::VectorXd::Constant(EV.rows(),0.0);
    
    viewer.data.clear();
    viewer.data.set_mesh(bigV,bigT);
    viewer.data.set_colors(bigTC);
    viewer.data.compute_normals();
    viewer.data.set_edges(V,EV,OrigEdgeColors);
    return true;
}

bool mouse_move(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y)
{
    if (!Editing)
        return false;
    
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    Eigen::RowVector3f NewPos=igl::unproject<float>(Eigen::Vector3f(x,y,CurrWinZ),
                                                    viewer.core.view * viewer.core.model,
                                                    viewer.core.proj,
                                                    viewer.core.viewport);
    
    HandlePoses[HandlePoses.size()-1]=NewPos.cast<double>();
    Eigen::RowVector3d Diff=HandlePoses[HandlePoses.size()-1]-VOrig.row(Handles[HandlePoses.size()-1]);
    
    Eigen::MatrixXd bc(Handles.size(),V.cols());
    for (int i=0;i<Handles.size();i++)
        bc.row(i)=HandlePoses[i].transpose();
    
   
    UpdateCurrentView(viewer);
    return true;

}


bool mouse_up(igl::viewer::Viewer& viewer, int button, int modifier)
{
    if (((igl::viewer::Viewer::MouseButton)button==igl::viewer::Viewer::MouseButton::Left))
        return false;
    
    Editing=false;
 
    return true;
}

bool mouse_down(igl::viewer::Viewer& viewer, int button, int modifier)
{
    if (((igl::viewer::Viewer::MouseButton)button==igl::viewer::Viewer::MouseButton::Left))
        return false;
    int vid, fid;
    Eigen::Vector3f bc;
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    if (!ChoosingHandleMode){
        Editing=true;
        return false;
    }
    if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view * viewer.core.model,
                                viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
    {
        //add the closest vertex to the handles
        Eigen::MatrixXf::Index maxRow, maxCol;
        bc.maxCoeff(&maxRow);
        int CurrVertex=F(fid, maxRow);
        bool Found=false;
        for (int i=0;i<Handles.size();i++)
            if (Handles[i]==CurrVertex){
                CurrVertex=Handles[i];
                Found=true;
            }
        
        if (!Found){
            Handles.push_back(CurrVertex);
            HandlePoses.push_back(V.row(CurrVertex));
        }
    
        Eigen::Vector3f WinCoords=igl::project<float>(V.row(CurrVertex).cast<float>(),
                                               viewer.core.view * viewer.core.model,
                                               viewer.core.proj,
                                               viewer.core.viewport);

        CurrWinZ=WinCoords(2);
        std::cout<<"Choosing Vertex :"<<CurrVertex<<std::endl;

        Eigen::VectorXi b(Handles.size());
        for (int i=0;i<Handles.size();i++)
            b(i)=Handles[i];
        
        //hedra::shapeup_precompute(V,D, F,SDFaceRegular,SFaceRegular, b,1.0,100.0, sudata);

        UpdateCurrentView(viewer);
    }
    return true;
}

bool key_up(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
    switch(key)
    {
            
        case '1': ChoosingHandleMode=false;
            break;
        case '2': ViewingMode=PLANARITY; break;
        case '3': ViewingMode=CONCYCLITY; break;
        case '4': ViewingMode=REGULARITY; break;

    }
    UpdateCurrentView(viewer);
    return false;
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
    switch(key)
    {
        case '1': ChoosingHandleMode=true;
            break;
    }
    return false;
}


int main(int argc, char *argv[])
{
    
    // Load a mesh in OFF format
    using namespace std;
    using namespace Eigen;
    
    std::cout<<R"(
    1 choose handles (with right mouse button)
    2 Show planarity between [0,1]
    3 Show concyclity between [0,5]
    4 Show regularity between [0,1]
    )";
    
    hedra::polygonal_read_OFF(DATA_PATH "/intersection.off", V, D, F);
    hedra::polygonal_edge_topology(D, F, EV, FE, EF,EFi,FEs,innerEdges);
    hedra::triangulate_mesh(D, F, T, TF);

    
    
    spans=V.colwise().maxCoeff()-V.colwise().minCoeff();
    
    VOrig=V;
    ViewingMode=REGULARITY;
    igl::viewer::Viewer viewer;
    viewer.callback_mouse_down = &mouse_down;
    viewer.callback_mouse_move = &mouse_move;
    viewer.callback_mouse_up=&mouse_up;
    viewer.callback_key_down=&key_down;
    viewer.callback_key_up=&key_up;
    viewer.core.background_color<<0.75,0.75,0.75,1.0;
    UpdateCurrentView(viewer);
    viewer.launch();
    
    cout<<"press 1+right button to select new handles"<<endl;
    cout<<"press the right button and drag the edit the mesh"<<endl;
}
