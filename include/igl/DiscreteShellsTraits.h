// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_DISCRETE_SHELLS_TRAITS_H
#define HEDRA_DISCRETE_SHELLS_TRAITS_H
#include <igl/igl_inline.h>
#include <igl/harmonic.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>
#include <set>


namespace hedra { namespace optimization {
        
        //this class is a traits class for optimization of discrete shells deformation by given positional constraints. It is an implementation on [Froehlich and Botsch 2012] for general polyhedral meshes, using a triangulation of them
    
        //the solution vector is assumed to be arranged as xyzxyzxyz... where each triplet is a coordinate of the free vertices.
    
        class DiscreteShellsTraits{
        public:
            
            //Requisites of the traits class
            Eigen::VectorXi JRows, JCols;  //rows and column indices for the jacobian matrix
            Eigen::VectorXd JVals;         //values for the jacobian matrix.
            int xSize;                  //size of the solution
            Eigen::VectorXd EVec;          //energy vector
            
            //These are for the the full Jacobian matrix without removing the handles
            Eigen::VectorXi fullJRows, fullJCols;
            Eigen::VectorXd fullJVals;
            Eigen::MatrixXi flapVertexIndices;  //vertices (i,j,k,l) of a flap on edge e=(k,i) where the triangles are f=(i,j,k) and g=(i,k,l)
            Eigen::MatrixXi EV;
            Eigen::VectorXi h;              //list of handles
            Eigen::MatrixXd qh;             //#h by 3 positions
            Eigen::VectorXi a2x;            //from the entire set of variables "a" to the free variables in the optimization "x".
            Eigen::VectorXi colMap;         //raw map of a2x into the columns from fullJCols into JCols
            Eigen::VectorXd origLengths;    //the original edge lengths.
            Eigen::VectorXd origDihedrals;  //Original dihedral angles
            Eigen::MatrixXd VOrig;          //original positions
            Eigen::MatrixXi T;              //triangles
            Eigen::VectorXd Wd, Wl;         //geometric weights for the lengths and the dihedral angles.
            double bendCoeff, lengthCoeff;  //coefficients for the different terms
            
            //Eigen::VectorXd x0;                 //the initial solution to the optimization
            Eigen::MatrixXd fullSolution;       //The final solution of the last optimization
            
            void init(const Eigen::MatrixXd& _VOrig,
                      const Eigen::MatrixXi& _T,
                      const Eigen::VectorXi& _h,
                      const Eigen::MatrixXi& _EV,
                      const Eigen::MatrixXi& ET,
                      const Eigen::MatrixXi& ETi,
                      const Eigen::VectorXi& innerEdges){
                
                using namespace std;
                using namespace Eigen;
                
                std::set<pair<int, int> > edgeIndicesList;
                
                VOrig=_VOrig;
                T=_T;
                h=_h;
                EV=_EV;
                lengthCoeff=10.0;
                bendCoeff=1.0;
                
                //lengths of edges and diagonals
                flapVertexIndices.resize(innerEdges.size(),4);
                
                //dihedral angles
                for (int i=0;i<innerEdges.size();i++){
                    int f=ET(innerEdges(i),0);
                    int g=ET(innerEdges(i),1);
                 
                    int vi=EV(innerEdges(i), 1);
                    int vj=T(f,(ETi(innerEdges(i),0)+2)%3);
                    int vk=EV(innerEdges(i),0);
                    int vl=T(g,(ETi(innerEdges(i),1)+2)%3);
                    flapVertexIndices.row(i)<<vi,vj,vk,vl;
                 }
                
                //across edges
                EVec.resize(EV.rows()+flapVertexIndices.rows());
                origLengths.resize(EV.rows());
                origDihedrals.resize(flapVertexIndices.rows());
                Wl.resize(EV.rows());
                Wd.resize(flapVertexIndices.rows());
                
                
                for (int i=0;i<EV.rows();i++){
                    origLengths(i)=(VOrig.row(EV(i,1))-VOrig.row(EV(i,0))).norm();
                    Wl(i)=1.0/origLengths(i);
                }
                
                for (int i=0;i<flapVertexIndices.rows();i++){
                    RowVector3d eji=VOrig.row(flapVertexIndices(i,0))-VOrig.row(flapVertexIndices(i,1));
                    RowVector3d ejk=VOrig.row(flapVertexIndices(i,2))-VOrig.row(flapVertexIndices(i,1));
                    RowVector3d eli=VOrig.row(flapVertexIndices(i,0))-VOrig.row(flapVertexIndices(i,3));
                    RowVector3d elk=VOrig.row(flapVertexIndices(i,2))-VOrig.row(flapVertexIndices(i,3));
                    RowVector3d eki=VOrig.row(flapVertexIndices(i,0))-VOrig.row(flapVertexIndices(i,2));
                    
                    RowVector3d n1 = (ejk.cross(eji));
                    RowVector3d n2 = (eli.cross(elk));
                    double sign=((n1.cross(n2)).dot(eki) >= 0 ? 1.0 : -1.0);
                    double sinHalf=sign*sqrt((1.0-n1.normalized().dot(n2.normalized()))/2.0);
                    origDihedrals(i)=2.0*asin(sinHalf);
                    double areaSum=(n1.norm()+n2.norm())/2.0;
                    Wd(i)=eki.norm()/sqrt(areaSum);
                }
                
                
                //creating the Jacobian pattern
                xSize=3*(VOrig.rows()-h.size());
                fullJRows.resize(6*EV.rows()+12*flapVertexIndices.rows());
                fullJCols.resize(6*EV.rows()+12*flapVertexIndices.rows());
                fullJVals.resize(6*EV.rows()+12*flapVertexIndices.rows());
                
                
                //Jacobian indices for edge lengths
                for (int i=0;i<EV.rows();i++){
                    fullJRows.segment(6*i,6).setConstant(i);
                    fullJCols.segment(6*i,3)<<3*EV(i,0),3*EV(i,0)+1,3*EV(i,0)+2;
                    fullJCols.segment(6*i+3,3)<<3*EV(i,1),3*EV(i,1)+1,3*EV(i,1)+2;
                }
                
                //Jacobian indices for dihedral angles
                for (int i=0;i<flapVertexIndices.rows();i++){
                    fullJRows.segment(6*EV.rows()+12*i,12).setConstant(EV.rows()+i);
                    for (int k=0;k<4;k++)
                        fullJCols.segment(6*EV.rows()+12*i+3*k,3)<<3*flapVertexIndices(i,k),3*flapVertexIndices(i,k)+1,3*flapVertexIndices(i,k)+2;
                }
                
                
                a2x=Eigen::VectorXi::Zero(VOrig.rows());
                int CurrIndex=0;
                for (int i=0;i<h.size();i++)
                    a2x(h(i))=-1;
                
                for (int i=0;i<VOrig.rows();i++)
                    if (a2x(i)!=-1)
                        a2x(i)=CurrIndex++;
                
                colMap.resize(3*VOrig.rows());
                for (int i=0;i<VOrig.rows();i++)
                    if (a2x(i)!=-1)
                        colMap.segment(3*i,3)<<3*a2x(i),3*a2x(i)+1,3*a2x(i)+2;
                    else
                        colMap.segment(3*i,3)<<-1,-1,-1;
                
                //setting up the Jacobian rows and columns
                int actualGradCounter=0;
                for (int i=0;i<fullJCols.size();i++){
                    if (colMap(fullJCols(i))!=-1)  //not a removed variable
                        actualGradCounter++;
                }
                
                JRows.resize(actualGradCounter);
                JCols.resize(actualGradCounter);
                JVals.resize(actualGradCounter);
                
                
                actualGradCounter=0;
                for (int i=0;i<fullJCols.size();i++){
                    if (colMap(fullJCols(i))!=-1){  //not a removed variable
                        JRows(actualGradCounter)=fullJRows(i);
                        JCols(actualGradCounter++)=colMap(fullJCols(i));
                    }
                }
            }
            
            //provide the initial solution to the solver
            void initial_solution(Eigen::VectorXd& x0){
                //using biharmonic deformation fields
                Eigen::MatrixXd origHandlePoses(qh.rows(),3);
                for (int i=0;i<qh.rows();i++)
                    origHandlePoses.row(i)=VOrig.row(h(i));
                Eigen::MatrixXd D;
                Eigen::MatrixXd D_bc = qh - origHandlePoses;
                igl::harmonic(VOrig,T,h,D_bc,2,D);
                Eigen::MatrixXd V0 = VOrig+D;
                for (int i=0;i<VOrig.rows();i++)
                    if (a2x(i)!=-1)
                        x0.segment(3*a2x(i),3)<<V0.row(i).transpose();
                
            }
            
            void pre_iteration(const Eigen::VectorXd& prevx){}
            bool post_iteration(const Eigen::VectorXd& x){return false;  /*never stop after an iteration*/}
            
            
            //updating the energy vector for a given current solution
            void update_energy(const Eigen::VectorXd& x){
                
                using namespace std;
                using namespace Eigen;
                
                MatrixXd fullx(xSize+h.size(),3);
                for (int i=0;i<a2x.size();i++)
                    if (a2x(i)!=-1)
                        fullx.row(i)<<x.segment(3*a2x(i),3).transpose();
                
                for (int i=0;i<h.size();i++)
                    fullx.row(h(i))=qh.row(i);
                
                fullJVals.setZero();
                for (int i=0;i<EV.rows();i++)
                    EVec(i)=((fullx.row(EV(i,1))-fullx.row(EV(i,0))).norm()-origLengths(i))*Wl(i)*lengthCoeff;
                
                
                for (int i=0;i<flapVertexIndices.rows();i++){
                    RowVector3d eji=fullx.row(flapVertexIndices(i,0))-fullx.row(flapVertexIndices(i,1));
                    RowVector3d ejk=fullx.row(flapVertexIndices(i,2))-fullx.row(flapVertexIndices(i,1));
                    RowVector3d eli=fullx.row(flapVertexIndices(i,0))-fullx.row(flapVertexIndices(i,3));
                    RowVector3d elk=fullx.row(flapVertexIndices(i,2))-fullx.row(flapVertexIndices(i,3));
                    RowVector3d eki=fullx.row(flapVertexIndices(i,0))-fullx.row(flapVertexIndices(i,2));
                    
                    RowVector3d n1 = (ejk.cross(eji));
                    RowVector3d n2 = (eli.cross(elk));
                    double sign=((n1.cross(n2)).dot(eki) >= 0 ? 1.0 : -1.0);
                    double dotn1n2=1.0-n1.normalized().dot(n2.normalized());
                    if (dotn1n2<0.0) dotn1n2=0.0;  //sanitizing
                    double sinHalf=sign*sqrt(dotn1n2/2.0);
                    if (sinHalf>1.0) sinHalf=1.0; if (sinHalf<-1.0) sinHalf=-1.0;
                    EVec(EV.rows()+i)=(2.0*asin(sinHalf)-origDihedrals(i))*Wd(i)*bendCoeff;
                }
                
                for (int i=0;i<EVec.size();i++)
                    if (isnan(EVec(i)))
                        cout<<"nan in EVec("<<i<<")"<<endl;
            }
            
            
            //update the jacobian values for a given current solution
            void update_jacobian(const Eigen::VectorXd& x){
                using namespace std;
                using namespace Eigen;
                
                MatrixXd fullx(xSize+h.size(),3);
                for (int i=0;i<a2x.size();i++)
                    if (a2x(i)!=-1)
                        fullx.row(i)<<x.segment(3*a2x(i),3).transpose();
                
                for (int i=0;i<h.size();i++)
                    fullx.row(h(i))=qh.row(i);
                
                fullJVals.setZero();
                for (int i=0;i<EV.rows();i++){
                   
                    RowVector3d normedEdgeVector=(fullx.row(EV(i,1))-fullx.row(EV(i,0))).normalized();
                    fullJVals.segment(6*i,3)<<-normedEdgeVector.transpose()*Wl(i)*lengthCoeff;
                    fullJVals.segment(6*i+3,3)<<normedEdgeVector.transpose()*Wl(i)*lengthCoeff;
                }
                
                for (int i=0;i<flapVertexIndices.rows();i++){
                    RowVector3d eji=fullx.row(flapVertexIndices(i,0))-fullx.row(flapVertexIndices(i,1));
                    RowVector3d ejk=fullx.row(flapVertexIndices(i,2))-fullx.row(flapVertexIndices(i,1));
                    RowVector3d eli=fullx.row(flapVertexIndices(i,0))-fullx.row(flapVertexIndices(i,3));
                    RowVector3d elk=fullx.row(flapVertexIndices(i,2))-fullx.row(flapVertexIndices(i,3));
                    RowVector3d eki=fullx.row(flapVertexIndices(i,0))-fullx.row(flapVertexIndices(i,2));
                    
                    RowVector3d n1 = (ejk.cross(eji));
                    RowVector3d n2 = (eli.cross(elk));
                    double sign=((n1.cross(n2)).dot(eki) >= 0 ? 1.0 : -1.0);
                    double dotn1n2=1.0-n1.normalized().dot(n2.normalized());
                    if (dotn1n2<0.0) dotn1n2=0.0;  //sanitizing
                    double sinHalf=sign*sqrt(dotn1n2/2.0);
                    if (sinHalf>1.0) sinHalf=1.0; if (sinHalf<-1.0) sinHalf=-1.0;
                    
                    fullJVals.segment(6*EV.rows()+12*i,3)<<(Wd(i)*((ejk.dot(-eki)/(n1.squaredNorm()*eki.norm()))*n1+(elk.dot(-eki)/(n2.squaredNorm()*eki.norm()))*n2)).transpose()*bendCoeff;
                    fullJVals.segment(6*EV.rows()+12*i+3,3)<<(Wd(i)*(-eki.norm()/n1.squaredNorm())*n1).transpose()*bendCoeff;
                    fullJVals.segment(6*EV.rows()+12*i+6,3)<<(Wd(i)*((eji.dot(eki)/(n1.squaredNorm()*eki.norm()))*n1+(eli.dot(eki)/(n2.squaredNorm()*eki.norm()))*n2)).transpose()*bendCoeff;
                    fullJVals.segment(6*EV.rows()+12*i+9,3)<<(Wd(i)*(-eki.norm()/n2.squaredNorm())*n2).transpose()*bendCoeff;
                    
                }
                

                int actualGradCounter=0;
                for (int i=0;i<fullJCols.size();i++)
                    if (colMap(fullJCols(i))!=-1)  //not a removed variable
                        JVals(actualGradCounter++)=fullJVals(i);
                
                for (int i=0;i<JVals.size();i++)
                    if (isnan(JVals(i)))
                        cout<<"nan in JVals("<<i<<")"<<endl;
            }

            
            
            bool post_optimization(const Eigen::VectorXd& x){
                fullSolution.conservativeResize(a2x.size(),3);
                for (int i=0;i<a2x.size();i++)
                    if (a2x(i)!=-1)
                        fullSolution.row(i)<<x.segment(3*a2x(i),3).transpose();
                
                for (int i=0;i<h.size();i++)
                    fullSolution.row(h(i))=qh.row(i);
                
                return true;  //stop optimization after this
            }
            
            DiscreteShellsTraits(){}
            ~DiscreteShellsTraits(){}
        };
        
        
    } }


#endif
