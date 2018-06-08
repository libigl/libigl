#include "fast_winding_number.h"
#include <vector>
#include <iostream>
#include <igl/parallel_for.h>

void fast_winding_number_precompute(const Eigen::MatrixXd & P,
                                        const Eigen::MatrixXd & N,
                                        const Eigen::VectorXd & A,
                                        const std::vector<std::vector<int> > & point_indices,
                                        const std::vector<Eigen::Matrix<int,8,1>, Eigen::aligned_allocator<Eigen::Matrix<int,8,1>>> & children,
                                        const int expansion_order,
                                        Eigen::MatrixXd & CM,
                                        Eigen::VectorXd & R,
                                        Eigen::MatrixXd & EC
                                    ){
    int m = children.size();
    int num_terms;
    if(expansion_order == 0){
        num_terms = 3;
    } else if(expansion_order ==1){
        num_terms = 3 + 9;
    } else if(expansion_order == 2){
        num_terms = 3 + 9 + 27;
    } else {
        assert(false);
    }
    
    R.resize(m);
    CM.resize(m,3);
    EC.resize(m,num_terms);
    EC = Eigen::MatrixXd::Zero(m,num_terms);
    std::function< void(const int) > helper;
    helper = [&helper,
              &P,&N,&A,&expansion_order,&point_indices,&children,&EC,&R,&CM]
    (const int index)-> void
    {
        double sum_area = 0;
        
        Eigen::RowVector3d masscenter = Eigen::RowVector3d::Zero();
        Eigen::RowVector3d zeroth_expansion = Eigen::RowVector3d::Zero();
        double areatotal = 0.0;
        for(int j = 0; j < point_indices.at(index).size(); j++){
            int curr_point_index = point_indices.at(index).at(j);
            
            areatotal += A(curr_point_index);
            masscenter += A(curr_point_index)*P.row(curr_point_index);
            zeroth_expansion += A(curr_point_index)*N.row(curr_point_index);
        }
        
        masscenter = masscenter/areatotal;
        CM.row(index) = masscenter;
        EC.block<1,3>(index,0) = zeroth_expansion;
        
        double max_norm = 0;
        double curr_norm;
        
        for(int i = 0; i < point_indices.at(index).size(); i++){
            //Get max distance from center of mass:
            int curr_point_index = point_indices.at(index).at(i);
            Eigen::RowVector3d point = P.row(curr_point_index)-masscenter;
            curr_norm = point.norm();
            if(curr_norm > max_norm){
                max_norm = curr_norm;
            }
            
            //Calculate higher order terms if necessary
            Eigen::Matrix3d TempCoeffs;
            if(EC.cols() >= (3+9)){
                TempCoeffs = A(curr_point_index) * point.transpose() * N.row(curr_point_index);
                EC.block<1,9>(index,3) += Eigen::Map<Eigen::RowVectorXd>(TempCoeffs.data(), TempCoeffs.size());
            }
            
            if(EC.cols() == (3+9+27)){
                for(int k = 0; k < 3; k++){
                    TempCoeffs = 0.5 * point(k) * (A(curr_point_index) * point.transpose() * N.row(curr_point_index)) ;
                    EC.block<1,9>(index,12+9*k) += Eigen::Map<Eigen::RowVectorXd>(TempCoeffs.data(), TempCoeffs.size());
                }
            }
        }
        
        R(index) = max_norm;
        if(children.at(index)(0) != -1)
        {
            for(int i = 0; i < 8; i++){
                int child = children.at(index)(i);
                helper(child);
            }
        }
    };
    helper(0);
}

double direct_eval(const Eigen::RowVector3d & loc, const Eigen::RowVector3d & anorm){
    double wn = (loc(0)*anorm(0) + loc(1)*anorm(1) + loc(2)*anorm(2))/(4.0*M_PI*std::pow(loc.norm(),3));
    if(std::isnan(wn)){
        return 0.5;
    }else{
        return wn;
    }
}

double expansion_eval(const Eigen::RowVector3d & loc, const Eigen::RowVectorXd & EC){
    double wn = direct_eval(loc,EC.head<3>());
    double r = loc.norm();
    if(EC.size()>3){
        Eigen::Matrix3d SecondDerivative = Eigen::Matrix3d::Identity()/(4.0*M_PI*std::pow(r,3));
        SecondDerivative += -3.0*loc.transpose()*loc/(4.0*M_PI*std::pow(r,5));
        Eigen::RowVectorXd derivative_vector = Eigen::Map<Eigen::RowVectorXd>(SecondDerivative.data(), SecondDerivative.size());
        wn += derivative_vector.cwiseProduct(EC.segment<9>(3)).sum();
    }
    if(EC.size()>3+9){
        Eigen::Matrix3d ThirdDerivative;
        for(int i = 0; i < 3; i++){
            ThirdDerivative = 15.0*loc(i)*loc.transpose()*loc/(4.0*M_PI*std::pow(r,7));
            Eigen::Matrix3d Diagonal;
            Diagonal << loc(i), 0, 0,
                        0, loc(i), 0,
                        0, 0, loc(i);
            Eigen::Matrix3d RowCol = Eigen::Matrix3d::Zero();
            RowCol.row(i) = loc;
            RowCol = RowCol + RowCol.transpose();
            ThirdDerivative += -3.0/(4.0*M_PI*std::pow(r,5)) * (RowCol + Diagonal);
            Eigen::RowVectorXd derivative_vector = Eigen::Map<Eigen::RowVectorXd>(ThirdDerivative.data(), ThirdDerivative.size());
            wn += derivative_vector.cwiseProduct(EC.segment<9>(12 + i*9)).sum();
        }
    }
    return wn;
}

void fast_winding_number(const Eigen::MatrixXd & P,
                         const Eigen::MatrixXd & N,
                         const Eigen::VectorXd & A,
                         const std::vector<std::vector<int> > & point_indices,
                         const std::vector<Eigen::Matrix<int,8,1>,               Eigen::aligned_allocator<Eigen::Matrix<int,8,1>>> & children,
                         const Eigen::MatrixXd & CM,
                         const Eigen::VectorXd & R,
                         const Eigen::MatrixXd & EC,
                         const Eigen::MatrixXd & Q,
                         const double & beta,
                         Eigen::VectorXd & WN
                         ){
    int m = Q.rows();
    WN.resize(m);
    
    std::function< double(const Eigen::RowVector3d, const std::vector<int>) > helper;
    helper = [&helper,
              &P,&N,&A,
              &point_indices,&children,
              &CM,&R,&EC,&beta]
    (const Eigen::RowVector3d query, const std::vector<int> near_indices)-> double
    {
        std::vector<int> new_near_indices;
        double wn = 0;
        for(int i = 0; i < near_indices.size(); i++){
            int index = near_indices.at(i);
            //Leaf Case, Brute force
            if(children.at(index)(0) == -1){
                for(int j = 0; j < point_indices.at(index).size(); j++){
                    int curr_row = point_indices.at(index).at(j);
                    wn += direct_eval(P.row(curr_row)-query,N.row(curr_row)*A(curr_row));
                }
            }
            //Non-Leaf Case
            else {
                for(int child = 0; child < 8; child++){
                    int child_index = children.at(index)(child);
                    if(point_indices.at(child_index).size() > 0){
                        if((CM.row(child_index)-query).norm() > beta*R(child_index)){
                            if(children.at(child_index)(0) == -1){
                                for(int j = 0; j < point_indices.at(child_index).size(); j++){
                                    int curr_row = point_indices.at(child_index).at(j);
                                    wn += direct_eval(P.row(curr_row)-query,N.row(curr_row)*A(curr_row));
                                }
                            }else{
                                wn += expansion_eval(CM.row(child_index)-query,EC.row(child_index));
                            }
                        }else {
                            new_near_indices.emplace_back(child_index);
                        }
                    }
                }
            }
        }
        if(new_near_indices.size() > 0){
            wn += helper(query,new_near_indices);
        }
        return wn;
    };
    
    
    if(beta >= 0){
        std::vector<int> near_indices_start = {0};
        igl::parallel_for(m,[&](int iter){
            WN(iter) = helper(Q.row(iter),near_indices_start);
        },1000);
    } else {
        igl::parallel_for(m,[&](int iter){
            double wn = 0;
            for(int j = 0; j <P.rows(); j++){
                wn += direct_eval(P.row(j)-Q.row(iter),N.row(j)*A(j));
            }
            WN(iter) = wn;
        },1000);
    }
}























