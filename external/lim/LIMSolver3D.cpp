// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "LIMSolver3D.h"
#include "TetrahedronMesh.h"

LIMSolver3D::LIMSolver3D()
{
  dim = 3;
  mesh = NULL;
}

LIMSolver3D::~LIMSolver3D()
{
}

void LIMSolver3D::Init(DeformableMesh* mesh)
{
  Init(static_cast<TetrahedronMesh*>(mesh));
}

void LIMSolver3D::Init(TetrahedronMesh* mesh)
{
  this->mesh = mesh;
  LIMSolver::Init(mesh);
}

void LIMSolver3D::prepareNMProblemData()
{
  std::vector<Eigen::Triplet<double> > triplets;
  const int numNodes =  mesh->InitalVertices->rows();
  const int numTets = mesh->Tetrahedra->rows();

  TetrahedronVertexIdx.resize(12,numTets);
  for(int t=0;t<numTets;t++)
  {
    Vector4i nodes = mesh->Tetrahedra->row(t);

    // create tets vertex indicies
    for(int v=0;v<4;v++)
    {
      for(int i=0;i<3;i++)
      {
        TetrahedronVertexIdx(v*3+i,t) = nodes[v]*3+i;
      }
    }
    
    for(int r=0;r<12;r++)
    {
      for(int c=r;c<12;c++)
      {
        int row = TetrahedronVertexIdx(r,t);
        int col = TetrahedronVertexIdx(c,t);
      
        // std::sort for upper triangule matrix				
        if(row > col)
        {
          int temp = col;
          col = row;
          row = temp;
        }
      
        triplets.push_back(Eigen::Triplet<double>(row,col,1));
      }
    }
  }

  // get positional constraint matrix structure of problem
  for (int k=0;k<linearConstraintsMatrix2.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(linearConstraintsMatrix2,k);it;++it)
    {
      int row = it.row();
      int col = it.col();

      // std::sort for upper triangule matrix
      if(row < col)
        triplets.push_back(Eigen::Triplet<double>(row,col,1));
    }
  }

  // get hessian matrix structure of problem
  std::vector<int> hessRowIdx;
  std::vector<int> hessColIdx;
  prepareProblemData(hessRowIdx,hessColIdx);
  assert(hessRowIdx.size() == hessColIdx.size());
  
  for(int i=0;i<hessRowIdx.size();i++)
  {
    // std::sort for upper triangle matrix
    if(hessRowIdx[i] <= hessColIdx[i])
      triplets.push_back(Eigen::Triplet<double>(hessRowIdx[i],hessColIdx[i],1));
  }
  hessian.setFromTriplets(triplets.begin(),triplets.end());
  numNonZeroHessian = hessian.nonZeros();
  
  // Init vertex indices link vector
  VertexPositionIndices.resize(numVariables);
  for(int n=0;n<numNodes;n++)
  {
    for(int i=0;i<3;i++)
    {
      VertexPositionIndices[n*3+i] = n+numNodes*i;
    }
  }

  initialNodes.resize(numVariables);
  for(int n=0;n<numVariables;n++)
  {
    initialNodes[n] = mesh->InitalVertices->coeff(VertexPositionIndices[n]);
  }

    // set link indices of non zero elements of constraint hessian matrix
  problemHessianCoeffs.resize(1,hessRowIdx.size());
  denseHessianCoeffs.resize(78,numTets);
  diagHessianCoeffs.resize(numVariables);
  posConstraintsHessianCoeffs.resize(linearConstraintsMatrix2.nonZeros());

  for(int i=0;i<hessRowIdx.size();i++)
  {
    if(hessRowIdx[i] <= hessColIdx[i])
      problemHessianCoeffs[i] = &hessian.coeffRef(hessRowIdx[i],hessColIdx[i]);
  }

  for(int t=0;t<numTets;t++)
  {
    int numElem = 0;
    for(int r=0;r<12;r++)
    {
      for(int c=r;c<12;c++)
      {
        int row = TetrahedronVertexIdx(r,t);
        int col = TetrahedronVertexIdx(c,t);

        // std::sort for upper triangule matrix				
        if(row > col)
        {
          int temp = col;
          col = row;
          row = temp;
        }

        denseHessianCoeffs(numElem,t) = &hessian.coeffRef(row,col);

        if(col == row)
          diagHessianCoeffs[row] = denseHessianCoeffs(numElem,t);

        numElem++;
      }
    }
  }
  
  // non-flip constraints
  nonFlipHessianCoeffs.resize(36,numTets);
  for(int t=0;t<numTets;t++)
  {
    for(int i=0;i<36;i++)
    {
      int row = TetrahedronVertexIdx(NonFlipHessian3DIdx[i][0],t);
      int col = TetrahedronVertexIdx(NonFlipHessian3DIdx[i][1],t);
    
      // std::sort for upper triangule matrix				
      if(row > col)
      {
        int temp = col;
        col = row;
        row = temp;
      }
      
      nonFlipHessianCoeffs(i,t) = &hessian.coeffRef(row,col);
    }
  }

  // get positional constraint hessian entries
  int count = 0;
  for (int k=0;k<linearConstraintsMatrix2.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(linearConstraintsMatrix2,k);it;++it)
    {
      int row = it.row();
      int col = it.col();

      // std::sort for upper triangule matrix
      if(row <= col)
        posConstraintsHessianCoeffs[count++] = &hessian.coeffRef(row,col);
    }
  }
}

void LIMSolver3D::computeRestPoseFunctionParameters()
{ 
  const int numTets = mesh->Tetrahedra->rows();
  
  initalSize.resize(numTets);
  barrierParam1.resize(numTets);
  barrierParam2.resize(numTets);
  
  for(int t=0;t<numTets;t++)
  {
    Eigen::VectorXi indices = TetrahedronVertexIdx.col(t);
    
    Eigen::Vector3d A(initialNodes[indices[0]],initialNodes[indices[1]],initialNodes[indices[2]]);
    Eigen::Vector3d B(initialNodes[indices[3]],initialNodes[indices[4]],initialNodes[indices[5]]);
    Eigen::Vector3d C(initialNodes[indices[6]],initialNodes[indices[7]],initialNodes[indices[8]]);
    Eigen::Vector3d D(initialNodes[indices[9]],initialNodes[indices[10]],initialNodes[indices[11]]);

    Eigen::Vector3d a = A-D;
    Eigen::Vector3d b = B-D;
    Eigen::Vector3d c = C-D;
    
    double det = a.dot(c.cross(b));
    double detEPS = det-mesh->EPS3;

    initalSize[t] = detEPS;

    double quot = detEPS*CompensationExp*pow(det,CompensationExp-1);

    // must be multiplied by beta later
    barrierParam1[t] = 1/quot;
    barrierParam2[t] = log(detEPS) - pow(det,CompensationExp)/quot;
  }
}

double LIMSolver3D::computeNMFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{  
  const int numTets = mesh->Tetrahedra->rows();
  
  double dObj;
  
  // position constraints
  double diff = (*dmesh->ConstraintMatrix * x - subStepConstraints).squaredNorm();

  // non-flip constraint energy
  double ic = 0.0;
  if(EnableBarriers)
  {
    if(EnableLogBarriers)
    {
      for(int t=0;t<numTets;t++)
      {
        Eigen::Vector3d a(x[TetrahedronVertexIdx.coeff(0,t)]-x[TetrahedronVertexIdx.coeff(9,t)],x[TetrahedronVertexIdx.coeff(1,t)]-x[TetrahedronVertexIdx.coeff(10,t)],x[TetrahedronVertexIdx.coeff(2,t)]-x[TetrahedronVertexIdx.coeff(11,t)]);
        Eigen::Vector3d b(x[TetrahedronVertexIdx.coeff(3,t)]-x[TetrahedronVertexIdx.coeff(9,t)],x[TetrahedronVertexIdx.coeff(4,t)]-x[TetrahedronVertexIdx.coeff(10,t)],x[TetrahedronVertexIdx.coeff(5,t)]-x[TetrahedronVertexIdx.coeff(11,t)]);
        Eigen::Vector3d c(x[TetrahedronVertexIdx.coeff(6,t)]-x[TetrahedronVertexIdx.coeff(9,t)],x[TetrahedronVertexIdx.coeff(7,t)]-x[TetrahedronVertexIdx.coeff(10,t)],x[TetrahedronVertexIdx.coeff(8,t)]-x[TetrahedronVertexIdx.coeff(11,t)]);

        double det = a.dot(c.cross(b));
        double detEPS = det-mesh->EPS3;
        if(detEPS > 0.0)
        {
          ic += -log(detEPS); // barrier function
          if(EnableBarrierCompensation) ic += barrierParam1[t]*pow((double)det,CompensationExp)+barrierParam2[t]; // barrier compensation function
        }
        else
        {
          ic = std::numeric_limits<double>::infinity();
          break;
        }
      }
    }
    else
    {
      for(int t=0;t<numTets;t++)
      {
        Eigen::Vector3d a(x[TetrahedronVertexIdx.coeff(0,t)]-x[TetrahedronVertexIdx.coeff(9,t)],x[TetrahedronVertexIdx.coeff(1,t)]-x[TetrahedronVertexIdx.coeff(10,t)],x[TetrahedronVertexIdx.coeff(2,t)]-x[TetrahedronVertexIdx.coeff(11,t)]);
        Eigen::Vector3d b(x[TetrahedronVertexIdx.coeff(3,t)]-x[TetrahedronVertexIdx.coeff(9,t)],x[TetrahedronVertexIdx.coeff(4,t)]-x[TetrahedronVertexIdx.coeff(10,t)],x[TetrahedronVertexIdx.coeff(5,t)]-x[TetrahedronVertexIdx.coeff(11,t)]);
        Eigen::Vector3d c(x[TetrahedronVertexIdx.coeff(6,t)]-x[TetrahedronVertexIdx.coeff(9,t)],x[TetrahedronVertexIdx.coeff(7,t)]-x[TetrahedronVertexIdx.coeff(10,t)],x[TetrahedronVertexIdx.coeff(8,t)]-x[TetrahedronVertexIdx.coeff(11,t)]);

        double det = a.dot(c.cross(b));
        double detEPS = det-mesh->EPS3;

        double minDet = initalSize[t]*Gamma;
        if(detEPS < minDet)
        {
          double coeffA =  1/pow(minDet,3);
          double coeffB = -3/pow(minDet,2);
          double coeffC =  3/minDet;

          double detEPS2 = detEPS*detEPS;
          double detEPS3 = detEPS2*detEPS;

          if(detEPS > 0.0)
          {
            ic += 1.0/(coeffA*detEPS3 + coeffB*detEPS2 + coeffC*detEPS) - 1.0; // barrier function
          }
          else
          {
            ic = std::numeric_limits<double>::infinity();
            break;
          }
        }
      }	
    }
  }

  if(ic == std::numeric_limits<double>::infinity())
    dObj = std::numeric_limits<double>::infinity();
  else
  {
    CurrentPositionalSubStepEnergy = Alpha*diff;
    CurrentConstraintEnergy = Beta*ic;
    CurrentDeformationEnergy = computeFunction(x);
    dObj = CurrentPositionalSubStepEnergy + CurrentConstraintEnergy + CurrentDeformationEnergy;
  }

    return dObj ;
}

void LIMSolver3D::computeNMGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad)
{
  const int numTets = mesh->Tetrahedra->rows();
  
  grad.setZero();
  
  // call subclass function
  computeGradient(x,grad);

  // position constraints
  Eigen::Matrix<double,Eigen::Dynamic,1> temp1 = linearConstraintsMatrix2 * x;
  Eigen::Matrix<double,Eigen::Dynamic,1> temp2 = subStepConstraints.transpose() * *dmesh->ConstraintMatrix;
  grad += 2.0*Alpha*(temp1 - temp2);

  // non-flip constraints
    if(EnableBarriers)
  {
    if(EnableLogBarriers)
    {
      for(int t=0;t<numTets;t++)
      {
        Eigen::VectorXi indices = TetrahedronVertexIdx.col(t);
    
        Eigen::Vector3d A(x[indices[0]],x[indices[1]],x[indices[2]]);
        Eigen::Vector3d B(x[indices[3]],x[indices[4]],x[indices[5]]);
        Eigen::Vector3d C(x[indices[6]],x[indices[7]],x[indices[8]]);
        Eigen::Vector3d D(x[indices[9]],x[indices[10]],x[indices[11]]);

        Eigen::Vector3d a = A-D;
        Eigen::Vector3d b = B-D;
        Eigen::Vector3d c = C-D;
    
        double det = a.dot(c.cross(b));
        double detEPS = det-mesh->EPS3;

        double term = -1/detEPS; // barrier function
        if(detEPS > 0)
        {
          term += barrierParam1[t]*CompensationExp*pow((double)det,CompensationExp-1); // barrier compensation function - barrier function
          term *= Beta;
        }
        else
          term = std::numeric_limits<double>::infinity();


        // partial derivatives in respect to a
        grad[indices[0]]  += term*(c[1]* b[2] - c[2]* b[1]);
        grad[indices[1]]  += term*(c[2]* b[0] - c[0]* b[2]);
        grad[indices[2]]  += term*(c[0]* b[1] - c[1]* b[0]);
    
        // partial derivatives in respect to b
        grad[indices[3]]  += term*(a[1]* c[2] + a[2]*-c[1]);
        grad[indices[4]]  += term*(a[0]*-c[2] + a[2]* c[0]);
        grad[indices[5]]  += term*(a[0]* c[1] + a[1]*-c[0]);
    
        // partial derivatives in respect to c
        grad[indices[6]]  += term*(a[1]*-b[2] + a[2]* b[1]);
        grad[indices[7]]  += term*(a[0]* b[2] + a[2]*-b[0]);
        grad[indices[8]]  += term*(a[0]*-b[1] + a[1]* b[0]);
    
        // partial derivatives in respect to d
        grad[indices[9]]  += term*(-c[1]*b[2] + c[2]* b[1] + a[1]*(-C[2]+B[2]) + a[2]*(-B[1]+C[1]));
        grad[indices[10]] += term*(a[0]*(-B[2]+C[2]) - c[2]*b[0] + c[0]*b[2] + a[2]*(-C[0]+B[0]));
        grad[indices[11]] += term*(a[0]*(-C[1]+B[1]) + a[1]*(-B[0]+C[0]) - c[0]*b[1] + c[1]*b[0]);
      }
    }
    else
    {
      for(int t=0;t<numTets;t++)
      {
        Eigen::VectorXi indices = TetrahedronVertexIdx.col(t);
    
        Eigen::Vector3d A(x[indices[0]],x[indices[1]],x[indices[2]]);
        Eigen::Vector3d B(x[indices[3]],x[indices[4]],x[indices[5]]);
        Eigen::Vector3d C(x[indices[6]],x[indices[7]],x[indices[8]]);
        Eigen::Vector3d D(x[indices[9]],x[indices[10]],x[indices[11]]);

        Eigen::Vector3d a = A-D;
        Eigen::Vector3d b = B-D;
        Eigen::Vector3d c = C-D;
    
        double det = a.dot(c.cross(b));
        double detEPS = det-mesh->EPS3;

        double minDet = initalSize[t]*Gamma;
        if(detEPS < minDet)
        {
          double coeffA =  1/pow(minDet,3);
          double coeffB = -3/pow(minDet,2);
          double coeffC =  3/minDet;

          double detEPS2 = detEPS*detEPS;
          double detEPS3 = detEPS2*detEPS;

          double term = -Beta*(3*coeffA*detEPS2 + 2*coeffB*detEPS + coeffC) / pow(coeffA*detEPS3 + coeffB*detEPS2 + coeffC*detEPS,2);

          // partial derivatives in respect to a
          grad[indices[0]]  += term*(c[1]* b[2] - c[2]* b[1]);
          grad[indices[1]]  += term*(c[2]* b[0] - c[0]* b[2]);
          grad[indices[2]]  += term*(c[0]* b[1] - c[1]* b[0]);
    
          // partial derivatives in respect to b
          grad[indices[3]]  += term*(a[1]* c[2] + a[2]*-c[1]);
          grad[indices[4]]  += term*(a[0]*-c[2] + a[2]* c[0]);
          grad[indices[5]]  += term*(a[0]* c[1] + a[1]*-c[0]);
    
          // partial derivatives in respect to c
          grad[indices[6]]  += term*(a[1]*-b[2] + a[2]* b[1]);
          grad[indices[7]]  += term*(a[0]* b[2] + a[2]*-b[0]);
          grad[indices[8]]  += term*(a[0]*-b[1] + a[1]* b[0]);
    
          // partial derivatives in respect to d
          grad[indices[9]]  += term*(-c[1]*b[2] + c[2]* b[1] + a[1]*(-C[2]+B[2]) + a[2]*(-B[1]+C[1]));
          grad[indices[10]] += term*(a[0]*(-B[2]+C[2]) - c[2]*b[0] + c[0]*b[2] + a[2]*(-C[0]+B[0]));
          grad[indices[11]] += term*(a[0]*(-C[1]+B[1]) + a[1]*(-B[0]+C[0]) - c[0]*b[1] + c[1]*b[0]);
        }
      }
    }
  }
}

void LIMSolver3D::computeNMHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  const int numTets = mesh->Tetrahedra->rows();

  hessian *= 0;

  // compute problem function hessian
  computeHessian(x,problemHessianCoeffs);

  // position constraints: 2*alpha
  int count = 0;
  for (int k=0;k<linearConstraintsMatrix2.outerSize();++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(linearConstraintsMatrix2,k);it;++it)
    {
      if(it.row() <= it.col())
        *posConstraintsHessianCoeffs[count++] += 2*Alpha*it.value();
    }
  }
  
  // non-flip constraints
  if(EnableBarriers)
  {
    if(EnableLogBarriers)
    {
      for(int t=0;t<numTets;t++)
      {
        Eigen::VectorXi indices = TetrahedronVertexIdx.col(t);
    
        Eigen::Vector3d A(x[indices[0]],x[indices[1]],x[indices[2]]);
        Eigen::Vector3d B(x[indices[3]],x[indices[4]],x[indices[5]]);
        Eigen::Vector3d C(x[indices[6]],x[indices[7]],x[indices[8]]);
        Eigen::Vector3d D(x[indices[9]],x[indices[10]],x[indices[11]]);

        Eigen::Vector3d a = A-D;
        Eigen::Vector3d b = B-D;
        Eigen::Vector3d c = C-D;

        double det = a.dot(c.cross(b));
        double detEPS = det-mesh->EPS3;
        double detEPSInv = 1/detEPS;
        double term = -detEPSInv;
        if(detEPS > 0)
        {
          if(EnableBarrierCompensation) term += barrierParam1[t]*CompensationExp*pow((double)det,CompensationExp-1); // barrier compensation function - barrier function
          term *= Beta;
        }
        else
          term = std::numeric_limits<double>::infinity();

        Eigen::Matrix<double*,Eigen::Dynamic,1> elems = nonFlipHessianCoeffs.col(t);
    
        *elems[0]   += term * ( c[2]);
        *elems[1]   += term * (-c[1]);
        *elems[2]   += term * (-c[2]);
        *elems[3]   += term * ( c[0]);
        *elems[4]   += term * ( c[1]);
        *elems[5]   += term * (-c[0]);
        *elems[6]   += term * (-b[2]);
        *elems[7]   += term * ( b[1]); 
        *elems[8]   += term * ( a[2]);
        *elems[9]   += term * (-a[1]);
        *elems[10]  += term * ( b[2]);
        *elems[11]  += term * (-b[0]);
        *elems[12]  += term * (-a[2]);
        *elems[13]  += term * ( a[0]); 
        *elems[14]  += term * (-b[1]);
        *elems[15]  += term * ( b[0]);
        *elems[16]  += term * ( a[1]);
        *elems[17]  += term * (-a[0]);
        *elems[18]  += term * (-C[2]+B[2]);
        *elems[19]  += term * (-B[1]+C[1]);
        *elems[20]  += term * ( C[2]-A[2]);
        *elems[21]  += term * ( A[1]-C[1]);
        *elems[22]  += term * (-B[2]+A[2]);
        *elems[23]  += term * (-A[1]+B[1]);
        *elems[24]  += term * (-B[2]+C[2]);
        *elems[25]  += term * (-C[0]+B[0]);
        *elems[26]  += term * (-C[2]+A[2]);
        *elems[27]  += term * (-A[0]+C[0]);
        *elems[28]  += term * ( B[2]-A[2]);
        *elems[29]  += term * ( A[0]-B[0]);
        *elems[30]  += term * (-C[1]+B[1]);
        *elems[31]  += term * (-B[0]+C[0]);
        *elems[32]  += term * ( C[1]-A[1]);
        *elems[33]  += term * ( A[0]-C[0]);
        *elems[34]  += term * (-B[1]+A[1]);
        *elems[35]  += term * (-A[0]+B[0]);
 
        Eigen::Matrix<double,Eigen::Dynamic,1> g(12);

        // partial derivatives in respect to a
        g[0] = (c[1]* b[2] - c[2]* b[1]);
        g[1] = (c[2]* b[0] - c[0]* b[2]);
        g[2] = (c[0]* b[1] - c[1]* b[0]);
    
        // partial derivatives in respect to b
        g[3] = (a[1]* c[2] + a[2]*-c[1]);
        g[4] = (a[0]*-c[2] + a[2]* c[0]);
        g[5] = (a[0]* c[1] + a[1]*-c[0]);
    
        // partial derivatives in respect to c
        g[6] = (a[1]*-b[2] + a[2]* b[1]);
        g[7] = (a[0]* b[2] + a[2]*-b[0]);
        g[8] = (a[0]*-b[1] + a[1]* b[0]);
    
        // partial derivatives in respect to d
        g[9] =  (-c[1]*b[2] + c[2]* b[1] + a[1]*(-C[2]+B[2]) + a[2]*(-B[1]+C[1]));
        g[10] = (a[0]*(-B[2]+C[2]) - c[2]*b[0] + c[0]*b[2] + a[2]*(-C[0]+B[0]));
        g[11] = (a[0]*(-C[1]+B[1]) + a[1]*(-B[0]+C[0]) - c[0]*b[1] + c[1]*b[0]);

        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> g2(12,12);
        g2 = g*g.transpose();
        g2 *= Beta*(barrierParam1[t]*CompensationExp*(CompensationExp-1)*pow((double)det,CompensationExp-2) + detEPSInv*detEPSInv); // barrier compensation function - barrier function

        int numElem = 0;
        for(int r=0;r<12;r++)
        {
          for(int c=r;c<12;c++)
          {
            *denseHessianCoeffs(numElem++,t) += g2(r,c);
          }
        }
      }
    }
    else
    {
      for(int t=0;t<numTets;t++)
      {			
        Eigen::VectorXi indices = TetrahedronVertexIdx.col(t);
    
        Eigen::Vector3d A(x[indices[0]],x[indices[1]],x[indices[2]]);
        Eigen::Vector3d B(x[indices[3]],x[indices[4]],x[indices[5]]);
        Eigen::Vector3d C(x[indices[6]],x[indices[7]],x[indices[8]]);
        Eigen::Vector3d D(x[indices[9]],x[indices[10]],x[indices[11]]);

        Eigen::Vector3d a = A-D;
        Eigen::Vector3d b = B-D;
        Eigen::Vector3d c = C-D;

        double det = a.dot(c.cross(b));
        double detEPS = det-mesh->EPS3;

        double minDet = initalSize[t]*Gamma;
        if(detEPS < minDet)
        {
          double coeffA =  1/pow(minDet,3);
          double coeffB = -3/pow(minDet,2);
          double coeffC =  3/minDet;

          double detEPS2 = detEPS*detEPS;
          double detEPS3 = detEPS2*detEPS;

          double divTerm  = (coeffA*detEPS3 + coeffB*detEPS2 + coeffC*detEPS);
          double divTerm2 = divTerm*divTerm;
          double divTerm3 = divTerm2*divTerm;

          double factTerm = (3*coeffA*detEPS2 + 2*coeffB*detEPS + coeffC);
        
          // barrier function
          double term = -Beta*factTerm/divTerm2;

          Eigen::Matrix<double*,Eigen::Dynamic,1> elems = nonFlipHessianCoeffs.col(t);
    
          *elems[0]   += term * ( c[2]);
          *elems[1]   += term * (-c[1]);
          *elems[2]   += term * (-c[2]);
          *elems[3]   += term * ( c[0]);
          *elems[4]   += term * ( c[1]);
          *elems[5]   += term * (-c[0]);
          *elems[6]   += term * (-b[2]);
          *elems[7]   += term * ( b[1]); 
          *elems[8]   += term * ( a[2]);
          *elems[9]   += term * (-a[1]);
          *elems[10]  += term * ( b[2]);
          *elems[11]  += term * (-b[0]);
          *elems[12]  += term * (-a[2]);
          *elems[13]  += term * ( a[0]); 
          *elems[14]  += term * (-b[1]);
          *elems[15]  += term * ( b[0]);
          *elems[16]  += term * ( a[1]);
          *elems[17]  += term * (-a[0]);
          *elems[18]  += term * (-C[2]+B[2]);
          *elems[19]  += term * (-B[1]+C[1]);
          *elems[20]  += term * ( C[2]-A[2]);
          *elems[21]  += term * ( A[1]-C[1]);
          *elems[22]  += term * (-B[2]+A[2]);
          *elems[23]  += term * (-A[1]+B[1]);
          *elems[24]  += term * (-B[2]+C[2]);
          *elems[25]  += term * (-C[0]+B[0]);
          *elems[26]  += term * (-C[2]+A[2]);
          *elems[27]  += term * (-A[0]+C[0]);
          *elems[28]  += term * ( B[2]-A[2]);
          *elems[29]  += term * ( A[0]-B[0]);
          *elems[30]  += term * (-C[1]+B[1]);
          *elems[31]  += term * (-B[0]+C[0]);
          *elems[32]  += term * ( C[1]-A[1]);
          *elems[33]  += term * ( A[0]-C[0]);
          *elems[34]  += term * (-B[1]+A[1]);
          *elems[35]  += term * (-A[0]+B[0]);
 
          Eigen::Matrix<double,Eigen::Dynamic,1> g(12);

          // partial derivatives in respect to a
          g[0] = (c[1]* b[2] - c[2]* b[1]);
          g[1] = (c[2]* b[0] - c[0]* b[2]);
          g[2] = (c[0]* b[1] - c[1]* b[0]);
    
          // partial derivatives in respect to b
          g[3] = (a[1]* c[2] + a[2]*-c[1]);
          g[4] = (a[0]*-c[2] + a[2]* c[0]);
          g[5] = (a[0]* c[1] + a[1]*-c[0]);
    
          // partial derivatives in respect to c
          g[6] = (a[1]*-b[2] + a[2]* b[1]);
          g[7] = (a[0]* b[2] + a[2]*-b[0]);
          g[8] = (a[0]*-b[1] + a[1]* b[0]);
    
          // partial derivatives in respect to d
          g[9] =  (-c[1]*b[2] + c[2]* b[1] + a[1]*(-C[2]+B[2]) + a[2]*(-B[1]+C[1]));
          g[10] = (a[0]*(-B[2]+C[2]) - c[2]*b[0] + c[0]*b[2] + a[2]*(-C[0]+B[0]));
          g[11] = (a[0]*(-C[1]+B[1]) + a[1]*(-B[0]+C[0]) - c[0]*b[1] + c[1]*b[0]);

          Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> g2(12,12);
          g2 = g*g.transpose();

          // barrier function
          term = Beta*(2*factTerm*factTerm - divTerm*(6*coeffA*detEPS + 2*coeffB)) / divTerm3;
          g2 *= term; 

          int numElem = 0;
          for(int r=0;r<12;r++)
          {
            for(int c=r;c<12;c++)
            {
              *denseHessianCoeffs(numElem++,t) += g2(r,c);
            }
          }
        }
      }
    }
  }
}
