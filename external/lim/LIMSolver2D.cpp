// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "LIMSolver2D.h"
#include "TriangleMesh.h"

LIMSolver2D::LIMSolver2D()
{
  dim = 2;
  mesh = NULL;
}

LIMSolver2D::~LIMSolver2D()
{
}

void LIMSolver2D::Init(DeformableMesh* mesh)
{
  Init(static_cast<TriangleMesh*>(mesh));
}

void LIMSolver2D::Init(TriangleMesh* mesh)
{
  this->mesh = mesh;
  LIMSolver::Init(mesh);
}

void LIMSolver2D::prepareNMProblemData()
{
  std::vector<Eigen::Triplet<double> > triplets;
  const int numNodes = mesh->InitalVertices->rows();

  TriangleVertexIdx.resize(6,mesh->Triangles->rows());
  for(int t=0;t<mesh->Triangles->rows();t++)
  {
    Eigen::Vector3i nodes = mesh->Triangles->row(t);

    // create tets vertex indicies
    for(int v=0;v<3;v++)
    {
      for(int i=0;i<2;i++)
      {
        TriangleVertexIdx(v*2+i,t) = nodes[v]*2+i;
      }
    }

    for(int r=0;r<6;r++)
    {
      for(int c=r;c<6;c++)
      {
        int row = TriangleVertexIdx(r,t);
        int col = TriangleVertexIdx(c,t);
      
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
    for(int i=0;i<2;i++)
    {
      VertexPositionIndices[n*2+i] = n+numNodes*i;
    }
  }

  initialNodes.resize(numVariables);
  for(int n=0;n<numVariables;n++)
  {
    initialNodes[n] = mesh->InitalVertices->coeff(VertexPositionIndices[n]);
  }
  
  // set link indices of non zero elements of constraint hessian matrix
  problemHessianCoeffs.resize(1,hessRowIdx.size());
  denseHessianCoeffs.resize(21,mesh->Triangles->rows());
  diagHessianCoeffs.resize(numVariables);
  posConstraintsHessianCoeffs.resize(linearConstraintsMatrix2.nonZeros());
  
  for(int i=0;i<hessRowIdx.size();i++)
  {
    if(hessRowIdx[i] <= hessColIdx[i])
      problemHessianCoeffs[i] = &hessian.coeffRef(hessRowIdx[i],hessColIdx[i]);
  }
  
  for(int t=0;t<mesh->Triangles->rows();t++)
  {
    int numElem = 0;
    for(int r=0;r<6;r++)
    {
      for(int c=r;c<6;c++)
      {
        int row = TriangleVertexIdx(r,t);
        int col = TriangleVertexIdx(c,t);

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
  nonFlipHessianCoeffs.resize(6,mesh->Triangles->rows());
  for(int t=0;t<mesh->Triangles->rows();t++)
  {
    for(int i=0;i<6;i++)
    {
      int row = TriangleVertexIdx(NonFlipHessian2DIdx[i][0],t);
      int col = TriangleVertexIdx(NonFlipHessian2DIdx[i][1],t);
    
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

void LIMSolver2D::computeRestPoseFunctionParameters()
{
  int numTets = mesh->Triangles->rows();
  
  initalSize.resize(numTets);
  barrierParam1.resize(numTets);
  barrierParam2.resize(numTets);
  
  for(int t=0;t<numTets;t++)
  {
    Eigen::VectorXi indices = TriangleVertexIdx.col(t);
    
    Eigen::Vector3d A = mesh->InitalVertices->row(mesh->Triangles->coeff(t,0));
    Eigen::Vector3d B = mesh->InitalVertices->row(mesh->Triangles->coeff(t,1));
    Eigen::Vector3d C = mesh->InitalVertices->row(mesh->Triangles->coeff(t,2));

    double area = ((A-C).cross(B-C)).norm();
    
    double areaEPS = area - mesh->EPS3;
    initalSize[t] = areaEPS;

    double quot = areaEPS*CompensationExp*pow(area,CompensationExp-1);

    // must be later multiplied by beta
    barrierParam1[t] = 1/quot;
    barrierParam2[t] = log(areaEPS) - pow(area,CompensationExp)/quot;
  }
}

double LIMSolver2D::computeNMFunction(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  double dObj;

  // position constraints
  double diff = (*dmesh->ConstraintMatrix * x - subStepConstraints).squaredNorm();

  // non-flip constraint energy
  double ic = 0.0;
  if(EnableBarriers)
  {
    if(EnableNeoHookeanBarriers)
    {
      for(int t=0;t<mesh->Triangles->rows();t++)
      {
        Eigen::Vector2d a(x[TriangleVertexIdx.coeff(0,t)]-x[TriangleVertexIdx.coeff(4,t)],x[TriangleVertexIdx.coeff(1,t)]-x[TriangleVertexIdx.coeff(5,t)]);
        Eigen::Vector2d b(x[TriangleVertexIdx.coeff(2,t)]-x[TriangleVertexIdx.coeff(4,t)],x[TriangleVertexIdx.coeff(3,t)]-x[TriangleVertexIdx.coeff(5,t)]);
        
        Eigen::Matrix<double,2,2> D;
        D.col(0) = a;
        D.col(1) = b;
      
        double det = D.determinant();
        double detEPS = det - mesh->EPS3;
        
        if(detEPS > 0.0)
        {
          double s = 1/(initalSize[t]-mesh->EPS3);
          double logdet = log(s*detEPS);
          ic += logdet*logdet; // barrier function
        }
        else
        {
          ic = std::numeric_limits<double>::infinity();
          break;
        }
      }
    }
    else if(EnableLogBarriers)
    {
      for(int t=0;t<mesh->Triangles->rows();t++)
      {
        Eigen::Vector2d a(x[TriangleVertexIdx.coeff(0,t)]-x[TriangleVertexIdx.coeff(4,t)],x[TriangleVertexIdx.coeff(1,t)]-x[TriangleVertexIdx.coeff(5,t)]);
        Eigen::Vector2d b(x[TriangleVertexIdx.coeff(2,t)]-x[TriangleVertexIdx.coeff(4,t)],x[TriangleVertexIdx.coeff(3,t)]-x[TriangleVertexIdx.coeff(5,t)]);
        
        Eigen::Matrix<double,2,2> D;
        D.col(0) = a;
        D.col(1) = b;
      
        double det = D.determinant();
        double detEPS = det - mesh->EPS3;

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
      for(int t=0;t<mesh->Triangles->rows();t++)
      {
        Eigen::Vector2d a(x[TriangleVertexIdx.coeff(0,t)]-x[TriangleVertexIdx.coeff(4,t)],x[TriangleVertexIdx.coeff(1,t)]-x[TriangleVertexIdx.coeff(5,t)]);
        Eigen::Vector2d b(x[TriangleVertexIdx.coeff(2,t)]-x[TriangleVertexIdx.coeff(4,t)],x[TriangleVertexIdx.coeff(3,t)]-x[TriangleVertexIdx.coeff(5,t)]);
        
        Eigen::Matrix<double,2,2> D;
        D.col(0) = a;
        D.col(1) = b;
      
        double det = D.determinant();
        double detEPS = det - mesh->EPS3;

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

void LIMSolver2D::computeNMGradient(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, Eigen::Matrix<double,Eigen::Dynamic,1>& grad)
{
  grad.setZero();

  // compute problem function value
  computeGradient(x,grad);
  
  Eigen::Matrix<double,Eigen::Dynamic,1> temp1 = linearConstraintsMatrix2 * x;
  Eigen::Matrix<double,Eigen::Dynamic,1> temp2 = subStepConstraints.transpose() * *dmesh->ConstraintMatrix;
  grad += 2.0*Alpha*(temp1 - temp2);

  // non-flip constraints
  if(EnableBarriers)
  {
    if(EnableNeoHookeanBarriers)
    {
      for(int t=0;t<mesh->Triangles->rows();t++)
      {
        Eigen::Matrix<int,6,1> indices = TriangleVertexIdx.col(t);
    
        Eigen::Vector2d A(x[indices[0]],x[indices[1]]);
        Eigen::Vector2d B(x[indices[2]],x[indices[3]]);
        Eigen::Vector2d C(x[indices[4]],x[indices[5]]);

        Eigen::Vector2d a = A-C;
        Eigen::Vector2d b = B-C;

        Eigen::Matrix<double,2,2> D;
        D.col(0) = a;
        D.col(1) = b;
    
        double det = D.determinant();
        double detEPS = det - mesh->EPS3;
        double s = 1/(initalSize[t]-mesh->EPS3);

        double term = 2*log(s*detEPS)/detEPS; // barrier function
        term *= Beta;

        // partial derivatives in respect to a
        grad[indices[0]] += term *  b[1];
        grad[indices[1]] += term * -b[0];
    
        // partial derivatives in respect to b
        grad[indices[2]] += term * -a[1];
        grad[indices[3]] += term *  a[0];
    
        // partial derivatives in respect to c
        grad[indices[4]] += term * (A[1]-B[1]);
        grad[indices[5]] += term * (B[0]-A[0]);
      }
    }
    else if(EnableLogBarriers)
    {
      for(int t=0;t<mesh->Triangles->rows();t++)
      {
        Eigen::Matrix<int,6,1> indices = TriangleVertexIdx.col(t);
    
        Eigen::Vector2d A(x[indices[0]],x[indices[1]]);
        Eigen::Vector2d B(x[indices[2]],x[indices[3]]);
        Eigen::Vector2d C(x[indices[4]],x[indices[5]]);

        Eigen::Vector2d a = A-C;
        Eigen::Vector2d b = B-C;

        Eigen::Matrix<double,2,2> D;
        D.col(0) = a;
        D.col(1) = b;
    
        double det = D.determinant();
        double detEPS = det-mesh->EPS3;

        double term = -1/detEPS; // barrier function
        if(EnableBarrierCompensation) term += barrierParam1[t]*CompensationExp*pow((double)det,CompensationExp-1); // barrier compensation function
        term *= Beta;

        // partial derivatives in respect to a
        grad[indices[0]] += term *  b[1];
        grad[indices[1]] += term * -b[0];
    
        // partial derivatives in respect to b
        grad[indices[2]] += term * -a[1];
        grad[indices[3]] += term *  a[0];
    
        // partial derivatives in respect to c
        grad[indices[4]] += term * (A[1]-B[1]);
        grad[indices[5]] += term * (B[0]-A[0]);
      }
    }
    else
    {
      for(int t=0;t<mesh->Triangles->rows();t++)
      {
        Eigen::Matrix<int,6,1> indices = TriangleVertexIdx.col(t);
    
        Eigen::Vector2d A(x[indices[0]],x[indices[1]]);
        Eigen::Vector2d B(x[indices[2]],x[indices[3]]);
        Eigen::Vector2d C(x[indices[4]],x[indices[5]]);

        Eigen::Vector2d a = A-C;
        Eigen::Vector2d b = B-C;

        Eigen::Matrix<double,2,2> D;
        D.col(0) = a;
        D.col(1) = b;
    
        double det = D.determinant();
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
          grad[indices[0]] += term *  b[1];
          grad[indices[1]] += term * -b[0];
    
          // partial derivatives in respect to b
          grad[indices[2]] += term * -a[1];
          grad[indices[3]] += term *  a[0];
    
          // partial derivatives in respect to c
          grad[indices[4]] += term * (A[1]-B[1]);
          grad[indices[5]] += term * (B[0]-A[0]);
        }
      }	
    }
  }
}

void LIMSolver2D::computeNMHessian(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  hessian *= 0;

  // compute problem function hessian
  computeHessian(x,problemHessianCoeffs);

  // position constraints:
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
    if(EnableNeoHookeanBarriers)
    {
      for(int t=0;t<mesh->Triangles->rows();t++)
      {
        Eigen::Matrix<int,6,1> indices = TriangleVertexIdx.col(t);
    
        Eigen::Vector2d A(x[indices[0]],x[indices[1]]);
        Eigen::Vector2d B(x[indices[2]],x[indices[3]]);
        Eigen::Vector2d C(x[indices[4]],x[indices[5]]);

        Eigen::Vector2d a = A-C;
        Eigen::Vector2d b = B-C;

        Eigen::Matrix<double,2,2> D;
        D.col(0) = a;
        D.col(1) = b;
    
        double det = D.determinant();
        double detEPS = det - mesh->EPS3;
        double s = 1/(initalSize[t]-mesh->EPS3);
        double logdet = 2*log(s*detEPS);
        double term = logdet/detEPS; // barrier function
        term *= Beta;

        Eigen::Matrix<double*,Eigen::Dynamic,1> elems = nonFlipHessianCoeffs.col(t);
    
        *elems[0] += -term;
        *elems[1] +=  term;
        *elems[2] +=  term;
        *elems[3] += -term;
        *elems[4] += -term;
        *elems[5] +=  term;
 
        Eigen::Matrix<double,Eigen::Dynamic,1> g(6);

        // partial derivatives in respect to a
        g[0] =  b[1];
        g[1] = -b[0];
    
        // partial derivatives in respect to b
        g[2] = -a[1];
        g[3] =  a[0];
    
        // partial derivatives in respect to c
        g[4] = (A[1]-B[1]);
        g[5] = (B[0]-A[0]);

        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> g2(6,6);
        g2 = g*g.transpose();
        
        term = (2+logdet)/(detEPS*detEPS); // barrier function
        g2 *= Beta*term; 

        int numElem = 0;
        for(int r=0;r<6;r++)
        {
          for(int c=r;c<6;c++)
          {
            *denseHessianCoeffs(numElem++,t) += g2(r,c);
          }
        }
      }
    }
    else if(EnableLogBarriers)
    {
      for(int t=0;t<mesh->Triangles->rows();t++)
      {
        Eigen::Matrix<int,6,1> indices = TriangleVertexIdx.col(t);
    
        Eigen::Vector2d A(x[indices[0]],x[indices[1]]);
        Eigen::Vector2d B(x[indices[2]],x[indices[3]]);
        Eigen::Vector2d C(x[indices[4]],x[indices[5]]);

        Eigen::Vector2d a = A-C;
        Eigen::Vector2d b = B-C;

        Eigen::Matrix<double,2,2> D;
        D.col(0) = a;
        D.col(1) = b;
    
        double det = D.determinant();
        double detEPS = det-mesh->EPS3;
        double detEPSInv = 1/detEPS;
        double term = -detEPSInv; // barrier function
        if(EnableBarrierCompensation) term += barrierParam1[t]*CompensationExp*pow((double)det,CompensationExp-1); // barrier compensation function
        term *= Beta;

        Eigen::Matrix<double*,Eigen::Dynamic,1> elems = nonFlipHessianCoeffs.col(t);
    
        *elems[0] += -term;
        *elems[1] +=  term;
        *elems[2] +=  term;
        *elems[3] += -term;
        *elems[4] += -term;
        *elems[5] +=  term;
 
        Eigen::Matrix<double,Eigen::Dynamic,1> g(6);

        // partial derivatives in respect to a
        g[0] =  b[1];
        g[1] = -b[0];
    
        // partial derivatives in respect to b
        g[2] = -a[1];
        g[3] =  a[0];
    
        // partial derivatives in respect to c
        g[4] = (A[1]-B[1]);
        g[5] = (B[0]-A[0]);

        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> g2(6,6);
        g2 = g*g.transpose();
        term = detEPSInv*detEPSInv; // barrier function
        if(EnableBarrierCompensation) term += barrierParam1[t]*CompensationExp*(CompensationExp-1)*pow((double)det,CompensationExp-2); // barrier compensation function 
        g2 *= Beta*term; 

        int numElem = 0;
        for(int r=0;r<6;r++)
        {
          for(int c=r;c<6;c++)
          {
            *denseHessianCoeffs(numElem++,t) += g2(r,c);
          }
        }
      }
    }
    else
    {
      for(int t=0;t<mesh->Triangles->rows();t++)
      {
        Eigen::Matrix<int,6,1> indices = TriangleVertexIdx.col(t);
    
        Eigen::Vector2d A(x[indices[0]],x[indices[1]]);
        Eigen::Vector2d B(x[indices[2]],x[indices[3]]);
        Eigen::Vector2d C(x[indices[4]],x[indices[5]]);

        Eigen::Vector2d a = A-C;
        Eigen::Vector2d b = B-C;

        Eigen::Matrix<double,2,2> D;
        D.col(0) = a;
        D.col(1) = b;
    
        double det = D.determinant();
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
    
          *elems[0] += -term;
          *elems[1] +=  term;
          *elems[2] +=  term;
          *elems[3] += -term;
          *elems[4] += -term;
          *elems[5] +=  term;
 
          Eigen::Matrix<double,Eigen::Dynamic,1> g(6);

          // partial derivatives in respect to a
          g[0] =  b[1];
          g[1] = -b[0];
    
          // partial derivatives in respect to b
          g[2] = -a[1];
          g[3] =  a[0];
    
          // partial derivatives in respect to c
          g[4] = (A[1]-B[1]);
          g[5] = (B[0]-A[0]);

          Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> g2(6,6);
          g2 = g*g.transpose();
          
          // barrier function
          term = Beta*(2*factTerm*factTerm - divTerm*(6*coeffA*detEPS + 2*coeffB)) / divTerm3;
          g2 *= term; 

          int numElem = 0;
          for(int r=0;r<6;r++)
          {
            for(int c=r;c<6;c++)
            {
              *denseHessianCoeffs(numElem++,t) += g2(r,c);
            }
          }
        }
      }
    }
  }
}
