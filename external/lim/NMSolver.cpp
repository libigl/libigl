// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich

#include "NMSolver.h"

NMSolver::NMSolver()
{
}

int NMSolver::init()
{
  // state of last step
  CurrentLambda = 0;
  CurrentStepSize = 1;
  CurrentFV = 0;

  // newton step parameters
  maxStepSize = 1;
  minStepSize = 1e-21;
  stepSize = 1.0;
  
  // hessian correction parameters
  lambda = 1e-1;
  minLambda = 1e-5;
  maxLambda = 10e16;
  
  // Wolfe condition step size parameters
  EnableWolfeConditions = false;
  wc1 = 0.001;
  wc2 = 0.1;
  wAlphaMax = 100;

  isPrevHessianSingular= true;
  
  // Init problem
  getNMProblemSize();

  solution.resize(numVariables);
  prevSolution.resize(numVariables);
  gradient.resize(numVariables);
  hessian.resize(numVariables,numVariables);
  step.resize(numVariables);
  
  prepareNMProblemData();

  lltSolver.analyzePattern(hessian);
  
  return 1;
}

int NMSolver::solve()
{
  // compute function value, gradient and hessian at current point
  functionValue = computeNMFunction(solution);
  computeNMGradient(solution,gradient);
  computeNMHessian(solution);
  //validateWithFD(solution);

  subSolve();

  if(CurrentFV == std::numeric_limits<double>::infinity())
    return -1;
  else
    return 1;

  return true;
}

bool NMSolver::subSolve()
{
  // factorize hessian
  lltSolver.factorize(hessian);
  int status = lltSolver.info();

  if(status == Eigen::Success)
  {
    lambda = minLambda;
    CurrentLambda = 0;
  }
  else
  {
    int run = 0;
    while(status != Eigen::Success && lambda < maxLambda)
    {
      // singular hessian -> trust region lambda approach
      for(int i=0;i<numVariables;i++)
        *diagHessianCoeffs[i] += lambda;

      // try to factorize hessian
      lltSolver.factorize(hessian);
      status = lltSolver.info();
    
      if(status != Eigen::Success)
      {
        if(lambda < maxLambda)
          lambda *= 10;
      }
      else 
      {
        if(CurrentLambda != lambda)
          stepSize = 1;
        
        CurrentLambda = lambda;
        if(run == 0 && lambda > minLambda)
        {
          lambda *= 0.1;
        }
        
      }
      
      run++;
    }
  }

  // call functions depending on factorization
  afterHessianFactorization();
  
  if(lambda < maxLambda)
  {
    // compute newton step direction
    step = lltSolver.solve(-gradient);			
  }
  else
  {
    // gradient descent
    CurrentLambda = -1;
    lambda *= 0.1;
    step = -gradient;
  }

  // find step size
  if(EnableWolfeConditions)
  {
    // wolfe condtions
    double sS = computeStepSizeWolfe(solution,functionValue,step);
    solution += sS*step;
    
    CurrentFV = computeNMFunction(solution);
    CurrentStepSize = sS;
  }
  else
  {
    // adaptive backtracking
    bool repeat = true;
    int run = 0;
    while(repeat)
    {
      prevSolution = solution + stepSize*step;
      double newfval = computeNMFunction(prevSolution);

      if(newfval < functionValue || stepSize < minStepSize)
      {
        CurrentFV = newfval;
        CurrentStepSize = stepSize;
        repeat = false;
        if(run == 0 && stepSize < maxStepSize)
          stepSize *= 2;
      }
      else
      {
        stepSize *= 0.5f;
      }
      run++;
    }
    solution = prevSolution;
  }

  return true;
}

double NMSolver::computeStepSizeWolfe(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, double of_x, const Eigen::Matrix<double,Eigen::Dynamic,1>& dir)
{
  double alpha = 1;
  double alpha_0 = 0;
  double alpha_1 = alpha;

  double of_0 = of_x;
  
  Eigen::Matrix<double,Eigen::Dynamic,1> gt(numVariables);
  Eigen::Matrix<double,Eigen::Dynamic,1> xc(numVariables);
  
  computeNMGradient(x,gt);
  double slope0 = gt.dot(dir);

  int iter = 0;
  while(true)
  {
    xc = x+alpha_1*dir;
    double of = computeNMFunction(xc);
        
    // check if current iterate violates sufficient decrease
    if((of > of_0 + wc1*alpha_1*slope0) || ((of >= of_x ) && (iter > 0)))
    {
      // there has to be an acceptable point between alpha_0 and alpha_1
      // (because c1 > c2)
      alpha = nocZoomWolfe(x,dir,slope0,alpha_0,alpha_1,of_0,of_x,wc1,wc2);
      break; 
    }

    computeNMGradient(xc,gt);
    double slopec = gt.dot(dir);

    // current iterate has sufficient decrease, but are we too close?
    if(abs(slopec) <= -wc2*slope0)
    {
      // strong wolfe fullfilled, quit
      alpha = alpha_1;
      break;
    }

    // are we behind the minimum?
    if (slopec >= 0)
    {
      // there has to be an acceptable point between alpha_0 and alpha_1
      alpha = nocZoomWolfe(x,dir,slope0,alpha_1,alpha_0,of_0,of,wc1,wc2);
      break;
    }

    alpha_0 = alpha_1;
    alpha_1 = std::min(wAlphaMax,alpha_1*3);
    of_x = of;

    iter++;
  }

  return alpha;
}

double NMSolver::nocZoomWolfe(const Eigen::Matrix<double,Eigen::Dynamic,1>& x, const Eigen::Matrix<double,Eigen::Dynamic,1>& dir, double slope0, double alphaLo, double alphaHi, double of_0, double ofLo, double c1, double c2)
{
  double alpha;
  double of = 0;
  for(int i=0;i<10 || of == std::numeric_limits<double>::infinity();i++)
  {
    alpha = (alphaLo+alphaHi)/2;
    Eigen::Matrix<double,Eigen::Dynamic,1> xc = x + alpha*dir;
    of = computeNMFunction(xc);
    
    if((of > of_0 + c1*alpha*slope0) || (of >= ofLo))
    {
      // if we do not observe sufficient decrease in point alpha, we set
      // the maximum of the feasible interval to alpha
      alphaHi = alpha;
    }
    else
    {
      Eigen::Matrix<double,Eigen::Dynamic,1> gt(numVariables);
      computeNMGradient(xc,gt);
      double slopec = gt.dot(dir);
      
      // strong wolfe fullfilled?
      if(abs(slopec) <= -c2*slope0)
        break;
      
      // if slope positive and alphaHi > alphaLo
      if(slopec*(alphaHi-alphaLo) >= 0)
        alphaHi = alphaLo;
      
      alphaLo = alpha;
      ofLo    = of;
    }
  }

  if(computeNMFunction(x + alpha*dir) == std::numeric_limits<double>::infinity())
    int test = 0;

  return alpha;
}

void NMSolver::validateWithFD(const Eigen::Matrix<double,Eigen::Dynamic,1>& x)
{
  double eps  = 1e-3;
  double eps2 = eps*eps;
  Eigen::Matrix<double,Eigen::Dynamic,1> fx(x);

  // check gradient
  std::cout << "Check gradient with finite differencies\n";
  double fv = computeNMFunction(fx);
  for(int i=0;i<numVariables;i++)
  {
    fx[i] = x[i] + eps;
    double fv0 = computeNMFunction(fx);		
    double g = (fv0 - fv) / eps;

    if(fv0 == std::numeric_limits<double>::infinity())
    {
      fx[i] = x[i] - eps;
      fv0 = computeNMFunction(fx);
      g = (fv0 - fv) / -eps;
    }
    fx[i] = x[i];

    double max = abs(fv) > abs(fv0) ? fv : fv0;
    double diff = gradient[i] - g;
    if(max > 1) diff /= abs(max);

    std::cout << i << ": " << diff << std::endl;
  }
}

void NMSolver::afterHessianFactorization()
{
}
