#ifndef IGL_SOLVER_STATUS_H
#define IGL_SOLVER_STATUS_H
namespace igl
{
  enum SolverStatus
  {
    // Good
    SOLVER_STATUS_CONVERGED = 0,
    // OK
    SOLVER_STATUS_MAX_ITER = 1,
    // Bad
    SOLVER_STATUS_ERROR = 2,
    NUM_SOLVER_STATUSES = 3,
  };
};
#endif
