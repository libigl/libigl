py::enum_<igl::SolverStatus>(m, "SolverStatus")
    .value("SOLVER_STATUS_CONVERGED", igl::SOLVER_STATUS_CONVERGED)
    .value("SOLVER_STATUS_MAX_ITER", igl::SOLVER_STATUS_MAX_ITER)
    .value("SOLVER_STATUS_ERROR", igl::SOLVER_STATUS_ERROR)
    .value("NUM_SOLVER_STATUSES", igl::NUM_SOLVER_STATUSES)
    .export_values();
