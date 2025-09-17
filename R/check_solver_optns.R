check_solver_optns <- function(solver, solver_optns) {
  # If no options are provided use default values
  if (is.null(solver_optns)) {
    solver_optns <- switch (solver,
      "sinkhorn" = list(numIterations = 1e3,
                        epsilon = 0.01,
                        maxErr = 1e-9),
      "sinkhorn_log" = list(numIterations = 1e3,
                        epsilon = 0.01,
                        maxErr = 1e-9,
                        tau = 1e4),
      "transport" = list(fullreturn = TRUE)
    )

    return(solver_optns)
  }

  # If options are provided, check correctness for sinkhorn solvers
  if (solver == "sinkhorn" || solver == "sinkhorn_log") {
    stopifnot(all(names(solver_optns) %in% c("numIterations", "epsilon", "maxErr", "tau")))

    if (!exists("numIterations", solver_optns)) {
      solver_optns[["numIterations"]] <- 1e3
    }
    if (!exists("epsilon", solver_optns)) {
      solver_optns[["epsilon"]] <- 0.01
    }
    if (!exists("maxErr", solver_optns)) {
      solver_optns[["maxErr"]] <- 1e-9
    }
    if (!exists("tau", solver_optns) && solver == "sinkhorn_log") {
      solver_optns[["tau"]] <- 1e5
    }

    return(solver_optns)
  }

  # Last but not least, check correctness for wasserstein solver
  stopifnot(all(names(solver_optns) %in% c("method", "control", "threads")))
  solver_optns[["fullreturn"]] <- TRUE

  return(solver_optns)
}
