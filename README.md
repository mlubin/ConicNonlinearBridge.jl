# ConicNonlinearBridge

[![Build Status](https://travis-ci.org/mlubin/ConicNonlinearBridge.jl.svg?branch=master)](https://travis-ci.org/mlubin/ConicNonlinearBridge.jl) [![ConicNonlinearBridge](http://pkg.julialang.org/badges/ConicNonlinearBridge_0.5.svg)](http://pkg.julialang.org/?pkg=ConicNonlinearBridge)

This package implements a wrapper to allow derivative-based nonlinear solvers to function as [conic solvers](http://mathprogbasejl.readthedocs.org/en/latest/conic.html), in particular to be used from Convex.jl.

Example:
    
    using Convex, ConicNonlinearBridge, Ipopt
    x = Variable(1)
    problem = maximize(x,exp(x) <= 4)
    solve!(problem, ConicNLPWrapper(nlp_solver=IpoptSolver()))

    evaluate(x) # roughly log(4)

You may replace ``IpoptSolver`` above with your favorite NLP solver (e.g., Knitro, Mosek). You may also pass options to the solver, e.g., ``IpoptSolver(print_level=0)``.

This wrapper is experimental. If you are experiencing convergence troubles with existing conic solvers, this wrapper *may* be helpful. In general, however, specialized conic solvers are more reliable than derivative-based nonlinear solvers, especially for detection of infeasibility and unboundedness. If you find this wrapper useful, please let us know.

TODO:

  1. Phase I for better infeasibility detection
