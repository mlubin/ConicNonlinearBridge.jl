using ConicNonlinearBridge
using MathProgBase
using Ipopt
using Base.Test

solver = ConicNLPWrapper(nlp_solver=IpoptSolver(print_level=0))

include(Pkg.dir("MathProgBase", "test", "conicinterface.jl"))

for testfun in [coniclineartest, conicSOCtest, conicSOCRotatedtest, conicEXPtest]
    testfun(solver; duals=false, tol=1e-6)
end
