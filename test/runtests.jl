using ConicNonlinearBridge
using Convex
using Ipopt
using FactCheck

tests = ["test_utilities.jl",
         "test_affine.jl",
         "test_lp.jl"]
tests_socp = ["test_socp.jl"]
tests_exp = ["test_exp.jl"]

println("Running tests:")

# The following syntax can be used to solve it using other solvers
set_default_solver(ConicNLPWrapper(IpoptSolver(print_level=0)))


for curtest in tests
    info(" Test: $(curtest)")
    include(Pkg.dir("Convex","test",curtest))
end

for curtest in tests_socp
    info(" Test: $(curtest)")
    include(Pkg.dir("Convex","test",curtest))
end

for curtest in tests_exp
    info(" Test: $(curtest)")
    include(Pkg.dir("Convex","test",curtest))
end
