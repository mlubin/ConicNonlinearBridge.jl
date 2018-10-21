using ConicNonlinearBridge
using MathProgBase
using Ipopt
using Base.Test

include(Pkg.dir("MathProgBase", "test", "conicinterface.jl"))

@testset "linear and exponential cone tests (remove_single_rows is $remove_single_rows)" for remove_single_rows in [true, false]
    solver = ConicNLPWrapper(nlp_solver=IpoptSolver(print_level=0), remove_single_rows=remove_single_rows)
    coniclineartest(solver; duals=false, tol=1e-6)
    conicEXPtest(solver; duals=false, tol=1e-6)
end

@testset "second-order cone tests (soc_as_quadratic is $soc_as_quadratic)" for soc_as_quadratic in [true, false]
    solver = ConicNLPWrapper(nlp_solver=IpoptSolver(print_level=0), soc_as_quadratic=soc_as_quadratic, disaggregate_soc=false)
    conicSOCRotatedtest(solver; duals=false, tol=1e-6)

    @testset "disaggregation tests (disaggregate_soc is $disaggregate_soc)" for disaggregate_soc in [true, false]
        solver = ConicNLPWrapper(nlp_solver=IpoptSolver(print_level=0), soc_as_quadratic=soc_as_quadratic, disaggregate_soc=disaggregate_soc)
        conicSOCtest(solver; duals=false, tol=1e-6)
    end
end

return nothing
