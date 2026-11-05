# Reference values for the derivative rules of the error functions and the exponential and trigonometric integrals in
# fastmath.calculus.autodiff (phase 5b): first and second derivatives, gradients and hessians of erf with two arguments.
#
# Tool: Julia ForwardDiff.jl 1.3.3 with DiffRules.jl 1.15.1 and SpecialFunctions.jl 2.7.2, Float64 precision, run in the
#       ~/.julia/environments/v1.12 project:
#
#   julia --project=$HOME/.julia/environments/v1.12 utils/fastmath/dev/generate_autodiff_erf_expint_reference.jl
#
# Output: test/resources/calculus/autodiff_erf_expint_reference.edn
# {:unary {<id> [[x d1 d2] ..]}
#  :En    [[n x d1 d2] ..]                                    derivatives in x, n constant
#  :erf2  [{:x [x1 x2] :gradient [..] :hessian [[..] ..]} ..]}
# The ids are the names of fastmath.calculus.autodiff (E1 = expint, Si = sinint, Ci = cosint). Lambert W, E0, Ein, Cin,
# si, li and Li have no Julia rule, and the ForwardDiff derivative of expinti (Ei) is off by 1e-13 at x = 7 and 6e-13
# at x = 30 (exp(x)/x and mpmath agree): these are in generate_autodiff_erf_expint_mpmath_reference.py.

using ForwardDiff, SpecialFunctions

d1(f, x) = ForwardDiff.derivative(f, x)
d2(f, x) = ForwardDiff.derivative(y -> ForwardDiff.derivative(f, y), x)

const general = [-3.7, -1.2, -0.3, 0.2, 0.9, 2.5, 7.0]
const positive = [0.01, 0.3, 1.0, 2.5, 7.0, 30.0]

const unary = [
    ("erf", erf, general), ("erfc", erfc, [-3.7, -1.2, -0.3, 0.2, 0.9, 2.5, 7.0, 20.0]),
    ("inv-erf", erfinv, [-0.999, -0.7, -0.1, 0.3, 0.9, 0.99999]),
    ("inv-erfc", erfcinv, [1.0e-10, 0.01, 0.4, 1.0, 1.6, 1.999]),
    ("E1", expint, positive),
    ("Si", sinint, [-30.0, -3.7, -0.3, 1.0e-3, 0.2, 2.5, 30.0]), ("Ci", cosint, positive)]

fmt(v::AbstractVector) = "[" * join(map(repr, v), " ") * "]"
fmt(m::AbstractMatrix) = "[" * join([fmt(m[i, :]) for i in 1:size(m, 1)], " ") * "]"

open(joinpath(@__DIR__, "..", "..", "..", "test", "resources", "calculus", "autodiff_erf_expint_reference.edn"), "w") do io
    println(io, "{:unary")
    println(io, " {")
    for (name, f, xs) in unary
        print(io, "  :", name, " [")
        for x in xs
            print(io, "\n   ", fmt([x, d1(f, x), d2(f, x)]))
        end
        println(io, "]")
    end
    println(io, " }")
    print(io, " :En [")
    for n in (0.0, 1.0, 2.0, 3.5, 7.0), x in (0.3, 1.0, 2.5, 10.0)
        f = y -> expint(n, y)
        print(io, "\n  ", fmt([n, x, d1(f, x), d2(f, x)]))
    end
    println(io, "]")
    print(io, " :erf2 [")
    for x in ([0.1, 0.3], [-1.0, 2.0], [3.0, 3.5], [-0.5, -0.4])
        g(v) = erf(v[1], v[2])
        print(io, "\n  {:x ", fmt(x), " :gradient ", fmt(ForwardDiff.gradient(g, x)), " :hessian ", fmt(ForwardDiff.hessian(g, x)), "}")
    end
    println(io, "]}")
end
