# Reference values for the derivative rules of fastmath.calculus.autodiff: elementary functions and the gamma/beta
# functions with closed-form derivatives, first and second derivatives.
#
# Tool: Julia ForwardDiff.jl 1.3.3 with DiffRules.jl 1.15.1 and SpecialFunctions.jl 2.7.2, Float64 precision, run in the
#       ~/.julia/environments/v1.12 project:
#
#   julia --project=$HOME/.julia/environments/v1.12 utils/fastmath/dev/generate_autodiff_rules_reference.jl
#
# Output: test/resources/calculus/autodiff_rules_reference.edn
# {:unary   {<id> [[x d1 d2] ..]}                      first and second derivative
#  :polygamma [[m x d1 d2] ..]                         derivatives in x, m constant
#  :binary  {<id> [{:x [..] :gradient [..] :hessian [[..] ..]} ..]}
#  :ternary {<id> [{:x [..] :gradient [..] :hessian [[..] ..]} ..]}}
# The ids are the names of fastmath.calculus.autodiff. The points are inside the domains (ForwardDiff throws a
# DomainError outside), away from the kinks and the poles (tested by hand).

using ForwardDiff, SpecialFunctions

d1(f, x) = ForwardDiff.derivative(f, x)
d2(f, x) = ForwardDiff.derivative(y -> ForwardDiff.derivative(f, y), x)

const general = [-3.7, -1.2, -0.3, 0.2, 0.9, 2.5, 7.0]
const unit = [-0.95, -0.5, -0.1, 0.3, 0.8]
const positive = [0.05, 0.5, 1.5, 4.0, 30.0]
const gammapts = [-2.5, -1.5, -0.5, 0.1, 0.5, 1.0, 2.5, 7.0, 25.0]

const unary = [
    ("sqrt", sqrt, positive),
    ("cbrt", cbrt, [-8.0, -0.3, 0.2, 2.5, 27.0]),
    ("sq", x -> x * x, general),
    ("tan", tan, [-1.2, -0.3, 0.2, 0.9, 1.5]),
    ("asin", asin, unit),
    ("acos", acos, unit),
    ("atan", atan, general),
    ("sinh", sinh, general),
    ("cosh", cosh, general),
    ("tanh", tanh, general),
    ("asinh", asinh, general),
    ("acosh", acosh, [1.05, 1.5, 3.0, 20.0]),
    ("atanh", atanh, unit),
    ("expm1", expm1, general),
    ("log1p", log1p, [-0.9, -0.3, 0.2, 2.5, 30.0]),
    ("log2", log2, positive),
    ("log10", log10, positive),
    ("abs", abs, [-2.5, -0.3, 0.2, 7.0]),
    ("gamma", gamma, gammapts),
    ("log-gamma", x -> logabsgamma(x)[1], gammapts),
    ("log-gamma-1p", x -> logabsgamma(1 + x)[1], [-0.9, -0.4, 0.1, 0.5, 1.4, 3.0]),
    ("inv-gamma-1pm1", x -> 1 / gamma(1 + x) - 1, [-0.9, -0.4, 0.1, 0.5, 1.4, 3.0]),
    ("digamma", digamma, gammapts),
    ("trigamma", trigamma, gammapts)]

const binary = [
    ("atan2", (y, x) -> atan(y, x), [[1.0, 2.0], [-0.5, 0.3], [2.0, -1.0], [-1.5, -0.7]]),
    ("hypot", hypot, [[3.0, 4.0], [-0.5, 0.3], [1.0e-3, 2.0], [-1.5, -0.7]]),
    ("max", max, [[1.0, 2.0], [2.5, -1.0]]),
    ("min", min, [[1.0, 2.0], [2.5, -1.0]]),
    ("beta", beta, [[0.5, 0.5], [1.5, 2.0], [3.0, 0.7], [10.0, 4.5]]),
    ("log-beta", logbeta, [[0.5, 0.5], [1.5, 2.0], [3.0, 0.7], [10.0, 4.5]])]

const ternary = [
    ("hypot", (x, y, z) -> hypot(x, y, z), [[1.0, 2.0, 2.0], [-0.5, 0.3, 1.2], [1.0e-3, -2.0, 0.5]])]

fmt(v::AbstractVector) = "[" * join(map(repr, v), " ") * "]"
fmt(m::AbstractMatrix) = "[" * join([fmt(m[i, :]) for i in 1:size(m, 1)], " ") * "]"

function multivariate(io, name, f, points)
    g(v) = f(v...)
    print(io, "  :", name, " [")
    for x in points
        print(io, "\n   {:x ", fmt(x), " :gradient ", fmt(ForwardDiff.gradient(g, x)),
              " :hessian ", fmt(ForwardDiff.hessian(g, x)), "}")
    end
    println(io, "]")
end

open(joinpath(@__DIR__, "..", "..", "..", "test", "resources", "calculus", "autodiff_rules_reference.edn"), "w") do io
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
    print(io, " :polygamma [")
    for m in (0, 1, 2, 3), x in (-1.5, 0.5, 2.5, 10.0)
        f = y -> polygamma(m, y)
        print(io, "\n  ", fmt([Float64(m), x, d1(f, x), d2(f, x)]))
    end
    println(io, "]")
    println(io, " :binary")
    println(io, " {")
    for (name, f, xs) in binary
        multivariate(io, name, f, xs)
    end
    println(io, " }")
    println(io, " :ternary")
    println(io, " {")
    for (name, f, xs) in ternary
        multivariate(io, name, f, xs)
    end
    println(io, " }}")
end
