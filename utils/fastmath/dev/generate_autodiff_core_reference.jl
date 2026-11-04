# Reference values for the derivative rules of the fastmath.core functions in fastmath.calculus.autodiff (phase 5a):
# first and second derivatives of the functions of one argument, gradients and hessians of the functions of several
# arguments and of sequences.
#
# Tool: Julia ForwardDiff.jl 1.3.3 with DiffRules.jl 1.15.1 and LogExpFunctions.jl, Float64 precision, run in the
#       ~/.julia/environments/v1.12 project:
#
#   julia --project=$HOME/.julia/environments/v1.12 utils/fastmath/dev/generate_autodiff_core_reference.jl
#
# Output: test/resources/calculus/autodiff_core_reference.edn
# {:unary    {<id> [[x d1 d2] ..]}                               first and second derivative
#  :multi    {<id> [{:x [..] :gradient [..] :hessian [[..] ..]} ..]}  functions of 2, 3 and 5 arguments
#  :sequence {<id> [{:xs [..] :alpha a :gradient [..] :hessian [[..] ..]} ..]}}  in the elements, alpha constant
# The ids are the names of fastmath.calculus.autodiff (smooth-max as smooth-max-<family>). Where Julia has no rule
# (tanpi, cotpi, secpi, cscpi, exprel, expexp, cexpexp, cloglog, loglog, frac, sfrac, spow, safe-sqrt, lerp,
# smoothstep, norm, smooth-max) the defining formula is differentiated. The points are inside the domains, away from
# the kinks, jumps and removable singularities (tested by hand and against mpmath).

using ForwardDiff, LogExpFunctions

d1(f, x) = ForwardDiff.derivative(f, x)
d2(f, x) = ForwardDiff.derivative(y -> ForwardDiff.derivative(f, y), x)

const general = [-3.7, -1.2, -0.3, 0.2, 0.9, 2.5, 7.0]
const nonzero = [-3.7, -1.2, -0.3, 0.2, 0.9, 2.5]
const outside_unit = [-3.7, -1.2, 1.05, 2.5, 30.0]
const positive = [0.05, 0.5, 1.5, 4.0, 30.0]
const probability = [0.01, 0.2, 0.5, 0.7, 0.99]
const negative = [-30.0, -2.0, -0.7, -0.1, -1.0e-3]
const halfturns = [-0.7, -0.3, 0.1, 0.2, 0.45, 1.3]
const noninteger = [-3.7, -1.2, -0.3, 0.2, 0.9, 2.5]

const unary = [
    ("sec", sec, nonzero), ("csc", csc, nonzero), ("cot", cot, nonzero),
    ("sech", sech, general), ("csch", csch, nonzero), ("coth", coth, nonzero),
    ("asec", asec, outside_unit), ("acsc", acsc, outside_unit), ("acot", acot, general),
    ("asech", asech, [0.05, 0.3, 0.6, 0.9]), ("acsch", acsch, nonzero), ("acoth", acoth, outside_unit),
    ("sinpi", sinpi, general), ("cospi", cospi, general),
    ("tanpi", x -> tan(pi * x), halfturns), ("cotpi", x -> cot(pi * x), halfturns),
    ("secpi", x -> sec(pi * x), halfturns), ("cscpi", x -> csc(pi * x), halfturns),
    ("sinc", sinc, [-3.7, -1.2, -0.75, 0.6, 0.9, 2.5, 7.3]),
    ("exp2", exp2, general), ("exp10", exp10, [-3.7, -1.2, -0.3, 0.2, 0.9, 2.5]),
    ("ln", log, positive),
    ("xlogx", xlogx, positive), ("log1pexp", log1pexp, [-40.0, -3.7, -0.3, 0.2, 2.5, 20.0, 40.0]),
    ("logexpm1", logexpm1, positive), ("log1mexp", log1mexp, negative),
    ("log1pmx", log1pmx, [-0.9, -0.5, -0.3, 0.2, 0.5, 0.95, 2.5, 30.0]), ("log1psq", log1psq, general),
    ("logcosh", logcosh, general), ("exprel", x -> expm1(x) / x, [-30.0, -3.7, -1.2, 0.6, 0.9, 2.5, 30.0]),
    ("expexp", x -> exp(-exp(-x)), general), ("cexpexp", x -> -expm1(-exp(x)), general),
    ("cloglog", x -> log(-log1p(-x)), probability), ("loglog", x -> -log(-log(x)), probability),
    ("sigmoid", logistic, general), ("logit", logit, probability),
    ("trunc", trunc, noninteger), ("rint", round, noninteger), ("round", x -> floor(x + 0.5), noninteger),
    ("frac", x -> abs(x - trunc(x)), noninteger), ("sfrac", x -> x - trunc(x), noninteger),
    ("cb", x -> x^3, general), ("pow2", x -> x^2, general), ("pow3", x -> x^3, general),
    ("pow10", x -> x^10, general), ("safe-sqrt", x -> x < 0 ? zero(x) : sqrt(x), [-2.0, 0.05, 0.5, 4.0])]

const multi = [
    ("logb", (b, x) -> log(b, x), [[2.0, 3.0], [0.5, 7.0], [10.0, 0.2], [3.0, 1.0]]),
    ("xlogy", xlogy, [[2.0, 3.0], [-0.5, 7.0], [1.5, 0.2]]),
    ("xlog1py", xlog1py, [[2.0, 3.0], [-0.5, 7.0], [1.5, -0.2]]),
    ("logaddexp", logaddexp, [[1.0, 2.0], [-3.0, 0.5], [30.0, -30.0], [2.0, 2.0]]),
    ("logsubexp", logsubexp, [[1.0, 2.0], [-3.0, 0.5], [30.0, -30.0], [2.0, -1.0]]),
    ("mod", mod, [[5.5, 2.0], [-5.5, 2.0], [5.5, -2.0], [0.3, 0.7], [7.2, 1.5]]),
    ("rem", rem, [[5.5, 2.0], [-5.5, 2.0], [5.5, -2.0], [0.3, 0.7], [7.2, 1.5]]),
    ("copy-sign", copysign, [[2.0, -1.0], [-2.0, 3.0], [1.5, 0.5], [-0.5, -4.0]]),
    ("spow", (x, e) -> sign(x) * abs(x)^e, [[2.0, 1.5], [-2.0, 1.5], [0.5, 3.0], [-0.7, 0.4], [3.0, -1.2]]),
    ("fma", fma, [[1.5, -2.0, 0.3], [0.2, 4.0, -7.0]]),
    ("muladd", muladd, [[1.5, -2.0, 0.3], [0.2, 4.0, -7.0]]),
    ("lerp", (a, b, t) -> a + t * (b - a), [[1.0, 3.0, 0.25], [-2.0, 5.0, 1.5], [0.5, -0.5, -0.2]]),
    ("smoothstep", (e0, e1, x) -> (t = clamp((x - e0) / (e1 - e0), 0, 1); t * t * (3 - 2 * t)),
     [[0.0, 1.0, 0.3], [-1.0, 2.0, 1.4], [2.0, -1.0, 0.5], [0.0, 1.0, 1.5], [0.0, 1.0, -0.5]]),
    ("norm", (v, s, e) -> (v - s) / (e - s), [[0.3, 0.0, 1.0], [5.0, -2.0, 3.0], [-1.0, 4.0, 2.0]]),
    ("norm5", (v, s1, e1, s2, e2) -> s2 + (e2 - s2) * ((v - s1) / (e1 - s1)),
     [[0.3, 0.0, 1.0, 10.0, 20.0], [5.0, -2.0, 3.0, 1.0, -1.0], [-1.0, 4.0, 2.0, 0.5, 0.25]])]

function smu(xs, alpha)
    epsilon = abs(1 / alpha)
    reduce((a, b) -> alpha > 0 ? 0.5 * (a + b + sqrt((a - b)^2 + epsilon)) : 0.5 * (a + b - sqrt((a - b)^2 + epsilon)), xs)
end

const sequence = [
    ("logsumexp", (xs, alpha) -> logsumexp(xs)),
    ("smooth-max-lse", (xs, alpha) -> logsumexp(alpha .* xs) / alpha),
    ("smooth-max-mellowmax", (xs, alpha) -> (logsumexp(alpha .* xs) - log(length(xs))) / alpha),
    ("smooth-max-boltzmann", (xs, alpha) -> sum(xs .* exp.(alpha .* xs)) / sum(exp.(alpha .* xs))),
    ("smooth-max-p-norm", (xs, alpha) -> sum(abs.(xs) .^ alpha)^(1 / alpha)),
    ("smooth-max-smu", smu)]

const sequence_points = [([1.0, 2.0, 0.5], 2.0), ([-3.0, 0.5, 2.5, -0.2], 0.7), ([0.3, -1.2], -1.5),
                         ([4.0, 3.9, -2.0, 1.0, 0.1], 5.0), ([2.5], 2.0)]

fmt(v::AbstractVector) = "[" * join(map(repr, v), " ") * "]"
fmt(m::AbstractMatrix) = "[" * join([fmt(m[i, :]) for i in 1:size(m, 1)], " ") * "]"

open(joinpath(@__DIR__, "..", "..", "..", "test", "resources", "calculus", "autodiff_core_reference.edn"), "w") do io
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
    println(io, " :multi")
    println(io, " {")
    for (name, f, points) in multi
        g(v) = f(v...)
        print(io, "  :", name, " [")
        for x in points
            print(io, "\n   {:x ", fmt(x), " :gradient ", fmt(ForwardDiff.gradient(g, x)),
                  " :hessian ", fmt(ForwardDiff.hessian(g, x)), "}")
        end
        println(io, "]")
    end
    println(io, " }")
    println(io, " :sequence")
    println(io, " {")
    for (name, f) in sequence
        print(io, "  :", name, " [")
        for (xs, alpha) in sequence_points
            (name == "smooth-max-p-norm" && alpha < 1) && continue
            g(v) = f(v, alpha)
            print(io, "\n   {:xs ", fmt(xs), " :alpha ", repr(alpha), " :gradient ", fmt(ForwardDiff.gradient(g, xs)),
                  " :hessian ", fmt(ForwardDiff.hessian(g, xs)), "}")
        end
        println(io, "]")
    end
    println(io, " }}")
end
