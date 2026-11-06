# Reference values for the gradient, jacobian and hessian of fastmath.calculus.autodiff.
#
# Tool: Julia ForwardDiff.jl 1.3.3 (ForwardDiff.gradient, jacobian, hessian with default chunking;
#       the Jacobian-vector product as derivative(t -> f(x + t v), 0) and the Hessian-vector product as
#       derivative(t -> gradient(f, x + t v), 0), ForwardDiff has no jvp/hvp functions),
#       Float64 precision, run in the ~/.julia/environments/v1.12 project:
#
#   julia --project=$HOME/.julia/environments/v1.12 utils/fastmath/dev/generate_autodiff_multivariate_reference.jl
#
# Output: test/resources/calculus/autodiff_multivariate_reference.edn
# {:scalar {<id> [{:x [..] :value v :gradient [..] :hessian [[..] ..] :v [..] :jvp p :hvp [..]} ..]}
#  :vector {<id> [{:x [..] :value [..] :jacobian [[..] ..] :v [..] :jvp [..]} ..]}
#  :univariate-vector {<id> [{:x x :value [..] :d1 [..] :d2 [..] :d3 [..]} ..]}
#  :reduction {<sum|prod> [{:x [..] :value v :gradient [..] :hessian [[..] ..]} ..]}
#  :sum-negative-zero v}
# The univariate vector functions R -> R^m are differentiated with nested ForwardDiff.derivative (orders 1-3);
# the reductions are Base sum and prod, bit-identical to a left fold below 16 elements (checked on random vectors,
# Float64 and duals; from 16 elements Julia reorders), all points here are shorter; :sum-negative-zero is sum([-0.0]).
# The functions are written identically in test/fastmath/calculus/autodiff_test.clj (powers as Float64 exponents).

using ForwardDiff
using SpecialFunctions

# chained Rosenbrock: sparse, each term depends on two inputs
rosen(x) = sum((100.0 * (x[i+1] - x[i]^2.0)^2.0 + (1.0 - x[i])^2.0 for i in 1:length(x)-1); init = 0.0)

# dense: every term depends on all inputs
function dense(x)
    s = sum(x; init = 0.0)
    q = sum((a * a for a in x); init = 0.0)
    sum((sin(a * s) * exp(q / 100.0) for a in x); init = 0.0)
end

# 3 -> 2 and 2 -> 3 vector functions
vec32(x) = [x[1] * x[2] * sin(x[3]), exp(x[1]) / (1.0 + x[2]^2.0 + x[3]^2.0)]
vec23(x) = [x[1] / x[2], log(x[1] + x[2]^2.0), x[1]^x[2]]

# R -> R^m vector functions: trigonometric (the ForwardDiff.jl docs example), constants mixed with duals, special
uvtrig(y) = [sin(y), cos(y), tan(y)]
uvmixed(y) = [1.0, y^2.0, exp(y) * y]
uvspecial(y) = [gamma(y), log(y), sqrt(y)]
uvrational(y) = [y / (1.0 + y^2.0), atan(y)]

point(n) = [0.3 * i - 0.8 for i in 1:n]
direction(n) = [0.7 - 0.25 * i for i in 1:n]

jvp(f, x, v) = ForwardDiff.derivative(t -> f(x .+ t .* v), 0.0)
hvp(f, x, v) = ForwardDiff.derivative(t -> ForwardDiff.gradient(f, x .+ t .* v), 0.0)

const scalars = [("rosen", rosen, [point(n) for n in (1, 2, 5, 10)]),
                 ("dense", dense, [point(n) for n in (1, 2, 5, 10)])]
const vectors = [("vec32", vec32, [[0.5, -1.2, 2.0], [1.5, 0.3, -0.7]]),
                 ("vec23", vec23, [[0.5, 1.2], [2.0, 0.7]])]

const univariate_vectors = [("uvtrig", uvtrig), ("uvmixed", uvmixed), ("uvspecial", uvspecial), ("uvrational", uvrational)]
const univariate_points = [0.5, 1.3, 2.7]
const reductions = [("sum", sum), ("prod", prod)]
const reduction_points = [[1.0, 2.0, 3.0], [0.0, 2.0, 3.0], [0.0, 0.0, 3.0], [0.5, -1.2, 2.0, 0.7]]

d1(f, y) = ForwardDiff.derivative(f, y)
d2(f, y) = ForwardDiff.derivative(t -> d1(f, t), y)
d3(f, y) = ForwardDiff.derivative(t -> d2(f, t), y)

fmt(v::AbstractVector) = "[" * join(map(repr, v), " ") * "]"
fmt(m::AbstractMatrix) = "[" * join([fmt(m[i, :]) for i in 1:size(m, 1)], " ") * "]"

open(joinpath(@__DIR__, "..", "..", "..", "test", "resources", "calculus", "autodiff_multivariate_reference.edn"), "w") do io
    println(io, "{:scalar")
    println(io, " {")
    for (name, f, xs) in scalars
        print(io, "  :", name, " [")
        for x in xs
            print(io, "\n   {:x ", fmt(x), " :value ", repr(f(x)),
                  " :gradient ", fmt(ForwardDiff.gradient(f, x)),
                  " :hessian ", fmt(ForwardDiff.hessian(f, x)),
                  " :v ", fmt(direction(length(x))),
                  " :jvp ", repr(jvp(f, x, direction(length(x)))),
                  " :hvp ", fmt(hvp(f, x, direction(length(x)))), "}")
        end
        println(io, "]")
    end
    println(io, " }")
    println(io, " :vector")
    println(io, " {")
    for (name, f, xs) in vectors
        print(io, "  :", name, " [")
        for x in xs
            print(io, "\n   {:x ", fmt(x), " :value ", fmt(f(x)),
                  " :jacobian ", fmt(ForwardDiff.jacobian(f, x)),
                  " :v ", fmt(direction(length(x))),
                  " :jvp ", fmt(jvp(f, x, direction(length(x)))), "}")
        end
        println(io, "]")
    end
    println(io, " }")
    println(io, " :univariate-vector")
    println(io, " {")
    for (name, f) in univariate_vectors
        print(io, "  :", name, " [")
        for y in univariate_points
            print(io, "\n   {:x ", repr(y), " :value ", fmt(f(y)),
                  " :d1 ", fmt(d1(f, y)), " :d2 ", fmt(d2(f, y)), " :d3 ", fmt(d3(f, y)), "}")
        end
        println(io, "]")
    end
    println(io, " }")
    println(io, " :reduction")
    println(io, " {")
    for (name, f) in reductions
        print(io, "  :", name, " [")
        for x in reduction_points
            print(io, "\n   {:x ", fmt(x), " :value ", repr(f(x)),
                  " :gradient ", fmt(ForwardDiff.gradient(f, x)),
                  " :hessian ", fmt(ForwardDiff.hessian(f, x)), "}")
        end
        println(io, "]")
    end
    println(io, " }")
    println(io, " :sum-negative-zero ", repr(sum([-0.0])), "}")
end
