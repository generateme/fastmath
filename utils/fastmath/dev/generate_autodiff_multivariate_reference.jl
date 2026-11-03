# Reference values for the gradient, jacobian and hessian of fastmath.calculus.autodiff.
#
# Tool: Julia ForwardDiff.jl 1.3.3 (ForwardDiff.gradient, jacobian, hessian with default chunking),
#       Float64 precision, run in the ~/.julia/environments/v1.12 project:
#
#   julia --project=$HOME/.julia/environments/v1.12 utils/fastmath/dev/generate_autodiff_multivariate_reference.jl
#
# Output: test/resources/calculus/autodiff_multivariate_reference.edn
# {:scalar {<id> [{:x [..] :value v :gradient [..] :hessian [[..] ..]} ..]}
#  :vector {<id> [{:x [..] :value [..] :jacobian [[..] ..]} ..]}}
# The functions are written identically in test/fastmath/calculus/autodiff_test.clj (powers as Float64 exponents).

using ForwardDiff

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

point(n) = [0.3 * i - 0.8 for i in 1:n]

const scalars = [("rosen", rosen, [point(n) for n in (1, 2, 5, 10)]),
                 ("dense", dense, [point(n) for n in (1, 2, 5, 10)])]
const vectors = [("vec32", vec32, [[0.5, -1.2, 2.0], [1.5, 0.3, -0.7]]),
                 ("vec23", vec23, [[0.5, 1.2], [2.0, 0.7]])]

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
                  " :hessian ", fmt(ForwardDiff.hessian(f, x)), "}")
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
                  " :jacobian ", fmt(ForwardDiff.jacobian(f, x)), "}")
        end
        println(io, "]")
    end
    println(io, " }}")
end
