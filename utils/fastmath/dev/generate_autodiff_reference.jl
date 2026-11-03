# Reference values for fastmath.calculus.autodiff (forward-mode AD with tagged duals).
#
# Tool: Julia ForwardDiff.jl 1.3.3 (nested ForwardDiff.derivative for orders 2 and 3),
#       Float64 precision, run in the ~/.julia/environments/v1.12 project:
#
#   julia --project=$HOME/.julia/environments/v1.12 utils/fastmath/dev/generate_autodiff_reference.jl
#
# Output: test/resources/calculus/autodiff_reference.edn
# Map keyed by function id, each with parallel vectors :x, :d0 (value), :d1, :d2, :d3.
# The functions are written identically in test/fastmath/calculus/autodiff_test.clj.
# Domains avoid points where a value or a derivative is not finite; those are tested by hand.

using ForwardDiff
const D = ForwardDiff.derivative

d(f, n) = n == 0 ? f : d(x -> D(f, x), n - 1)

const all_x = [-1.5, -0.3, 0.0, 0.7, 2.0]
const pos_x = [0.3, 1.0, 2.5]

const functions = [
    ("sin",        x -> sin(x),                                  all_x),
    ("cos",        x -> cos(x),                                  all_x),
    ("exp",        x -> exp(x),                                  all_x),
    ("exp-sin",    x -> exp(x) * sin(x),                         all_x),
    ("poly",       x -> 3.0 * x^4 - 2.0 * x^3 + x - 5.0,         all_x),
    ("rational",   x -> (x^2 + 1.0) / (cos(x) + 2.0),            all_x),
    ("log",        x -> log(x),                                  pos_x),
    ("log-pow",    x -> log(x^2.5 + 1.0),                        pos_x),
    ("const-pow",  x -> 2.0^x,                                   all_x),
    ("pow-pow",    x -> x^x,                                     pos_x),
    ("pow-sin",    x -> (sin(x) + 2.0)^(cos(x)),                 all_x),
    ("reciprocal", x -> 1.0 / x - x / (x^2 + 3.0),               pos_x),
]

fmt(v) = "[" * join(map(repr, v), " ") * "]"

open(joinpath(@__DIR__, "..", "..", "..", "test", "resources", "calculus", "autodiff_reference.edn"), "w") do io
    println(io, "{")
    for (name, f, xs) in functions
        print(io, " :", name, " {:x ", fmt(xs))
        for n in 0:3
            print(io, " :d", n, " ", fmt([Float64(d(f, n)(x)) for x in xs]))
        end
        println(io, "}")
    end
    println(io, "}")
end
