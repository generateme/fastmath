# Generator of reference values for `fastmath.calculus/extrapolate` (Richardson extrapolation).
#
# Oracle: Richardson.jl (Julia), `extrapolate(f, h; x0, contract, power, atol, rtol, maxeval, breaktol)`.
# The Clojure function is a port of it; the one difference is the stopping rule: fastmath tests the error estimate of
# the previous iteration against the tolerance (one iteration later than Richardson.jl), so it may evaluate `g` once
# more and return a slightly different value. The two agree to the tolerance of the extrapolation, which is
# what the test asserts (absolute 1e-8, measured: at most 2.8e-13 in 15 of the 16 cases and 2.6e-9 in the
# noisy one).
#
# Domain: central and one-sided finite differences (power 2 and 1), a removable singularity (sin(x)/x), a limit
# at infinity of a sequence ((1 + h)^(1/h) -> e), an integrand with a square root, a function with
# a high frequency noise term, absolute and relative tolerances and a limit of the evaluations.
# Contraction 0.5, `breaktol` 2.0 (the fastmath defaults).
#
# Produces `test/resources/calculus/richardson_reference.edn`: a vector of maps
#   {:id :x :init-h :power :abs :rel :max-evals :value}
# The functions `g(x, h)` are identified by `:id` and defined identically in `test/fastmath/calculus/finite_test.clj`.
#
# Run from anywhere (deterministic; the output is byte-identical on rerun):
#     /home/ts/.juliaup/bin/julia /home/ts/clojure/fastmath/utils/fastmath/dev/generate_richardson_reference.jl
using Richardson

const GS = Dict{String,Function}(
    "central-sin" => (x, h) -> (sin(x + h) - sin(x - h)) / (2h),
    "forward-exp" => (x, h) -> (exp(x + h) - exp(x)) / h,
    "sinc" => (x, h) -> sin(x + h) / (x + h),
    "euler" => (x, h) -> (1 + h)^(1 / h),
    "central-atan" => (x, h) -> (atan(x + h) - atan(x - h)) / (2h),
    "central-sqrt" => (x, h) -> (sqrt(x + h) - sqrt(x - h)) / (2h),
    "noisy" => (x, h) -> 1.0 + h^2 + 1e-9 * sin(1 / h),
)

# (id, x, init-h, power, abs, rel, max-evals)
const CASES = [
    ("central-sin", 1.0, 0.5, 2.0, 0.0, 1.0e-10, typemax(Int32)),
    ("central-sin", 1.0, 0.5, 2.0, 0.0, 1.0e-14, typemax(Int32)),
    ("central-sin", 1.0, 0.5, 2.0, 1.0e-12, 0.0, typemax(Int32)),
    ("central-sin", 1.0, 0.5, 2.0, 0.0, 1.0e-10, 3),
    ("central-sin", 1.0, 0.5, 2.0, 0.0, 1.0e-10, 4),
    ("central-sin", 1.0, 1.0, 2.0, 0.0, 1.0e-10, 6),
    ("forward-exp", 0.3, 0.1, 1.0, 0.0, 1.0e-10, typemax(Int32)),
    ("forward-exp", 0.3, 0.1, 1.0, 1.0e-13, 0.0, typemax(Int32)),
    ("forward-exp", 0.3, 0.5, 1.0, 0.0, 1.0e-8, 5),
    ("sinc", 0.0, 1.0, 2.0, 0.0, 1.0e-10, typemax(Int32)),
    ("sinc", 0.0, 1.0, 2.0, 0.0, 1.0e-14, typemax(Int32)),
    ("euler", 0.0, 0.1, 1.0, 0.0, 1.0e-10, typemax(Int32)),
    ("central-atan", 2.0, 0.5, 2.0, 0.0, 1.0e-12, typemax(Int32)),
    ("central-sqrt", 1.0, 0.5, 2.0, 0.0, 1.0e-12, typemax(Int32)),
    ("noisy", 0.0, 0.5, 2.0, 0.0, 1.0e-12, typemax(Int32)),
    ("noisy", 0.0, 0.5, 2.0, 1.0e-7, 0.0, typemax(Int32)),
]

function run_case(id, x, h0, power, atol, rtol, maxeval)
    g = GS[id]
    v, _ = extrapolate(h -> g(x, h), h0; x0=0.0, contract=0.5, power=power, atol=atol, rtol=rtol,
        maxeval=maxeval, breaktol=2.0)
    return v
end

const OUT = normpath(joinpath(@__DIR__, "..", "..", "..", "test", "resources", "calculus", "richardson_reference.edn"))
mkpath(dirname(OUT))
open(OUT, "w") do io
    println(io, "[")
    for (id, x, h0, power, atol, rtol, maxeval) in CASES
        v = run_case(id, x, h0, power, atol, rtol, maxeval)
        println(io, " {:id \"", id, "\" :x ", repr(x), " :init-h ", repr(h0), " :power ", repr(power), " :abs ", repr(atol),
            " :rel ", repr(rtol), " :max-evals ", maxeval, " :value ", repr(v), "}")
    end
    println(io, "]")
end
println("wrote ", OUT)
