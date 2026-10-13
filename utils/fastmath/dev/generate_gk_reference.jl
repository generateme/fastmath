# Generator of reference values for `fastmath.calculus.quadrature/gk-quadrature` (and `fastmath.calculus/integrate`
# with `:integrator :gauss-kronrod`).
#
# Oracle: QuadGK.jl (Julia), `quadgk(f, segments...; rtol, atol, order, maxevals)`. Both implementations are
# adaptive Gauss-Kronrod with a max-error priority queue, the same rule (`order` Gauss points, 2*order+1 Kronrod
# points) and the same substitutions for infinite bounds, so `result`, `error` and the number of integrand
# evaluations are expected to agree (value and error to rounding, evaluations exactly).
#
# Domain: integrands with endpoint singularities, kinks, oscillations, many zero-error segments (polynomials),
# a narrow peak (needs `atol = 0`), improper intervals, orders 1..21, user break points and equal initial pieces
# (`initdiv`, emulated here with equally spaced break points). Cases never converge by an exact tie of
# segment errors that could be resolved differently by two queues.
#
# Produces `test/resources/calculus/gk_reference.edn`: a vector of maps
#   {:id :lower :upper :rel :abs :order :points :initdiv :result :error :evaluations}
# The integrands are identified by `:id` and defined identically in `test/fastmath/calculus/quadrature_test.clj`.
#
# Also produces `test/resources/calculus/kronrod_reference.edn`: `QuadGK.kronrod(n)` (nodes `:x`, Kronrod weights
# `:w`, embedded Gauss weights `:gw` on the half interval) for selected orders 1..50.
#
# Run from anywhere (deterministic; the output is byte-identical on rerun):
#     /home/ts/.juliaup/bin/julia /home/ts/clojure/fastmath/utils/fastmath/dev/generate_gk_reference.jl
using QuadGK

const FS = Dict{String,Function}(
    "poly" => x -> x^2,
    "const" => x -> 1.0,
    "sqrtabs" => x -> sqrt(abs(x)),
    "abssin20" => x -> abs(sin(20x)),
    "invsqrt" => x -> 1 / sqrt(x),
    "logx" => x -> log(x),
    "expinf" => x -> exp(-x),
    "halfexp" => x -> exp(-x) / sqrt(x),
    "cauchy" => x -> 1 / (1 + x^2),
    "osc" => x -> sin(50x),
    "runge" => x -> 1 / (1 + 25x^2),
    "gauss" => x -> exp(-x^2),
    "peak" => x -> exp(-1.0e6 * x^2),
)

# (id, lower, upper, rel, abs, order, points, initdiv)
const CASES = [
    ("poly", -1.0, 1.0, 1.0e-10, 0.0, 7, Float64[], 1),
    ("poly", 0.0, 3.0, 1.0e-10, 0.0, 7, [0.5, 1.0, 2.0], 1),
    ("const", 0.0, 1.0, 1.0e-10, 0.0, 7, [0.25, 0.5, 0.75], 1),
    ("sqrtabs", -1.0, 1.0, 1.0e-10, 0.0, 7, Float64[], 1),
    ("sqrtabs", -1.0, 1.0, 1.0e-10, 0.0, 7, [0.0], 1),
    ("sqrtabs", -1.0, 1.0, 1.0e-8, 0.0, 3, Float64[], 1),
    ("abssin20", -1.0, 1.0, 1.0e-8, 0.0, 7, Float64[], 1),
    ("abssin20", -1.0, 1.0, 1.0e-8, 0.0, 15, Float64[], 1),
    ("abssin20", -1.0, 1.0, 1.0e-8, 0.0, 7, Float64[], 5),
    ("invsqrt", 0.0, 1.0, 1.0e-10, 0.0, 7, Float64[], 1),
    ("invsqrt", 0.0, 1.0, 1.0e-10, 0.0, 3, Float64[], 1),
    ("invsqrt", 0.0, 1.0, 1.0e-10, 0.0, 10, Float64[], 1),
    ("logx", 0.0, 1.0, 1.0e-10, 0.0, 7, Float64[], 1),
    ("expinf", 0.0, Inf, 1.0e-10, 0.0, 7, Float64[], 1),
    ("expinf", 0.0, Inf, 1.0e-10, 0.0, 21, Float64[], 1),
    ("expinf", 0.0, Inf, 1.0e-6, 1.0e-12, 7, Float64[], 1),
    ("halfexp", 0.0, Inf, 1.0e-10, 0.0, 7, Float64[], 1),
    ("cauchy", -Inf, Inf, 1.0e-10, 0.0, 7, Float64[], 1),
    ("cauchy", -Inf, Inf, 1.0e-10, 0.0, 15, Float64[], 1),
    ("cauchy", -Inf, 0.0, 1.0e-10, 0.0, 7, Float64[], 1),
    ("osc", 0.0, 1.0, 1.0e-10, 0.0, 7, Float64[], 1),
    ("osc", 0.0, 1.0, 1.0e-10, 0.0, 15, Float64[], 1),
    ("osc", 0.0, 1.0, 1.0e-10, 0.0, 4, Float64[], 1),
    ("runge", -1.0, 1.0, 1.0e-12, 0.0, 7, Float64[], 1),
    ("gauss", -5.0, 5.0, 1.0e-12, 0.0, 7, Float64[], 1),
    ("gauss", -5.0, 5.0, 1.0e-12, 0.0, 3, Float64[], 1),
    ("gauss", -5.0, 5.0, 1.0e-12, 0.0, 21, Float64[], 1),
    ("gauss", -5.0, 5.0, 1.0e-12, 0.0, 7, Float64[], 4),
    ("peak", -1.0, 1.0, 1.0e-8, 0.0, 7, Float64[], 1),
    ("peak", -1.0, 1.0, 1.0e-8, 0.0, 7, [0.0], 1),
]

fmt(x::Float64) = isnan(x) ? "##NaN" : isinf(x) ? (x > 0 ? "##Inf" : "##-Inf") : repr(x)

function run_case(id, a, b, rel, abs, order, points, initdiv)
    f = FS[id]
    n = Ref(0)
    g = x -> (n[] += 1; f(x))
    pts = isempty(points) && initdiv > 1 ? [a + (b - a) * k / initdiv for k in 1:initdiv-1] : points
    I, E = quadgk(g, a, pts..., b; rtol=rel, atol=abs, order=order, maxevals=10^9)
    return I, E, n[]
end

const OUT = normpath(joinpath(@__DIR__, "..", "..", "..", "test", "resources", "calculus", "gk_reference.edn"))
mkpath(dirname(OUT))
open(OUT, "w") do io
    println(io, "[")
    for (id, a, b, rel, abs, order, points, initdiv) in CASES
        I, E, n = run_case(id, a, b, rel, abs, order, points, initdiv)
        println(io, " {:id \"", id, "\" :lower ", fmt(a), " :upper ", fmt(b), " :rel ", fmt(rel), " :abs ", fmt(abs),
            " :order ", order, " :points [", join(fmt.(points), " "), "] :initdiv ", initdiv,
            " :result ", fmt(I), " :error ", fmt(E), " :evaluations ", n, "}")
    end
    println(io, "]")
end
println("wrote ", OUT)

const KRONROD_ORDERS = [1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 16, 20, 21, 25, 30, 40, 50]
const KOUT = normpath(joinpath(@__DIR__, "..", "..", "..", "test", "resources", "calculus", "kronrod_reference.edn"))
open(KOUT, "w") do io
    println(io, "[")
    for n in KRONROD_ORDERS
        x, w, gw = QuadGK.kronrod(n)
        println(io, " {:order ", n, " :x [", join(fmt.(x), " "), "] :w [", join(fmt.(w), " "), "] :gw [", join(fmt.(gw), " "), "]}")
    end
    println(io, "]")
end
println("wrote ", KOUT)
