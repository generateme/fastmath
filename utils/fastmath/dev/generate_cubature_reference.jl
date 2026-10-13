# Generator of reference values for `fastmath.calculus.cubature/cubature` (h-adaptive Genz-Malik cubature).
#
# Oracle: HCubature.jl, `hcubature(f, a, b; rtol, atol, initdiv, maxevals)`. Both implementations use the same
# degree 7/5 Genz-Malik rule, the same choice of the dimension to subdivide, the same max-error priority queue and the
# same stopping rule `E <= max(rtol * |I|, atol)`, so `result`, `error` and the number of integrand evaluations are
# expected to agree (value and error to rounding, evaluations exactly).
#
# Domain: smooth integrands (exp, sin, cos), a peak that needs many subdivisions, an integrand with a kink
# (abs), integrands which depend on one coordinate only through a tiny amplitude (the choice of the dimension to
# split decides the cost), an anisotropic integrand, a non-unit box, a 5-d integrand, `initdiv` 1..3 and dimensions 2..5.
# Integrands with discontinuities are not included: exact ties between the errors of boxes can be resolved
# differently by two queues.
#
# Produces `test/resources/calculus/cubature_reference.edn`: a vector of maps
#   {:id :lower :upper :rel :abs :initdiv :result :error :evaluations}
# The integrands are identified by `:id` and defined identically in `test/fastmath/calculus/cubature_test.clj`.
#
# Also produces `test/resources/calculus/cubature_axis_reference.edn`: the dimension (0-based) which one
# application of the Genz-Malik rule chooses to subdivide next (`HCubature.GenzMalik`, third result), for
# integrands `amp * sin(freq * x_k)` (small and large amplitudes, two frequencies), a constant integrand on boxes of
# different widths and a sum of two equal terms.
#
# Run from anywhere (deterministic; the output is byte-identical on rerun):
#     /home/ts/.juliaup/bin/julia /home/ts/clojure/fastmath/utils/fastmath/dev/generate_cubature_reference.jl
using HCubature, StaticArrays

const FS = Dict{String,Function}(
    "sincosexp3" => v -> sin(v[1]) * cos(v[2]) * exp(v[3]),
    "gauss2" => v -> exp(-(v[1]^2 + v[2]^2)),
    "peak2" => v -> 1 / (1e-2 + (v[1] - 0.3)^2 + (v[2] - 0.6)^2),
    "abs3" => v -> abs(v[1]^2 + v[2]^2 + v[3]^2 - 0.125),
    "smalldim3" => v -> 1e-3 * sin(9 * v[3]) + 1.0,
    "anis4" => v -> exp(10 * v[1]) * (1 + 0.01 * v[2] * v[3] * v[4]),
    "poly2" => v -> v[1]^2 * v[2],
    "sincos2" => v -> sin(v[1]) * cos(v[2]),
    "expsum5" => v -> exp(-(v[1]^2 + v[2]^2 + v[3]^2 + v[4]^2 + v[5]^2)),
)

# (id, lower, upper, rel, abs, initdiv)
const CASES = [
    ("poly2", [0.0, 0.0], [1.0, 1.0], 1.0e-7, 0.0, 1),
    ("poly2", [0.0, 0.0], [1.0, 1.0], 1.0e-7, 0.0, 3),
    ("gauss2", [0.0, 0.0], [1.0, 1.0], 1.0e-7, 0.0, 1),
    ("gauss2", [0.0, 0.0], [1.0, 1.0], 1.0e-7, 0.0, 2),
    ("gauss2", [0.0, 0.0], [1.0, 1.0], 1.0e-10, 0.0, 1),
    ("sincosexp3", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 1.0e-7, 0.0, 1),
    ("sincosexp3", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 1.0e-7, 0.0, 3),
    ("sincos2", [0.0, -1.0], [2.0, 3.0], 1.0e-9, 0.0, 1),
    ("peak2", [0.0, 0.0], [1.0, 1.0], 1.0e-6, 0.0, 1),
    ("peak2", [0.0, 0.0], [1.0, 1.0], 1.0e-4, 0.0, 2),
    ("abs3", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 1.0e-5, 0.0, 1),
    ("smalldim3", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 1.0e-9, 0.0, 1),
    ("anis4", [0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0], 1.0e-6, 0.0, 1),
    ("expsum5", [0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0, 1.0], 1.0e-5, 0.0, 1),
    ("gauss2", [0.0, 0.0], [1.0, 1.0], 1.0e-3, 1.0e-2, 1),
]

fmt(x::Float64) = isnan(x) ? "##NaN" : isinf(x) ? (x > 0 ? "##Inf" : "##-Inf") : repr(x)
fmtv(v) = string("[", join(fmt.(v), " "), "]")

function run_case(id, lower, upper, rel, abs, initdiv)
    f = FS[id]
    n = Ref(0)
    g = v -> (n[] += 1; f(v))
    I, E = hcubature(g, lower, upper; rtol=rel, atol=abs, initdiv=initdiv, maxevals=10^9)
    return I, E, n[]
end

const OUT = normpath(joinpath(@__DIR__, "..", "..", "..", "test", "resources", "calculus", "cubature_reference.edn"))
mkpath(dirname(OUT))
open(OUT, "w") do io
    println(io, "[")
    for (id, lower, upper, rel, abs, initdiv) in CASES
        I, E, n = run_case(id, lower, upper, rel, abs, initdiv)
        println(io, " {:id \"", id, "\" :lower ", fmtv(lower), " :upper ", fmtv(upper), " :rel ", fmt(rel), " :abs ", fmt(abs),
            " :initdiv ", initdiv, " :result ", fmt(I), " :error ", fmt(E), " :evaluations ", n, "}")
    end
    println(io, "]")
end
println("wrote ", OUT)

# (kind, n, lower, upper, k, amp, freq): :sin = amp * sin(freq * x_k), :const = amp, :sum = amp * (sin(freq * x_1) + sin(freq * x_2))
function axis_cases()
    cases = []
    for n in (2, 3, 4), k in 0:n-1, amp in (1.0e-3, 1.0, 1.0e3), freq in (4.0, 9.0)
        push!(cases, (:sin, n, zeros(n), ones(n), k, amp, freq))
    end
    push!(cases, (:const, 3, [0.0, 0.0, 0.0], [1.0, 3.0, 2.0], 0, 1.0, 0.0))
    push!(cases, (:const, 3, [0.0, 0.0, 0.0], [2.0, 1.0, 0.5], 0, 1.0, 0.0))
    push!(cases, (:const, 2, [0.0, 0.0], [1.0, 1.0], 0, 1.0, 0.0))
    push!(cases, (:sum, 3, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0, 1.0, 9.0))
    push!(cases, (:sum, 3, [0.0, 0.0, 0.0], [1.0, 2.0, 1.0], 0, 1.0, 9.0))
    return cases
end

function axis_fn(kind, k, amp, freq)
    kind == :sin && return v -> amp * sin(freq * v[k+1])
    kind == :const && return v -> amp
    return v -> amp * (sin(freq * v[1]) + sin(freq * v[2]))
end

const AOUT = normpath(joinpath(@__DIR__, "..", "..", "..", "test", "resources", "calculus", "cubature_axis_reference.edn"))
open(AOUT, "w") do io
    println(io, "[")
    for (kind, n, lower, upper, k, amp, freq) in axis_cases()
        g = HCubature.GenzMalik(Val(n))
        _, _, kd = g(axis_fn(kind, k, amp, freq), SVector{n}(lower), SVector{n}(upper))
        println(io, " {:kind :", kind, " :lower ", fmtv(lower), " :upper ", fmtv(upper), " :k ", k, " :amp ", fmt(amp),
            " :freq ", fmt(freq), " :axis ", kd - 1, "}")
    end
    println(io, "]")
end
println("wrote ", AOUT)
