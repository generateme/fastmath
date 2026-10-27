# Reference values for fastmath.quaternion (test/resources/quaternion/quaternion_reference.edn).
#
# Tool: Julia 1.13.3 with Quaternions.jl 0.7.7 (/home/ts/.juliaup/bin/julia), Float64.
# Domain: quaternions a+bi+cj+dk with a in {-2, -0.5, 0.3, 1.5}, b, c, d in {-1.2, 0, 0.7, 2} and a non zero
# vector part (the choice of the axis for a real argument is a convention, tested separately), a few points
# with a tiny and with a huge vector part, and a few points with a small scalar part.
# Output: one map, :points [[a b c d] ...] and for every function a vector of [a b c d] (nil when
# Quaternions.jl has no method); :mult-args / :mult and :div, :pow (quaternion exponent) for pairs.
# Run: julia generate_quaternion_reference.jl > quaternion_reference.edn
# Regeneration is byte-identical (fixed order, repr of doubles).
using Quaternions

as = [-2.0, -0.5, 0.3, 1.5]
vs = [-1.2, 0.0, 0.7, 2.0]
points = [(a, b, c, d) for a in as for b in vs for c in vs for d in vs if !(b == 0.0 && c == 0.0 && d == 0.0)]
append!(points, [(0.5, 1e-9, 2e-9, -1e-9), (1.0, 1e5, -2e5, 3e5), (0.0, 1.0, 2.0, 3.0), (1e-9, 0.3, -0.4, 0.5), (3.0, 1e-4, 0.0, 0.0)])
fns = [:exp, :log, :sqrt, :sin, :cos, :tan, :sinh, :cosh, :tanh, :asin, :acos, :atan, :asinh, :acosh, :atanh,
       :sec, :csc, :cot, :sech, :csch, :coth, :inv, :conj]

fmt(x) = isnan(x) ? "##NaN" : isinf(x) ? (x > 0 ? "##Inf" : "##-Inf") : repr(x)
row(q) = "[" * join(fmt.([q.s, q.v1, q.v2, q.v3]), " ") * "]"

println("{")
println(":points [", join(["[" * join(fmt.(p), " ") * "]" for p in points], " "), "]")
for f in fns
    vals = String[]
    for p in points
        q = Quaternion(p...)
        r = try
            w = @eval($f($q))
            row(w)
        catch
            "nil"
        end
        push!(vals, r)
    end
    println(":", f, " [", join(vals, " "), "]")
end

pairs = [(Quaternion(a...), Quaternion(b...)) for a in points[1:7:end] for b in points[3:37:end]]
println(":pair-args [", join(["[[" * join(fmt.([a.s, a.v1, a.v2, a.v3]), " ") * "] [" * join(fmt.([b.s, b.v1, b.v2, b.v3]), " ") * "]]" for (a, b) in pairs], " "), "]")
println(":mult [", join([row(a * b) for (a, b) in pairs], " "), "]")
println(":div [", join([row(a / b) for (a, b) in pairs], " "), "]")
println(":pow [", join([let r = try row(exp(log(a) * b)) catch; "nil" end; r end for (a, b) in pairs], " "), "]")
println("}")
