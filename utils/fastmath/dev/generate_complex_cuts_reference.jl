# Reference values of fastmath.complex on the branch cuts and the axes, with signed zeros
# (test/resources/complex/complex_cuts_reference.edn).
#
# Tool: Julia 1.13.3, ComplexF64 (/home/ts/.juliaup/bin/julia), C99 conventions: the sign of a zero part selects
# the side of a branch cut. mpmath has no signed zero, so this file complements
# generate_complex_reference.py, which covers the points off the axes.
# Domain: a and b from a list of values including -0.0 and 0.0, with a == 0 or b == 0 and not both
# (the pole at 0 is excluded: Julia returns NaN where fastmath returns the infinities of the C99 limits).
# Output: one map, function => vector of [a b re im]. Infinities are written as ##Inf.
# Run: julia generate_complex_cuts_reference.jl > complex_cuts_reference.edn
# Regeneration is byte-identical (fixed order, repr of doubles).

fs = [:sqrt, :exp, :log, :sin, :cos, :tan, :sec, :csc, :cot, :sinh, :cosh, :tanh, :sech, :csch, :coth,
      :asin, :acos, :atan, :asec, :acsc, :acot, :asinh, :acosh, :atanh, :asech, :acsch, :acoth]
xs = [-100.0, -5.0, -2.0, -1.0, -0.5, -1e-10, -0.0, 0.0, 1e-10, 0.5, 1.0, 2.0, 5.0, 100.0]

function fmt(x)
    isnan(x) && return "##NaN"
    isinf(x) && return x > 0 ? "##Inf" : "##-Inf"
    return repr(x)
end

println("{")
for f in fs
    print(":", f, " [")
    first = true
    for a in xs, b in xs
        (a == 0.0 || b == 0.0) || continue
        (a == 0.0 && b == 0.0) && continue
        z = ComplexF64(a, b)
        w = @eval($f($z))
        first || print(" ")
        first = false
        print("[", fmt(a), " ", fmt(b), " ", fmt(real(w)), " ", fmt(imag(w)), "]")
    end
    println("]")
end
println("}")
