(ns fastmath.core-examples
  (:refer-clojure
   :exclude [* + - / > < >= <= == rem quot mod bit-or bit-and bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd?])
  (:require [metadoc.examples :refer :all]
            [fastmath.core :refer :all]))

(add-examples atan2
  (example-session "atan2 values"
    (degrees (atan2 1 1))
    (degrees (atan2 0 1))
    (degrees (atan2 0 -1))))

(add-examples bessel-j
  (example-session "bessel-j values"
    (bessel-j 1 1)
    (bessel-j 0 1)
    (bessel-j 0 10)))

(add-examples between?
  (example-session "Examples"
    {:test-values [true false false true
                   true false false true]}
    (between? [1 4] 2)
    (between? [1 4] 5)
    (between? [1 4] -2)
    (between? [1 4] 1)
    (between? 1 4 2)
    (between? 1 4 5)
    (between? 1 4 -2)
    (between? 1 4 1)))

(add-examples co-intervals
  (example (co-intervals [1 2 3 1 2 3 4 5 6 1 2 3 4 5 6 7 1 1 1 1 -1]))
  (example "Higher overlap " (co-intervals [1 2 3 1 2 3 4 5 6 1 2 3 4 5 6 7 1 1 1 1 -1] 5 0.95)))

(add-examples group-by-intervals
  (example-session "Examples"
    (group-by-intervals [1 2 3 4 1 2 3 4 5 1 -2])
    (group-by-intervals [[1 2] [3 4]] [1 2 3 4 1 2 3 4 5 1 -2])))

(add-examples rem
  (example-session "Remainder (compared to `clojure.core` version)."
    (rem 3.123 0.333)
    (rem -3.123 0.333)
    (rem -3.123 -0.333)
    (rem 3.123 -0.333)
    (clojure.core/rem 3.123 0.333)
    (clojure.core/rem -3.123 0.333)
    (clojure.core/rem -3.123 -0.333)
    (clojure.core/rem 3.123 -0.333)))

(add-examples quot
  (example-session "Quotient (compared to `clojure.core` version)."
    (quot 3.123 0.333)
    (quot -3.123 0.333)
    (quot -3.123 -0.333)
    (quot 3.123 -0.333)
    (clojure.core/quot 3.123 0.333)
    (clojure.core/quot -3.123 0.333)
    (clojure.core/quot -3.123 -0.333)
    (clojure.core/quot 3.123 -0.333)))

(add-examples mod
  (example-session "Modulus (compared to `clojure.core` version)."
    (mod 3.123 0.333)
    (mod -3.123 0.333)
    (mod -3.123 -0.333)
    (mod 3.123 -0.333)
    (clojure.core/mod 3.123 0.333)
    (clojure.core/mod -3.123 0.333)
    (clojure.core/mod -3.123 -0.333)
    (clojure.core/mod 3.123 -0.333)))

(add-examples qsin
  (example-session "Compare [[sin]] and [[qsin]]." (sin 1.123) (qsin 1.123)))

(add-examples qcos
  (example-session "Compare [[cos]] and [[qcos]]." (cos 1.123) (qcos 1.123)))

(add-examples qexp
  (example-session "Compare [[exp]] and [[qexp]]." (exp 1.123) (qexp 1.123)))

(add-examples radians
  (example "Let's convert 180 degrees to radians." (radians 180))) 

(add-examples degrees
  (example "Let's convert \\\\(\\pi\\\\) radians to degrees." (degrees PI)))

(add-examples qlog
  (example-session "Compare [[log]] and [[qlog]]." (log 23.123) (qlog 23.123)))

(add-examples qpow
  (example-session "Compare [[pow]] and [[qpow]]." (pow 1.23 43.3) (qpow 1.23 43.3)))

(add-examples fpow
  (example-session "Example" (fpow 1.23 4) (fpow 1.23 4.123) (fpow 1.23 4.999) (fpow 1.23 5) (fpow 1.23 -2)))

(add-examples qsqrt
  (example-session "Compare [[sqrt]] and [[qsqrt]]." (sqrt 23.123) (qsqrt 23.123)))

(add-examples dist
  (example "Distance between two points." (dist 1 3 -2 10)))

(add-examples qdist
  (example "Distance between two points (quick version)." (qdist 1 3 -2 10))
  (example "Distance between two points (accurate version)." (dist 1 3 -2 10)))

(add-examples round
  (example "Round to long." (round PI)))

(add-examples rint
  (example "Round to double." (rint PI)))

(add-examples remainder
  (example-session "Compare with [[rem]] and [[mod]]."
    (remainder 3.123 0.2)
    (rem 3.123 0.2)
    (mod 3.123 0.2)))

(add-examples trunc
  (example-session "Examples" (trunc 1.234) (trunc -1.544)))

(add-examples itrunc
  (example-session "Examples" (itrunc 1.234) (itrunc -1.544)))

(add-examples approx
  (example "Default rounding (2 digits)." (approx 1.232323))
  (example "Rounding up to 4 digits." (approx 1.232323 4)))

(add-examples approx-eq
  (example "Default rounding (2 digits)." (approx-eq 1.232323 1.231999))
  (example "Rounding up to 4 digits." (approx-eq 1.232323 1.23231999 4))
  (example "Keep an eye on rounding" (approx-eq 1.2349 1.2350)))

(add-examples frac
  (example-session "Examples" (frac 0.555) (frac -0.555)))

(add-examples sfrac
  (example-session "Examples" (sfrac 0.555) (sfrac -0.555)))

(add-examples low-2-exp
  (example "Result 4 means, that \\\\(2^4=16\\\\) is lower than 23.11. Next exponent (5) gives greater value (32)." (low-2-exp 23.11))
  (example "For `x` less than 1.0 gives negative exponent." (low-2-exp 0.11)))

(add-examples high-2-exp
  (example "Result 5 means, that \\\\(2^5=32\\\\) is greater than 23.11. Lower exponent (4) gives lower value (16)." (high-2-exp 23.11))
  (example "For `x` less than 1.0 gives negative exponent." (high-2-exp 0.11)))

(add-examples low-exp
  (example "Result `1` means, that \\\\(9^1=9\\\\) is lower than `23.11`. Next exponent `2` gives greater value `82`." (low-exp 9 23.11))
  (example "For `x` less than `1.0` gives negative exponent." (low-exp 10 0.011)))

(add-examples high-exp
  (example "Result `2` means, that \\\\(9^2=81\\\\) is greater than `23.11`. Lower exponent `1` gives lower value `9`." (high-exp 9 23.11))
  (example "For `x` less than 1.0 gives negative exponent." (high-exp 10 0.011)))

(add-examples round-up-pow2
  (example-session "Examples" (round-up-pow2 1023) (round-up-pow2 1024) (round-up-pow2 1025)))

(add-examples next-double
  (example "Next double." (next-double 1234.56789))
  (example "Next double with delta." (next-double 1234.56789 1000)))

(add-examples prev-double
  (example "Prev double." (prev-double 1234.56789))
  (example "Prev double with delta." (prev-double 1234.56789 1000)))

(add-examples constrain
  (example-session "Examples" (constrain 0.5 1 2) (constrain 1.5 1 2) (constrain 2.5 1 2)))

(add-examples norm
  (example "Normalize from [1,-1] to [0,1]" (norm 0.234 -1.0 1.0))
  (example "Normalize from [-1,1] to [0,1]" (norm 0.234 1.0 -1.0))
  (example "Normalize cos() to [0,255]" (norm (cos HALF_PI) -1.0 1.0 0.0 255.0))
  (example "Normalize cos() to [255,0]" (norm (cos HALF_PI) -1.0 1.0 255.0 0.0)))

(add-examples mnorm
  (example "Normalize from [1,-1] to [0,1]" (mnorm 0.234 -1.0 1.0))
  (example "Normalize from [-1,1] to [0,1]" (mnorm 0.234 1.0 -1.0))
  (example "Normalize cos() to [0,255]" (mnorm (cos HALF_PI) -1.0 1.0 0.0 255.0))
  (example "Normalize cos() to [255,0]" (mnorm (cos HALF_PI) -1.0 1.0 255.0 0.0)))


(add-examples make-norm
  (example "Make cos() normalizer from [-1.0,1.0] to [0.0, 1.0]." (let [norm-cos (make-norm -1.0 1.0 0.0 1.0)]
                                                                    (norm-cos (cos 2.0))))
  (example "Make normalizer from [0,255] to any range." (let [norm-0-255 (make-norm 0 255)]
                                                          [(norm-0-255 123 -10 -20)
                                                           (norm-0-255 123 20 10)])))

(add-examples cnorm
  (example-session "Constrain result of norm." (cnorm 1.5 0 1 100 200) (cnorm 555 200 500)))

(add-examples lerp
  (example-session "Examples" (lerp 0.0 1.0 0.123) (lerp 0.0 100.0 0.123) (lerp 100 200 0.5))
  (example "Interpolate outside given range." (lerp -1.0 1.0 1.5)))

(add-examples mlerp
  (example-session "Examples" (mlerp 0.0 1.0 0.123) (mlerp 0.0 100.0 0.123) (mlerp 100 200 0.5))
  (example "Interpolate outside given range." (mlerp -1.0 1.0 1.5)))

(add-examples cos-interpolation
  (example "Example" (cos-interpolation 0.0 1.0 0.123)))

(add-examples smooth-interpolation
  (example "Example" (smooth-interpolation 0.0 1.0 0.123)))

(add-examples quad-interpolation
  (example "Example" (quad-interpolation 0.0 1.0 0.123)))

(add-examples smoothstep
  (example "x from range." (smoothstep 100 200 120))
  (example "corner case (< x edge0)" (smoothstep 100 200 50))
  (example "corner case (> x edge1)" (smoothstep 100 200 250)))

(add-examples wrap
  (example "Example 1" (wrap 0 -1 1))
  (example "Example 2 (value outside range)" (wrap -1.1 -1 1.5))
  (example "Example 3 (reversed range)" (wrap 0.7 0.5 1.0)))

(add-examples gcd
  (example-session "Usage"
    (gcd 340 440)
    (gcd (* 123 5544331) (* 123 123))
    (gcd -234 -432)))

(add-examples lcm
  (example-session "Usage"
    (lcm 340 440)
    (lcm (* 123 331) (* 123 123))
    (lcm 331 (* 123 123))
    (lcm -234 -432)))

(add-examples sample
  (example-session "Usage"
    (sample identity 10)
    (sample identity -11 22 5)
    (sample sq 1 5 5)
    (sample sq 1 5 5 true)))

(add-examples double-array->seq
  (example "Convert" (double-array->seq (seq->double-array [4 3 2]))))

(add-examples seq->double-array
  (example-session "Convert"
    (seq->double-array [1 2 3])
    (seq (seq->double-array [1 2 3]))
    (double-array->seq (seq->double-array [1 2 3])))
  (example "Also works on number (treated as one element list)." (seq (seq->double-array 1))))

(add-examples double-double-array->seq
  (example "Convert" (double-double-array->seq (seq->double-double-array [[4 3 2] (double-array [1 2 3])]))))

(add-examples seq->double-double-array
  (example-session "Convert"
    (seq->double-double-array [[1 2] [3 4]])
    (double-double-array->seq (seq->double-double-array [[1 2] [3 4]])))
  (example "Also works on seq of numbers" (seq (second (seq->double-double-array [1 2 3])))))

(add-examples haversine
  (example-session "Usage"
    (haversine 10)
    (haversine [10 20] [30 10])
    (haversine 10 20 30 10)))

(add-examples haversine-dist
  (example-session "Usage"
    (haversine-dist [10 20] [30 10])
    (haversine-dist 10 20 30 10)))

(add-examples inf?
  (example-session "Usage"
    {:test-values [true true false true false]}
    (inf? ##Inf)
    (inf? (/ 0))
    (inf? 1)
    (inf? ##-Inf)
    (inf? ##NaN)))

(add-examples neg-inf?
  (example-session "Usage"
    {:test-values [false false false true false]}
    (neg-inf? ##Inf)
    (neg-inf? (/ 0))
    (neg-inf? 1)
    (neg-inf? ##-Inf)
    (neg-inf? ##NaN)))

(add-examples pos-inf?
  (example-session "Usage"
    {:test-values [true true false false false]}
    (pos-inf? ##Inf)
    (pos-inf? (/ 0))
    (pos-inf? 1)
    (pos-inf? ##-Inf)
    (pos-inf? ##NaN)))

(add-examples invalid-double?
  (example-session "Usage"
    {:test-values [true true false true true]}
    (invalid-double? ##Inf)
    (invalid-double? (/ 0))
    (invalid-double? 1)
    (invalid-double? ##-Inf)
    (invalid-double? ##NaN)))

(add-examples valid-double?
  (example-session "Usage"
    {:test-values [false false true false false]}
    (valid-double? ##Inf)
    (valid-double? (/ 0))
    (valid-double? 1)
    (valid-double? ##-Inf)
    (valid-double? ##NaN)))

(add-examples nan?
  (example-session "Usage"
    {:test-values [false false false false true]}
    (nan? ##Inf)
    (nan? (/ 0))
    (nan? 1)
    (nan? ##-Inf)
    (nan? ##NaN)))

(add-examples order
  (example "Usage" {:test-value '(0 1 4 5 3 6 2 7 9 8)} (order [1 1 3 2 1 1 2 3 4 3]))
  (example "Reverse order" {:test-value '(8 2 7 9 3 6 0 1 4 5)} (order [1 1 3 2 1 1 2 3 4 3] true)))

(add-examples rank
  (example-session "Usage" {:test-values ['(3.5 0.5 5.0 0.5 7.0 10.0 2.0 9.0 7.0 3.5 7.0)
                                          '(3.5 0.5 5.0 0.5 7.0 10.0 2.0 9.0 7.0 3.5 7.0)
                                          [3 0 5 1 6 10 2 9 7 4 8]
                                          [4 1 5 0 8 10 2 9 7 3 6]
                                          '(3 0 5 0 6 10 2 9 6 3 6)
                                          '(4 1 5 1 8 10 2 9 8 4 8)
                                          '(2 0 3 0 4 6 1 5 4 2 4)]}
    (rank [3 1 4 1 5 9 2 6 5 3 5])
    (rank [3 1 4 1 5 9 2 6 5 3 5] :average)
    (rank [3 1 4 1 5 9 2 6 5 3 5] :first)
    (rank [3 1 4 1 5 9 2 6 5 3 5] :last)
    (rank [3 1 4 1 5 9 2 6 5 3 5] :min)
    (rank [3 1 4 1 5 9 2 6 5 3 5] :max)
    (rank [3 1 4 1 5 9 2 6 5 3 5] :dense))
  (example-session "Random ties"
    (rank [3 1 4 1 5 9 2 6 5 3 5] :random)
    (rank [3 1 4 1 5 9 2 6 5 3 5] :random)
    (rank [3 1 4 1 5 9 2 6 5 3 5] :random)
    (rank [3 1 4 1 5 9 2 6 5 3 5] :random)))

(def single-list `(sin cos tan cot sec csc asin acos atan acot asec acsc haversine
                     sinh cosh tanh coth sech csch asinh acosh atanh acoth asech acsch
                     qsin qcos exp log log10 ln log1p expm1 sqrt cbrt qexp qsqrt rqsqrt
                     erf erfc inv-erf inv-erfc sinc log2 qlog log1pexp
                     sq cb pow2 pow3 safe-sqrt floor ceil round rint abs iabs trunc itrunc
                     frac sfrac low-2-exp high-2-exp round-up-pow2
                     signum sgn sigmoid logit
                     gamma log-gamma digamma log-gamma-1p trigamma inv-gamma-1pm1))

(add-examples floor
  (example-session "Scaled floor" {:test-values [230.0 234.3]} (floor 234.312 10) (floor 234.312 0.1)))

(add-examples ceil
  (example-session "Scaled ceil" {:test-values [240.0 234.4]} (ceil 234.312 10) (ceil 234.312 0.1)))

(def double-list `(bessel-j atan2 hypot hypot-sqrt log-beta))

(def interp-list `(quad-interpolation smooth-interpolation wrap lerp cos-interpolation))

(def fn-list (concat single-list interp-list double-list))

(defmacro ^:private add-image-examples
  []
  `(do
     ~@(for [x fn-list]
         `(add-examples ~x
            (example-image ~(str "Plot of " (name x)) ~(str "images/m/" (name x) ".png"))))))

(defmacro ^:private add-call-examples
  []
  `(do
     ~@(for [x single-list]
         `(add-examples ~x
            (example (~x 1.0))))))

(add-call-examples)
(add-image-examples)

(add-examples erf
  (example-image "Plort of 2d erf" "images/m/erf2.png"))

