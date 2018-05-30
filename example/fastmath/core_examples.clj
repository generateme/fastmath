(ns fastmath.core-examples
  (:refer-clojure
   :exclude [* + - / > < >= <= == rem quot mod bit-or bit-and bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd?])
  (:require [metadoc.examples :refer :all]
            [fastmath.core :refer :all]))


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

(add-examples sample
  (example-session "Usage" (sample identity 10) (sample identity -11 22 5) (sample sq 1 5 5)))

(add-examples double-array->seq
  (example "Convert" (double-array->seq (seq->double-array [4 3 2]))))

(add-examples seq->double-array
  (example-session "Convert"
    (seq->double-array [1 2 3])
    (seq (seq->double-array [1 2 3]))
    (double-array->seq (seq->double-array [1 2 3]))))

(add-examples double-double-array->seq
  (example "Convert" (double-double-array->seq (seq->double-double-array [[4 3 2] (double-array [1 2 3])]))))

(add-examples seq->double-double-array
  (example-session "Convert"
    (seq->double-double-array [[1 2] [3 4]])
    (double-double-array->seq (seq->double-double-array [[1 2] [3 4]]))))


(def single-list `(sin cos tan cot sec csc asin acos atan acot asec acsc
                       sinh cosh tanh coth sech csch asinh acosh atanh acoth asech acsch
                       qsin qcos exp log log10 ln log1p sqrt cbrt qexp qsqrt rqsqrt
                       erf erfc inv-erf inv-erfc sinc log2 qlog
                       sq pow2 pow3 safe-sqrt floor ceil round rint abs iabs trunc
                       frac sfrac low-2-exp high-2-exp round-up-pow2 next-double prev-double
                       signum sgn sigmoid
                       gamma log-gamma digamma log-gamma-1p trigamma inv-gamma-1pm1))

(def interp-list `(quad-interpolation smooth-interpolation wrap lerp cos-interpolation))

(def fn-list (concat single-list interp-list))

(defmacro ^:private add-image-examples
  []
  `(do
     ~@(for [x fn-list]
         `(add-examples ~x
            (example-image ~(str "Plot of " (name x)) ~(str "images/m/" (name x) ".png"))))))

(add-image-examples)
