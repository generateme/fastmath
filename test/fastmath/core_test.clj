(ns fastmath.core-test
  (:require [fastmath.core :as m]
            [clojure.test :as t]))

;; primitive ops

(let [numbers (repeatedly 10000 #(m/* 100.0 (m/- (double (rand)) 0.5)))]
  (t/deftest add
    (t/is (int? (m/+ 1 2 3)))
    (t/is (double? (m/+ 1 2 3.0)))
    (t/is (m/== 6 (m/+ 1 2 3)))
    (t/is (m/== 6.0 (m/+ 1 2 3.0)))
    (t/is (m/== 6 (m/long-add 1 2 3)))
    (t/is (m/== 6 (m/long-add 1 2 3.0)))
    (t/is (m/== (reduce + numbers)
                (reduce m/+ numbers)))
    (t/is (m/== (apply + numbers)
                (apply m/+ numbers)))
    (t/is (m/== (reduce + (map long numbers))
                (reduce m/long-add numbers)))))

(let [numbers (repeatedly 10000 #(m/* 100.0 (m/- (double (rand)) 0.5)))]
  (t/deftest sub
    (t/is (int? (m/- 1)))
    (t/is (double? (m/- 1.0)))
    (t/is (m/== -1.0 (m/- 1.0)))
    (t/is (m/== -1 (m/- 1)))
    (t/is (int? (m/- 1 2 3)))
    (t/is (double? (m/- 1 2 3.0)))
    (t/is (m/== -4 (m/- 1 2 3)))
    (t/is (m/== -4.0 (m/- 1 2 3.0)))
    (t/is (m/== -4 (m/long-sub 1 2 3)))
    (t/is (m/== -4 (m/long-sub 1 2 3.0)))
    (t/is (m/== (reduce - numbers)
                (reduce m/- numbers)))
    (t/is (m/== (apply - numbers)
                (apply m/- numbers)))
    (t/is (m/== (reduce - (map long numbers))
                (reduce m/long-sub numbers)))))

;; Reference: * and / are the independent reference
;; implementation for the double-typed n-ary arithmetic ops (exact arithmetic,
;; no floating-point-sensitive reduction order difference at these arities).
(let [numbers (repeatedly 10000 #(m/inc (double (rand 100.0))))]
  (t/deftest mult
    (t/is (int? (m/* 1 2 3)))
    (t/is (double? (m/* 1 2 3.0)))
    (t/is (m/== 6 (m/* 1 2 3)))
    (t/is (m/== 24 (m/* 1 2 3 4)))
    (t/is (m/== 24 (apply m/* [1 2 3 4])) "apply bypasses inlining, exercises the actual fn body")
    (t/is (m/== (reduce * numbers)
                (reduce m/* numbers)))
    (t/is (m/== (apply * numbers)
                (apply m/* numbers)))))

;; Reference: hand-computed. Regression test for a confirmed bug where the
;; non-inlined 3-arity body computed (a*b)+c instead of (a*b)*c; only visible
;; via `apply`/`reduce` (the compiler inlines direct calls, masking the bug).
(t/deftest long-mult
  (t/is (m/== 24 (m/long-mult 2 3 4)) "inlined path (direct call)")
  (t/is (m/== 24 (apply m/long-mult [2 3 4])) "non-inlined path (apply) -- was 10 before fix")
  (t/is (m/== 120 (apply m/long-mult [2 3 4 5])))
  (t/is (m/== 1 (m/long-mult)))
  (t/is (m/== 7 (m/long-mult 7)))
  (t/is (m/== (reduce * (map long (range 1 20)))
              (reduce m/long-mult (map long (range 1 20))))))

(let [numbers (repeatedly 10000 #(m/inc (double (rand 100.0))))]
  (t/deftest div
    (t/is (m/== 0.5 (m// 2)))
    (t/is (m/== 5 (m// 100 5 4)))
    (t/is (m/== 5 (apply m// [100 5 4])))
    (t/is (m/== ##Inf (m// 1.0 0.0)) "IEEE 754 semantics, no exception")
    (t/is (m/== (reduce / numbers)
                (reduce m// numbers)))))

;; Reference: hand-computed, per the 1-arg-returns-double exception documented
;; in long-div's docstring (PrimitiveMath has no long-returning reciprocal).
(t/deftest long-div
  (t/is (double? (m/long-div 5)))
  (t/is (m/delta-eq 0.2 (m/long-div 5)))
  (t/is (m/== 5 (m/long-div 100 5 4)))
  (t/is (m/== 5 (apply m/long-div [100 5 4])))
  (t/is (thrown? ArithmeticException (m/long-div 1 0))))

(let [numbers (repeatedly 1000 #(rand 100.0))]
  (t/deftest min-max
    (t/is (m/== (reduce min numbers) (reduce m/min numbers)))
    (t/is (m/== (reduce max numbers) (reduce m/max numbers)))
    (t/is (m/== (apply min numbers) (apply m/min numbers)))
    (t/is (m/== (apply max numbers) (apply m/max numbers)))
    (t/is (m/== 1 (m/min 5 3 1 4)))
    (t/is (m/== 5 (m/max 5 3 1 4)))
    (t/is (m/== 1 (apply m/long-min [5 3 1 4])))
    (t/is (m/== 5 (apply m/long-max [5 3 1 4])))))

;; Reference: ==, </>/<=/>=, hand-computed. Regression tests for a
;; confirmed bug where the non-inlined 3+-arity body of ==, <, >, <=, >= never
;; compared the 2nd argument against the 3rd (only a-vs-b and then a chain
;; starting at the 3rd arg); only visible via `apply`/`reduce` (the compiler
;; inlines direct calls, masking the bug). not== was unaffected (different,
;; correct implementation via a single sorted reduce over all arguments).
(t/deftest chained-comparison-non-inlined-path
  (t/is (false? (apply m/== [1 1 5])) "was true before fix")
  (t/is (false? (apply m/== [1 1 5 5])) "was true before fix")
  (t/is (true? (apply m/== [3 3 3 3])))
  (t/is (false? (apply m/< [1 2 0])) "was true before fix")
  (t/is (true? (apply m/< [1 2 3 4])))
  (t/is (false? (apply m/> [3 2 5])) "was true before fix")
  (t/is (true? (apply m/> [5 4 3 2])))
  (t/is (false? (apply m/<= [1 2 0])) "was true before fix")
  (t/is (true? (apply m/<= [1 2 2 3])))
  (t/is (false? (apply m/>= [3 2 5])) "was true before fix")
  (t/is (true? (apply m/>= [5 5 4 2])))
  (t/is (false? (apply m/eq [1 1 5])) "alias for ==, was true before fix")
  (t/is (false? (apply m/not== [1 2 3 2])))
  (t/is (true? (apply m/not== [1 2 3 4]))))

(t/deftest chained-comparison-inlined-path
  (t/is (false? (m/== 1 1 5)))
  (t/is (false? (m/== 1 1 5 5)))
  (t/is (true? (m/== 3 3 3 3)))
  (t/is (false? (m/< 1 2 0)))
  (t/is (true? (m/< 1 2 3 4)))
  (t/is (false? (m/> 3 2 5)))
  (t/is (true? (m/> 5 4 3 2)))
  (t/is (false? (m/<= 1 2 0)))
  (t/is (true? (m/<= 1 2 2 3)))
  (t/is (false? (m/>= 3 2 5)))
  (t/is (true? (m/>= 5 5 4 2)))
  (t/is (false? (m/eq 1 1 5)))
  (t/is (false? (m/not== 1 2 3 2)))
  (t/is (true? (m/not== 1 2 3 4))))

(let [numbers (vec (repeatedly 500 #(rand 100.0)))
      sorted-asc (vec (sort numbers))
      sorted-desc (vec (sort m/> numbers))]
  (t/deftest chained-comparison-cross-check
    (t/is (= (apply == [1.0 1.0 5.0]) (apply m/== [1.0 1.0 5.0])))
    (t/is (= (apply < sorted-asc) (apply m/< sorted-asc)))
    (t/is (true? (apply m/< sorted-asc)))
    (t/is (= (apply > sorted-desc) (apply m/> sorted-desc)))
    (t/is (true? (apply m/> sorted-desc)))
    (t/is (= (apply <= (conj sorted-asc (last sorted-asc)))
             (apply m/<= (conj sorted-asc (last sorted-asc)))))
    (t/is (= (apply >= (conj sorted-desc (last sorted-desc)))
             (apply m/>= (conj sorted-desc (last sorted-desc)))))))

;; Reference: clojure.core/bit-and, bit-or, bit-xor (associative, n-ary
;; semantics match directly) for the 3 operators clojure.core also defines;
;; hand-computed left-to-right pairwise folds (matching the library's own
;; documented semantics) for bit-nand/bit-nor/bit-xnor/bit-and-not, which have
;; no direct n-ary clojure.core equivalent. Fixture: 12=0b1100, 10=0b1010, 9=0b1001.
(t/deftest bitwise-n-ary
  (t/is (= 8 (m/bit-and 12 10 9) (apply m/bit-and [12 10 9]) (apply clojure.core/bit-and [12 10 9])))
  (t/is (= 15 (m/bit-or 12 10 9) (apply m/bit-or [12 10 9]) (apply clojure.core/bit-or [12 10 9])))
  (t/is (= 15 (m/bit-xor 12 10 9) (apply m/bit-xor [12 10 9]) (apply clojure.core/bit-xor [12 10 9])))
  (t/is (= -2 (m/bit-nand 12 10 9) (apply m/bit-nand [12 10 9])
            (clojure.core/bit-not (clojure.core/bit-and (clojure.core/bit-not (clojure.core/bit-and 12 10)) 9))))
  (t/is (= 6 (m/bit-nor 12 10 9) (apply m/bit-nor [12 10 9])
            (clojure.core/bit-not (clojure.core/bit-or (clojure.core/bit-not (clojure.core/bit-or 12 10)) 9))))
  (t/is (= 15 (m/bit-xnor 12 10 9) (apply m/bit-xnor [12 10 9])
             (clojure.core/bit-not (clojure.core/bit-xor (clojure.core/bit-not (clojure.core/bit-xor 12 10)) 9))))
  (t/is (= 4 (m/bit-and-not 12 10 9) (apply m/bit-and-not [12 10 9])
            (clojure.core/bit-and-not (clojure.core/bit-and-not 12 10) 9)))
  ;; single-argument identity
  (t/is (= 12 (m/bit-and 12) (m/bit-or 12) (m/bit-xor 12) (m/bit-nand 12) (m/bit-nor 12) (m/bit-xnor 12) (m/bit-and-not 12))))

;; Reference: hand-computed. xor/bool-xor are identical implementations
;; (associative, so pairwise-fold n-ary is correct); negative-zero?/integer?
;; verified against IEEE 754 bit-pattern facts; identity-double/identity-long
;; are trivially correct by inspection (single-line, no PrimitiveMath call).
(t/deftest boolean-and-identity-utilities
  (t/is (true? (m/xor true false)))
  (t/is (false? (m/xor true true)))
  (t/is (false? (apply m/xor [true false true])) "chained: t^f=t, t^t=f")
  (t/is (false? (apply m/xor [true true true true])))
  (t/is (true? (apply m/bool-xor [true false false])) "chained: t^f=t, t^f=t")
  (t/is (= (apply m/xor [true false true]) (apply m/bool-xor [true false true])) "identical implementations")
  (t/is (true? (m/negative-zero? -0.0)))
  (t/is (false? (m/negative-zero? 0.0)))
  (t/is (false? (m/negative-zero? 1.0)))
  (t/is (false? (m/negative-zero? -1.0)))
  (t/is (m/== 0.0 -0.0) "sanity: normal equality treats -0.0 and 0.0 as equal, unlike negative-zero?")
  (t/is (= 5.5 (m/identity-double 5.5)))
  (t/is (= 7 (m/identity-long 7)))
  (t/is (true? (m/integer? 5.0)))
  (t/is (true? (m/integer? -3.0)))
  (t/is (true? (m/integer? 0.0)))
  (t/is (false? (m/integer? 5.5)))
  (t/is (false? (m/integer? -0.1))))

;; Reference: fractions.Fraction (Python, exact rational arithmetic on the
;; IEEE 754 double bit patterns) for the catastrophic-cancellation case;
;; hand-computed for the basic cases. Fixture is the classic Kahan
;; two-product test case (a*b and c*d agree to ~8 significant digits).
(t/deftest fma-kahan-products
  (t/is (m/== 10.0 (m/muladd 2.0 3.0 4.0)))
  (t/is (m/== 10.0 (m/fma 2.0 3.0 4.0)))
  (t/is (m/== -2.0 (m/negmuladd 2.0 3.0 4.0)))
  (t/is (m/== 5.0 (m/difference-of-products 2.0 3.0 1.0 1.0)))
  (t/is (m/== 7.0 (m/sum-of-products 2.0 3.0 1.0 1.0)))
  (let [a 33962.035, b -30438.8, c 41563.4, d -24871.3
        exact -27800.538000075645
        naive (- (* a b) (* c d))
        dop (m/difference-of-products a b c d)]
    (t/is (m/delta-eq naive -27800.53800010681 1.0e-9) "naive computation loses precision to cancellation")
    (t/is (m/== exact dop) "Kahan algorithm recovers full precision")
    (t/is (< (Math/abs (- dop exact)) (Math/abs (- naive exact))) "strictly more accurate than naive"))
  (let [a 33962.035, b -30438.8, c 41563.4, d2 24871.3
        exact -27800.538000075645
        naive (+ (* a b) (* c d2))
        sop (m/sum-of-products a b c d2)]
    (t/is (m/delta-eq naive -27800.53800010681 1.0e-9) "naive computation loses precision to cancellation")
    (t/is (m/== exact sop) "Kahan algorithm recovers full precision")
    (t/is (< (Math/abs (- sop exact)) (Math/abs (- naive exact))) "strictly more accurate than naive")))

;; Reference: Python mpmath (50 digits precision) for pi/e/euler-gamma/catalan;
;; Apache Commons Math's own Gamma/GAMMA and Gamma/LANCZOS_G cross-checked
;; against mpmath.euler and the standard Lanczos g=607/128 parameter
;; respectively (fastmath.core's GAMMA/LANCZOS_G are direct aliases of these,
;; so verifying the upstream constants verifies fastmath's).
(t/deftest math-constants-block-1
  (let [pi 3.141592653589793
        e 2.718281828459045
        euler 0.5772156649015329
        catalan 0.915965594177219
        pi2 9.869604401089358
        macheps 1.1102230246251565e-16]
    (t/is (m/== m/PI pi))
    (t/is (m/== m/HALF_PI (/ pi 2.0)))
    (t/is (m/== m/THIRD_PI (/ pi 3.0)))
    (t/is (m/== m/QUARTER_PI (/ pi 4.0)))
    (t/is (m/== m/TWO_PI (* 2.0 pi)))
    (t/is (m/== m/TAU (* 2.0 pi)))
    (t/is (m/== m/E e))
    (t/is (m/== m/-PI (- pi)))
    (t/is (m/== m/-HALF_PI (- (/ pi 2.0))))
    (t/is (m/== m/-THIRD_PI (- (/ pi 3.0))) "regression: was +pi/3 before fix (double-negative sign bug)")
    (t/is (m/== m/-QUARTER_PI (- (/ pi 4.0))))
    (t/is (m/== m/-TWO_PI (- (* 2.0 pi))))
    (t/is (m/== m/-TAU (- (* 2.0 pi))))
    (t/is (m/== m/-E (- e)))
    (t/is (m/== m/INV_PI (/ 1.0 pi)))
    (t/is (m/== m/TWO_INV_PI (/ 2.0 pi)))
    (t/is (m/== m/FOUR_INV_PI (/ 4.0 pi)))
    (t/is (m/== m/INV_TWO_PI (/ 1.0 (* 2.0 pi))))
    (t/is (m/== m/INV_FOUR_PI (/ 1.0 (* 4.0 pi))))
    (t/is (m/== m/EPSILON 1.0e-10))
    (t/is (m/== m/GAMMA euler))
    (t/is (m/== m/LANCZOS_G (/ 607.0 128.0)))
    (t/is (m/== m/CATALAN_G catalan))
    (t/is (m/== m/PI2 pi2))
    (t/is (m/== m/MACHINE-EPSILON macheps))
    (t/is (m/== m/MACHINE-EPSILON10 (* 10.0 macheps)))
    (t/is (m/== m/THIRD (/ 1.0 3.0)))
    (t/is (m/== m/ONE_THIRD (/ 1.0 3.0)))
    (t/is (m/== m/TWO_THIRD (/ 2.0 3.0)))
    (t/is (m/== m/TWO_THIRDS (/ 2.0 3.0)))
    (t/is (m/== m/SIXTH (/ 1.0 6.0)))
    (t/is (m/== m/ONE_SIXTH (/ 1.0 6.0)))))

;; Reference: hand-computed sign semantics, incl. -0.0 edge case (neither
;; signum nor sgn treat -0.0 as negative, matching IEEE 754 numeric equality).
(t/deftest signum-sgn
  (t/is (= 1.0 (m/signum 5.0)))
  (t/is (= 0.0 (m/signum 0.0)))
  (t/is (= 0.0 (m/signum -0.0)))
  (t/is (= -1.0 (m/signum -3.0)))
  (t/is (= 1.0 (m/sgn 5.0)))
  (t/is (= 1.0 (m/sgn 0.0)))
  (t/is (= 1.0 (m/sgn -0.0)))
  (t/is (= -1.0 (m/sgn -3.0))))

(t/deftest agm
  (t/is (m/delta-eq 13.4581714817256154207668 (m/agm 24 6 1.0e-16) 1.0e-16)))

(t/deftest angles
  (t/is (m/delta-eq 180.0 (m/degrees m/PI)))
  (t/is (m/delta-eq m/PI (m/radians 180.0))))

(t/deftest frac
  (t/is (m/== 0.342225 (m/frac 3.342225)))
  (t/is (m/== 0.342225 (m/frac -3.342225)))
  (t/is (m/== 0.342225 (m/sfrac 3.342225)))
  (t/is (m/== -0.342225 (m/sfrac -3.342225))))

(t/deftest round-up-down
  (t/is (m/> (m/next-double 4.44) 4.44))
  (t/is (m/< (m/prev-double 4.44) 4.44))
  (t/is (m/== (m/round-up-pow2 1023) 1024))
  (t/is (m/== (m/round-up-pow2 1024) 1024))
  (t/is (m/== (m/round-up-pow2 1025) 2048)))

(t/deftest sgn
  (t/is (m/== -1.0 (m/signum -2)))
  (t/is (m/== 0.0 (m/signum 0)))
  (t/is (m/== 1.0 (m/signum 2)))
  (t/is (m/== -1.0 (m/sgn -2)))
  (t/is (m/== 1.0 (m/sgn 0)))
  (t/is (m/== 1.0 (m/sgn 2))))

(t/deftest norm
  (t/is (m/== (m/constrain -2 -1 1) -1))
  (t/is (m/== (m/constrain 0 -1 1) 0))
  (t/is (m/== (m/constrain 2 -1 1) 1))
  (t/is (m/== (m/norm 2 0 10) 0.2))
  (t/is (m/== (m/norm 2 0 10 0 100) 20.0)))

(t/deftest floating-points
  (t/is (m/== (m/double-exponent 2.0) 1))
  (t/is (m/== (m/double-exponent 0.5) -1))
  (t/is (m/== (m/double-significand 3.0) (Long/parseLong "1000000000000000000000000000000000000000000000000000" 2)))
  (t/is (m/== (m/double-significand 7.0) (Long/parseLong "1100000000000000000000000000000000000000000000000000" 2))))
;;


(t/deftest bernoulli
  (t/are [n res] (m/delta-eq (m/bernoulli n) (double res))
    0 1
    1 1/2
    2 1/6
    4 -1/30
    6 1/42
    8 -1/30
    10 5/66
    12 -691/2730
    14 7/6
    16 -3617/510
    18 43867/798
    20 -174611/330)
  (t/are [n] (m/zero? (m/bernoulli n))
    3 5 7 9 11 13 15 17 19 21))
