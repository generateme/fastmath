(ns fastmath.polynomials-test
  (:require [fastmath.polynomials :as sut]
            [clojure.test :as t]
            [clojure.java.io :as io]
            [clojure.edn :as edn]
            [fastmath.core :as m]
            [fastmath.complex :as cplx]))

(defn lp5 [^double x] (/ (sut/mevalpoly x 120.0 -600.0 600.0 -200.0 25.0 -1.0) 120.0))

(t/deftest laguerre
  (t/is (== 1.0 (sut/eval-laguerre-L 0 1)))
  (t/is (== 1.0 (sut/eval-laguerre-L 0 2 1)))
  (t/is (== -1.0 (sut/eval-laguerre-L 1 2)))
  (t/is (== -0.5 (sut/eval-laguerre-L 1 0.5 2)))
  (t/is (== (lp5 1) (sut/eval-laguerre-L 5 1)))
  (t/is (== (lp5 5) (sut/eval-laguerre-L 5 5)))
  (t/is (m/delta-eq (lp5 -1.234) (sut/eval-laguerre-L 5 -1.234))))

;; Scalar evaluators: evalpoly, mevalpoly, makepoly and their complex variants.
;;
;; Reference values: `test/resources/polynomials/evalpoly_reference.edn`, computed
;; in `mpmath` (50 digits) by `utils/fastmath/dev/generate_polynomials_reference.py`.
;; Cross-checked once against R: `polynom::predict` (real; 810 points; max error
;; 4.6e-16 of the error scale), `pracma::polyval` (complex coefficients; 392 points;
;; 3.6e-16) and `polynom::predict` with complex `z` (real coefficients; 630 points;
;; 5.5e-16).
;;
;; Tolerance: Horner's forward error is bounded by `c * n * eps * sum |c_i||x|^i`,
;; so every comparison uses that scale instead of a fixed absolute/relative number
;; (the polynomials are evaluated near their roots too). Largest observed ratio
;; error / (n * eps * scale) on the reference grid: real 0.62, complex 0.61,
;; real coefficients with complex z 3.2. Factors used: 2, 4 and 8.

(def ^:private evalpoly-reference
  (delay (edn/read-string (slurp (io/resource "polynomials/evalpoly_reference.edn")))))

(def ^:private EPS (Math/pow 2.0 -53))

(defn- exact->double
  "The nearest double of an exact rational (`double` of a ratio keeps only 16 decimal digits)."
  ^double [r] (@#'sut/rational->double r))

(defn- error-scale
  "Sum of `|c_i| r^i`, the scale of Horner's forward error. `moduli` are `|c_i|`, `r` is `|x|`."
  ^double [moduli ^double r]
  (reduce + (map-indexed (fn [^long i ^double c] (* c (Math/pow r i))) moduli)))

(defn- within-horner-bound?
  [^double error ^double scale ^long n ^double factor]
  (<= error (* factor n EPS scale)))

(defn- real-ok?
  [got expected cs x factor]
  (within-horner-bound? (m/abs (- (double got) (double expected)))
                        (error-scale (map m/abs cs) (m/abs (double x)))
                        (count cs) factor))

(defn- complex-ok?
  [got [re im] moduli r n factor]
  (within-horner-bound? (Math/hypot (- (cplx/re got) re) (- (cplx/im got) im))
                        (error-scale moduli (double r))
                        n factor))

(defn- cplx-modulus ^double [[re im]] (Math/hypot re im))

(defn- bits ^long [^double v] (Double/doubleToLongBits v))

(defn- cplx-form [[re im]] (list `cplx/complex re im))

;; function value built from the macro for fixed coefficients (the macros need the
;; coefficients at compile time)
(defn- mevalpoly-fn [cs] (eval `(fn [x#] (sut/mevalpoly x# ~@cs))))
(defn- mevalpoly-complex-fn [cs] (eval `(fn [z#] (sut/mevalpoly-complex z# ~@(map cplx-form cs)))))
(defn- mevalpoly-scalar-complex-fn [cs] (eval `(fn [z#] (sut/mevalpoly-scalar-complex z# ~@cs))))

(t/deftest evalpoly-reference-real
  (let [{:keys [polys x ref]} (:real @evalpoly-reference)]
    (doseq [[cs row] (map vector polys ref)
            :let [p (sut/makepoly cs)
                  mp (mevalpoly-fn cs)]
            [xx expected] (map vector x row)]
      (t/is (real-ok? (apply sut/evalpoly xx cs) expected cs xx 2.0) (str "evalpoly " cs " " xx))
      (t/is (real-ok? (p xx) expected cs xx 2.0) (str "makepoly " cs " " xx))
      (t/is (real-ok? (mp xx) expected cs xx 2.0) (str "mevalpoly " cs " " xx)))))

(t/deftest evalpoly-reference-complex
  (let [{:keys [polys z ref]} (:complex @evalpoly-reference)]
    (doseq [[cs row] (map vector polys ref)
            :let [ccs (map (fn [[re im]] (cplx/complex re im)) cs)
                  moduli (map cplx-modulus cs)
                  p (sut/makepoly-complex ccs)
                  mp (mevalpoly-complex-fn cs)]
            [zz expected] (map vector z row)
            :let [zc (cplx/complex (first zz) (second zz))
                  r (cplx-modulus zz)]]
      (t/is (complex-ok? (apply sut/evalpoly-complex zc ccs) expected moduli r (count cs) 4.0)
            (str "evalpoly-complex " cs " " zz))
      (t/is (complex-ok? (p zc) expected moduli r (count cs) 4.0) (str "makepoly-complex " cs " " zz))
      (t/is (complex-ok? (mp zc) expected moduli r (count cs) 4.0) (str "mevalpoly-complex " cs " " zz)))))

(t/deftest evalpoly-reference-scalar-complex
  (let [{:keys [polys z ref]} (:scalar-complex @evalpoly-reference)]
    (doseq [[cs row] (map vector polys ref)
            :let [moduli (map m/abs cs)
                  p (sut/makepoly-scalar-complex cs)
                  mp (mevalpoly-scalar-complex-fn cs)]
            [zz expected] (map vector z row)
            :let [zc (cplx/complex (first zz) (second zz))
                  r (cplx-modulus zz)]]
      (t/is (complex-ok? (apply sut/evalpoly-scalar-complex zc cs) expected moduli r (count cs) 8.0)
            (str "evalpoly-scalar-complex " cs " " zz))
      (t/is (complex-ok? (p zc) expected moduli r (count cs) 8.0)
            (str "makepoly-scalar-complex " cs " " zz))
      (t/is (complex-ok? (mp zc) expected moduli r (count cs) 8.0)
            (str "mevalpoly-scalar-complex " cs " " zz)))))

(t/deftest evalpoly-r-spot-check
  ;; R: polynom::predict(polynomial(c(1,-2,3,0.5,-4)), c(0.7,-1.3,2.5))
  (doseq [[x expected] [[0.7 0.28109999999999991] [-1.3 -3.85290000000000088] [2.5 -133.6875]]]
    (t/is (m/delta-eq expected (sut/evalpoly x 1.0 -2.0 3.0 0.5 -4.0) 1.0e-13))
    (t/is (m/delta-eq expected (sut/mevalpoly x 1.0 -2.0 3.0 0.5 -4.0) 1.0e-13))
    (t/is (m/delta-eq expected ((sut/makepoly [1.0 -2.0 3.0 0.5 -4.0]) x) 1.0e-13)))
  ;; R: polynom::predict(polynomial(c(1,-2,3,0.5)), complex(real=0.7, imaginary=-1.3))
  (let [z (cplx/complex 0.7 -1.3)
        expected (cplx/complex -5.6030000000000015 -2.7169999999999996)]
    (t/is (cplx/delta-eq expected (sut/evalpoly-scalar-complex z 1.0 -2.0 3.0 0.5) 1.0e-13))
    (t/is (cplx/delta-eq expected (sut/evalpoly-complex z 1.0 -2.0 3.0 0.5) 1.0e-13))
    (t/is (cplx/delta-eq expected (sut/mevalpoly-scalar-complex z 1.0 -2.0 3.0 0.5) 1.0e-13)))
  ;; R: pracma::polyval(rev(c(1+2i, 1i, 3-1i)), 0.7-1.3i)
  (let [z (cplx/complex 0.7 -1.3)
        expected (cplx/complex -3.12 -1.56)]
    (t/is (cplx/delta-eq expected (sut/evalpoly-complex z (cplx/complex 1 2) (cplx/complex 0 1) (cplx/complex 3 -1)) 1.0e-13))
    (t/is (cplx/delta-eq expected (sut/mevalpoly-complex z (cplx/complex 1 2) (cplx/complex 0 1) (cplx/complex 3 -1)) 1.0e-13))))

;; no coefficients, one coefficient, many; zero coefficients; special arguments

(t/deftest evalpoly-real-edge-cases
  ;; empty: zero polynomial, as a double
  (t/is (identical? Double (class (sut/evalpoly 2.0))))
  (t/is (identical? Double (class (apply sut/evalpoly 2.0 []))))
  (t/is (identical? Double (class (sut/mevalpoly 2.0))))
  (t/is (== 0.0 (sut/evalpoly 2.0) (apply sut/evalpoly 2.0 []) (sut/mevalpoly 2.0)))
  (t/is (== 0.0 ((sut/makepoly []) 2.0) ((sut/makepoly nil) 2.0)))
  ;; one coefficient: a double in every form, whatever the type of the coefficient
  (doseq [c [5 5.0 5N 10/2]]
    (t/is (identical? Double (class (sut/evalpoly 2.0 c))) (str "evalpoly " (class c)))
    (t/is (identical? Double (class (apply sut/evalpoly 2.0 [c]))) (str "apply evalpoly " (class c)))
    (t/is (identical? Double (class (sut/mevalpoly 2.0 c))) (str "mevalpoly " (class c)))
    (t/is (identical? Double (class ((sut/makepoly [c]) 2.0))) (str "makepoly " (class c)))
    (t/is (== 5.0 (sut/evalpoly 2.0 c) (apply sut/evalpoly 2.0 [c]) (sut/mevalpoly 2.0 c) ((sut/makepoly [c]) 2.0))))
  (let [c 7]
    (t/is (identical? Double (class (sut/mevalpoly 2.0 c))) "symbol coefficient")
    (t/is (== 7.0 (sut/mevalpoly 2.0 c))))
  ;; two and more coefficients
  (t/is (== 17.0 (sut/evalpoly 2.0 1 2 3) (sut/mevalpoly 2.0 1 2 3) ((sut/makepoly [1 2 3]) 2.0)))
  (t/is (== 17.0 (sut/evalpoly 2 1 2 3)) "integer argument")
  (t/is (== 2.75 (sut/evalpoly 1/2 1 2 3)) "ratio argument")
  (t/is (== 3.0 (sut/evalpoly 0.0 3 4 5) (sut/evalpoly -0.0 3 4 5) (sut/mevalpoly 0.0 3 4 5)) "x=0 gives c0")
  (t/is (== 12.0 (sut/mevalpoly 2.0 0 0 3) (sut/evalpoly 2.0 0 0 3)) "zero constant and linear coefficients")
  (t/is (== 4.0 (sut/evalpoly 2.0 0 0 1) (sut/mevalpoly 2.0 0 0 1)))
  (t/is (== 1.0 (sut/evalpoly 2.0 1 0 0 0) (sut/mevalpoly 2.0 1 0 0 0)) "trailing zero coefficients")
  ;; special arguments
  (t/is (== ##Inf (sut/evalpoly ##Inf 1 2 3) (sut/mevalpoly ##Inf 1 2 3) (sut/evalpoly ##Inf 0 0 1)
            (sut/mevalpoly ##Inf 0 0 1) (sut/evalpoly 1.0e308 0 0 1) (sut/mevalpoly 1.0e308 0 0 1)))
  (t/is (== ##-Inf (sut/evalpoly ##-Inf 1 2) (sut/mevalpoly ##-Inf 1 2)))
  (t/is (== ##Inf (sut/evalpoly ##-Inf 1 2 3) (sut/mevalpoly ##-Inf 1 2 3)))
  (t/is (m/nan? (sut/evalpoly ##NaN 1 2)))
  (t/is (m/nan? (sut/mevalpoly ##NaN 1 2)))
  (t/is (m/nan? (sut/mevalpoly ##NaN 0 0 1)))
  (t/is (m/nan? ((sut/makepoly [1 2 3]) ##NaN)))
  ;; sign of a zero result is the same in all forms (IEEE Horner: the constant coefficient is always added)
  (doseq [x [0.0 -0.0 1.0 -1.0 -2.5 -1.0e308]]
    (t/is (= (bits 0.0) (bits (sut/mevalpoly x 0 0 0)) (bits (apply sut/evalpoly x [0 0 0])) (bits ((sut/makepoly [0 0 0]) x)))
          (str "zero polynomial at " x))))

(t/deftest evalpoly-macro-and-function-forms-agree
  ;; every coefficient count 0..8 with leading/trailing/all zeros and signed zeros, against a set
  ;; of special arguments; results must be identical bit for bit
  (let [xs [0.0 -0.0 1.0 -1.0 0.5 -2.5 3.0 ##Inf ##-Inf ##NaN 1.0e308 -1.0e308 1.0e-300]
        bases [[1 2 3 4 5 6 7 8] [0 2 3 4 5 6 7 8] [0.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]
               [1.5 0 0 0 0 0 0 0] [0 0 0 0 0 0 0 0] [-0.0 1 2 3 4 5 6 7]
               [1 2 3 4 5 6 7 0] [0 0 0 1 0 0 0 0]]
        sets (distinct (for [n (range 0 9) b bases] (vec (take n b))))]
    (doseq [cs sets
            :let [mf (mevalpoly-fn cs)
                  p (sut/makepoly cs)]
            x xs
            :let [f (double (apply sut/evalpoly x cs))]]
      (t/is (= (bits f) (bits (double (mf x))) (bits (double (p x))))
            (str "coefficients " cs " at " x)))))

(t/deftest evalpoly-complex-edge-cases
  (let [z (cplx/complex 0.7 -1.3)
        one (cplx/complex 5 0)]
    ;; empty
    (t/is (= cplx/ZERO (sut/evalpoly-complex z) (apply sut/evalpoly-complex z [])
             (sut/evalpoly-scalar-complex z) (apply sut/evalpoly-scalar-complex z [])
             (sut/mevalpoly-complex z) (sut/mevalpoly-scalar-complex z)
             ((sut/makepoly-complex []) z) ((sut/makepoly-complex nil) z)
             ((sut/makepoly-scalar-complex []) z) ((sut/makepoly-scalar-complex nil) z)))
    ;; one coefficient: the coefficient as a complex number, independent of z
    (doseq [c [5 5.0 10/2]]
      (t/is (= one (sut/evalpoly-complex z c)) (str "evalpoly-complex " (class c)))
      (t/is (= one (sut/evalpoly-scalar-complex z c)) (str "evalpoly-scalar-complex " (class c)))
      (t/is (= one (sut/mevalpoly-scalar-complex z c)) (str "mevalpoly-scalar-complex " (class c)))
      (t/is (= one ((sut/makepoly-complex [c]) z)) (str "makepoly-complex " (class c)))
      (t/is (= one ((sut/makepoly-scalar-complex [c]) z)) (str "makepoly-scalar-complex " (class c))))
    (t/is (= one (sut/evalpoly-complex z one) (sut/mevalpoly-complex z one) ((sut/makepoly-complex [one]) z)))
    (let [c (cplx/complex 1 2)]
      (t/is (= c (sut/evalpoly-complex z c) (sut/mevalpoly-complex z c) ((sut/makepoly-complex [c]) z)))
      (t/is (= cplx/ZERO (sut/evalpoly-complex z cplx/ZERO) ((sut/makepoly-complex [cplx/ZERO]) z))))
    ;; z given as a real number
    (t/is (= (cplx/complex 17 0) (sut/evalpoly-complex 2.0 1 2 3) (sut/evalpoly-scalar-complex 2.0 1 2 3)
             ((sut/makepoly-complex [1 2 3]) 2) ((sut/makepoly-scalar-complex [1 2 3]) 2)))
    ;; z = 0 gives the constant coefficient
    (t/is (= (cplx/complex 3 0) (sut/evalpoly-scalar-complex cplx/ZERO 3 4 5) (sut/evalpoly-complex cplx/ZERO 3 4 5)))
    ;; (1+i)^2 = 2i
    (let [w (cplx/complex 1 1)]
      (t/is (cplx/delta-eq (cplx/complex 0 2) (sut/evalpoly-complex w 0 0 1) 1.0e-15))
      (t/is (cplx/delta-eq (cplx/complex 0 2) (sut/evalpoly-scalar-complex w 0 0 1) 1.0e-15))
      (t/is (cplx/delta-eq (cplx/complex 0 2) (sut/mevalpoly-scalar-complex w 0 0 1) 1.0e-15))
      (t/is (cplx/delta-eq (cplx/complex 0 2) (sut/mevalpoly-complex w cplx/ZERO cplx/ZERO cplx/ONE) 1.0e-15)))
    ;; real coefficients: p(conj z) = conj p(z)
    (doseq [zz [(cplx/complex 0.3 0.4) (cplx/complex -2.0 5.0) (cplx/complex 1.0e-3 -7.0)]
            :let [zc (cplx/conjugate zz)]]
      (t/is (cplx/delta-eq (cplx/conjugate (sut/evalpoly-scalar-complex zz 1 -2 3 0.5 4))
                           (sut/evalpoly-scalar-complex zc 1 -2 3 0.5 4) 1.0e-9)))))

(t/deftest evalpoly-scalar-complex-large-and-small-arguments
  ;; |z|^2 overflows (|z| > ~1.3e154) while the polynomial value is still finite: p(z) = 1 + 2z
  ;; R: predict(polynomial(c(1,2)), complex(real=1e200, imaginary=1)) = 2e200+2i
  (let [z (cplx/complex 1.0e200 1.0)]
    (t/is (m/delta-eq 2.0 (cplx/im (sut/evalpoly-scalar-complex z 1 2)) 1.0e-12))
    (t/is (m/delta-eq 1.0 (/ (cplx/re (sut/evalpoly-scalar-complex z 1 2)) 2.0e200) 1.0e-12))
    (t/is (m/delta-eq 1.0 (/ (cplx/re ((sut/makepoly-scalar-complex [1 2]) z)) 2.0e200) 1.0e-12))
    (t/is (m/delta-eq 1.0 (/ (cplx/re (sut/mevalpoly-scalar-complex z 1 2)) 2.0e200) 1.0e-12)))
  ;; a real overflow of the result gives Inf, not NaN (real part of 3 z^2 + 2 z + 1 at z = 1e200)
  (let [z (cplx/complex 1.0e200 0.0)]
    (t/is (== ##Inf (cplx/re (sut/evalpoly-scalar-complex z 1 2 3))))
    (t/is (== ##Inf (cplx/re (sut/evalpoly-complex z 1 2 3)))))
  (t/is (== ##Inf (cplx/re (sut/evalpoly-scalar-complex (cplx/complex ##Inf 0.0) 1 2))))
  ;; tiny arguments: |z|^2 underflows, polynomial value is c0 + c1 z
  (let [z (cplx/complex 1.0e-200 1.0e-200)
        v (sut/evalpoly-scalar-complex z 1 2 3)]
    (t/is (== 1.0 (cplx/re v)))
    (t/is (cplx/delta-eq (sut/evalpoly-complex z 1 2 3) v 1.0e-300))))

;; Polynomial objects: `Polynomial` (double coefficients, from `polynomial`) and
;; `PolynomialR` (exact rational coefficients, from `ratio-polynomial`).
;;
;; Reference values: `test/resources/polynomials/polynomial_ops_reference.edn`, exact
;; rational arithmetic (Python `fractions`) from decimal coefficients, generated by
;; `utils/fastmath/dev/generate_polynomial_ops_reference.py`. `PolynomialR` results must
;; equal them exactly. `Polynomial` results must be within a few roundoff units of the error
;; scale of each coefficient (the sum of the absolute values of the terms that form it).
;; Cross-checked once against R `polynom` (`+`, `-`, `*`, `deriv`, `predict`, `poly.calc`):
;; max relative differences sum 8.3e-15, difference 9.1e-16, product 2.9e-15 (decimals
;; cancel, so these are input-rounding amplification), deriv 2.7e-16, predict 2.8e-16 of the
;; evaluation scale, poly.calc 9.5e-15.
;;
;; Tolerances (roundoff units, eps = 2^-53, of the per-coefficient error scale) against the
;; largest observed values: add/sub 4 (observed 1.9), scale 4 (2.1), mult 4 per 1 + the shorter
;; length (0.73), derivative 4 per 1 + length (0.38), evaluate 4 per 1 + length (0.66).

(def ^:private ops-reference
  (delay (edn/read-string (slurp (io/resource "polynomials/polynomial_ops_reference.edn")))))

(defn- ratio [[n d]] (/ (bigint n) (bigint d)))

(defn- attempt
  "Result of `(apply f args)` or the thrown exception, so a crash fails one assertion, not the whole test."
  [f & args]
  (try (apply f args) (catch Throwable e e)))

(defn- failed? [result] (instance? Throwable result))

(defn- bits-or-failure
  "Bit pattern of the double returned by `(apply f args)`, or -1 (never a bit pattern of a double result) if it threw."
  [f & args]
  (let [result (apply attempt f args)]
    (if (failed? result) -1 (bits (double result)))))

(defn- type-name [p] (.getSimpleName (class p)))

(defn- abs-vec [cs] (mapv #(m/abs (double %)) cs))

(defn- padded [v n] (into v (repeat (- n (count v)) 0.0)))

(defn- sum-scales [a b]
  (let [n (max (count a) (count b))]
    (mapv + (padded (abs-vec a) n) (padded (abs-vec b) n))))

(defn- product-scales [a b]
  (let [aa (abs-vec a) bb (abs-vec b)]
    (vec (for [k (range (+ (count a) (count b) -1))]
           (reduce + 0.0 (for [i (range (count a)) :let [j (- k i)] :when (< -1 j (count bb))]
                           (* (aa i) (bb j))))))))

(defn- derivative-scales [cs order]
  (let [aa (abs-vec cs)]
    (vec (for [j (range (- (count cs) order))]
           (* (aa (+ j order)) (double (reduce * 1.0 (range (inc j) (+ j order 1)))))))))

(defn- coefficients-within?
  "True when every double coefficient in `got` is within `units` roundoff units of the exact
  coefficient in `expected`, scaled by the matching entry of `scales`."
  [got expected scales units]
  (and (not (failed? got))
       (= (count got) (count expected) (count scales))
       (every? true? (map (fn [g e s] (<= (m/abs (- (double g) (double e))) (* units EPS (double s))))
                          got expected scales))))

(t/deftest polynomial-add-sub-mult-reference
  (doseq [{:keys [a b sum diff product]} (:pairs @ops-reference)
          [op-name op expected scales units] [["add" sut/add sum (sum-scales a b) 4.0]
                                              ["sub" sut/sub diff (sum-scales a b) 4.0]
                                              ["mult" sut/mult product (product-scales a b)
                                               (* 4.0 (inc (min (count a) (count b))))]]
          :let [exact (mapv ratio expected)
                p (attempt op (sut/polynomial a) (sut/polynomial b))
                r (attempt op (sut/ratio-polynomial a) (sut/ratio-polynomial b))]]
    (t/is (and (not (failed? p))
               (coefficients-within? (vec (sut/coeffs p)) exact scales units)
               (= (dec (count exact)) (sut/degree p)))
          (str "Polynomial " op-name " " a " " b " -> " p))
    (t/is (and (not (failed? r))
               (= exact (vec (sut/coeffs r)))
               (= (dec (count exact)) (sut/degree r)))
          (str "PolynomialR " op-name " " a " " b " -> " r))))

(t/deftest polynomial-scale-reference
  (let [{:keys [polys scalars ref]} (:scale @ops-reference)]
    (doseq [[cs row] (map vector polys ref)
            [s expected] (map vector scalars row)
            :let [exact (mapv ratio expected)
                  scales (mapv #(* % (m/abs s)) (abs-vec cs))
                  p (attempt sut/scale (sut/polynomial cs) s)
                  r (attempt sut/scale (sut/ratio-polynomial cs) s)]]
      (t/is (and (not (failed? p))
                 (coefficients-within? (vec (sut/coeffs p)) exact scales 4.0)
                 (= (dec (count cs)) (sut/degree p)))
            (str "Polynomial scale " cs " " s))
      (t/is (and (not (failed? r))
                 (= exact (vec (sut/coeffs r)))
                 (= (dec (count cs)) (sut/degree r)))
            (str "PolynomialR scale " cs " " s)))))

(t/deftest polynomial-derivative-reference
  ;; orders 0 to degree+2: orders above the degree give the zero polynomial (degree 0)
  (let [{:keys [polys ref]} (:derivative @ops-reference)]
    (doseq [[cs row] (map vector polys ref)
            [order expected] (map-indexed vector row)
            :let [exact (mapv ratio expected)
                  scales (if (> order (dec (count cs))) [0.0] (derivative-scales cs order))
                  p (attempt sut/derivative (sut/polynomial cs) order)
                  r (attempt sut/derivative (sut/ratio-polynomial cs) order)]]
      (t/is (and (not (failed? p))
                 (coefficients-within? (vec (sut/coeffs p)) exact scales (* 4.0 (inc (count cs))))
                 (= (dec (count exact)) (sut/degree p)))
            (str "Polynomial derivative " cs " order " order " -> " p))
      (t/is (and (not (failed? r))
                 (= exact (vec (sut/coeffs r)))
                 (= (dec (count exact)) (sut/degree r)))
            (str "PolynomialR derivative " cs " order " order " -> " r))))
  (let [cs [1 2 3 4]]
    (t/is (= (vec (sut/coeffs (sut/derivative (sut/ratio-polynomial cs) 1)))
             (vec (sut/coeffs (sut/derivative (sut/ratio-polynomial cs)))))
          "default order is 1")
    (doseq [make [sut/polynomial sut/ratio-polynomial]]
      (t/is (= (vec (sut/coeffs (make cs))) (vec (sut/coeffs (sut/derivative (make cs) 0)))) "order 0")
      (t/is (thrown? IllegalArgumentException (sut/derivative (make cs) -1)) (str (type-name (make cs)) " order -1"))
      (t/is (thrown? IllegalArgumentException (sut/derivative (make cs) -5))))))

(t/deftest ratio-polynomial-derivative-is-exact-at-high-order
  ;; independent formula: coefficient j of the order k derivative is c(j+k) * (j+1)(j+2)...(j+k)
  (let [cs (mapv (fn [i] (/ (inc i) (+ 3 i))) (range 41))
        r (sut/ratio-polynomial cs)]
    (doseq [k [1 10 20 22 23 25 30 39 40]
            :let [expected (vec (for [j (range (- 41 k))]
                                  (*' (cs (+ j k)) (reduce *' 1 (range (inc j) (+ j k 1))))))
                  d (attempt sut/derivative r k)]]
      (t/is (and (not (failed? d)) (= expected (vec (sut/coeffs d))) (= (- 40 k) (sut/degree d)))
            (str "order " k))))
  ;; order beyond the range of a double factorial (170!): still exact, no exception
  (let [r (sut/ratio-polynomial (range 1 202))
        d (attempt sut/derivative r 180)]
    (t/is (and (not (failed? d))
               (= (*' 181 (reduce *' 1 (range 1 181))) (first (sut/coeffs d)))
               (= 20 (sut/degree d))))))

(t/deftest polynomial-evaluate-reference
  (let [{:keys [polys x ref]} (:evaluate @ops-reference)]
    (doseq [[cs row] (map vector polys ref)
            [xx expected] (map vector x row)
            :let [exact (ratio expected)
                  bound (* 4.0 (inc (count cs)) EPS (error-scale (abs-vec cs) (m/abs xx)))
                  p (sut/polynomial cs)
                  r (sut/ratio-polynomial cs)]]
      (t/is (<= (m/abs (- (sut/evaluate p xx) (exact->double exact))) bound) (str "Polynomial evaluate " cs " at " xx))
      (t/is (<= (m/abs (- (double (p xx)) (exact->double exact))) bound) (str "Polynomial call " cs " at " xx))
      (t/is (= exact (r xx)) (str "PolynomialR call " cs " at " xx))
      (t/is (== (exact->double exact) (sut/evaluate r xx)) (str "PolynomialR evaluate " cs " at " xx)))))

(t/deftest zero-polynomial
  ;; the zero polynomial has degree 0 and a single zero coefficient
  (doseq [[label z] [["polynomial []" (sut/polynomial [])]
                     ["polynomial nil" (sut/polynomial nil)]
                     ["coeffs->polynomial" (sut/coeffs->polynomial)]
                     ["ratio-polynomial []" (sut/ratio-polynomial [])]
                     ["ratio-polynomial nil" (sut/ratio-polynomial nil)]
                     ["coeffs->ratio-polynomial" (sut/coeffs->ratio-polynomial)]]]
    (t/is (= 0 (sut/degree z)) label)
    (t/is (= [0.0] (mapv double (sut/coeffs z))) label)
    (t/is (every? #(== 0.0 %) (map #(sut/evaluate z %) [-1.0 0.0 2.5 ##Inf])) label)
    (t/is (== 0.0 (z 3.0)) label))
  ;; derivatives beyond the degree, and of constants and of the zero polynomial itself
  (doseq [make [sut/polynomial sut/ratio-polynomial]
          [cs order] [[[1 2 3] 3] [[1 2 3] 4] [[1 2 3] 10] [[7] 1] [[7] 5] [[0] 1] [[0 0] 3]]
          :let [d (attempt sut/derivative (make cs) order)]]
    (t/is (and (not (failed? d))
               (= 0 (sut/degree d))
               (= [0.0] (mapv double (sut/coeffs d)))
               (== 0.0 (sut/evaluate d 1.5)))
          (str (type-name (make cs)) " derivative of " cs " order " order " -> " d))))

(t/deftest polynomial-equality-and-hash
  (let [a (sut/polynomial [1 2 3]) b (sut/polynomial [1.0 2.0 3.0]) c (sut/polynomial [1 2 4])
        ra (sut/ratio-polynomial [1 2 3]) rb (sut/ratio-polynomial [1.0 2.0 3.0]) rc (sut/ratio-polynomial [1 2 4])]
    (t/is (= a b))
    (t/is (= (hash a) (hash b)))
    (t/is (not= a c))
    (t/is (= ra rb))
    (t/is (= (hash ra) (hash rb)))
    (t/is (not= ra rc))
    (t/is (not= a ra))
    (t/is (not= ra a))
    (t/is (not= a nil))
    (t/is (not= ra nil))
    (t/is (not= ra 5))
    (t/is (not= a "polynomial"))
    ;; usable as keys of hash maps and elements of hash sets
    (t/is (contains? #{ra} rb))
    (t/is (= :x ({ra :x} rb)))
    (t/is (contains? #{a} b))
    (t/is (= 2 (count (set [ra rb rc]))))
    ;; the same rational coefficients built in different numeric types are equal, with equal hashes
    (let [product (sut/mult (sut/ratio-polynomial [1 1]) (sut/ratio-polynomial [1 -1]))
          expected (sut/ratio-polynomial [1 0 -1])]
      (t/is (= expected product))
      (t/is (= (hash expected) (hash product))))
    (let [sum (sut/add (sut/ratio-polynomial [1/2 1/3]) (sut/ratio-polynomial [1/2 -1/3]))]
      (t/is (= (sut/ratio-polynomial [1 0]) sum))
      (t/is (= (hash (sut/ratio-polynomial [1 0])) (hash sum))))
    ;; double coefficients are compared by bit pattern: 0.0 and -0.0 differ
    (t/is (not= (sut/polynomial [0.0]) (sut/polynomial [-0.0])))))

(t/deftest evaluate-non-finite-arguments
  ;; PolynomialR returns what Polynomial returns for the same coefficients
  (doseq [cs [[1 2 3] [5] [0 0 -1] [1 -2 0.5 4] [0] [3 0]]
          x [##NaN ##Inf ##-Inf]
          :let [d (sut/polynomial cs)
                r (sut/ratio-polynomial cs)
                expected (bits (sut/evaluate d x))]]
    (t/is (= expected (bits-or-failure sut/evaluate r x)) (str "evaluate " cs " at " x))
    (t/is (= expected (bits-or-failure r x)) (str "call " cs " at " x))))

(t/deftest mixed-polynomial-types
  (let [p (sut/polynomial [1 2]) r (sut/ratio-polynomial [1 2])]
    (doseq [[label f] [["add" sut/add] ["sub" sut/sub] ["mult" sut/mult]]
            [a b] [[p r] [r p]]]
      (t/is (thrown-with-msg? IllegalArgumentException #"cannot be combined" (f a b))
            (str label " " (type-name a) " " (type-name b))))
    (doseq [[label f] [["add" sut/add] ["mult" sut/mult]]
            other [5 nil "x" [1 2]]
            poly [p r]]
      (t/is (thrown? IllegalArgumentException (f poly other)) (str label " " (type-name poly) " with " (pr-str other))))))

(t/deftest scale-and-unary-operations
  (let [p (sut/polynomial [1 2 3]) r (sut/ratio-polynomial [1 2 3])]
    (t/is (every? m/nan? (sut/coeffs (sut/scale p ##NaN))))
    (t/is (every? #(== ##Inf %) (sut/coeffs (sut/scale p ##Inf))))
    (t/is (thrown? IllegalArgumentException (sut/scale r ##NaN)))
    (t/is (thrown? IllegalArgumentException (sut/scale r ##Inf)))
    (t/is (= [0 0 0] (vec (sut/coeffs (sut/scale r 0)))))
    (t/is (= 2 (sut/degree (sut/scale p 0))))
    (t/is (= [1/2 1 3/2] (vec (sut/coeffs (sut/scale r 1/2)))))
    ;; one-argument forms
    (t/is (identical? p (sut/add p)))
    (t/is (identical? r (sut/mult r)))
    (t/is (= [-1.0 -2.0 -3.0] (vec (sut/coeffs (sut/sub p)))))
    (t/is (= [-1 -2 -3] (vec (sut/coeffs (sut/sub r)))))
    (t/is (= (sut/derivative p 1) (sut/derivative p)))))

(t/deftest polynomial-string-form
  (t/is (= "#polynomial{2}(x) = 1+2x+3x^2" (str (sut/polynomial [1 2 3]))))
  (t/is (= "#polynomial{2}(x) = 1+2x+3x^2" (str (sut/ratio-polynomial [1 2 3]))))
  (t/is (= "#polynomial{2}(x) = 1-2x+0.5x^2" (str (sut/polynomial [1 -2 0.5]))))
  (t/is (= "#polynomial{1}(x) = 0.3333+0.5x" (str (sut/ratio-polynomial [1/3 1/2]))))
  (t/is (= "#polynomial{0}(x) = 0" (str (sut/polynomial []))))
  (t/is (= "#polynomial{2}(x) = 1+2x+3x^2" (pr-str (sut/polynomial [1 2 3]))))
  ;; no digit grouping
  (t/is (= "#polynomial{2}(x) = 1234567.891+1x+10000000000x^2" (str (sut/polynomial [1234567.891 1 10000000000]))))
  ;; a coefficient that rounds to zero is not printed as -0
  (t/is (= "#polynomial{1}(x) = 0+1x" (str (sut/polynomial [-0.00001 1]))))
  (t/is (= "#polynomial{2}(x) = 1+0x+1x^2" (str (sut/polynomial [1 -0.00001 1]))))
  ;; non-finite coefficients are printed, not thrown on
  (t/is (re-find #"NaN" (str (sut/polynomial [##NaN 1]))))
  (t/is (string? (str (sut/polynomial [##Inf ##-Inf 1]))))
  ;; long polynomials are abbreviated
  (t/is (.endsWith ^String (str (sut/polynomial (range 13))) "+..."))
  ;; independent of the default locale
  (let [old (java.util.Locale/getDefault)]
    (try
      (java.util.Locale/setDefault (java.util.Locale/forLanguageTag "pl-PL"))
      (t/is (= "#polynomial{2}(x) = 1+0.5x+1234567.891x^2" (str (sut/polynomial [1 0.5 1234567.891]))))
      (finally (java.util.Locale/setDefault old)))))

(t/deftest polynomial-constructors
  (t/is (= "Polynomial" (type-name (sut/polynomial [1 2]))))
  (t/is (= "PolynomialR" (type-name (sut/ratio-polynomial [1 2]))))
  (t/is (= "PolynomialR" (type-name (sut/coeffs->ratio-polynomial 1 2))))
  (t/is (= "Polynomial" (type-name (sut/coeffs->polynomial 1 2))))
  ;; accepted collections
  (doseq [cs [[0 1 2] '(0 1 2) (range 3) (double-array [0 1 2]) (map identity [0 1 2])]]
    (t/is (= [0.0 1.0 2.0] (vec (sut/coeffs (sut/polynomial cs)))) (str (class cs)))
    (t/is (= [0 1 2] (vec (sut/coeffs (sut/ratio-polynomial (vec (seq cs)))))) (str (class cs))))
  (t/is (= (sut/polynomial [1 2 3]) (sut/coeffs->polynomial 1 2 3)))
  (t/is (= (sut/ratio-polynomial [1 2 3]) (sut/coeffs->ratio-polynomial 1 2 3)))
  (t/is (every? double? (sut/coeffs (sut/polynomial [1 2 3]))))
  ;; ratio coefficients are the exact decimal values of the numbers given
  (t/is (= [1/10 1/4 1/3] (vec (sut/coeffs (sut/ratio-polynomial [0.1 0.25 1/3])))))
  (t/is (thrown? IllegalArgumentException (sut/ratio-polynomial [##NaN])))
  (t/is (thrown? IllegalArgumentException (sut/ratio-polynomial [1 ##Inf])))
  ;; the degree is nominal: trailing zero coefficients are kept
  (t/is (= 3 (sut/degree (sut/polynomial [1 2 0 0]))))
  (t/is (= 3 (sut/degree (sut/ratio-polynomial [1 2 0 0])))))

(t/deftest polynomial-fitted-constructors
  (doseq [{:keys [xs ys ref]} (:fit @ops-reference)
          :let [exact (mapv ratio ref)
                p (attempt sut/polynomial xs ys)
                r (attempt sut/ratio-polynomial xs ys)
                ;; observed max error of the fitted coefficients: 4.4e-14 relative to 1 + |exact|
                ;; (nodes spread over [-1, 1], up to 8 points: the Vandermonde system is ill-conditioned)
                tolerance (fn [e] (* 1.0e-12 (+ 1.0 (m/abs (double e)))))]]
    (t/is (and (not (failed? p)) (= (count xs) (inc (sut/degree p)))
               (every? true? (map (fn [g e] (<= (m/abs (- g (double e))) (tolerance e))) (sut/coeffs p) exact)))
          (str "Polynomial fit " xs))
    (t/is (and (not (failed? r)) (= "PolynomialR" (type-name r)) (= (count xs) (inc (sut/degree r)))
               (every? true? (map (fn [g e] (<= (m/abs (- (double g) (double e))) (tolerance e))) (sut/coeffs r) exact)))
          (str "PolynomialR fit " xs))
    ;; the fitted polynomial passes through its points
    (t/is (every? true? (map (fn [x y] (<= (m/abs (- (sut/evaluate p x) y)) 1.0e-9)) xs ys))
          (str "fit passes through " xs)))
  (doseq [make [sut/polynomial sut/ratio-polynomial]]
    (t/is (thrown? IllegalArgumentException (make [0 1 1] [1 2 3])) "duplicate abscissae")
    (t/is (thrown? IllegalArgumentException (make [0 1 2] [1 2])) "different lengths")
    (t/is (thrown? IllegalArgumentException (make [2] [5])) "a single point")))

(t/deftest nominal-degree
  (doseq [make [sut/polynomial sut/ratio-polynomial]
          :let [p (make [1 2 3]) q (make [1 2 3 4 5])]]
    (t/is (= 2 (sut/degree (sut/sub p p))) "cancellation keeps the degree")
    (t/is (= 4 (sut/degree (sut/add p q))))
    (t/is (= 6 (sut/degree (sut/mult p q))))
    (t/is (= 6 (sut/degree (sut/mult (make [0 0 0]) q))) "zero coefficients keep the degree")
    (t/is (= 2 (sut/degree (sut/mult p (make [0])))))
    (t/is (= 2 (sut/degree (sut/derivative q 2))))
    (t/is (= 0 (sut/degree (make [5]))))))

(t/deftest polynomial-call-forms
  (let [p (sut/polynomial [1 2 3]) r (sut/ratio-polynomial [1 2 3])]
    (t/is (== 17.0 (p 2) (p 2.0) (sut/evaluate p 2) (sut/evaluate p 2.0)))
    (t/is (== 2.75 (p 1/2)))
    (t/is (= 17 (r 2)))
    (t/is (= 11/4 (r 1/2)))
    (t/is (= [6.0 17.0] (mapv p [1 2])))
    (t/is (= [6 17] (mapv r [1 2])))
    (t/is (== 17.0 (sut/evaluate r 2)))
    ;; `apply` with exactly one argument works; other argument counts are an arity error
    (t/is (== 17.0 (apply p [2])))
    (t/is (= 17 (apply r [2])))
    (t/is (= [6.0 17.0] (map #(apply p [%]) [1 2])))
    (doseq [poly [p r]
            args [[] [1 2] [1 2 3]]]
      (t/is (thrown? clojure.lang.ArityException (apply poly args)) (str (type-name poly) " with " (count args) " arguments")))))

;; Orthogonal polynomial families: shared checks.
;;
;; Each family has three implementations: `eval-*` (a recurrence or a closed form in doubles),
;; `*-ratio` (exact rational coefficients, `PolynomialR`) and the polynomial object (`Polynomial`,
;; double coefficients). `three-forms-disagreements` evaluates all three on the same points and
;; compares each with the exact value of the rational form, so a wrong implementation of one form
;; is found even when the other two are right.

(defn- exact-ratio
  "The exact rational value of a finite double."
  [x]
  (rationalize (java.math.BigDecimal. (double x))))

(defn- three-forms-disagreements
  "Points where one of the three forms of a family differs from the exact value of its rational form.

  Options: `:eval-fn` `(fn [args x])`, `:ratio-fn` and `:object-fn` `(fn [args])`, `:cases` (a
  collection of `args`), `:xs` (doubles), `:eval-units` and `:object-units` (allowed roundoff units).

  The `eval-*` function may differ by `eval-units * n * eps * (max(1, |value|) + |x * derivative|)`
  (backward stable: the exact value for an argument changed by a few roundoff units); the object, by
  `object-units * n * eps * sum |c_i| |x|^i` (Horner's bound). Values above 1e290 are skipped.
  Returns a sequence of maps, empty when the three forms agree."
  [{:keys [eval-fn ratio-fn object-fn cases xs eval-units object-units]
    :or {eval-units 8.0 object-units 4.0}}]
  (for [args cases
        :let [ratio-poly (ratio-fn args)
              object-poly (object-fn args)
              derivative-poly (sut/derivative ratio-poly 1)
              abs-cfs (abs-vec (sut/coeffs ratio-poly))
              n (count abs-cfs)]
        x xs
        :let [rx (exact-ratio x)
              exact (double (ratio-poly rx))
              exact-derivative (double (derivative-poly rx))]
        :when (< (m/abs exact) 1.0e290)
        :let [eval-value (eval-fn args x)
              object-value (object-poly x)
              eval-bound (* eval-units n EPS (+ (max 1.0 (m/abs exact)) (m/abs (* x exact-derivative))))
              object-bound (* object-units n EPS (error-scale abs-cfs (m/abs x)))]
        :when (not (and (<= (m/abs (- eval-value exact)) eval-bound)
                        (<= (m/abs (- object-value exact)) object-bound)))]
    {:args args :x x :exact exact :eval eval-value :object object-value
     :eval-bound eval-bound :object-bound object-bound}))

(t/deftest three-forms-helper-flags-disagreement
  (let [good {:eval-fn (fn [[n] x] (sut/eval-chebyshev-T n x))
              :ratio-fn (fn [[n]] (sut/chebyshev-T-ratio n))
              :object-fn (fn [[n]] (sut/chebyshev-T n))
              :cases [[3] [8]]
              :xs [0.3 -0.7 1.5]}]
    (t/is (empty? (three-forms-disagreements good)))
    ;; a recurrence that uses the wrong kind (the shape of the eval-bessel-y bug)
    (t/is (seq (three-forms-disagreements (assoc good :eval-fn (fn [[n] x] (sut/eval-chebyshev-U n x))))))
    ;; a result that is off by a relative 1e-9
    (t/is (seq (three-forms-disagreements
                (assoc good :eval-fn (fn [[n] x] (* (sut/eval-chebyshev-T n x) (+ 1.0 1.0e-9)))))))
    ;; a wrong polynomial object
    (t/is (seq (three-forms-disagreements (assoc good :object-fn (fn [[n]] (sut/chebyshev-U n))))))
    ;; a wrong rational form: the other two forms are then reported
    (t/is (seq (three-forms-disagreements (assoc good :ratio-fn (fn [[n]] (sut/chebyshev-U-ratio n))))))))

;; Chebyshev polynomials of the first to fourth kind: T, U, V, W.
;;
;; Reference values: `test/resources/polynomials/chebyshev_reference.edn`, exact rational
;; arithmetic (`utils/fastmath/dev/generate_chebyshev_reference.py`): integer coefficients of degree
;; 0 to 40; values and first derivatives at the exact binary value of 133 arguments (0, +-1,
;; +-(1 -+ delta) for many delta down to the neighbouring doubles of 1, moderate and large |x|) for
;; degrees up to 500; and large arguments with finite values.
;; `scipy.special.eval_chebyt/u` against the reference: largest difference 2.0e-12 (T) and 6.7e-13 (U)
;; relative to max(1, |value|), both at degree 200 next to x = -1 (the trigonometric form of scipy
;; amplifies the error of its argument); V and W are the exact identities U(n) -+ U(n-1), checked
;; against `mpmath` sums with difference 0.

(def ^:private chebyshev-reference
  (delay (edn/read-string (slurp (io/resource "polynomials/chebyshev_reference.edn")))))

(def ^:private chebyshev-forms
  {:T {:eval sut/eval-chebyshev-T :ratio sut/chebyshev-T-ratio :object sut/chebyshev-T}
   :U {:eval sut/eval-chebyshev-U :ratio sut/chebyshev-U-ratio :object sut/chebyshev-U}
   :V {:eval sut/eval-chebyshev-V :ratio sut/chebyshev-V-ratio :object sut/chebyshev-V}
   :W {:eval sut/eval-chebyshev-W :ratio sut/chebyshev-W-ratio :object sut/chebyshev-W}})

;; Roundoff units allowed per kind; the largest observed values on the reference grid are T 2.0, U 2.6,
;; V 7.7 and W 7.7 (V and W are the difference and the sum of two values of U; the difference of two
;; nearly equal values loses digits next to x = 1).
(def ^:private chebyshev-eval-units {:T 4.0 :U 4.0 :V 12.0 :W 12.0})

(defn- eval-error-bound
  "Allowed absolute error of an `eval-*` value `v` with the derivative `d` at `x` for degree `n`, in units of
  eps. The result is the exact value for an argument changed by a few units (the term with `x * d`) and
  the trigonometric forms add an error of about `n` units of the angle, relative to max(1, |v|), or to
  |v| (1 + log(2|x|)) outside [-1, 1], where the argument enters through an exponential."
  [n x v d units]
  (let [ax (m/abs x)
        scale (* (max 1.0 (m/abs v)) (if (<= ax 1.0) 1.0 (+ 1.0 (Math/log (* 2.0 ax)))))]
    (* units EPS (+ (* (inc n) scale) (m/abs (* x d))))))

(defn- reference-failures
  "Reference rows `[n x value derivative]` of `block` for which the kind's `eval-*` function is outside the
  bound; rows without a derivative (the large arguments) use zero."
  [kind block]
  (let [f (get-in chebyshev-forms [kind :eval])
        units (chebyshev-eval-units kind)]
    (for [[n x v d] (get-in @chebyshev-reference [block kind])
          :let [got (attempt f n x)]
          :when (not (and (not (failed? got))
                          (<= (m/abs (- got v)) (eval-error-bound n x v (or d 0.0) units))))]
      [n x v got])))

(t/deftest chebyshev-eval-reference
  (doseq [kind [:T :U :V :W]
          :let [failures (reference-failures kind :grid)]]
    (t/is (empty? failures)
          (str kind ": " (count failures) " grid points outside the bound; first (n x expected got): "
               (pr-str (take 4 failures))))))

(t/deftest chebyshev-large-arguments
  ;; values that are finite while an intermediate quantity of a naive formula overflows
  (doseq [kind [:T :U :V :W]
          :let [failures (reference-failures kind :large)]]
    (t/is (seq (get-in @chebyshev-reference [:large kind])))
    (t/is (empty? failures)
          (str kind ": " (count failures) " large arguments outside the bound; first (n x expected got): "
               (pr-str (take 4 failures))))))

(t/deftest chebyshev-endpoints-are-exact
  ;; T(n, +-1) = (+-1)^n, U(n, +-1) = (+-1)^n (n+1), V(n, 1) = 1, V(n, -1) = (-1)^n (2n+1),
  ;; W(n, 1) = 2n+1, W(n, -1) = (-1)^n
  (doseq [n [0 1 2 3 4 5 6 7 8 11 20 51 100 500]
          :let [s (if (even? n) 1.0 -1.0)]]
    (t/is (== 1.0 (sut/eval-chebyshev-T n 1.0)) (str "T " n))
    (t/is (== s (sut/eval-chebyshev-T n -1.0)) (str "T " n))
    (t/is (== (+ n 1.0) (sut/eval-chebyshev-U n 1.0)) (str "U " n))
    (t/is (== (* s (+ n 1.0)) (sut/eval-chebyshev-U n -1.0)) (str "U " n))
    (t/is (== 1.0 (sut/eval-chebyshev-V n 1.0)) (str "V " n))
    (t/is (== (* s (+ (* 2.0 n) 1.0)) (sut/eval-chebyshev-V n -1.0)) (str "V " n))
    (t/is (== (+ (* 2.0 n) 1.0) (sut/eval-chebyshev-W n 1.0)) (str "W " n))
    (t/is (== s (sut/eval-chebyshev-W n -1.0)) (str "W " n))))

(t/deftest chebyshev-degree-zero-is-exactly-one
  (doseq [f [sut/eval-chebyshev-T sut/eval-chebyshev-U sut/eval-chebyshev-V sut/eval-chebyshev-W]
          x [0.3 -0.9 0.0 1.0 -1.0 2.5 -7.0 ##Inf ##-Inf ##NaN]]
    (t/is (= (bits 1.0) (bits (f 0 x))) (str "degree 0 at " x))))

(t/deftest chebyshev-non-finite-arguments
  (doseq [[kind {f :eval}] chebyshev-forms
          n [1 2 3 4 5 6 7 20]
          :let [odd-degree? (odd? n)]]
    (t/is (m/nan? (f n ##NaN)) (str kind " " n " at NaN"))
    ;; the leading coefficients of all four kinds are positive: the sign of the infinity follows the parity
    (t/is (== ##Inf (f n ##Inf)) (str kind " " n " at +Inf"))
    (t/is (== (if odd-degree? ##-Inf ##Inf) (f n ##-Inf)) (str kind " " n " at -Inf"))))

(t/deftest chebyshev-three-forms-agree
  (doseq [[kind {:keys [eval ratio object]}] chebyshev-forms
          :let [disagreements (three-forms-disagreements
                               {:eval-fn (fn [[n] x] (eval n x))
                                :ratio-fn (fn [[n]] (ratio n))
                                :object-fn (fn [[n]] (object n))
                                :cases (map vector [0 1 2 3 4 5 6 7 8 9 10 12 15 18])
                                :xs [-3.0 -1.5 -1.0 -0.999999 -0.9 -0.5 -0.1 0.0 0.1 0.5 0.9 0.999999 1.0 1.5 3.0]})]]
    (t/is (empty? disagreements)
          (str kind ": " (count disagreements) " disagreements; first: " (pr-str (take 2 disagreements))))))

(t/deftest chebyshev-exact-coefficients
  (doseq [[kind {:keys [ratio object]}] chebyshev-forms
          [n exact] (map-indexed vector (get-in @chebyshev-reference [:coefficients kind]))
          :let [r (attempt ratio n)
                o (attempt object n)]]
    (t/is (and (not (failed? r)) (= exact (vec (sut/coeffs r))) (= n (sut/degree r)))
          (str kind " ratio, degree " n))
    (t/is (and (not (failed? o)) (= (mapv exact->double exact) (vec (sut/coeffs o))) (= n (sut/degree o)))
          (str kind " object, degree " n))))

(t/deftest chebyshev-negative-degree
  (doseq [[kind {:keys [eval ratio object]}] chebyshev-forms
          n [-1 -2 -100]]
    (t/is (thrown? IllegalArgumentException (eval n 0.3)) (str kind " eval " n))
    (t/is (thrown? IllegalArgumentException (ratio n)) (str kind " ratio " n))
    (t/is (thrown? IllegalArgumentException (object n)) (str kind " object " n))))

;; Legendre, Gegenbauer (ultraspherical) and Jacobi polynomials.
;;
;; Reference values: `test/resources/polynomials/legendre_gegenbauer_jacobi_reference.edn`, exact
;; rational arithmetic (`utils/fastmath/dev/generate_legendre_gegenbauer_jacobi_reference.py`) at the
;; exact binary values of the parameters and arguments: 17 values of the Gegenbauer parameter
;; (including 0, negative values and the neighbours 1 +- 2^-30, 0.5 +- 2^-30 of the two shortcuts) and 20
;; Jacobi parameter pairs (including those for which the three term recurrence divides by zero or is
;; close to doing so), 27 arguments around +-1, degrees up to 30 (Legendre up to 500); and the exact
;; coefficients. Jacobi values come from the explicit sum, which holds for every real parameters.
;; `scipy.special.eval_legendre/gegenbauer/jacobi` against the reference: largest difference 2.6e-11
;; (Legendre, degree 500 outside [-1, 1]), 1.8e-14 (Gegenbauer) and 3.7e-9 (Jacobi next to a degenerate
;; pair), relative to max(1, |value|); scipy returns NaN for 942 Gegenbauer and 2088 Jacobi rows of
;; negative parameters, and 0 for degree 0 and parameter 0 where the polynomial is 1.

(def ^:private lgj-reference
  (delay (edn/read-string (slurp (io/resource "polynomials/legendre_gegenbauer_jacobi_reference.edn")))))

;; Observed largest error over the grid (in the units of `eval-error-bound`): Legendre 0.7, Gegenbauer 1.4,
;; Jacobi 1.1 apart from the pair below.
(def ^:private lgj-eval-units {:legendre 4.0 :gegenbauer 4.0 :jacobi 4.0})


(defn- grid-failures
  "Rows `[n x value derivative]` for which `(f n x)` is outside the bound of `eval-error-bound`."
  [f rows units]
  (for [[n x v d] rows
        :let [got (attempt f n x)]
        :when (not (and (not (failed? got)) (<= (m/abs (- got v)) (eval-error-bound n x v d units))))]
    [n x v got]))

(defn- failures-message [label failures]
  (str label ": " (count failures) " grid points outside the bound; first (n x expected got): "
       (pr-str (take 4 failures))))

(t/deftest legendre-eval-reference
  (let [failures (grid-failures sut/eval-legendre-P (get-in @lgj-reference [:legendre :grid]) (:legendre lgj-eval-units))]
    (t/is (empty? failures) (failures-message "Legendre" failures))))

(t/deftest gegenbauer-eval-reference
  (doseq [{:keys [alpha grid]} (:gegenbauer @lgj-reference)
          :let [failures (grid-failures (fn [n x] (sut/eval-gegenbauer-C n alpha x)) grid (:gegenbauer lgj-eval-units))]]
    (t/is (empty? failures) (failures-message (str "Gegenbauer alpha " alpha) failures))))

(t/deftest jacobi-eval-reference
  (doseq [{:keys [alpha beta grid]} (:jacobi @lgj-reference)
          :let [failures (grid-failures (fn [n x] (sut/eval-jacobi-P n alpha beta x)) grid (:jacobi lgj-eval-units))]]
    (t/is (empty? failures) (failures-message (str "Jacobi alpha " alpha " beta " beta) failures))))

(defn- exact-coefficients [pairs] (mapv ratio pairs))

(t/deftest legendre-exact-coefficients
  (doseq [[n exact] (map-indexed vector (get-in @lgj-reference [:legendre :coefficients]))
          :let [r (attempt sut/legendre-P-ratio n)
                o (attempt sut/legendre-P n)
                expected (exact-coefficients exact)]]
    (t/is (and (not (failed? r)) (= expected (vec (sut/coeffs r))) (= n (sut/degree r))) (str "ratio, degree " n))
    (t/is (and (not (failed? o)) (= (mapv exact->double expected) (vec (sut/coeffs o))) (= n (sut/degree o))) (str "object, degree " n))))

(t/deftest gegenbauer-exact-coefficients
  (doseq [{:keys [alpha decimal-exact? coefficients]} (:gegenbauer @lgj-reference)
          :when decimal-exact?
          [n exact] (map-indexed vector coefficients)
          :let [r (attempt sut/gegenbauer-C-ratio n alpha)
                o (attempt sut/gegenbauer-C n alpha)
                expected (exact-coefficients exact)]]
    (t/is (and (not (failed? r)) (= expected (vec (sut/coeffs r))) (= n (sut/degree r)))
          (str "ratio, alpha " alpha " degree " n))
    (t/is (and (not (failed? o)) (= (mapv exact->double expected) (vec (sut/coeffs o))) (= n (sut/degree o)))
          (str "object, alpha " alpha " degree " n))))

(t/deftest jacobi-exact-coefficients
  (doseq [{:keys [alpha beta decimal-exact? coefficients]} (:jacobi @lgj-reference)
          :when decimal-exact?
          [n exact] (map-indexed vector coefficients)
          :let [r (attempt sut/jacobi-P-ratio n alpha beta)
                o (attempt sut/jacobi-P n alpha beta)
                expected (exact-coefficients exact)]]
    (t/is (and (not (failed? r)) (= expected (vec (sut/coeffs r))) (= n (sut/degree r)))
          (str "ratio, alpha " alpha " beta " beta " degree " n))
    (t/is (and (not (failed? o)) (= (mapv exact->double expected) (vec (sut/coeffs o))) (= n (sut/degree o)))
          (str "object, alpha " alpha " beta " beta " degree " n))))

(defn- pochhammer [a k] (reduce *' 1 (map #(+ a %) (range k))))

(defn- factorial-exact [n] (reduce *' 1 (range 1 (inc n))))

(defn- gegenbauer-explicit-coefficients
  "Coefficients, ascending, of the Gegenbauer polynomial of a rational parameter from the explicit sum
  `sum_k (-1)^k (alpha)_(n-k) / (k! (n-2k)!) (2x)^(n-2k)`."
  [n alpha]
  (reduce (fn [cs k]
            (let [m (- n (* 2 k))]
              (assoc cs m (*' (if (even? k) 1 -1) (pochhammer alpha (- n k))
                              (/ (reduce *' 1 (repeat m 2)) (*' (factorial-exact k) (factorial-exact m)))))))
          (vec (repeat (inc n) 0))
          (range (inc (quot n 2)))))

(defn- jacobi-explicit-ratio
  "The Jacobi polynomial of rational parameters as the explicit sum
  `sum_s C(n+a, n-s) C(n+b, s) ((x-1)/2)^s ((x+1)/2)^(n-s)`, built with the exact polynomial operations."
  [n alpha beta]
  (let [binomial (fn [top k] (/ (reduce *' 1 (map #(- top %) (range k))) (factorial-exact k)))
        powers (fn [base] (vec (take (inc n) (iterate #(sut/mult % base) (sut/ratio-polynomial [1])))))
        lower (powers (sut/ratio-polynomial [-1/2 1/2]))
        upper (powers (sut/ratio-polynomial [1/2 1/2]))]
    (reduce sut/add
            (for [s (range (inc n))]
              (sut/scale (sut/mult (lower s) (upper (- n s)))
                         (* (binomial (+ n alpha) (- n s)) (binomial (+ n beta) s)))))))

(t/deftest ratio-forms-are-exact-for-decimal-parameters
  ;; the parameter is converted with rationalize (0.3 is 3/10): the coefficients are then exact
  (doseq [alpha [0.1 0.3 1.7 0.123 -0.6] n (range 0 13)]
    (t/is (= (gegenbauer-explicit-coefficients n (rationalize alpha)) (vec (sut/coeffs (sut/gegenbauer-C-ratio n alpha))))
          (str "Gegenbauer alpha " alpha " degree " n)))
  (doseq [[alpha beta] [[0.1 0.2] [1.7 -0.3] [-0.6 2.4] [0.3 0.3]] n (range 0 9)]
    (t/is (= (vec (sut/coeffs (jacobi-explicit-ratio n (rationalize alpha) (rationalize beta))))
             (vec (sut/coeffs (sut/jacobi-P-ratio n alpha beta))))
          (str "Jacobi alpha " alpha " beta " beta " degree " n))))

(t/deftest legendre-gegenbauer-jacobi-three-forms-agree
  (let [xs [-3.0 -1.5 -1.0 -0.9 -0.5 0.0 0.5 0.9 1.0 1.5 3.0]
        run (fn [eval-fn ratio-fn object-fn cases]
              (three-forms-disagreements {:eval-fn eval-fn :ratio-fn ratio-fn :object-fn object-fn
                                          :cases cases :xs xs}))
        ns [0 1 2 3 4 5 6 8 10 15]]
    (let [d (run (fn [[n] x] (sut/eval-legendre-P n x)) (fn [[n]] (sut/legendre-P-ratio n)) (fn [[n]] (sut/legendre-P n))
                 (map vector ns))]
      (t/is (empty? d) (str "Legendre: " (pr-str (take 2 d)))))
    (let [d (run (fn [[n a] x] (sut/eval-gegenbauer-C n a x)) (fn [[n a]] (sut/gegenbauer-C-ratio n a))
                 (fn [[n a]] (sut/gegenbauer-C n a))
                 (for [n ns a [0.25 0.5 0.75 1.0 1.5 2.0 5.5 -0.25 -0.5 -1.5 0.3 1.7]] [n a]))]
      (t/is (empty? d) (str "Gegenbauer: " (count d) " disagreements; first: " (pr-str (take 2 d)))))
    (let [d (run (fn [[n a b] x] (sut/eval-jacobi-P n a b x)) (fn [[n a b]] (sut/jacobi-P-ratio n a b))
                 (fn [[n a b]] (sut/jacobi-P n a b))
                 (for [n (remove #{15} ns) [a b] [[0.0 0.0] [0.5 -0.5] [1.5 2.5] [-0.5 0.25] [0.1 0.2] [1.7 -0.3]
                                                  [-1.0 -1.0] [-1.0 -2.0] [-2.0 0.0] [-0.5 -1.5] [-1.0 3.0]]]
                   [n a b]))]
      (t/is (empty? d) (str "Jacobi: " (count d) " disagreements; first: " (pr-str (take 2 d)))))))

(t/deftest gegenbauer-parameter-conventions
  ;; parameter 0: the polynomial of degree 0 is 1, those of higher degree are identically 0 (as in scipy
  ;; and mpmath)
  (doseq [x [-3.0 -0.5 0.0 0.3 1.0 2.0]]
    (t/is (== 1.0 (sut/eval-gegenbauer-C 0 0.0 x)))
    (doseq [n [1 2 3 4 7 20]]
      (t/is (== 0.0 (sut/eval-gegenbauer-C n 0.0 x)) (str "degree " n " at " x))))
  (doseq [n [1 2 5]]
    (t/is (every? zero? (sut/coeffs (sut/gegenbauer-C-ratio n 0.0))) "zero coefficients"))
  (t/is (= [1] (vec (sut/coeffs (sut/gegenbauer-C-ratio 0 0.0)))))
  ;; the two shortcuts are the same polynomials as the Chebyshev U and the Legendre ones
  (doseq [n [0 1 2 3 4 5 6 10 31] x [-2.0 -0.9 0.0 0.3 0.999999 1.0 2.5]]
    (t/is (= (bits (sut/eval-chebyshev-U n x)) (bits (sut/eval-gegenbauer-C n 1.0 x)) (bits (sut/eval-gegenbauer-C n x))))
    (t/is (= (bits (sut/eval-legendre-P n x)) (bits (sut/eval-gegenbauer-C n 0.5 x)))))
  (t/is (= (sut/chebyshev-U-ratio 6) (sut/gegenbauer-C-ratio 6 1.0)))
  (t/is (= (sut/legendre-P-ratio 6) (sut/gegenbauer-C-ratio 6 0.5)))
  ;; a NaN parameter
  (t/is (m/nan? (sut/eval-gegenbauer-C 3 ##NaN 0.3)))
  (t/is (== 1.0 (sut/eval-gegenbauer-C 0 ##NaN 0.3)))
  (t/is (thrown? IllegalArgumentException (sut/gegenbauer-C-ratio 3 ##NaN)))
  (t/is (thrown? IllegalArgumentException (sut/gegenbauer-C-ratio 3 ##Inf))))

(t/deftest jacobi-degenerate-parameters
  ;; the three term recurrence divides by zero where alpha + beta is a negative integer of at most -2;
  ;; the values are those of the explicit sum
  (t/is (m/delta-eq -0.2275 (sut/eval-jacobi-P 2 -1.0 -1.0 0.3) 1.0e-15) "scipy: eval_jacobi(2, -1, -1, 0.3)")
  (t/is (m/delta-eq -0.2275 (sut/evaluate (sut/jacobi-P 2 -1.0 -1.0) 0.3) 1.0e-15))
  (t/is (= -91/400 ((sut/jacobi-P-ratio 2 -1.0 -1.0) 0.3)))
  (t/is (== 0.0 (sut/eval-jacobi-P 5 -3.0 -4.0 0.3)) "the explicit sum has no non-zero term")
  (t/is (every? zero? (sut/coeffs (sut/jacobi-P-ratio 5 -3.0 -4.0))))
  (doseq [make [(fn [n] (sut/jacobi-P-ratio n -1.0 -1.0)) (fn [n] (sut/jacobi-P n -1.0 -1.0))]
          n [2 3 4 7 12]]
    (t/is (not-any? #(or (m/nan? %) (m/inf? %)) (map double (sut/coeffs (make n)))) (str "degree " n)))
  ;; close to a degenerate pair the three term recurrence would lose digits (relative error about
  ;; 3e-16 / delta); the grid contains alpha + beta = -2 + 2^-10, -2 + 2^-20 and -2 + 2^-30
  (doseq [delta [1.0e-2 1.0e-4 1.0e-6 1.0e-9 1.0e-12]
          :let [alpha (+ -1.0 (* 0.5 delta))]
          n [2 3 4 6 8]
          x [-0.5 0.3 0.9]]
    (let [exact (double ((jacobi-explicit-ratio n (exact-ratio alpha) (exact-ratio alpha)) (exact-ratio x)))
          got (sut/eval-jacobi-P n alpha alpha x)]
      (t/is (<= (m/abs (- got exact)) (* 16.0 (inc n) EPS (max 1.0 (m/abs exact))))
            (str "delta " delta " degree " n " at " x ": got " got " expected " exact))))
  (t/is (m/nan? (sut/eval-jacobi-P 3 ##NaN 0.5 0.3)))
  (t/is (m/nan? (sut/eval-jacobi-P 3 0.5 ##NaN 0.3)))
  (t/is (m/nan? (sut/eval-jacobi-P 3 -1.0 -1.0 ##NaN)))
  (t/is (thrown? IllegalArgumentException (sut/jacobi-P-ratio 3 ##NaN 0.5)))
  (t/is (thrown? IllegalArgumentException (sut/jacobi-P-ratio 3 0.5 ##Inf))))

(defn- value-at-infinity
  "The value of the polynomial with the ascending exact coefficients `cs` at an infinite `x`: the limit of its
  leading non-zero term (a finite constant for a constant polynomial, 0 for the zero polynomial)."
  [cs x]
  (let [k (last (keep-indexed (fn [i c] (when-not (zero? c) i)) cs))]
    (cond (nil? k) 0.0
          (zero? k) (double (cs 0))
          :else (* (if (pos? (cs k)) 1.0 -1.0)
                   (if (and (neg? x) (odd? k)) -1.0 1.0)
                   ##Inf))))

(t/deftest legendre-gegenbauer-jacobi-non-finite-arguments
  (let [reference @lgj-reference]
    ;; infinite argument: the limit of the leading non-zero term of the exact coefficients
    (doseq [x [##Inf ##-Inf]
            [n exact] (map-indexed vector (get-in reference [:legendre :coefficients]))
            :when (<= n 12)]
      (t/is (== (value-at-infinity (exact-coefficients exact) x) (sut/eval-legendre-P n x)) (str "Legendre " n " at " x)))
    (doseq [x [##Inf ##-Inf]
            {:keys [alpha decimal-exact? coefficients]} (:gegenbauer reference)
            :when decimal-exact?
            [n exact] (map-indexed vector coefficients)
            :when (<= n 12)]
      (t/is (== (value-at-infinity (exact-coefficients exact) x) (sut/eval-gegenbauer-C n alpha x))
            (str "Gegenbauer alpha " alpha " degree " n " at " x)))
    (doseq [x [##Inf ##-Inf]
            {:keys [alpha beta decimal-exact? coefficients]} (:jacobi reference)
            :when decimal-exact?
            [n exact] (map-indexed vector coefficients)
            :when (<= n 12)]
      (t/is (== (value-at-infinity (exact-coefficients exact) x) (sut/eval-jacobi-P n alpha beta x))
            (str "Jacobi alpha " alpha " beta " beta " degree " n " at " x))))
  ;; NaN argument
  (doseq [n [1 2 3 7]]
    (t/is (m/nan? (sut/eval-legendre-P n ##NaN)))
    (t/is (m/nan? (sut/eval-gegenbauer-C n 2.5 ##NaN)))
    (t/is (m/nan? (sut/eval-jacobi-P n 0.5 1.5 ##NaN))))
  ;; degree 0 is 1 for every argument, degree 1 is linear
  (doseq [x [-2.0 0.3 ##Inf ##NaN]]
    (t/is (== 1.0 (sut/eval-legendre-P 0 x) (sut/eval-gegenbauer-C 0 2.5 x) (sut/eval-jacobi-P 0 0.5 1.5 x))))
  (t/is (== 0.3 (sut/eval-legendre-P 1 0.3)))
  (t/is (== 1.5 (sut/eval-gegenbauer-C 1 2.5 0.3)))
  (t/is (== 1.5 (sut/eval-jacobi-P 1 0.5 1.5 1.0)) "(alpha + 1) + (alpha + beta + 2)/2 (x - 1) at x = 1")
  (t/is (== -0.5 (sut/eval-jacobi-P 1 0.5 1.5 0.0)))
  ;; huge argument with a representable value
  (t/is (m/delta-eq 1.0 (/ (sut/eval-legendre-P 5 1.0e60) 7.875e300) 1.0e-14))
  (t/is (m/delta-eq 1.0 (/ (sut/eval-legendre-P 10 1.0e30) 1.8042578125e302) 1.0e-14)))

(t/deftest legendre-gegenbauer-jacobi-negative-degree
  (doseq [n [-1 -2 -100]]
    (doseq [[label f] {"eval-legendre-P" #(sut/eval-legendre-P % 0.3)
                       "legendre-P-ratio" sut/legendre-P-ratio
                       "legendre-P" sut/legendre-P
                       "eval-gegenbauer-C, default parameter" #(sut/eval-gegenbauer-C % 0.3)
                       "eval-gegenbauer-C, parameter 1" #(sut/eval-gegenbauer-C % 1.0 0.3)
                       "eval-gegenbauer-C, parameter 0.5" #(sut/eval-gegenbauer-C % 0.5 0.3)
                       "eval-gegenbauer-C, parameter 2.5" #(sut/eval-gegenbauer-C % 2.5 0.3)
                       "gegenbauer-C-ratio, parameter 1" #(sut/gegenbauer-C-ratio % 1.0)
                       "gegenbauer-C-ratio, parameter 0.5" #(sut/gegenbauer-C-ratio % 0.5)
                       "gegenbauer-C-ratio, parameter 2.5" #(sut/gegenbauer-C-ratio % 2.5)
                       "gegenbauer-C, default parameter" sut/gegenbauer-C
                       "gegenbauer-C, parameter 2.5" #(sut/gegenbauer-C % 2.5)
                       "eval-jacobi-P" #(sut/eval-jacobi-P % 0.5 1.5 0.3)
                       "jacobi-P-ratio" #(sut/jacobi-P-ratio % 0.5 1.5)
                       "jacobi-P" #(sut/jacobi-P % 0.5 1.5)}]
      (t/is (thrown? IllegalArgumentException (f n)) (str label " " n)))))

;; Gegenbauer and Jacobi outside the range of the three term recurrence, and at infinity.
;;
;; The recurrence is accurate for a Gegenbauer order above -1 and for Jacobi `alpha > -1`, `beta > -1`,
;; `alpha + beta > -1.9` (measured: at most 3.9 and 2.5 units on random parameters). Below that it lost up to
;; 1e14 units (an identically zero polynomial gave noise of 2.7e6); there the evaluators sum the explicit
;; formula in decimal arithmetic of adaptive precision, so the result is the nearest double up to about 2^-60.

(defn- horner-exact [cs rx] (reduce (fn [acc c] (+ (* acc rx) c)) 0 (reverse cs)))

(defn- random-x
  "A random argument: wide, inside [-1, 1], next to ±1 on both sides, or tiny."
  [^java.util.Random rng]
  (let [sign (if (.nextBoolean rng) 1.0 -1.0)
        exponent (fn [lo span] (m/pow 10.0 (- (+ lo (* span (.nextDouble rng))))))]
    (case (.nextInt rng 5)
      0 (- (* 6.0 (.nextDouble rng)) 3.0)
      1 (- (* 2.0 (.nextDouble rng)) 1.0)
      2 (* sign (- 1.0 (exponent 1.0 14.0)))
      3 (* sign (exponent 1.0 11.0))
      4 (* sign (+ 1.0 (exponent 1.0 11.0))))))

(defn- correct-to-rounding?
  "True when `got` is the nearest double of `exact` up to 4 units of roundoff, relative to `exact` (exactly 0
  for an exact 0)."
  [got exact]
  (<= (m/abs (- got exact)) (* 4.0 EPS (m/abs exact))))

(t/deftest gegenbauer-negative-integer-order-is-identically-zero
  ;; C_n^(-k) = 0 for n > 2k (the exact coefficients are all zero); the recurrence gave noise for |x| > 1
  (doseq [[k n x] [[2 16 1.3] [2 23 5.856421080721988] [3 20 2.1586225252597906] [3 36 2.1586225252597906]
                   [5 28 5.7] [1 12 5.7] [4 40 0.3] [6 60 -7.5]]]
    (t/is (every? zero? (sut/coeffs (sut/gegenbauer-C-ratio n (double (- k))))))
    (t/is (== 0.0 (sut/eval-gegenbauer-C n (double (- k)) x)) (str "order " (- k) " degree " n " at " x)))
  ;; degree 2k is the constant 1 (also in the limit)
  (doseq [k [1 2 3 5] x [-4.5 0.25 7.0 ##Inf ##-Inf]]
    (t/is (== 1.0 (sut/eval-gegenbauer-C (* 2 k) (double (- k)) x)) (str "order " (- k) " at " x))))

(t/deftest gegenbauer-order-below-minus-one-is-correct-to-rounding
  (let [rng (java.util.Random. 20261010)]
    (dotimes [_ 300]
      (let [n (+ 2 (.nextInt rng 29))
            x (random-x rng)
            alpha (case (.nextInt rng 3)
                    0 (- (/ (+ 8 (.nextInt rng 24)) 8.0))
                    1 (- (double (inc (.nextInt rng 6))))
                    2 (+ (- (double (inc (.nextInt rng 6)))) (m/pow 2.0 (- (+ 5 (.nextInt rng 20))))))
            exact (exact->double (horner-exact (gegenbauer-explicit-coefficients n (exact-ratio alpha)) (exact-ratio x)))
            got (sut/eval-gegenbauer-C n alpha x)]
        (t/is (correct-to-rounding? got exact) (str "degree " n " order " alpha " at " x ": got " got " expected " exact))))))

(t/deftest jacobi-outside-the-recurrence-range-is-correct-to-rounding
  ;; the pairs of the recurrence range keep the looser bound relative to the largest value of degrees 0 to n
  (let [rng (java.util.Random. 20261011)
        recurrence-range? (fn [a b] (and (> a -1.0) (> b -1.0) (> (+ a b) -1.9)))]
    (dotimes [_ 250]
      (let [n (+ 2 (.nextInt rng 15))
            x (random-x rng)
            alpha (/ (- (.nextInt rng 33) 24) 8.0)
            ;; a quarter of the pairs have alpha + beta a negative integer of -2 and below (recurrence divides by 0)
            beta (if (zero? (.nextInt rng 4))
                   (- (- (double (+ 2 (.nextInt rng 8)))) alpha)
                   (/ (- (.nextInt rng 33) 24) 8.0))
            exact (exact->double (horner-exact (sut/coeffs (jacobi-explicit-ratio n (exact-ratio alpha) (exact-ratio beta))) (exact-ratio x)))
            got (sut/eval-jacobi-P n alpha beta x)]
        (if (recurrence-range? alpha beta)
          (let [envelope (reduce max 1.0 (map #(m/abs (sut/eval-jacobi-P % alpha beta x)) (range (inc n))))]
            (t/is (<= (m/abs (- got exact)) (* 16.0 (inc n) EPS (max 1.0 (m/abs exact) envelope)))
                  (str "recurrence range, degree " n " (" alpha ", " beta ") at " x)))
          (t/is (correct-to-rounding? got exact) (str "degree " n " (" alpha ", " beta ") at " x ": got " got " expected " exact)))))))

(t/deftest gegenbauer-small-order-keeps-its-relative-accuracy
  ;; the factors (a + i - 1) and (i - 2 + 2a) of the recurrence were formed as (a + i) - 1 and (i + 2a) - 2: a
  ;; relative error of eps / |a| (6e-9 for a = 1e-8) in every value, which is proportional to a
  (doseq [alpha [1.0e-8 -1.0e-8 1.0e-12 -1.0e-3 0.001 -0.5]
          [n x] [[2 0.3] [5 0.3] [12 0.3] [20 -2.5] [20 0.99] [30 0.7] [30 -1.0e-3] [8 4.5]]
          :let [exact (exact->double (horner-exact (gegenbauer-explicit-coefficients n (exact-ratio alpha)) (exact-ratio x)))
                got (sut/eval-gegenbauer-C n alpha x)]]
    (t/is (<= (m/abs (- got exact)) (* 16.0 (inc n) EPS (m/abs exact)))
          (str "degree " n " order " alpha " at " x ": got " got " expected " exact))))

(t/deftest jacobi-and-gegenbauer-below-minus-one-reproducers
  (doseq [[n a b x] [[25 0.375 -4.25 -1.0508027409467864] [18 -2.75 -1.375 1.087223071601323]
                     [23 -1.5 -2.625 -0.9999999294905596] [30 -2.75 -1.375 0.99999]]]
    (t/is (correct-to-rounding? (sut/eval-jacobi-P n a b x)
                                (exact->double (horner-exact (sut/coeffs (jacobi-explicit-ratio n (exact-ratio a) (exact-ratio b))) (exact-ratio x))))
          (str n " " a " " b " " x)))
  (doseq [[n a x] [[37 -9.125 -1.7388703311057259] [33 -5.999999940395355 -2.876693657987249]]]
    (t/is (correct-to-rounding? (sut/eval-gegenbauer-C n a x)
                                (exact->double (horner-exact (gegenbauer-explicit-coefficients n (exact-ratio a)) (exact-ratio x))))
          (str n " " a " " x))))

(t/deftest rational->double-rounds-to-nearest-and-ties-to-even
  ;; 1 + (2k+1)/2^53 lies exactly between the doubles 1 + k/2^52 and 1 + (k+1)/2^52: the even one wins
  (doseq [k (range 0 200)
          :let [tie (+ 1 (/ (+ (* 2 k) 1) 9007199254740992N))
                lower (+ 1.0 (* k (m/pow 2.0 -52)))
                upper (+ 1.0 (* (inc k) (m/pow 2.0 -52)))]]
    (t/is (== (if (even? k) lower upper) (exact->double tie)) (str "tie " k))
    (t/is (== (- (if (even? k) lower upper)) (exact->double (- tie))) (str "negative tie " k))
    ;; a hair above or below the tie goes to the nearer double
    (t/is (== upper (exact->double (+ tie (/ 1 (.pow (BigInteger/valueOf 2) 100))))) (str "above " k))
    (t/is (== lower (exact->double (- tie (/ 1 (.pow (BigInteger/valueOf 2) 100))))) (str "below " k)))
  ;; the subnormal range: the unit is 2^-1074 = MIN_VALUE
  (let [two-1075 (.pow (BigInteger/valueOf 2) 1075)]
    (t/is (== 0.0 (exact->double (/ 1 two-1075))) "0.5 MIN_VALUE: tie to 0")
    (t/is (== (* 2 Double/MIN_VALUE) (exact->double (/ 3 two-1075))) "1.5 MIN_VALUE: tie to 2 MIN_VALUE")
    (t/is (== (* 2 Double/MIN_VALUE) (exact->double (/ 5 two-1075))) "2.5 MIN_VALUE: tie to 2 MIN_VALUE")
    (t/is (== (* 4 Double/MIN_VALUE) (exact->double (/ 7 two-1075))) "3.5 MIN_VALUE: tie to 4 MIN_VALUE"))
  ;; range and exact values
  (t/is (== 0.3333333333333333 (exact->double 1/3)))
  (t/is (== -0.3333333333333333 (exact->double -1/3)))
  (t/is (== ##Inf (exact->double (.pow (BigInteger/valueOf 10) 400))))
  (t/is (== 0.0 (exact->double (/ 1 (.pow (BigInteger/valueOf 10) 400)))))
  (t/is (== Double/MAX_VALUE (exact->double (* (bigint 9007199254740991) (.pow (BigInteger/valueOf 2) 971))))))

(defmacro ^:private finishes-within
  "True when the form returns within `ms` milliseconds (a form that does not stop keeps its thread busy)."
  [ms form]
  `(not= ::timeout (deref (future ~form) ~ms ::timeout)))

(t/deftest high-precision-paths-do-not-depend-on-the-size-of-the-argument
  ;; the exact binary value of a tiny x has a denominator of 2^1074; the decimal sum does not care
  (doseq [x [1.0e-300 Double/MIN_VALUE -1.0e-100 1.0e300 Double/MAX_VALUE]]
    (t/is (finishes-within 5000 (sut/eval-jacobi-P 100 -1.0 -1.0 x)) (str "Jacobi at " x))
    (t/is (finishes-within 5000 (sut/eval-gegenbauer-C 100 -3.5 x)) (str "Gegenbauer at " x)))
  (t/is (finishes-within 5000 (sut/eval-jacobi-P 300 -1.0 -1.0 Double/MIN_VALUE)))
  (t/is (finishes-within 5000 (sut/eval-jacobi-P 1000 -1.0 -1.0 0.3)))
  ;; the values themselves at tiny x: P_n(0) and neighbours are finite and equal the exact ones
  (doseq [x [1.0e-300 Double/MIN_VALUE]]
    (let [exact (exact->double (horner-exact (sut/coeffs (jacobi-explicit-ratio 12 -1 -1)) (exact-ratio x)))]
      (t/is (correct-to-rounding? (sut/eval-jacobi-P 12 -1.0 -1.0 x) exact) (str "x = " x)))))

(t/deftest infinite-argument-is-found-from-the-signs
  (t/is (finishes-within 5000 (sut/eval-gegenbauer-C 100000 2.5 ##Inf)))
  (t/is (finishes-within 5000 (sut/eval-gegenbauer-C 2000 -2.5 ##-Inf)))
  (t/is (finishes-within 5000 (sut/eval-jacobi-P 100000 0.5 1.5 ##Inf)))
  (t/is (finishes-within 5000 (sut/eval-jacobi-P 100000 -3.0 0.5 ##-Inf)))
  ;; the limit of the leading non-zero term of the exact polynomial, for orders and parameters in steps of 1/4
  ;; including the negative integers where terms vanish
  (doseq [x [##Inf ##-Inf]
          n (range 1 15)]
    (doseq [k (range -32 49 3)
            :let [alpha (/ k 4.0)
                  cs (gegenbauer-explicit-coefficients n (exact-ratio alpha))
                  expected (value-at-infinity cs x)
                  got (sut/eval-gegenbauer-C n alpha x)]]
      (t/is (== expected got) (str "Gegenbauer order " alpha " degree " n " at " x)))
    (doseq [ka (range -12 9 2) kb (range -12 9 3)
            :let [alpha (/ ka 2.0) beta (/ kb 2.0) n (min n 10)
                  cs (sut/coeffs (jacobi-explicit-ratio n (exact-ratio alpha) (exact-ratio beta)))
                  expected (value-at-infinity cs x)
                  got (sut/eval-jacobi-P n alpha beta x)]]
      ;; a finite limit is a constant; `value-at-infinity` converts it with `double` (8 units of roundoff)
      (t/is (if (m/inf? expected) (== expected got) (<= (m/abs (- got expected)) (* 8.0 EPS (m/abs expected))))
            (str "Jacobi (" alpha ", " beta ") degree " n " at " x)))))

;; The evaluators dispatch on `(int degree)`: from 2^32 the degree wrapped around to a small one (degree
;; 2^32 gave the value of degree 0). A degree not below `Integer/MAX_VALUE` is rejected.
(t/deftest degree-must-be-below-int-max-value
  (doseq [n [Integer/MAX_VALUE (inc Integer/MAX_VALUE) 4294967296 4294967297 Long/MAX_VALUE]]
    (doseq [[kind {:keys [eval ratio object]}] chebyshev-forms]
      (t/is (thrown? IllegalArgumentException (eval n 0.3)) (str kind " eval " n))
      (t/is (thrown? IllegalArgumentException (ratio n)) (str kind " ratio " n))
      (t/is (thrown? IllegalArgumentException (object n)) (str kind " object " n)))
    (doseq [[label f] {"eval-legendre-P" #(sut/eval-legendre-P % 0.3)
                       "legendre-P-ratio" sut/legendre-P-ratio
                       "legendre-P" sut/legendre-P
                       "eval-gegenbauer-C, default parameter" #(sut/eval-gegenbauer-C % 0.3)
                       "eval-gegenbauer-C, parameter 2.5" #(sut/eval-gegenbauer-C % 2.5 0.3)
                       "gegenbauer-C-ratio" #(sut/gegenbauer-C-ratio % 2.5)
                       "gegenbauer-C" #(sut/gegenbauer-C % 2.5)
                       "eval-jacobi-P" #(sut/eval-jacobi-P % 0.5 1.5 0.3)
                       "jacobi-P-ratio" #(sut/jacobi-P-ratio % 0.5 1.5)
                       "jacobi-P" #(sut/jacobi-P % 0.5 1.5)}]
      (t/is (thrown? IllegalArgumentException (f n)) (str label " " n))))
  ;; the message names the limit
  (t/is (thrown-with-msg? IllegalArgumentException #"below 2147483647" (sut/eval-legendre-P 4294967296 0.3)))
  ;; small degrees are unaffected (a degree just below the limit would need billions of steps)
  (t/is (== 1.0 (sut/eval-legendre-P 0 0.3)))
  (t/is (== 26.0 (sut/eval-chebyshev-T 3 2.0))))

;; Generalized Laguerre and Hermite (physicists' H and probabilists' He) polynomials.
;;
;; Reference values: `test/resources/polynomials/laguerre_hermite_reference.edn`, exact rational
;; arithmetic from the explicit sums (`utils/fastmath/dev/generate_laguerre_hermite_reference.py`) at the
;; exact binary values of the order and of the arguments: 13 orders (including 0, negative integers and
;; fractions, 0.3 and 1e-9 which are not binary exact, and 1 + 2^-30), 26 or 24 arguments from 1e-9 to 700 of
;; both signs, degrees up to 100 (Laguerre) and 200 (Hermite); and the exact coefficients.
;; `scipy.special.eval_genlaguerre/eval_hermite/eval_hermitenorm` against the reference: largest difference
;; 8.9e-14 (Laguerre, order 7, degree 100), 1.1e-13 (H, degree 200) and 1.4e-14 (He), relative to
;; max(1, |value|); scipy returns NaN for the 1300 Laguerre rows of orders -1 and below.

(def ^:private lh-reference
  (delay (edn/read-string (slurp (io/resource "polynomials/laguerre_hermite_reference.edn")))))

;; Hermite: largest observed ratio to `eval-error-bound` 0.35 (H) and 0.3 (He); limit 4.
(def ^:private lh-eval-units {:hermite-h 4.0 :hermite-he 4.0})

(defn- laguerre-grid-failures
  "Rows `[n x value]` for which `eval-laguerre-L` is outside `units * eps * (n+1)^2 * max(1, max_i |L_i(x)|)`.
  The forward recurrence loses accuracy like the square of the degree next to x = 0 (about 1e-13 relative at
  degree 100) and its error is relative to the largest of the intermediate values `L_0 ... L_n`, not to `L_n`."
  [order rows units]
  (for [[n x v] rows
        :let [got (attempt sut/eval-laguerre-L n order x)
              envelope (reduce max 1.0 (map #(m/abs (sut/eval-laguerre-L % order x)) (range (inc n))))]
        :when (not (and (not (failed? got))
                        (<= (m/abs (- got v)) (* units EPS (m/sq (inc n)) envelope))))]
    [n x v got]))

;; Largest observed ratio to this bound: 0.5 units (a rounding of the degree 1 value); limit 2.
(t/deftest laguerre-eval-reference
  (doseq [{:keys [order grid]} (:laguerre @lh-reference)
          :let [failures (laguerre-grid-failures order grid 2.0)]]
    (t/is (empty? failures) (failures-message (str "Laguerre order " order) failures))))

(t/deftest hermite-eval-reference
  (doseq [[label eval-fn key] [["H" sut/eval-hermite-H :hermite-h] ["He" sut/eval-hermite-He :hermite-he]]
          :let [failures (grid-failures eval-fn (get-in @lh-reference [key :grid]) (get lh-eval-units key))]]
    (t/is (empty? failures) (failures-message (str "Hermite " label) failures))))

(t/deftest laguerre-exact-coefficients
  (doseq [{:keys [order decimal-exact? coefficients]} (:laguerre @lh-reference)
          :when decimal-exact?
          [n exact] (map-indexed vector coefficients)
          :let [r (attempt sut/laguerre-L-ratio n order)
                o (attempt sut/laguerre-L n order)
                expected (exact-coefficients exact)]]
    (t/is (and (not (failed? r)) (= expected (vec (sut/coeffs r))) (= n (sut/degree r)))
          (str "ratio, order " order " degree " n))
    (t/is (and (not (failed? o)) (= (mapv exact->double expected) (vec (sut/coeffs o))) (= n (sut/degree o)))
          (str "object, order " order " degree " n))))

(t/deftest hermite-exact-coefficients
  (doseq [[label ratio-fn object-fn key] [["H" sut/hermite-H-ratio sut/hermite-H :hermite-h]
                                          ["He" sut/hermite-He-ratio sut/hermite-He :hermite-he]]
          [n exact] (map-indexed vector (get-in @lh-reference [key :coefficients]))
          :let [r (attempt ratio-fn n)
                o (attempt object-fn n)
                expected (exact-coefficients exact)]]
    (t/is (and (not (failed? r)) (= expected (vec (sut/coeffs r))) (= n (sut/degree r))) (str label " ratio, degree " n))
    (t/is (and (not (failed? o)) (= (mapv exact->double expected) (vec (sut/coeffs o))) (= n (sut/degree o)))
          (str label " object, degree " n))))

(defn- generalized-binomial [top k]
  (/ (reduce *' 1 (map #(- top %) (range k))) (factorial-exact k)))

(defn- laguerre-explicit-coefficients
  "Coefficients, ascending, of the generalized Laguerre polynomial of a rational order from the explicit sum
  `sum_k (-1)^k C(n+a, n-k) x^k / k!`."
  [n alpha]
  (mapv (fn [k] (/ (* (if (even? k) 1 -1) (generalized-binomial (+ n alpha) (- n k))) (factorial-exact k)))
        (range (inc n))))

(t/deftest laguerre-ratio-form-is-exact-for-decimal-orders
  ;; the order is converted with rationalize (0.3 is 3/10): the coefficients are then exact
  (doseq [alpha [0.1 0.3 1.7 0.123 -0.6 2.5 -3.0 1e-9] n (range 0 13)]
    (t/is (= (laguerre-explicit-coefficients n (rationalize alpha)) (vec (sut/coeffs (sut/laguerre-L-ratio n alpha))))
          (str "order " alpha " degree " n))))

(t/deftest laguerre-hermite-three-forms-agree
  (let [xs [-3.0 -1.5 -0.5 0.0 0.5 1.5 3.0]
        run (fn [eval-fn ratio-fn object-fn cases]
              (three-forms-disagreements {:eval-fn eval-fn :ratio-fn ratio-fn :object-fn object-fn
                                          :cases cases :xs xs}))
        ns [0 1 2 3 4 5 6 8 10 15]]
    (let [d (run (fn [[n a] x] (sut/eval-laguerre-L n a x)) (fn [[n a]] (sut/laguerre-L-ratio n a))
                 (fn [[n a]] (sut/laguerre-L n a))
                 (for [n ns a [0.0 0.5 1.0 2.5 7.0 0.3 -0.5 -1.0 -2.5]] [n a]))]
      (t/is (empty? d) (str "Laguerre: " (count d) " disagreements; first: " (pr-str (take 2 d)))))
    (let [d (run (fn [[n] x] (sut/eval-hermite-H n x)) (fn [[n]] (sut/hermite-H-ratio n)) (fn [[n]] (sut/hermite-H n))
                 (map vector ns))]
      (t/is (empty? d) (str "H: " (pr-str (take 2 d)))))
    (let [d (run (fn [[n] x] (sut/eval-hermite-He n x)) (fn [[n]] (sut/hermite-He-ratio n)) (fn [[n]] (sut/hermite-He n))
                 (map vector ns))]
      (t/is (empty? d) (str "He: " (pr-str (take 2 d)))))))

(t/deftest laguerre-order-conventions
  ;; the order defaults to 0
  (doseq [n [0 1 2 5 9] x [-2.0 0.0 0.5 3.0]]
    (t/is (= (bits (sut/eval-laguerre-L n x)) (bits (sut/eval-laguerre-L n 0.0 x))) (str "degree " n " at " x)))
  (t/is (= (vec (sut/coeffs (sut/laguerre-L 5))) (vec (sut/coeffs (sut/laguerre-L 5 0.0)))))
  ;; L_1 = 1 - x, L_2 = (x^2 - 4x + 2) / 2
  (t/is (= [1 -1] (vec (sut/coeffs (sut/laguerre-L-ratio 1 0.0)))))
  (t/is (= [1 -2 1/2] (vec (sut/coeffs (sut/laguerre-L-ratio 2 0.0)))))
  ;; a negative integer order -k with n >= k: L_n^(-k)(x) = (n-k)! / n! * (-x)^k * L_(n-k)^(k)(x)
  (doseq [k [1 2 3] n [3 4 6 9]
          :let [shifted (sut/laguerre-L-ratio (- n k) (double k))
                factor (/ (factorial-exact (- n k)) (factorial-exact n))
                monomial (sut/ratio-polynomial (concat (repeat k 0) [(if (even? k) 1 -1)]))]]
    (t/is (= (vec (sut/coeffs (sut/scale (sut/mult monomial shifted) factor)))
             (vec (sut/coeffs (sut/laguerre-L-ratio n (double (- k))))))
          (str "order " (- k) " degree " n)))
  ;; a NaN order: degree 0 is 1, the others NaN; no exact form for NaN or an infinite order
  (t/is (== 1.0 (sut/eval-laguerre-L 0 ##NaN 0.3)))
  (t/is (m/nan? (sut/eval-laguerre-L 3 ##NaN 0.3)))
  (t/is (= [1] (vec (sut/coeffs (sut/laguerre-L-ratio 0 ##NaN)))))
  (doseq [bad [##NaN ##Inf ##-Inf] n [1 3]]
    (t/is (thrown? IllegalArgumentException (sut/laguerre-L-ratio n bad)) (str "order " bad " degree " n))
    (t/is (thrown? IllegalArgumentException (sut/laguerre-L n bad)) (str "order " bad " degree " n))))

(t/deftest hermite-special-values
  ;; H_n(0) = (-1)^(n/2) n! / (n/2)!, He_n(0) = (-1)^(n/2) (n-1)!!, both 0 for an odd degree
  (doseq [[n h he] [[0 1.0 1.0] [2 -2.0 -1.0] [4 12.0 3.0] [6 -120.0 -15.0] [8 1680.0 105.0] [1 0.0 0.0] [3 0.0 0.0] [7 0.0 0.0]]]
    (t/is (== h (sut/eval-hermite-H n 0.0)) (str "H " n))
    (t/is (== he (sut/eval-hermite-He n 0.0)) (str "He " n)))
  ;; H_n(x) = 2^(n/2) He_n(sqrt(2) x)
  (doseq [n [2 3 6 9] x [-1.5 0.25 0.9 2.0]]
    (t/is (m/delta-eq (sut/eval-hermite-H n x) (* (m/pow 2.0 (/ n 2.0)) (sut/eval-hermite-He n (* m/SQRT2 x))) 1.0e-9)
          (str "degree " n " at " x)))
  ;; parity
  (doseq [n [1 2 5 8 13] x [0.3 1.7 4.0]]
    (let [sign (if (even? n) 1.0 -1.0)]
      (t/is (== (sut/eval-hermite-H n x) (* sign (sut/eval-hermite-H n (- x)))) (str "H " n))
      (t/is (== (sut/eval-hermite-He n x) (* sign (sut/eval-hermite-He n (- x)))) (str "He " n)))))

(t/deftest laguerre-hermite-non-finite-and-huge-arguments
  ;; an infinite argument gives the infinity of the leading term (sign of the coefficient and the parity);
  ;; the leading coefficients are (-1)^n / n!, 2^n and 1, never zero
  (doseq [n [0 1 2 3 4 5 6 9 10 15 40]
          :let [even (even? n)]]
    (doseq [order [0.0 0.5 2.5 -0.5 -1.0 -2.5 7.0 0.3 1.0e-9]]
      (t/is (== (cond (zero? n) 1.0 even ##Inf :else ##-Inf) (sut/eval-laguerre-L n order ##Inf)) (str "L " n " order " order " at +Inf"))
      (t/is (== (if (zero? n) 1.0 ##Inf) (sut/eval-laguerre-L n order ##-Inf)) (str "L " n " order " order " at -Inf")))
    (t/is (== (if (zero? n) 1.0 (if even ##Inf ##-Inf)) (sut/eval-laguerre-L n ##Inf)) (str "L " n " default order at +Inf"))
    (doseq [[label f] [["H" sut/eval-hermite-H] ["He" sut/eval-hermite-He]]]
      (t/is (== (if (zero? n) 1.0 ##Inf) (f n ##Inf)) (str label " " n " at +Inf"))
      (t/is (== (cond (zero? n) 1.0 even ##Inf :else ##-Inf) (f n ##-Inf)) (str label " " n " at -Inf"))))
  ;; an infinite argument with a NaN or infinite order is NaN
  (doseq [order [##NaN ##Inf ##-Inf] x [##Inf ##-Inf]]
    (t/is (m/nan? (sut/eval-laguerre-L 2 order x)) (str "order " order " at " x)))
  ;; NaN argument: degree 0 is 1, the others NaN
  (t/is (== 1.0 (sut/eval-laguerre-L 0 0.5 ##NaN) (sut/eval-hermite-H 0 ##NaN) (sut/eval-hermite-He 0 ##NaN)))
  (doseq [n [1 2 3 7]]
    (t/is (m/nan? (sut/eval-laguerre-L n 0.5 ##NaN)))
    (t/is (m/nan? (sut/eval-hermite-H n ##NaN)))
    (t/is (m/nan? (sut/eval-hermite-He n ##NaN))))
  ;; degree 0 and 1
  (t/is (== 1.0 (sut/eval-laguerre-L 0 2.5 ##Inf)))
  (t/is (== 1.5 (sut/eval-laguerre-L 1 2.5 2.0)) "1 + a - x")
  (t/is (== 0.6 (sut/eval-hermite-H 1 0.3)))
  (t/is (== 0.3 (sut/eval-hermite-He 1 0.3)))
  ;; huge arguments with a representable value
  (t/is (m/delta-eq 1.0 (/ (sut/eval-laguerre-L 2 0.0 1.0e100) 5.0e199) 1.0e-14) "x^2 / 2")
  (t/is (m/delta-eq 1.0 (/ (sut/eval-hermite-H 2 1.0e150) 4.0e300) 1.0e-14) "4 x^2 - 2")
  (t/is (m/delta-eq 1.0 (/ (sut/eval-hermite-He 2 1.0e150) 1.0e300) 1.0e-14) "x^2 - 1"))

(t/deftest laguerre-hermite-degree-limits
  (doseq [n [-1 -2 -100 Integer/MAX_VALUE 4294967296 Long/MAX_VALUE]]
    (doseq [[label f] {"eval-laguerre-L, default order" #(sut/eval-laguerre-L % 0.3)
                       "eval-laguerre-L" #(sut/eval-laguerre-L % 2.5 0.3)
                       "laguerre-L-ratio" #(sut/laguerre-L-ratio % 2.5)
                       "laguerre-L, default order" sut/laguerre-L
                       "laguerre-L" #(sut/laguerre-L % 2.5)
                       "eval-hermite-H" #(sut/eval-hermite-H % 0.3)
                       "hermite-H-ratio" sut/hermite-H-ratio
                       "hermite-H" sut/hermite-H
                       "eval-hermite-He" #(sut/eval-hermite-He % 0.3)
                       "hermite-He-ratio" sut/hermite-He-ratio
                       "hermite-He" sut/hermite-He}]
      (t/is (thrown? IllegalArgumentException (f n)) (str label " " n)))))

;; Bernstein basis polynomials and Bessel polynomials.
;;
;; Reference values: `test/resources/polynomials/bernstein_bessel_reference.edn`, exact rational arithmetic from
;; the explicit sums (`utils/fastmath/dev/generate_bernstein_bessel_reference.py`) at the exact binary values of
;; the arguments: Bernstein degrees 0 to 5000, orders from -1 to degree + 3, 19 arguments (including 0, 1 and
;; values outside [0, 1]) and the exact coefficients; the Bessel polynomials y and theta (degrees up to 100,
;; arguments of both signs) and their exact coefficients.
;; One-off checks of the reference: `scipy.stats.binom.pmf` for the 740 Bernstein rows in [0, 1]: largest
;; relative difference 3.0e-13 (degree 5000, scipy's own error); `mpmath` through the modified Bessel function
;; `y_n(x) = sqrt(2/(pi x)) exp(1/x) K_(n+1/2)(1/x)` and `theta_n(x) = x^n y_n(1/x)` for the 258 rows with
;; x > 0: 1.1e-16.

(def ^:private bb-reference
  (delay (edn/read-string (slurp (io/resource "polynomials/bernstein_bessel_reference.edn")))))

(defn- bernstein-failures
  "Rows `[degree order x value]` for which `eval-bernstein` is outside `units * eps * (n+1) * |value|`
  (a reference value of 0 must be returned exactly)."
  [rows units]
  (for [[n k x v] rows
        :let [got (attempt sut/eval-bernstein n k x)]
        :when (not (and (not (failed? got))
                        (if (zero? v) (zero? got) (<= (m/abs (- got v)) (* units EPS (inc n) (m/abs v))))))]
    [n k x v got]))

;; The product `C(n, k) x^k (1-x)^(n-k)` loses about (n+1) units of roundoff relative to the value (`m/combinations`
;; uses a log-beta formula from k = 30); largest observed ratio to this bound: 1.6 (degree 500); limit 4.
(t/deftest bernstein-eval-reference
  (let [failures (bernstein-failures (get-in @bb-reference [:bernstein :grid]) 4.0)]
    (t/is (empty? failures) (failures-message "Bernstein" failures))))

(t/deftest bernstein-exact-coefficients
  (doseq [[n k exact] (get-in @bb-reference [:bernstein :coefficients])
          :let [o (attempt sut/bernstein n k)
                expected (mapv double (exact-coefficients exact))]]
    (t/is (and (not (failed? o)) (= expected (vec (sut/coeffs o))) (= n (sut/degree o)))
          (str "degree " n " order " k))))

(t/deftest bernstein-basis-properties
  ;; the basis of one degree sums to 1, also for degrees where the binomial coefficient overflows
  (doseq [n [0 1 2 3 7 20 100 1000 2000 5000] x [0.0 0.1 0.37 0.5 0.9 1.0]]
    (t/is (m/delta-eq 1.0 (reduce + (map #(sut/eval-bernstein n % x) (range (inc n)))) 1.0e-9)
          (str "degree " n " at " x)))
  (doseq [n [1 2 3] x [-0.5 1.7]]
    (t/is (m/delta-eq 1.0 (reduce + (map #(sut/eval-bernstein n % x) (range (inc n)))) 1.0e-12) (str "degree " n " at " x)))
  ;; b(k, n)(x) = b(n-k, n)(1-x)
  (doseq [n [2 5 13 30] k (range (inc n)) x [0.05 0.3 0.5 0.8]
          :let [a (sut/eval-bernstein n k x) b (sut/eval-bernstein n (- n k) (- 1.0 x))]]
    (t/is (<= (m/abs (- a b)) (* 1.0e-12 (max a b))) (str n " " k " " x)))
  ;; the values at the ends of the interval are exact
  (doseq [n [1 2 5 40 3000]]
    (t/is (== 1.0 (sut/eval-bernstein n 0 0.0) (sut/eval-bernstein n n 1.0)))
    (doseq [k (range 1 (inc n)) :when (< k 4)]
      (t/is (== 0.0 (sut/eval-bernstein n k 0.0)))
      (t/is (== 0.0 (sut/eval-bernstein n (- n k) 1.0)))))
  ;; degrees 0 and 1
  (t/is (== 1.0 (sut/eval-bernstein 0 0 0.3) (sut/eval-bernstein 0 0 ##NaN) (sut/eval-bernstein 0 0 ##Inf)))
  (t/is (== 0.7 (sut/eval-bernstein 1 0 0.3)))
  (t/is (== 0.3 (sut/eval-bernstein 1 1 0.3)))
  (t/is (== 0.5 (sut/eval-bernstein 2 1 0.5)))
  ;; infinite and NaN arguments: the leading term; NaN for a positive degree
  (doseq [n [1 2 3 6] k (range (inc n))]
    (t/is (== (if (even? (- n k)) ##Inf ##-Inf) (sut/eval-bernstein n k ##Inf)) (str n " " k " at +Inf"))
    (t/is (== (if (even? k) ##Inf ##-Inf) (sut/eval-bernstein n k ##-Inf)) (str n " " k " at -Inf"))
    (t/is (m/nan? (sut/eval-bernstein n k ##NaN)) (str n " " k " at NaN"))))

(t/deftest bernstein-order-outside-the-degree-is-the-zero-function
  ;; b(k, n) = 0 for k < 0 and for k > n: for every argument, also the ends, NaN and the infinities
  (doseq [[n k] [[0 1] [0 -1] [1 5] [1 -1] [2 3] [5 6] [5 -2] [5 7] [10 -1] [3000 3001]]
          x [0.0 1.0 0.3 -2.0 1.5 ##NaN ##Inf ##-Inf]]
    (t/is (== 0.0 (sut/eval-bernstein n k x)) (str "degree " n " order " k " at " x)))
  (doseq [[n k] [[0 1] [0 -1] [1 5] [1 -1] [3 4] [5 -2] [40 41]]
          :let [o (sut/bernstein n k)]]
    (t/is (and (every? #(== 0.0 %) (sut/coeffs o)) (= n (sut/degree o))) (str "object, degree " n " order " k))
    (t/is (not-any? #(neg? (bits %)) (sut/coeffs o)) "positive zeros")))

;; Bessel y and theta, standard model of `eval-error-bound`: largest observed ratio 0.40 (y) and 0.44 (theta),
;; for x >= 0 from the recurrences and for x < 0 from the decimal explicit sum; limit 2. (The recurrences at
;; x < 0 gave relative errors from 1e-8 to total loss, theta_50(-30) with the wrong sign.)
(def ^:private bb-eval-units {:bessel-y 2.0 :bessel-t 2.0})

(t/deftest bessel-eval-reference
  (doseq [[label eval-fn key] [["y" sut/eval-bessel-y :bessel-y] ["theta" sut/eval-bessel-t :bessel-t]]
          :let [failures (grid-failures eval-fn (get-in @bb-reference [key :grid]) (get bb-eval-units key))]]
    (t/is (empty? failures) (failures-message (str "Bessel " label) failures))))

(t/deftest bessel-exact-coefficients
  (doseq [[label ratio-fn object-fn key] [["y" sut/bessel-y-ratio sut/bessel-y :bessel-y]
                                          ["theta" sut/bessel-t-ratio sut/bessel-t :bessel-t]]
          [n exact] (map-indexed vector (get-in @bb-reference [key :coefficients]))
          :let [r (attempt ratio-fn n)
                o (attempt object-fn n)
                expected (exact-coefficients exact)]]
    (t/is (and (not (failed? r)) (= expected (vec (sut/coeffs r))) (= n (sut/degree r))) (str label " ratio, degree " n))
    (t/is (and (not (failed? o)) (= (mapv exact->double expected) (vec (sut/coeffs o))) (= n (sut/degree o)))
          (str label " object, degree " n))))

(t/deftest bessel-three-forms-agree
  (let [xs [-3.0 -1.5 -0.5 0.0 0.5 1.5 3.0 10.0]
        ns [0 1 2 3 4 5 6 8 10 15]]
    (let [d (three-forms-disagreements {:eval-fn (fn [[n] x] (sut/eval-bessel-y n x)) :ratio-fn (fn [[n]] (sut/bessel-y-ratio n))
                                        :object-fn (fn [[n]] (sut/bessel-y n)) :cases (map vector ns) :xs xs})]
      (t/is (empty? d) (str "y: " (pr-str (take 2 d)))))
    (let [d (three-forms-disagreements {:eval-fn (fn [[n] x] (sut/eval-bessel-t n x)) :ratio-fn (fn [[n]] (sut/bessel-t-ratio n))
                                        :object-fn (fn [[n]] (sut/bessel-t n)) :cases (map vector ns) :xs xs})]
      (t/is (empty? d) (str "theta: " (pr-str (take 2 d)))))))

(t/deftest bessel-polynomials-satisfy-their-equations
  ;; independent of the recurrences: the exact residual of the differential equations is the zero polynomial
  (doseq [n (range 0 21)
          :let [y (sut/bessel-y-ratio n)
                theta (sut/bessel-t-ratio n)
                residual-y (sut/add (sut/add (sut/mult (sut/ratio-polynomial [0 0 1]) (sut/derivative y 2))
                                              (sut/mult (sut/ratio-polynomial [2 2]) (sut/derivative y 1)))
                                    (sut/scale y (- (* n (inc n)))))
                residual-theta (sut/add (sut/add (sut/mult (sut/ratio-polynomial [0 1]) (sut/derivative theta 2))
                                                 (sut/mult (sut/ratio-polynomial [(* -2 n) -2]) (sut/derivative theta 1)))
                                        (sut/scale theta (* 2 n)))]]
    ;; x^2 y'' + (2x + 2) y' - n (n+1) y = 0 and x theta'' - 2 (x + n) theta' + 2 n theta = 0
    (t/is (every? zero? (sut/coeffs residual-y)) (str "y, degree " n))
    (t/is (every? zero? (sut/coeffs residual-theta)) (str "theta, degree " n))
    ;; theta_n(x) = x^n y_n(1/x): the coefficients in the reverse order; y_n(0) = 1, theta_n is monic
    (t/is (= (vec (reverse (sut/coeffs y))) (vec (sut/coeffs theta))) (str "degree " n))
    (t/is (= 1 (first (sut/coeffs y)) (last (sut/coeffs theta))) (str "degree " n))
    ;; theta_n(0) = (2n-1)!!
    (t/is (= (reduce *' 1 (range 1 (inc (* 2 n)) 2)) (first (sut/coeffs theta))) (str "degree " n)))
  (t/is (= [1 3 3] (vec (sut/coeffs (sut/bessel-y-ratio 2)))))
  (t/is (= [1 6 15 15] (vec (sut/coeffs (sut/bessel-y-ratio 3)))))
  (t/is (= [15 15 6 1] (vec (sut/coeffs (sut/bessel-t-ratio 3)))))
  ;; the value that the recurrence with the wrong previous term gave: y_3(0.5) = 9.625, not 9.125
  (t/is (== 9.625 (sut/eval-bessel-y 3 0.5)))
  (t/is (== 36.9375 (sut/eval-bessel-y 4 0.5))))

(t/deftest bessel-non-finite-and-huge-arguments
  ;; the leading term is x^n in both: the sign follows the parity of the degree at -Inf
  (doseq [n [1 2 3 4 5 6 9 10 15]
          [label f] [["y" sut/eval-bessel-y] ["theta" sut/eval-bessel-t]]]
    (t/is (== ##Inf (f n ##Inf)) (str label " " n " at +Inf"))
    (t/is (== (if (even? n) ##Inf ##-Inf) (f n ##-Inf)) (str label " " n " at -Inf"))
    (t/is (m/nan? (f n ##NaN)) (str label " " n " at NaN")))
  (doseq [f [sut/eval-bessel-y sut/eval-bessel-t] x [##NaN ##Inf ##-Inf 0.3]]
    (t/is (== 1.0 (f 0 x))))
  (t/is (== 1.5 (sut/eval-bessel-y 1 0.5)))
  (t/is (== 1.5 (sut/eval-bessel-t 1 0.5)))
  ;; huge arguments with a representable value: y_2 = 3 x^2 + 3 x + 1, theta_2 = x^2 + 3 x + 3
  (t/is (m/delta-eq 1.0 (/ (sut/eval-bessel-y 2 1.0e150) 3.0e300) 1.0e-14))
  (t/is (m/delta-eq 1.0 (/ (sut/eval-bessel-t 2 1.0e150) 1.0e300) 1.0e-14))
  (t/is (m/delta-eq 1.0 (/ (sut/eval-bessel-t 2 -1.0e150) 1.0e300) 1.0e-14)))

(t/deftest bernstein-and-bessel-degree-limits
  (doseq [n [-1 -2 -100 Integer/MAX_VALUE 4294967296 Long/MAX_VALUE]]
    (doseq [[label f] {"eval-bernstein" #(sut/eval-bernstein % 0 0.3)
                       "eval-bernstein, order outside the degree" #(sut/eval-bernstein % 5 0.3)
                       "bernstein" #(sut/bernstein % 0)
                       "eval-bessel-y" #(sut/eval-bessel-y % 0.3)
                       "bessel-y-ratio" sut/bessel-y-ratio
                       "bessel-y" sut/bessel-y
                       "eval-bessel-t" #(sut/eval-bessel-t % 0.3)
                       "bessel-t-ratio" sut/bessel-t-ratio
                       "bessel-t" sut/bessel-t}]
      (t/is (thrown? IllegalArgumentException (f n)) (str label " " n)))))

;; Meixner-Pollaczek polynomials P_n^(lambda)(x; phi).
;;
;; Reference values: `test/resources/polynomials/meixner_pollaczek_reference.edn`, `mpmath` at 60 digits at the exact
;; binary values of lambda, phi and x (`utils/fastmath/dev/generate_meixner_pollaczek_reference.py`): 8 values of
;; lambda (0.25 to 7, and 0, -0.5, -1.5), 10 values of phi (including pi/2, pi, 0, negative and above pi), 21
;; arguments of both signs, degrees up to 50. For lambda > 0 the values come from the hypergeometric definition
;; `(2 lambda)_n / n! exp(i n phi) 2F1(-n, lambda + i x; 2 lambda; 1 - exp(-2 i phi))`; the generator asserts that
;; it equals the three term recurrence (largest relative difference 3e-60), which is used for the other values of
;; lambda.

(def ^:private mp-reference
  (delay (edn/read-string (slurp (io/resource "polynomials/meixner_pollaczek_reference.edn")))))

(defn- meixner-pollaczek-failures
  "Rows `[n x value derivative]` for which `eval-meixner-pollaczek-P` is outside
  `units * eps * (n+1)^2 * max(1, max_i |P_i(x)|)`: the forward recurrence loses accuracy like the square of the
  degree (as for the Laguerre polynomials), relative to the largest value of the degrees 0 to n."
  [lambda phi rows units]
  (for [[n x v] rows
        :let [got (attempt sut/eval-meixner-pollaczek-P n lambda phi x)
              envelope (reduce max 1.0 (map #(m/abs (sut/eval-meixner-pollaczek-P % lambda phi x)) (range (inc n))))]
        :when (not (and (not (failed? got)) (<= (m/abs (- got v)) (* units EPS (m/sq (inc n)) envelope))))]
    [n x v got]))

;; Largest observed ratio to this bound: 2.6 (a value that is a difference of two terms of size 7); limit 6.
(t/deftest meixner-pollaczek-eval-reference
  (doseq [{:keys [lambda phi grid]} @mp-reference
          :let [failures (meixner-pollaczek-failures lambda phi grid 6.0)]]
    (t/is (empty? failures) (failures-message (str "Meixner-Pollaczek lambda " lambda " phi " phi) failures))))

(defn- meixner-pollaczek-exact-value
  "P_n at the rational `x` by the three term recurrence in exact arithmetic, for the rational `lam`, `c` = cos(phi)
  and `s` = sin(phi)."
  [n lam c s x]
  (if (zero? n)
    1
    (loop [i 2 pprev 1 prev (* 2 (+ (* lam c) (* x s)))]
      (if (> i n)
        prev
        (recur (inc i) prev (/ (- (* 2 (+ (* x s) (* (+ lam (dec i)) c)) prev) (* (+ (- i 2) (* 2 lam)) pprev)) i))))))

(t/deftest meixner-pollaczek-ratio-form-is-exact
  ;; lambda, cos(phi) and sin(phi) are converted with rationalize (the decimal numbers the doubles print as)
  (doseq [lambda [0.5 2.5 0.3 0.0 -0.5 -1.5] phi [1.0 0.1 (/ Math/PI 2) 3.0 0.0 -1.0 7.0] n (range 0 11)
          :let [lam (rationalize lambda) c (rationalize (Math/cos phi)) s (rationalize (Math/sin phi))
                polynomial (sut/meixner-pollaczek-P-ratio n lambda phi)]]
    (t/is (= n (sut/degree polynomial)))
    (doseq [x [1/3 -2 7/5 0]]
      (t/is (= (meixner-pollaczek-exact-value n lam c s x) (polynomial x))
            (str "lambda " lambda " phi " phi " degree " n " at " x)))))

(t/deftest meixner-pollaczek-three-forms-agree
  (let [d (three-forms-disagreements
           {:eval-fn (fn [[n l p] x] (sut/eval-meixner-pollaczek-P n l p x))
            :ratio-fn (fn [[n l p]] (sut/meixner-pollaczek-P-ratio n l p))
            :object-fn (fn [[n l p]] (sut/meixner-pollaczek-P n l p))
            :cases (for [n [0 1 2 3 4 5 6 8 10] l [0.5 1.0 2.5 0.3 -0.5] p [0.5 (/ Math/PI 2) 2.0 3.0 -1.0]] [n l p])
            :xs [-3.0 -1.5 -0.5 0.0 0.5 1.5 3.0]
            ;; the error of the recurrence grows like (n+1)^2 units, not n: the default 8 units are for the linear family
            :eval-units 40.0})]
    (t/is (empty? d) (str (count d) " disagreements; first: " (pr-str (take 2 d))))))

(t/deftest meixner-pollaczek-closed-forms-and-symmetry
  ;; P_1 = 2 (lambda cos(phi) + x sin(phi)), and P_2 from the recurrence
  (t/is (m/delta-eq (* 2.0 (+ (* 0.5 (Math/cos 1.0)) (* 0.3 (Math/sin 1.0)))) (sut/eval-meixner-pollaczek-P 1 0.5 1.0 0.3) 1.0e-15))
  (let [lambda 2.5 phi 0.7 x -1.2
        p1 (* 2.0 (+ (* lambda (Math/cos phi)) (* x (Math/sin phi))))
        p2 (/ (- (* 2.0 (+ (* x (Math/sin phi)) (* (+ lambda 1.0) (Math/cos phi))) p1) (* 2.0 lambda)) 2.0)]
    (t/is (m/delta-eq p2 (sut/eval-meixner-pollaczek-P 2 lambda phi x) 1.0e-13)))
  ;; phi = 0: the polynomial does not depend on x, P_n = (2 lambda)_n / n!, for lambda = 1 it is n + 1, also at infinity
  (doseq [n [0 1 2 5 12] x [-3.0 0.0 2.5 ##Inf ##-Inf]]
    (t/is (== (inc n) (sut/eval-meixner-pollaczek-P n 1.0 0.0 x)) (str "degree " n " at " x)))
  (t/is (== -4.0 (sut/eval-meixner-pollaczek-P 3 1.0 Math/PI 0.3)) "phi = pi: (-1)^n (n + 1)")
  ;; P_n(x; pi - phi) = (-1)^n P_n(-x; phi)
  (doseq [n [1 2 3 6 9] lambda [0.5 2.5] phi [0.4 1.0 2.0] x [-2.0 0.3 1.7]
          :let [a (sut/eval-meixner-pollaczek-P n lambda (- Math/PI phi) x)
                b (* (if (even? n) 1.0 -1.0) (sut/eval-meixner-pollaczek-P n lambda phi (- x)))]]
    (t/is (<= (m/abs (- a b)) (* 1.0e-12 (max 1.0 (m/abs a)))) (str n " " lambda " " phi " " x)))
  ;; phi = pi/2: cos(phi) is 6e-17, not 0 (the double next to pi/2): P_1(0) is lambda times 1.2e-16, with its full
  ;; relative accuracy (`m/cos` had a relative error of 6e-11 there)
  (t/is (== (* 2.0 0.5 (Math/cos (/ Math/PI 2))) (sut/eval-meixner-pollaczek-P 1 0.5 (/ Math/PI 2) 0.0))))

(t/deftest meixner-pollaczek-non-finite-arguments
  ;; infinite x: the leading term (2 sin(phi) x)^n / n!; the sign follows sin(phi) x and the parity of the degree
  (doseq [n [1 2 3 4 7] lambda [0.5 2.5 -0.5 0.0] phi [0.5 2.0 3.0 -1.0 7.0]
          :let [sp (Math/sin phi)]
          x [##Inf ##-Inf]]
    (t/is (== (if (or (pos? (* sp x)) (even? n)) ##Inf ##-Inf) (sut/eval-meixner-pollaczek-P n lambda phi x))
          (str "degree " n " lambda " lambda " phi " phi " at " x)))
  ;; NaN: degree 0 is 1 for any arguments
  (t/is (== 1.0 (sut/eval-meixner-pollaczek-P 0 ##NaN ##NaN ##NaN) (sut/eval-meixner-pollaczek-P 0 0.5 1.0 ##Inf)))
  (doseq [n [1 2 3 6]]
    (t/is (m/nan? (sut/eval-meixner-pollaczek-P n 0.5 1.0 ##NaN)))
    (t/is (m/nan? (sut/eval-meixner-pollaczek-P n ##NaN 1.0 0.3)))
    (t/is (m/nan? (sut/eval-meixner-pollaczek-P n 0.5 ##NaN 0.3)))
    (t/is (m/nan? (sut/eval-meixner-pollaczek-P n 0.5 ##Inf 0.3)))
    (t/is (m/nan? (sut/eval-meixner-pollaczek-P n ##NaN 1.0 ##Inf)))
    (t/is (m/nan? (sut/eval-meixner-pollaczek-P n 0.5 ##NaN ##Inf)))
    (t/is (thrown? IllegalArgumentException (sut/meixner-pollaczek-P-ratio n ##NaN 1.0)))
    (t/is (thrown? IllegalArgumentException (sut/meixner-pollaczek-P-ratio n 0.5 ##Inf)))
    (t/is (thrown? IllegalArgumentException (sut/meixner-pollaczek-P n 0.5 ##NaN))))
  (t/is (= [1] (vec (sut/coeffs (sut/meixner-pollaczek-P-ratio 0 ##NaN ##NaN)))))
  ;; huge finite x with a representable value: P_2 ~ (2 sin(phi) x)^2 / 2
  (t/is (m/delta-eq 1.0 (/ (sut/eval-meixner-pollaczek-P 2 1.0 (/ Math/PI 2) 1.0e150) 2.0e300) 1.0e-14)))

(t/deftest meixner-pollaczek-degree-limits
  (doseq [n [-1 -2 -100 Integer/MAX_VALUE 4294967296 Long/MAX_VALUE]]
    (doseq [[label f] {"eval-meixner-pollaczek-P" #(sut/eval-meixner-pollaczek-P % 1.0 1.0 0.3)
                       "meixner-pollaczek-P-ratio" #(sut/meixner-pollaczek-P-ratio % 1.0 1.0)
                       "meixner-pollaczek-P" #(sut/meixner-pollaczek-P % 1.0 1.0)}]
      (t/is (thrown? IllegalArgumentException (f n)) (str label " " n)))))

;; Ince polynomials C_p^m and S_p^m, angular and radial.
;;
;; Reference: `test/resources/polynomials/ince_reference.edn`, `mpmath` at 60 digits
;; (`utils/fastmath/dev/generate_ince_reference.py`). The eigenproblem of Ince's equation
;; `w'' + e sin(2x) w' + (a - p e cos(2x)) w = 0` is built from the equation itself (not from the DLMF 28.31
;; coefficient recurrences the library uses) by projecting `L[b_r]` onto the Fourier basis; the eigenvalue of
;; the polynomial of degree `m` is the `m/2`-th, `(m-1)/2`-th, ... in increasing order. 406 entries: p from 0 to
;; 8 and 12, every valid m, e = 0, 0.3, 0.6, -0.6, 2, 5, 10. The generator asserts the equation at random points
;; to 1e-45 and `a = m^2` for `e = 0`. Normalization of the reference: Euclidean norm 1 of the coefficients,
;; C(0) > 0 and S'(0) > 0.

(def ^:private ince-reference
  (delay (edn/read-string (slurp (io/resource "polynomials/ince_reference.edn")))))

(defn- ince-frequencies
  "The frequencies k of the basis cos(k x) (kind :C) or sin(k x) (:S) of the Ince polynomial of order p."
  [kind p]
  (cond
    (and (= kind :C) (even? p)) (map #(* 2 %) (range (inc (quot p 2))))
    (= kind :C) (map #(inc (* 2 %)) (range (inc (quot (dec p) 2))))
    (even? p) (map #(* 2 (inc %)) (range (quot p 2)))
    :else (map #(inc (* 2 %)) (range (inc (quot (dec p) 2))))))

(defn- ince-series
  "[w w' w''] at `x` of the series with the `coefficients`; the trigonometric basis, or the hyperbolic one when
  `radial?`."
  [kind p coefficients x radial?]
  (reduce (fn [[w w1 w2] [c k]]
            (let [kx (* k x)
                  [f df] (cond (and (= kind :C) (not radial?)) [(Math/cos kx) (- (Math/sin kx))]
                               (= kind :C) [(Math/cosh kx) (Math/sinh kx)]
                               (not radial?) [(Math/sin kx) (Math/cos kx)]
                               :else [(Math/sinh kx) (Math/cosh kx)])
                  second-sign (if radial? 1.0 -1.0)]
              [(+ w (* c f)) (+ w1 (* c k df)) (+ w2 (* c k k second-sign f))]))
          [0.0 0.0 0.0]
          (map vector coefficients (ince-frequencies kind p))))

(defn- ince-residual
  "Left-hand side of Ince's equation (angular) or of the radial equation `R'' - e sinh(2x) R' - (a - p e cosh(2x)) R`
  for the series at `x`."
  [kind p e a coefficients x radial?]
  (let [[w w1 w2] (ince-series kind p coefficients x radial?)]
    (if radial?
      (- w2 (* e (Math/sinh (* 2 x)) w1) (* (- a (* p e (Math/cosh (* 2 x)))) w))
      (+ w2 (* e (Math/sin (* 2 x)) w1) (* (- a (* p e (Math/cos (* 2 x)))) w)))))

(def ^:private ince-xs [0.0 0.13 0.31 0.62 0.83 1.05 1.5 2.4 3.0])

;; Largest observed difference to the reference: 3.2e-15 (C 8 6 at e = 10); limit 1e-12.
(t/deftest ince-coefficients-match-the-reference
  (doseq [{:keys [kind p m e coefficients]} @ince-reference
          :let [got (attempt (if (= kind :C) sut/ince-C-coeffs sut/ince-S-coeffs) p m e :none)]]
    (t/is (and (not (failed? got))
               (= (count coefficients) (count got))
               (every? #(<= (m/abs %) 1.0e-12) (map - got coefficients)))
          (str (name kind) " p " p " m " m " e " e ": " (if (failed? got) got (vec got))))))

(t/deftest ince-coefficients-satisfy-the-equation-with-the-reference-eigenvalue
  ;; independent of the reference coefficients: the library's own coefficients and the reference `a`
  (doseq [{:keys [kind p m e a]} @ince-reference
          :when (pos? p)
          :let [coefficients (vec ((if (= kind :C) sut/ince-C-coeffs sut/ince-S-coeffs) p m e :none))]
          x ince-xs]
    ;; residual scaled by (1 + |a|)(p + 1): largest observed 3.5e-15; limit 1e-13
    (t/is (<= (m/abs (ince-residual kind p e a coefficients x false)) (* 1.0e-13 (+ 1.0 (m/abs a)) (inc p)))
          (str (name kind) " p " p " m " m " e " e " at " x))))

(t/deftest ince-functions-evaluate-the-series
  ;; the angular and the radial functions against the series of the reference coefficients
  (doseq [{:keys [kind p m e coefficients]} @ince-reference
          :let [angular ((if (= kind :C) sut/ince-C sut/ince-S) p m e)
                radial ((if (= kind :C) sut/ince-C-radial sut/ince-S-radial) p m e)]
          x [0.0 0.13 0.62 1.05]]
    (let [[w] (ince-series kind p coefficients x false)
          [r] (ince-series kind p coefficients x true)]
      (t/is (<= (m/abs (- (angular x) w)) 1.0e-12) (str "angular " (name kind) " p " p " m " m " e " e " at " x))
      (t/is (<= (m/abs (- (radial x) r)) (* 1.0e-12 (+ 1.0 (m/abs r))))
            (str "radial " (name kind) " p " p " m " m " e " e " at " x)))))

(t/deftest ince-radial-equation-with-the-same-eigenvalue
  ;; R'' - e sinh(2x) R' - (a - p e cosh(2x)) R = 0 with the same `a` as the angular equation
  (doseq [{:keys [kind p m e a coefficients]} @ince-reference
          :when (and (pos? p) (<= p 8))
          x [0.1 0.3 0.7 1.0]]
    ;; scaled likewise and by 1 + |R|: largest observed 1.0e-14; limit 1e-12
    (t/is (<= (m/abs (ince-residual kind p e a coefficients x true)) (* 1.0e-12 (+ 1.0 (m/abs a)) (inc p)
                                                                       (+ 1.0 (m/abs (first (ince-series kind p coefficients x true))))))
          (str "radial " (name kind) " p " p " m " m " e " e " at " x))))

(t/deftest ince-reduces-to-trigonometric-functions-at-zero-e
  (doseq [p (range 0 9) m (range 0 (inc p)) :when (even? (- p m))]
    (let [x 0.7]
      (t/is (m/delta-eq (Math/cos (* m x)) ((sut/ince-C p m 0.0) x) 1.0e-14) (str "C " p " " m))
      (t/is (m/delta-eq (Math/cosh (* m x)) ((sut/ince-C-radial p m 0.0) x) 1.0e-13) (str "C radial " p " " m))
      (when (pos? m)
        (t/is (m/delta-eq (Math/sin (* m x)) ((sut/ince-S p m 0.0) x) 1.0e-14) (str "S " p " " m))
        (t/is (m/delta-eq (Math/sinh (* m x)) ((sut/ince-S-radial p m 0.0) x) 1.0e-13) (str "S radial " p " " m))))))

(t/deftest ince-normalizations
  (doseq [{:keys [kind p m e]} @ince-reference
          :when (and (pos? p) (<= p 8) (<= (m/abs e) 2.0))
          :let [coeffs (if (= kind :C) sut/ince-C-coeffs sut/ince-S-coeffs)
                none (vec (coeffs p m e :none))
                trigonometric (vec (coeffs p m e :trigonometric))
                millers (vec (coeffs p m e :millers))
                function ((if (= kind :C) sut/ince-C sut/ince-S) p m e :trigonometric)
                ;; (1/pi) integral of w^2 over a period, exact for equal steps of a trigonometric polynomial
                steps 512
                mean-square (/ (reduce + (map #(let [w (function (* % (/ (* 2.0 Math/PI) steps)))] (* w w)) (range steps))) steps)]]
    ;; :none is the unit vector of the coefficients
    (t/is (m/delta-eq 1.0 (reduce + (map * none none)) 1.0e-13) (str "none " (name kind) " " p " " m " " e))
    ;; :trigonometric: (1/pi) integral over [0, 2 pi] of w^2 = 1
    (t/is (m/delta-eq 1.0 (* 2.0 mean-square) 1.0e-12) (str "trigonometric " (name kind) " " p " " m " " e))
    ;; the same as :none for S and for odd p
    (when (or (= kind :S) (odd? p))
      (t/is (every? #(<= (m/abs %) 1.0e-14) (map - none trigonometric)) (str "none = trigonometric " (name kind) " " p " " m)))
    ;; :millers is a positive multiple of :none (same direction and sign)
    (let [ratios (keep (fn [[a b]] (when (> (m/abs b) 1.0e-8) (/ a b))) (map vector millers none))]
      (t/is (and (every? m/valid-double? ratios) (pos? (first ratios))
                 (every? #(<= (m/abs (- % (first ratios))) (* 1.0e-10 (m/abs (first ratios)))) ratios))
            (str "millers " (name kind) " " p " " m " " e))))
  ;; the sign: C(0) > 0, S'(0) > 0
  (doseq [{:keys [kind p m e]} @ince-reference :when (pos? p)]
    (if (= kind :C)
      (t/is (pos? ((sut/ince-C p m e) 0.0)) (str "C(0) > 0 " p " " m " " e))
      (t/is (pos? ((sut/ince-S p m e) 1.0e-6)) (str "S'(0) > 0 " p " " m " " e)))))

(t/deftest ince-order-zero
  ;; C_0^0 is the constant: the single coefficient 1 (:none), 1/sqrt(2) (:trigonometric)
  (doseq [e [0.0 0.6 -2.0 10.0]]
    (t/is (= [1.0] (vec (sut/ince-C-coeffs 0 0 e :none))))
    (t/is (m/delta-eq (/ 1.0 (Math/sqrt 2.0)) (first (sut/ince-C-coeffs 0 0 e :trigonometric)) 1.0e-15))
    (t/is (== 1.0 ((sut/ince-C 0 0 e) 0.7) ((sut/ince-C 0 0 e) -3.0) ((sut/ince-C-radial 0 0 e) 1.5)))
    (t/is (pos? (first (sut/ince-C-coeffs 0 0 e :millers))))))

(t/deftest ince-invalid-arguments
  (let [entry-points {"ince-C-coeffs" #(sut/ince-C-coeffs %1 %2 %3 %4)
                      "ince-S-coeffs" #(sut/ince-S-coeffs %1 %2 %3 %4)
                      "ince-C" #(sut/ince-C %1 %2 %3 %4)
                      "ince-S" #(sut/ince-S %1 %2 %3 %4)
                      "ince-C-radial" #(sut/ince-C-radial %1 %2 %3 %4)
                      "ince-S-radial" #(sut/ince-S-radial %1 %2 %3 %4)}]
    (doseq [[label f] entry-points
            [p m e normalization why] [[4 6 0.6 :none "m above p"]
                                       [4 -2 0.6 :none "negative m"]
                                       [4 3 0.6 :none "parity"]
                                       [3 2 0.6 :none "parity"]
                                       [-2 0 0.6 :none "negative p"]
                                       [-1 -1 0.6 :none "negative p"]
                                       [4 2 ##NaN :none "NaN e"]
                                       [4 2 ##Inf :none "infinite e"]
                                       [4 2 ##-Inf :none "infinite e"]
                                       [4 2 0.6 :bogus "unknown normalization"]
                                       [4 2 0.6 "millers" "string instead of a keyword"]
                                       [4 2 0.6 nil "no normalization"]]]
      (t/is (thrown? IllegalArgumentException (f p m e normalization)) (str label ": " why)))
    ;; S starts at m = 1
    (doseq [f [(get entry-points "ince-S-coeffs") (get entry-points "ince-S") (get entry-points "ince-S-radial")]
            [p m] [[0 0] [2 0] [4 0] [3 0]]]
      (t/is (thrown? IllegalArgumentException (f p m 0.6 :none)) (str "S with m = " m " and p = " p)))
    ;; the edges of the valid range work
    (doseq [[label f] entry-points
            [p m] [[1 1] [2 2] [4 4] [5 5] [6 0] [7 1]]
            :when (or (pos? m) (#{"ince-C-coeffs" "ince-C" "ince-C-radial"} label))]
      (t/is (not (failed? (attempt f p m 0.6 :none))) (str label " " p " " m)))
    (t/is (not (failed? (attempt (get entry-points "ince-C") 0 0 0.6 :none))) "p = 0")))

;; ---------------------------------------------------------------------------------------------
;; Fixes after the #Test attack on groups 1-2 and 5-9 (vault: Test Review - Groups 1-2 and 5-9)

(t/deftest ratio-to-double-is-the-nearest-double                           ; T-10
  ;; `(double 1/6)` is 0.1666666666666667: Ratio.doubleValue keeps 16 decimal digits
  (t/is (== 0.16666666666666666 (sut/evaluate (sut/ratio-polynomial [1/6]) 0.0)))
  (t/is (== 0.16666666666666666 (first (sut/coeffs (sut/polynomial [1/6])))))
  (t/is (== 0.16666666666666666 (first (sut/coeffs (sut/coeffs->polynomial 1/6)))))
  (let [rng (java.util.Random. 3)]
    (dotimes [_ 600]
      (let [cs (repeatedly 4 #(- (.nextDouble rng) 0.5))
            x (- (* 4.0 (.nextDouble rng)) 2.0)
            ;; `evaluate` works at the decimal value of x
            exact (reduce (fn [acc c] (+ (* acc (rationalize x)) (rationalize c))) 0 (reverse cs))]
        (t/is (== (exact->double exact) (sut/evaluate (sut/ratio-polynomial cs) x)) (str cs " " x)))))
  ;; the object forms hold the nearest doubles of their exact coefficients
  (doseq [[label ratio-form object-form]
          [["laguerre-L 20" (sut/laguerre-L-ratio 20 0.0) (sut/laguerre-L 20 0.0)]
           ["laguerre-L 30 0.5" (sut/laguerre-L-ratio 30 0.5) (sut/laguerre-L 30 0.5)]
           ["legendre-P 30" (sut/legendre-P-ratio 30) (sut/legendre-P 30)]
           ["gegenbauer-C 25 0.75" (sut/gegenbauer-C-ratio 25 0.75) (sut/gegenbauer-C 25 0.75)]
           ["jacobi-P 25 0.5 1.5" (sut/jacobi-P-ratio 25 0.5 1.5) (sut/jacobi-P 25 0.5 1.5)]
           ["meixner-pollaczek-P 20" (sut/meixner-pollaczek-P-ratio 20 0.7 1.1) (sut/meixner-pollaczek-P 20 0.7 1.1)]
           ["bessel-t 30" (sut/bessel-t-ratio 30) (sut/bessel-t 30)]]]
    (t/is (= (map exact->double (sut/coeffs ratio-form)) (sut/coeffs object-form)) label)))

(t/deftest polynomial-printing-keeps-every-term-up-to-degree-10            ; T-19
  (let [s10 "#polynomial{10}(x) = 1+2x+3x^2+4x^3+5x^4+6x^5+7x^6+8x^7+9x^8+10x^9+11x^10"]
    (t/is (= s10 (str (sut/polynomial (range 1 12)))))
    (t/is (= s10 (str (sut/ratio-polynomial (range 1 12)))))
    (t/is (= s10 (pr-str (sut/polynomial (range 1 12)))))
    (t/is (= "#polynomial{11}(x) = 1+2x+3x^2+4x^3+5x^4+6x^5+7x^6+8x^7+9x^8+10x^9+11x^10+..."
             (str (sut/polynomial (range 1 13)))))))

(t/deftest derivative-of-a-high-order-keeps-the-factor-finite              ; T-20
  (let [cs (concat (repeat 171 0.0) [1e-30 0.0 0.0 0.0 0.0])
        got (vec (sut/coeffs (sut/derivative (sut/polynomial cs) 171)))
        exact (mapv exact->double (sut/coeffs (sut/derivative (sut/ratio-polynomial cs) 171)))]
    (t/is (== 5 (count got)))
    (t/is (<= (m/abs (- (first got) (first exact))) (* 1e-12 (first exact))))
    (t/is (every? #(== 0.0 %) (rest got)) "zero coefficients stay zero, not NaN"))
  ;; the true result overflows: an infinity, and zeros elsewhere
  (let [got (vec (sut/coeffs (sut/derivative (sut/polynomial (concat (repeat 200 0.0) [1.0 0.0 0.0])) 200)))]
    (t/is (= [##Inf 0.0 0.0] got)))
  ;; order 1000 of a degree 1200 polynomial: finite and equal to the exact value
  (let [cs (concat (repeat 1000 0.0) [1e-300] (repeat 199 0.0))
        got (vec (sut/coeffs (sut/derivative (sut/polynomial cs) 1000)))]
    (t/is (= 200 (count got)))
    (t/is (not-any? #(Double/isNaN %) got))))

(t/deftest inline-evalpoly-evaluates-the-coefficients-in-order             ; T-21
  (let [log (atom [])
        note (fn [k v] (swap! log conj k) v)]
    (t/is (== 6.0 (sut/evalpoly 1.0 (note :c0 1) (note :c1 2) (note :c2 3))))
    (t/is (= [:c0 :c1 :c2] @log))
    (reset! log [])
    (t/is (== 7.0 (sut/evalpoly (note :x 2.0) (note :c0 1) (note :c1 3))))
    (t/is (= [:x :c0 :c1] @log))
    (reset! log [])
    (t/is (== 5.0 (sut/evalpoly 1.0 (note :c0 5))))
    (t/is (= [:c0] @log))))

(t/deftest zero-sign-of-evalpoly                                           ; T-22
  (let [bits #(Double/doubleToRawLongBits (double %))]
    (t/is (= (bits -0.0) (bits (sut/evalpoly 1.0 -0.0 -0.0))))
    (t/is (= (bits -0.0) (bits (sut/evalpoly 1.0 -0.0))))
    (t/is (= (bits 0.0) (bits (sut/evalpoly 1.0 1.0 -1.0))) "an exact cancellation is +0.0")
    (t/is (= (bits 0.0) (bits (sut/mevalpoly 1.0 1.0 -1.0))))))

(t/deftest bernstein-with-a-subnormal-power                                ; T-14
  ;; 0.3^616 and 0.5935642792213479^1402 are subnormal while the products are about 1e-168
  (doseq [[n k x exact] [[794 616 0.3 2.321385782417405E-168]
                         [1546 1402 0.5935642792213479 4.095317600262778E-168]]]
    (let [v (sut/eval-bernstein n k x)]
      (t/is (<= (m/abs (- 1.0 (/ v exact))) (* 8.0 (inc n) EPS)) (str n " " k " " x " " v))))
  ;; exact values by BigDecimal for a range of cases next to the change of regime
  (let [mc (java.math.MathContext. 120)
        exact-value (fn [n k x]
                      (let [binomial (java.math.BigDecimal. (str (reduce (fn [^BigInteger a ^long j] (.divide (.multiply a (BigInteger/valueOf (- (inc n) j))) (BigInteger/valueOf j)))
                                                                         BigInteger/ONE (range 1 (inc k)))))
                            bx (java.math.BigDecimal. (double x))]
                        (.doubleValue (.multiply (.multiply binomial (.pow bx (int k) mc) mc)
                                                 (.pow (.subtract java.math.BigDecimal/ONE bx) (int (- n k)) mc) mc))))]
    (doseq [[n k x] [[700 560 0.3] [900 700 0.35] [1200 1000 0.4] [2000 1700 0.45] [1000 800 0.2] [800 780 0.7]]
            :let [exact (exact-value n k x)]
            :when (> (m/abs exact) 1e-290)]
      (t/is (<= (m/abs (- 1.0 (/ (sut/eval-bernstein n k x) exact))) (* 8.0 (inc n) EPS)) (str n " " k " " x)))))

;; T-15, T-16: a recurrence that leaves the double range gives the signed infinity, not NaN

(def ^:private max-exact-recurrence-degree
  "The exact ratio recurrences below are for small degrees only: the denominators (2^53 for a double) grow
  with the degree, and a degree in the hundreds takes minutes. The decimal references serve those."
  60)

(defn- check-exact-degree!
  [^long n]
  (when (> n max-exact-recurrence-degree)
    (throw (IllegalArgumentException.
            (str "The exact ratio reference is too slow for degree " n " (limit " max-exact-recurrence-degree
                 "): use the decimal reference")))))

(defn- exact-hermite-double
  "The nearest double of `H_n` (kind :H) or `He_n` (kind :He) at the double `x`, from the recurrence in ratios.
  Degree up to `max-exact-recurrence-degree`."
  [kind ^long n x]
  (check-exact-degree! n)
  (let [rx (exact-ratio x)
        h1 (if (= kind :H) (* 2 rx) rx)
        step (fn [i prev pprev] (let [t (- (* rx prev) (* (dec i) pprev))] (if (= kind :H) (* 2 t) t)))]
    (exact->double
     (cond (zero? n) 1
           (== n 1) h1
           :else (loop [i 2 pprev 1 prev h1]
                   (if (> i n) prev (recur (inc i) prev (step i prev pprev))))))))

(defn- exact-laguerre-double
  "The nearest double of the generalized Laguerre polynomial at the doubles `a` and `x`, from the recurrence in
  ratios. Degree up to `max-exact-recurrence-degree`."
  [^long n a x]
  (check-exact-degree! n)
  (let [ra (exact-ratio a)
        rx (exact-ratio x)
        l1 (- (+ 1 ra) rx)]
    (exact->double
     (cond (zero? n) 1
           (== n 1) l1
           :else (loop [i 2 pprev 1 prev l1]
                   (if (> i n)
                     prev
                     (recur (inc i) prev (/ (- (* (- (+ (dec (* 2 i)) ra) rx) prev) (* (+ (dec i) ra) pprev)) i))))))))

(defn- agrees-with-expected?
  "True when `got` is the infinity `expected`, or is finite and within 1e-6 relative (plus 1e-300) of the finite
  `expected`."
  [^double expected ^double got]
  (if (Double/isInfinite expected)
    (== expected got)
    (and (Double/isFinite got)
         (<= (m/abs (- got expected)) (+ (* 1e-6 (m/abs expected)) 1e-300)))))

(t/deftest overflow-gives-the-signed-infinity                              ; T-15
  (doseq [[n x e] [[6 1e100 ##Inf] [6 -1e100 ##Inf] [5 1e154 ##Inf] [5 -1e154 ##-Inf] [7 1e308 ##Inf] [7 -1e308 ##-Inf] [300 0.3 ##Inf]]]
    (t/is (= e (sut/eval-hermite-H n x)) (str "H " n " " x)))
  (doseq [[n x e] [[6 1e100 ##Inf] [5 -1e154 ##-Inf] [400 0.3 ##Inf]]]
    (t/is (= e (sut/eval-hermite-He n x)) (str "He " n " " x)))
  (doseq [[n a x e] [[6 0.0 1e100 ##Inf] [5 0.0 1e154 ##-Inf] [100 0.0 1e5 ##Inf] [4 2.0 -1e300 ##Inf] [3 2.0 -1e300 ##Inf]]]
    (t/is (= e (sut/eval-laguerre-L n a x)) (str "L " n " " a " " x)))
  (doseq [[n l phi x e] [[3 1e154 1.0 0.5 ##Inf] [5 1.0 1.0 1e154 ##Inf] [5 1.0 -1.0 1e154 ##-Inf] [4 1.0 1.0 -1e200 ##Inf]]]
    (t/is (= e (sut/eval-meixner-pollaczek-P n l phi x)) (str "MP " n " " l " " phi " " x)))
  ;; a NaN input still gives NaN
  (t/is (Double/isNaN (sut/eval-hermite-H 6 ##NaN)))
  (t/is (Double/isNaN (sut/eval-laguerre-L 6 ##NaN 1e100)))
  (t/is (Double/isNaN (sut/eval-laguerre-L 6 0.0 ##NaN)))
  (t/is (Double/isNaN (sut/eval-meixner-pollaczek-P 6 ##NaN 1.0 1e100))))

;; The references for degrees in the hundreds use 200 digit decimals (see `max-exact-recurrence-degree`).

(def ^:private reference-digits (java.math.MathContext. 200))

(defn- decimal-recurrence
  "Value at `n` of `P_0 = 1`, `P_1 = first-term`, `P_i = (step i prev pprev)`, in 200 digit decimals."
  [^long n first-term step]
  (cond (zero? n) java.math.BigDecimal/ONE
        (== n 1) first-term
        :else (loop [i 2 pprev java.math.BigDecimal/ONE prev first-term]
                (if (> i n) prev (recur (inc i) prev (step i prev pprev))))))

(defn- decimal-hermite-double
  "The nearest double of `H_n` (kind :H) or `He_n` (kind :He) at the double `x`, from 200 digit decimals; any
  degree."
  [kind ^long n x]
  (let [bx (java.math.BigDecimal. (double x))
        mc reference-digits
        two (java.math.BigDecimal. 2)
        t (fn [i ^java.math.BigDecimal prev ^java.math.BigDecimal pprev]
            (.subtract (.multiply bx prev mc) (.multiply (java.math.BigDecimal. (long (dec i))) pprev mc) mc))]
    (.doubleValue ^java.math.BigDecimal
     (if (= kind :H)
       (decimal-recurrence n (.multiply two bx mc) (fn [i prev pprev] (.multiply two ^java.math.BigDecimal (t i prev pprev) mc)))
       (decimal-recurrence n bx t)))))

(defn- decimal-laguerre-double
  "The nearest double of the generalized Laguerre polynomial at the doubles `a` and `x`, from 200 digit
  decimals; any degree."
  [^long n a x]
  (let [ba (java.math.BigDecimal. (double a))
        bx (java.math.BigDecimal. (double x))
        mc reference-digits]
    (.doubleValue ^java.math.BigDecimal
     (decimal-recurrence n (.subtract (.add java.math.BigDecimal/ONE ba) bx)
                         (fn [i ^java.math.BigDecimal prev ^java.math.BigDecimal pprev]
                           (let [factor (.subtract (.add (java.math.BigDecimal. (long (dec (* 2 i)))) ba) bx)
                                 weight (.add (java.math.BigDecimal. (long (dec i))) ba)]
                             (.divide (.subtract (.multiply factor prev mc) (.multiply weight pprev mc) mc)
                                      (java.math.BigDecimal. (long i)) mc)))))))

(t/deftest overflow-in-the-oscillatory-region-has-the-sign-of-the-exact-value  ; T-15
  (let [rng (java.util.Random. 11)]
    (dotimes [_ 40]
      (let [n (+ 150 (.nextInt rng 600))
            x (- (* 6.0 (.nextDouble rng)) 3.0)
            a (- (* 8.0 (.nextDouble rng)) 2.0)
            xl (* 40.0 (.nextDouble rng))]
        (t/is (agrees-with-expected? (decimal-hermite-double :H n x) (sut/eval-hermite-H n x)) (str "H " n " " x))
        (t/is (agrees-with-expected? (decimal-hermite-double :He n x) (sut/eval-hermite-He n x)) (str "He " n " " x))
        (t/is (agrees-with-expected? (decimal-laguerre-double n a xl) (sut/eval-laguerre-L n a xl)) (str "L " n " " a " " xl))))))

(t/deftest decimal-references-agree-with-the-exact-ones
  (doseq [n [5 17 40 60] x [0.3 -1.7 2.5]]
    (t/is (agrees-with-expected? (exact-hermite-double :H n x) (decimal-hermite-double :H n x)) (str "H " n " " x))
    (t/is (agrees-with-expected? (exact-hermite-double :He n x) (decimal-hermite-double :He n x)) (str "He " n " " x))
    (t/is (agrees-with-expected? (exact-laguerre-double n 1.5 x) (decimal-laguerre-double n 1.5 x)) (str "L " n " " x))))

(t/deftest exact-references-refuse-a-high-degree
  ;; the exact ratio recurrence takes minutes there: a clear error instead
  (t/is (thrown? IllegalArgumentException (exact-hermite-double :H 61 0.3)))
  (t/is (thrown? IllegalArgumentException (exact-hermite-double :He 750 0.3)))
  (t/is (thrown? IllegalArgumentException (exact-laguerre-double 600 0.5 2.5)))
  (t/is (number? (exact-laguerre-double 60 0.5 2.5))))

(t/deftest recurrences-near-the-end-of-the-double-range                    ; T-16
  ;; (a+1)(a+2)/2 = 1.125e308 is representable; the product in the recurrence is not
  (t/is (agrees-with-expected? (exact-laguerre-double 2 1.5e154 0.0) (sut/eval-laguerre-L 2 1.5e154 0.0)))
  (t/is (< 1.1e308 (sut/eval-laguerre-L 2 1.5e154 0.0) 1.2e308))
  (doseq [x [2.7e102 2.8e102 2.9e102 6e153 6.7e153 6.8e153]
          n [2 3]]
    (t/is (agrees-with-expected? (exact-hermite-double :H n x) (sut/eval-hermite-H n x)) (str "H " n " " x))
    (t/is (agrees-with-expected? (exact-hermite-double :He n x) (sut/eval-hermite-He n x)) (str "He " n " " x)))
  (doseq [a [1.0e100 1.0e150 1.4e154 1.5e154 1.9e154]]
    (t/is (agrees-with-expected? (exact-laguerre-double 2 a 0.0) (sut/eval-laguerre-L 2 a 0.0)) (str "L 2 " a))))

;; T-11, T-13: Ince polynomials for large and extreme e. The eigenproblem is solved in a symmetric form.
;; Reference: `ince_large_e_reference.edn` (`mpmath`, built from the differential equation; |e| up to 1e8,
;; p up to 30).

(def ^:private ince-large-e-reference
  (delay (edn/read-string (slurp (io/resource "polynomials/ince_large_e_reference.edn")))))

(defn- max-difference-up-to-sign
  "Smallest over the two signs of the largest coefficient difference."
  [as bs]
  [(apply max (map #(m/abs (- (double %1) (double %2))) as bs))
   (apply max (map #(m/abs (+ (double %1) (double %2))) as bs))])

(t/deftest ince-coefficients-for-large-e
  (let [entries @ince-large-e-reference
        worst (atom 0.0)]
    (t/is (== 560 (count entries)))
    (doseq [{:keys [kind p m e zero-value coefficients]} entries
            :let [got (vec (if (= kind :C) (sut/ince-C-coeffs p m e :none) (sut/ince-S-coeffs p m e :none)))
                  [same flipped] (max-difference-up-to-sign got coefficients)
                  ;; the sign follows C(0) > 0 or S'(0) > 0, which the sum of the coefficients decides: reliable
                  ;; only above rounding
                  difference (if (< (m/abs zero-value) 1e-10) (min same flipped) same)]]
      (swap! worst max difference)
      (t/is (< difference 1e-12) (str kind " p=" p " m=" m " e=" e " difference " difference)))
    ;; observed at most 3.5e-15
    (t/is (< @worst 1e-12))))

(t/deftest ince-for-extreme-finite-e                                       ; T-13
  (doseq [e [1e-100 1e-150 1e-200 1e-300 2.3e-308 1e-310 5e-324 -5e-324 -1e-300]
          [kind p m] [[:C 6 2] [:C 4 4] [:C 5 3] [:S 5 3] [:S 6 4] [:C 2 2]]
          :let [got (vec (if (= kind :C) (sut/ince-C-coeffs p m e :none) (sut/ince-S-coeffs p m e :none)))
                ;; the limit e -> 0: the unit vector of cos(m x) / sin(m x)
                index (case kind :C (if (even? p) (quot m 2) (quot (dec m) 2)) :S (if (even? p) (dec (quot m 2)) (quot (dec m) 2)))]]
    (t/is (every? #(Double/isFinite %) got) (str kind " " p " " m " " e))
    (t/is (< (m/abs (- 1.0 (m/abs (got index)))) 1e-12) (str kind " " p " " m " " e)))
  ;; huge finite e: the vectors stay unit vectors with finite entries and no exception
  (doseq [e [1e150 1e200 1e300 1.7e308 -1e300]
          [kind p m] [[:C 4 2] [:S 5 3] [:C 8 8] [:S 6 2]]
          :let [got (vec (if (= kind :C) (sut/ince-C-coeffs p m e :none) (sut/ince-S-coeffs p m e :none)))]]
    (t/is (every? #(Double/isFinite %) got) (str kind " " p " " m " " e))
    (t/is (< (m/abs (- 1.0 (Math/sqrt (reduce + (map #(* % %) got))))) 1e-12)))
  ;; huge e is continuous: the limit vector is reached
  (t/is (< (first (max-difference-up-to-sign (sut/ince-C-coeffs 4 2 1e50 :none) (sut/ince-C-coeffs 4 2 1e300 :none))) 1e-9)))

;; T-12: Miller's normalization for large p

(t/deftest ince-millers-normalization-for-large-p
  ;; the multiple of :none is the inverse of sqrt(sum (w_r a_r)^2): checked against the weights in exact form
  ;; for a p where the plain weights are still finite, and by continuity beyond
  (doseq [kind [:C :S]
          p [100 170 172 196 198 200 250 300 330]
          :let [m (if (= kind :C) 0 (if (even? p) 2 1))
                none (vec (if (= kind :C) (sut/ince-C-coeffs p m 0.5 :none) (sut/ince-S-coeffs p m 0.5 :none)))
                millers (vec (if (= kind :C) (sut/ince-C-coeffs p m 0.5 :millers) (sut/ince-S-coeffs p m 0.5 :millers)))
                ratios (keep (fn [[a b]] (when (and (not (zero? a)) (not (zero? b))) (/ b a))) (map vector none millers))]]
    (t/is (every? #(Double/isFinite %) millers) (str kind " " p))
    (t/is (pos? (first ratios)) (str kind " " p))
    (t/is (some pos? millers) (str kind " " p " not all zero"))
    ;; a positive multiple of :none (where the entries do not underflow)
    (t/is (every? #(<= (m/abs (- 1.0 (/ % (first ratios)))) 1e-9) (filter #(not (Double/isNaN %)) (take 5 ratios))) (str kind " " p)))
  ;; the multiple falls like 1/sqrt(p!): known values
  (t/is (m/delta-eq 2.9278788877268646E-65 (/ (first (sut/ince-C-coeffs 100 0 0.5 :millers)) (first (sut/ince-C-coeffs 100 0 0.5 :none))) 1e-77))
  (t/is (< 0.0 (first (sut/ince-C-coeffs 300 0 0.5 :millers)) 1e-260))
  ;; below the range of a double: a clear error, not zeros or NaN
  (doseq [p [400 500 1000]]
    (t/is (thrown? IllegalArgumentException (sut/ince-C-coeffs p 0 0.5 :millers)) (str "C " p))
    (t/is (thrown? IllegalArgumentException (sut/ince-S-coeffs (inc p) 1 0.5 :millers)) (str "S " (inc p))))
  ;; the other normalizations have no such limit
  (t/is (every? #(Double/isFinite %) (sut/ince-C-coeffs 1000 0 0.5 :trigonometric))))

;; T-23: the radial functions return the signed infinity where the series is out of the double range

(t/deftest ince-radial-overflow-is-an-infinity
  (let [fns {"C 4 2" (sut/ince-C-radial 4 2 0.5) "S 4 2" (sut/ince-S-radial 4 2 0.5) "C 4 0" (sut/ince-C-radial 4 0 0.5)
             "C 3 1" (sut/ince-C-radial 3 1 0.5) "S 3 3" (sut/ince-S-radial 3 3 0.5) "C 2 2" (sut/ince-C-radial 2 2 -0.5)}]
    (doseq [[label f] fns
            xi [400.0 710.0 1e5 1e300 ##Inf]]
      (let [v (f xi)]
        (t/is (and (not (Double/isNaN v)) (Double/isInfinite v)) (str label " at " xi " gave " v))))
    ;; the sign: the leading cosh / sinh term decides; for sinh the sign of xi too
    (let [coeffs (vec (sut/ince-S-coeffs 4 2 0.5 :none))
          last-coefficient (peek coeffs)]
      (t/is (= (Math/signum last-coefficient) (Math/signum ((fns "S 4 2") 400.0))))
      (t/is (= (- (Math/signum last-coefficient)) (Math/signum ((fns "S 4 2") -400.0))))
      (t/is (= (- (Math/signum last-coefficient)) (Math/signum ((fns "S 4 2") ##-Inf)))))
    (let [coeffs (vec (sut/ince-C-coeffs 4 2 0.5 :none))]
      (t/is (= (Math/signum (peek coeffs)) (Math/signum ((fns "C 4 2") 400.0))))
      (t/is (= (Math/signum (peek coeffs)) (Math/signum ((fns "C 4 2") -400.0)))))
    ;; before the range the value is finite and the series is unchanged
    (t/is (Double/isFinite ((fns "C 4 2") 100.0)))
    (t/is (Double/isNaN ((fns "C 4 2") ##NaN)))
    (t/is (m/delta-eq (/ 1.0 (Math/sqrt 2.0)) ((sut/ince-C-radial 0 0 0.5 :trigonometric) 1e300) 1e-15) "p = 0: no growth")))
