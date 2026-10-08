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
      (t/is (<= (m/abs (- (sut/evaluate p xx) (double exact))) bound) (str "Polynomial evaluate " cs " at " xx))
      (t/is (<= (m/abs (- (double (p xx)) (double exact))) bound) (str "Polynomial call " cs " at " xx))
      (t/is (= exact (r xx)) (str "PolynomialR call " cs " at " xx))
      (t/is (== (double exact) (sut/evaluate r xx)) (str "PolynomialR evaluate " cs " at " xx)))))

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
    (t/is (and (not (failed? o)) (= (mapv double exact) (vec (sut/coeffs o))) (= n (sut/degree o)))
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

;; Jacobi (-3, 0.5): the polynomials of degree 3 and more have a triple root at x = 1, so the values next to
;; it (down to 1e-25 at 1 - 1e-9) are far below the rounding error of the recurrence relative to the size
;; of the polynomial on the interval (1.5e-11 at degree 30; 4482 units of the bound).
(def ^:private jacobi-units-by-parameters {[-3.0 0.5] 5000.0})

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
          :let [failures (grid-failures (fn [n x] (sut/eval-jacobi-P n alpha beta x)) grid
                                        (get jacobi-units-by-parameters [alpha beta] (:jacobi lgj-eval-units)))]]
    (t/is (empty? failures) (failures-message (str "Jacobi alpha " alpha " beta " beta) failures))))

(defn- exact-coefficients [pairs] (mapv ratio pairs))

(t/deftest legendre-exact-coefficients
  (doseq [[n exact] (map-indexed vector (get-in @lgj-reference [:legendre :coefficients]))
          :let [r (attempt sut/legendre-P-ratio n)
                o (attempt sut/legendre-P n)
                expected (exact-coefficients exact)]]
    (t/is (and (not (failed? r)) (= expected (vec (sut/coeffs r))) (= n (sut/degree r))) (str "ratio, degree " n))
    (t/is (and (not (failed? o)) (= (mapv double expected) (vec (sut/coeffs o))) (= n (sut/degree o))) (str "object, degree " n))))

(t/deftest gegenbauer-exact-coefficients
  (doseq [{:keys [alpha decimal-exact? coefficients]} (:gegenbauer @lgj-reference)
          :when decimal-exact?
          [n exact] (map-indexed vector coefficients)
          :let [r (attempt sut/gegenbauer-C-ratio n alpha)
                o (attempt sut/gegenbauer-C n alpha)
                expected (exact-coefficients exact)]]
    (t/is (and (not (failed? r)) (= expected (vec (sut/coeffs r))) (= n (sut/degree r)))
          (str "ratio, alpha " alpha " degree " n))
    (t/is (and (not (failed? o)) (= (mapv double expected) (vec (sut/coeffs o))) (= n (sut/degree o)))
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
    (t/is (and (not (failed? o)) (= (mapv double expected) (vec (sut/coeffs o))) (= n (sut/degree o)))
          (str "object, alpha " alpha " beta " beta " degree " n))))

(defn- pochhammer [a k] (reduce *' 1 (map #(+ a %) (range k))))

(defn- factorial-exact [n] (reduce *' 1 (range 1 (inc n))))

(defn- gegenbauer-explicit-coefficients
  "Coefficients, ascending, of the Gegenbauer polynomial of a rational parameter from the explicit sum
  `sum_k (-1)^k (alpha)_(n-k) / (k! (n-2k)!) (2x)^(n-2k)`."
  [n alpha]
  (reduce (fn [cs k]
            (let [m (- n (* 2 k))]
              (assoc cs m (* (if (even? k) 1 -1) (pochhammer alpha (- n k))
                             (/ (reduce *' 1 (repeat m 2)) (* (factorial-exact k) (factorial-exact m)))))))
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
