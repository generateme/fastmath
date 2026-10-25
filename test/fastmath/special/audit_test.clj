(ns fastmath.special.audit-test
  "Tests added by the `fastmath.special` audit.

  Reference values come from `mpmath` (60 digits) and are stored in
  `special/audit_reference.edn` as `{function-name [[args... ref]...]}`; the generator is
  `utils/fastmath/dev/generate_special_audit_reference.py`.
  The tolerances are the measured worst relative error of the grid times about 10."
  (:require [fastmath.special :as sut]
            [fastmath.special.gamma :as gam]
            [fastmath.special.hypergeometric :as hg]
            [fastmath.core :as m]
            [fastmath.complex :as cplx]
            [clojure.test :as t]
            [clojure.java.io :as io]
            [clojure.edn :as edn]))

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "special/audit_reference.edn")))))

;; relative tolerances = measured worst error x ~10 (the comment gives the observed worst)
(def ^:private rtol
  {"erf" 1.0e-14 ;; 1.5e-15
   "erfc" 1.0e-12 ;; 7.2e-14 at x = 26 (Commons Math)
   "inv-erf" 1.0e-14 ;; 2.3e-16
   "inv-erfc" 1.0e-13 ;; 2.3e-16
   "gamma" 1.0e-12 ;; 7.0e-14, 747 points on [-185, 171.7]
   "log-gamma" 1.0e-13 ;; 2.2e-16
   "log-gamma-1p" 1.0e-14 ;; 1.8e-16
   "inv-gamma-1pm1" 1.0e-14 ;; 2.2e-16
   "lower-incomplete-gamma" 1.0e-12 ;; 7.0e-14
   "upper-incomplete-gamma" 1.0e-12 ;; 4.0e-14
   "regularized-gamma-p" 1.0e-12 ;; 3.2e-14
   "regularized-gamma-q" 1.0e-12 ;; 1.2e-13
   "beta" 1.0e-12 ;; 4.0e-14
   "log-beta" 1.0e-13 ;; 6.1e-15
   "regularized-beta" 1.0e-12 ;; 4.8e-14
   "incomplete-beta" 1.0e-12 ;; 3.5e-14
   "digamma" 1.0e-12 ;; 1.0e-14
   "trigamma" 1.0e-13 ;; 2.1e-16
   "polygamma" 1.0e-12 ;; 1.4e-14
   "zeta" 1.0e-10 ;; 1.3e-12, Hurwitz form degrades for very negative s (documented)
   "eta" 1.0e-11 ;; 1.1e-13
   "dirichlet-beta" 1.0e-12 ;; 3.8e-14
   "bessel-J0" 1.0e-11 ;; 2.8e-13
   "bessel-J1" 1.0e-13 ;; 4.8e-15
   "bessel-J" 3.0e-11 ;; 2.3e-12 (J_2(150))
   "bessel-Y0" 1.0e-12 ;; 1.3e-14
   "bessel-Y1" 3.0e-11 ;; 1.8e-12
   "bessel-Y" 1.0e-11 ;; 2.4e-13
   "bessel-I0" 1.0e-14 ;; 2.2e-16
   "bessel-I1" 1.0e-14 ;; 4.1e-16
   "bessel-I" 1.0e-12 ;; 9.0e-14
   "bessel-K0" 1.0e-14 ;; 2.3e-16
   "bessel-K1" 1.0e-14 ;; 1.5e-16
   "bessel-K" 1.0e-12 ;; 9.1e-14
   "spherical-bessel-j" 1.0e-12 ;; 1.5e-13
   "spherical-bessel-y" 1.0e-12 ;; 1.2e-13
   "spherical-bessel-1-i" 1.0e-13 ;; 1.6e-15
   "spherical-bessel-2-i" 1.0e-13 ;; 4.8e-16
   "spherical-bessel-k" 1.0e-13 ;; 5.2e-16
   "jinc" 1.0e-9 ;; 5.9e-11, at the zero x = 1.2197 (relative to a tiny value)
   "Cin" 1.0e-14 ;; 4.7e-16 (x = 1e-12 .. 50; the old gamma + ln x - Ci(x) form had 1e-2 at 1e-6)
   "harmonic-number" 1.0e-13 ;; 1.1e-15
   "minkowski" 1.0e-11 ;; 9.6e-13 (near the golden-ratio conjugate the reference is a truncated continued fraction)
   "lambert-W" 1.0e-12 ;; 1.3e-14 (ill-conditioned near -1/e)
   "lambert-W-1" 1.0e-12 ;; 2.9e-15
   "owens-t" 1.0e-13 ;; 3.9e-15
   "airy-Ai" 1.0e-11 ;; plus the amplitude term below
   "airy-Ai'" 1.0e-11
   "airy-Bi" 1.0e-11
   "airy-Bi'" 1.0e-11
   "elliptic-K" 1.0e-13 ;; 2.1e-16
   "elliptic-E" 1.0e-13 ;; 1.8e-15
   "elliptic-F" 1.0e-13 ;; 8.1e-16
   "elliptic-D" 1.0e-12 ;; 1.1e-13
   "elliptic-PI" 1.0e-13 ;; 1.1e-15
   "elliptic-Rc" 1.0e-14 ;; 2.1e-16
   "elliptic-Rd" 1.0e-14 ;; 2.1e-16
   "elliptic-Rf" 1.0e-14 ;; 3.4e-16
   "elliptic-Rg" 1.0e-13 ;; 5.8e-15
   "elliptic-Rj" 1.0e-14 ;; 4.2e-16
   "jacobi-am" 1.0e-13 ;; 1.4e-16
   "hypergeometric-2F1" 2.0e-9 ;; 1.0e-9 (cancellation in the c-a-b near-integer transformations); the variogram grid 2F1(1/2,-k/2;3/2;x), k <= 180, is within 1e-13
   "hypergeometric-0F1" 1.0e-11 ;; 1.1e-13
   "hypergeometric-0F2" 1.0e-11 ;; 1.0e-13, x >= -50
   "kummers-M" 5.0e-9 ;; 4.6e-10, x >= -10 (the alternating series for large negative x has only absolute accuracy)
   "whittaker-M" 1.0e-13 ;; 8.4e-15
   "whittaker-W" 1.0e-11 ;; 2.2e-13
   "tricomis-U" 1.0e-10}) ;; 2.9e-12 (a > 0, integer a, or 0 < b <= 8.5 for a < 0; the docstring lists the degradation for other a < 0)

(defn- jacobi? [fname] (.startsWith ^String fname "jacobi-"))

(defn- airy-atol
  "Absolute error allowed for the Airy functions at negative `x`: the phase `2/3 |x|^1.5` is
  known to the relative accuracy of one ulp times a few, so the error grows like phase * amplitude."
  ^double [^String fname ^double x]
  (if (m/>= x -10.0)
    0.0
    (let [ax (m/abs x)
          amp (m/* (m/pow ax (if (.endsWith fname "'") 0.25 -0.25)) m/INV_SQRTPI)]
      (m/* 20.0 m/MACHINE-EPSILON (m/* m/TWO_THIRD (m/pow ax 1.5)) amp))))

(defn- close?
  [fname v r tol args]
  (cond
    (m/nan? r) (m/nan? v)
    (m/inf? r) (== v r)
    :else (<= (m/abs (m/- v r))
              (m/+ (m/* tol (m/abs r))
                   (if (.startsWith ^String fname "airy-") (airy-atol fname (double (first args))) 0.0)))))

(defn- check-reference
  "Number of reference points of `fname` that `f` does not reproduce."
  [fname]
  (let [f (ns-resolve 'fastmath.special (symbol fname))
        tol (get rtol fname)]
    (count (remove (fn [c] (let [args (butlast c)
                                 r (double (last c))
                                 v (double (apply f args))]
                             (close? fname v r tol args)))
                   (get @reference fname)))))

(t/deftest reference-grids
  (doseq [fname (sort (keys rtol))]
    (t/testing fname
      (t/is (pos? (count (get @reference fname))) (str "reference for " fname))
      (t/is (zero? (check-reference fname)) (str fname " mismatches")))))

(t/deftest hankel-reference
  ;; principal branch, including the negative real axis
  (doseq [[nm order x re im] (get @reference "hankel")
          :let [^fastmath.vector.Vec2 v (if (= nm "hankel-1") (sut/hankel-1 order x) (sut/hankel-2 order x))
                tol (* 1.0e-11 (max 1.0 (m/abs re) (m/abs im)))]]
    (t/is (<= (m/abs (- (.x v) re)) tol) (str nm " re " order " " x))
    (t/is (<= (m/abs (- (.y v) im)) tol) (str nm " im " order " " x))))

(def ^:private complex-rtol
  {"gamma-complex" 1.0e-12 ;; 3.6e-14
   "log-gamma-complex" 1.0e-12 ;; 1.3e-14, the analytic branch (Im not wrapped to (-pi, pi]) as mpmath, Julia
   "bessel-K-half-odd-complex" 1.0e-12 ;; 1.4e-14
   "bessel-K-half-odd-scaled-complex" 1.0e-12 ;; 1.5e-14
   "hypergeometric-pFq-complex" 1.0e-12 ;; 2.4e-15
   "tricomis-U-complex" 1.0e-5}) ;; 1.1e-6 (a, b complex with Re z < 0: documented limit); typically 1e-12

(defn- cv [[re im]] (cplx/complex re im))

(t/deftest complex-reference
  (doseq [[nm args ref] (get @reference "complex")
          :let [f (ns-resolve 'fastmath.special (symbol nm))
                got (case nm
                      ("gamma-complex" "log-gamma-complex") (f (cv (first args)))
                      "tricomis-U-complex" (f (cv (nth args 0)) (cv (nth args 1)) (cv (nth args 2)))
                      ("bessel-K-half-odd-complex" "bessel-K-half-odd-scaled-complex") (f (long (first args)) (cv (second args)))
                      "hypergeometric-pFq-complex" (f (mapv cv (nth args 0)) (mapv cv (nth args 1)) (cv (nth args 2))))
                r (cv ref)
                n (cplx/abs r)]]
    (t/is (<= (cplx/abs (cplx/sub got r)) (* (get complex-rtol nm) (max n 1.0e-300))) (str nm " " args))))

(t/deftest log-gamma-complex-branch
  (t/testing "the imaginary part is the analytic continuation (not wrapped), also on the half-integer lines where the sign of the zero of sin(pi z) decides the branch of the logarithm"
    (doseq [[z im] [[(cplx/complex -8.5 -2.0) 23.862683602452897] [(cplx/complex -8.5 2.0) -23.862683602452897]
                    [(cplx/complex -8.5 3.5) -20.49800369634529] [(cplx/complex -1.2 0.3) -5.3263577371851625]
                    [(cplx/complex -3.0 0.5) -10.365857636927407]]]
      (t/is (m/delta-eq im (cplx/im (sut/log-gamma-complex z)) 0.0 1.0e-12) (str z)))
    (t/is (m/delta-eq 1.0 (cplx/abs (cplx/exp (cplx/sub (sut/log-gamma-complex (cplx/complex -8.5 -2.0))
                                                         (cplx/log (sut/gamma-complex (cplx/complex -8.5 -2.0)))))) 1.0e-12))))

(t/deftest bessel-k-half-odd-complex-cut
  (t/testing "on the negative real axis the value is the limit from above (Im z = +0), as mpmath"
    (t/is (m/delta-eq -1.7432789876671936E12 (cplx/im (sut/bessel-K-half-odd-complex 9 (cplx/complex -30.0 0.0))) 0.0 1.0e-12))
    (t/is (m/delta-eq -3732.3434432428962 (cplx/im (sut/bessel-K-half-odd-complex 13 (cplx/complex -1.2 0.0))) 0.0 1.0e-12))))

;; gamma family

(t/deftest gamma-edges
  (t/is (m/nan? (sut/gamma ##NaN)))
  (t/is (m/pos-inf? (sut/gamma ##Inf)))
  (t/is (m/nan? (sut/gamma ##-Inf)))
  (t/testing "poles"
    (t/is (m/nan? (sut/gamma 0.0)))
    (t/is (m/nan? (sut/gamma -0.0)))
    (t/is (m/nan? (sut/gamma -3.0)))
    (t/is (m/nan? (sut/gamma -150.0))))
  (t/testing "overflow at 171.62, not at 142 (Commons Math)"
    (t/is (m/delta-eq 3.80892263763057E260 (sut/gamma 150.0) 0.0 1.0e-12))
    (t/is (m/pos-inf? (sut/gamma 171.7)))
    (t/is (m/pos-inf? (sut/gamma 1.0e-310))))
  (t/testing "negative large non-integers are (sub)normal, then underflow to a signed zero"
    (t/is (m/delta-eq -3.3127395215386074E-308 (sut/gamma -170.5) 0.0 1.0e-12))
    (t/is (m/zero? (sut/gamma -200.5)))))

(t/deftest gamma-above-the-overflow-point
  (t/testing "the split power gives NaN above x = 1490 without the explicit overflow"
    (doseq [x [171.7 172.0 400.0 1489.0 1491.0 2000.0 1.0e5 1.0e300 ##Inf]]
      (t/is (m/pos-inf? (sut/gamma x)) (str x)))
    (t/is (m/delta-eq -1.0 (sut/inv-gamma-1pm1 2000.0) 0.0 1.0e-15)))
  (t/testing "negative arguments beyond the underflow point"
    (doseq [x [-200.5 -1500.5 -100000.5]]
      (t/is (m/zero? (sut/gamma x)) (str x)))
    (t/is (m/nan? (sut/gamma -1.0e300)))))   ;; an integer: a pole

(t/deftest log-gamma-edges
  (t/is (m/nan? (sut/log-gamma ##NaN)))
  (t/is (m/pos-inf? (sut/log-gamma 0.0)))
  (t/is (m/pos-inf? (sut/log-gamma -3.0)))
  (t/is (m/pos-inf? (sut/log-gamma ##Inf)))
  (t/is (m/nan? (sut/log-gamma ##-Inf)))
  (t/testing "log|Gamma| for negative non-integers"
    (t/is (m/delta-eq (m/log (m/abs (sut/gamma -2.5))) (sut/log-gamma -2.5) 1.0e-14))
    (t/is (m/delta-eq -0.05624371649767405 (sut/log-gamma -2.5) 1.0e-14))))

(t/deftest log-gamma-1p-domain
  (t/testing "outside [-0.5, 1.5] falls back to the general formula"
    (t/is (m/delta-eq (sut/log-gamma 0.4) (sut/log-gamma-1p -0.6) 1.0e-15))
    (t/is (m/delta-eq (m/log 2.0) (sut/log-gamma-1p 2.0) 1.0e-15))
    (t/is (m/delta-eq -0.5 (sut/inv-gamma-1pm1 2.0) 1.0e-15))
    (t/is (m/delta-eq (- (/ 1.0 (sut/gamma 0.4)) 1.0) (sut/inv-gamma-1pm1 -0.6) 1.0e-14))))

(t/deftest erf-edges
  (t/is (m/delta-eq 1.0 (sut/erf ##Inf) 0.0))
  (t/is (m/delta-eq -1.0 (sut/erf ##-Inf) 0.0))
  (t/is (m/nan? (sut/erf ##NaN)))
  (t/is (m/delta-eq 2.0 (sut/erfc ##-Inf) 0.0))
  (t/is (m/zero? (sut/erfc ##Inf)))
  (t/is (m/zero? (sut/erfc 30.0)))
  (t/testing "two arguments: erf(x2) - erf(x1), accurate for close arguments"
    (t/is (m/delta-eq (- (sut/erf 2.0) (sut/erf 1.0)) (sut/erf 1.0 2.0) 1.0e-15))
    (t/is (m/delta-eq 0.15729920705028488 (sut/erf 1.0 ##Inf) 1.0e-15)))
  (t/testing "inverse"
    (t/is (m/pos-inf? (sut/inv-erf 1.0)))
    (t/is (m/neg-inf? (sut/inv-erf -1.0)))
    (t/is (m/nan? (sut/inv-erf 1.5)))
    (t/is (m/pos-inf? (sut/inv-erfc 0.0)))
    (t/is (m/neg-inf? (sut/inv-erfc 2.0)))
    (t/is (m/zero? (sut/inv-erfc 1.0)))
    (t/is (m/nan? (sut/inv-erfc -0.1)))
    (t/is (m/nan? (sut/inv-erfc 2.1)))
    (t/is (m/nan? (sut/inv-erfc ##NaN))))
  (t/testing "round trip over the whole range of the tail (inv-erf(1-x) cannot be used there)"
    (doseq [x [1.0e-300 1.0e-200 1.0e-100 1.0e-50 1.0e-20 1.0e-10 1.0e-5 0.01 0.3 1.0 1.7 1.99999]]
      (t/is (m/delta-eq x (sut/erfc (sut/inv-erfc x)) 0.0 1.0e-9) (str "x=" x)))))

(t/deftest incomplete-gamma-edges
  (t/testing "lower"
    (t/is (m/zero? (sut/lower-incomplete-gamma 2.0 0.0)))
    (t/is (m/delta-eq 1.0 (sut/lower-incomplete-gamma 2.0 ##Inf) 1.0e-14))
    (t/is (m/nan? (sut/lower-incomplete-gamma 2.0 -1.0)))
    (t/is (m/nan? (sut/lower-incomplete-gamma -1.0 1.0)))
    (t/is (m/nan? (sut/lower-incomplete-gamma ##NaN 1.0))))
  (t/testing "no cancellation for small x (was -944 for s = 20, x = 1e-5)"
    (t/is (m/delta-eq 4.999952381179661E-102 (sut/lower-incomplete-gamma 20.0 1.0e-5) 0.0 1.0e-12))
    (t/is (m/delta-eq 5.0E-21 (sut/lower-incomplete-gamma 2.0 1.0e-10) 0.0 1.0e-9))
    (t/is (pos? (sut/lower-incomplete-gamma 5.0 1.0e-3))))
  (t/testing "upper"
    (t/is (m/delta-eq 1.0 (sut/upper-incomplete-gamma 2.0 0.0) 1.0e-15))
    (t/is (m/zero? (sut/upper-incomplete-gamma 2.0 ##Inf)))
    (t/is (m/nan? (sut/upper-incomplete-gamma 0.5 -1.0)))
    (t/is (m/delta-eq 0.21938393439552029 (sut/upper-incomplete-gamma 0.0 1.0) 1.0e-13))) ;; E1(1)
  (t/testing "regularized"
    (t/is (m/zero? (sut/regularized-gamma-p 2.0 0.0)))
    (t/is (m/delta-eq 1.0 (sut/regularized-gamma-q 2.0 0.0) 0.0))
    (t/is (m/delta-eq 1.0 (sut/regularized-gamma-p 2.0 ##Inf) 0.0))
    (t/is (m/zero? (sut/regularized-gamma-q 2.0 ##Inf)))
    (t/is (m/nan? (sut/regularized-gamma-p 2.0 -1.0)))
    (t/is (m/nan? (sut/regularized-gamma-q 2.0 ##NaN)))
    (t/testing "a non-positive integer: P = 1, Q = 0 (the limits)"
      (t/is (m/delta-eq 1.0 (sut/regularized-gamma-p 0.0 1.0) 0.0))
      (t/is (m/delta-eq 1.0 (sut/regularized-gamma-p -2.0 1.0) 0.0))
      (t/is (m/zero? (sut/regularized-gamma-q 0.0 1.0)))
      (t/is (m/zero? (sut/regularized-gamma-q -2.0 1.0))))
    (t/testing "P + Q = 1"
      (doseq [a [0.3 1.0 7.5 100.0] x [0.01 1.0 10.0 120.0]]
        (t/is (m/delta-eq 1.0 (+ (sut/regularized-gamma-p a x) (sut/regularized-gamma-q a x)) 1.0e-13))))))

;; beta, digamma, zeta

(t/deftest beta-pole-cancellation
  (t/testing "the pole of Gamma(q) cancels the pole of Gamma(p+q): (-1)^p (p-1)! k! / n!, q = -n, p+q = -k"
    (t/is (m/delta-eq -1.0 (sut/beta -1.0 1.0) 1.0e-14))
    (t/is (m/delta-eq -1.0 (sut/beta 1.0 -1.0) 1.0e-14))
    (t/is (m/delta-eq 0.05 (sut/beta 2.0 -5.0) 1.0e-14))
    (t/is (m/delta-eq 0.05 (sut/beta -5.0 2.0) 1.0e-14))
    (t/is (m/delta-eq (/ -1.0 12.0) (sut/beta 3.0 -4.0) 1.0e-14)))
  (t/testing "genuine poles stay NaN"
    (t/is (m/nan? (sut/beta -2.0 3.0)))
    (t/is (m/nan? (sut/beta 3.0 -1.0)))
    (t/is (m/nan? (sut/beta -2.0 -3.0)))
    (t/is (m/nan? (sut/log-beta -2.0 -3.0))))
  (t/testing "log-beta is log|B|"
    (t/is (m/delta-eq (m/log 0.05) (sut/log-beta 2.0 -5.0) 1.0e-13))
    (t/is (m/delta-eq 0.0 (sut/log-beta -1.0 1.0) 1.0e-14))
    (t/is (m/neg-inf? (sut/log-beta -0.25 -3.75)))))

(t/deftest poles-of-digamma-family
  (t/testing "digamma: NaN at negative integers, one-sided limits at the signed zeros"
    (t/is (m/nan? (sut/digamma -1.0)))
    (t/is (m/nan? (sut/digamma -42.0)))
    (t/is (m/neg-inf? (sut/digamma 0.0)))
    (t/is (m/pos-inf? (sut/digamma -0.0))))
  (t/testing "trigamma: double pole, +Inf from both sides"
    (t/is (m/pos-inf? (sut/trigamma 0.0)))
    (t/is (m/pos-inf? (sut/trigamma -1.0)))
    (t/is (m/pos-inf? (sut/trigamma -7.0)))
    (t/is (m/nan? (sut/trigamma ##-Inf))))
  (t/testing "polygamma: +Inf for odd order, NaN for even order at negative integers"
    (t/is (m/pos-inf? (sut/polygamma 3 -2.0)))
    (t/is (m/pos-inf? (sut/polygamma 5 -1.0)))
    (t/is (m/nan? (sut/polygamma 2 -1.0)))
    (t/is (m/nan? (sut/polygamma 4 -3.0)))
    (t/is (m/nan? (sut/polygamma 2 ##-Inf))))
  (t/testing "order must be a non-negative integer"
    (t/is (m/nan? (sut/polygamma 1.5 1.0)))
    (t/is (m/nan? (sut/polygamma -1 1.0)))
    (t/is (m/nan? (sut/polygamma ##NaN 1.0)))
    (t/is (m/delta-eq (sut/trigamma 2.5) (sut/polygamma 1.0 2.5) 1.0e-15)))
  (t/testing "polygamma at negative half-integers (cot(pi z) is exactly 0 there)"
    (t/is (m/delta-eq -156.229596593233 (sut/polygamma 10 -1.5) 0.0 1.0e-12))
    (t/is (m/delta-eq -42108.85876897549 (sut/polygamma 10 -0.5) 0.0 1.0e-12))))

(t/deftest hurwitz-zeta-negative-integer-s
  (t/is (m/delta-eq -1.0 (sut/zeta -50.0 2.0) 1.0e-12))
  (t/is (m/delta-eq (- (sut/zeta -3.0) 1.0 8.0 27.0) (sut/zeta -3.0 4.0) 1.0e-11))) ;; zeta(-3,4) = zeta(-3) - 1^3 - 2^3 - 3^3

(t/deftest harmonic-number-edges
  (t/is (m/zero? (sut/harmonic-number 0.0)))
  (t/is (m/nan? (sut/harmonic-number -1.0)))
  (t/is (m/nan? (sut/harmonic-number -5.0)))
  (t/is (m/delta-eq 1.5 (sut/harmonic-number 2.0) 1.0e-14))
  (t/is (m/delta-eq 5.0 (sut/harmonic-number 5.0 0.0) 1.0e-14))
  (t/is (m/delta-eq 1.3611111111111112 (sut/harmonic-number 3.0 2.0) 1.0e-13)))

;; Bessel family

(t/deftest bessel-y-overflow-and-reflection
  (t/testing "true overflow gives -Inf, not NaN (the upward recurrence overflows)"
    (t/is (m/neg-inf? (sut/bessel-Y 50.0 1.0e-6)))
    (t/is (m/neg-inf? (sut/bessel-Y 100.0 0.01))))
  (t/testing "negative non-integer order: Y_{-v} = cos(v pi) Y_v + sin(v pi) J_v with exact sin/cos"
    (t/is (m/delta-eq -2.659615202675952E-10 (sut/bessel-Y -1.5 1.0e-6) 0.0 1.0e-11))
    (t/is (m/delta-eq (- (sut/bessel-J 1.5 0.01)) (sut/bessel-Y -1.5 0.01) 0.0 1.0e-11)))
  (t/testing "edges"
    (t/is (m/neg-inf? (sut/bessel-Y 0.0 0.0)))
    (t/is (m/nan? (sut/bessel-Y 0.5 -1.0)))))

(t/deftest spherical-bessel-edges
  (t/testing "limits at 0 and infinity"
    (t/is (m/delta-eq 1.0 (sut/spherical-bessel-j 0.0 0.0) 0.0))
    (t/is (m/zero? (sut/spherical-bessel-j 2.0 0.0)))
    (t/is (m/zero? (sut/spherical-bessel-j 3.0 0.0)))
    (t/is (m/neg-inf? (sut/spherical-bessel-y 0.0 0.0)))
    (t/is (m/neg-inf? (sut/spherical-bessel-y 3.0 0.0)))
    (doseq [n [0.0 1.0 2.0 3.0 6.0]]
      (t/is (m/zero? (sut/spherical-bessel-j n ##Inf)) (str "j " n))
      (t/is (m/zero? (sut/spherical-bessel-y n ##Inf)) (str "y " n))
      (t/is (m/zero? (sut/spherical-bessel-k n ##Inf)) (str "k " n))
      (t/is (m/pos-inf? (sut/spherical-bessel-k n 0.0)) (str "k(0) " n))
      (t/is (m/pos-inf? (sut/spherical-bessel-1-i n ##Inf)) (str "1-i " n))
      (t/is (m/pos-inf? (sut/spherical-bessel-2-i n ##Inf)) (str "2-i " n))))
  (t/testing "small x: power series (the closed forms cancel)"
    (t/is (m/delta-eq 6.6666190477513225E-6 (sut/spherical-bessel-j 2.0 0.01) 0.0 1.0e-13))
    (t/is (m/delta-eq 6.666714285846561E-6 (sut/spherical-bessel-1-i 2.0 0.01) 0.0 1.0e-13))
    (t/is (m/delta-eq 0.0033333000001190475 (sut/spherical-bessel-j 1.0 0.01) 0.0 1.0e-13))
    (t/is (m/delta-eq 0.0033333666667857147 (sut/spherical-bessel-1-i 1.0 0.01) 0.0 1.0e-13)))
  (t/testing "series and closed form agree at the switch"
    (doseq [n [1.0 2.0]
            x [0.4999999 0.5000001 0.9999999 1.0000001]]
      (t/is (m/delta-eq (sut/spherical-bessel-j n x) (* (m/sqrt (/ m/HALF_PI x)) (sut/bessel-J (+ n 0.5) x)) 0.0 1.0e-12))))
  (t/testing "parity for negative x (integer orders); non-integer orders are NaN"
    (doseq [n [0.0 1.0 2.0 3.0 5.0] x [0.01 0.5 3.0 20.0]
            :let [s (if (even? (long n)) 1.0 -1.0)]]
      (t/is (m/delta-eq (* s (sut/spherical-bessel-j n x)) (sut/spherical-bessel-j n (- x)) 0.0 1.0e-12))
      (t/is (m/delta-eq (* -1.0 s (sut/spherical-bessel-y n x)) (sut/spherical-bessel-y n (- x)) 0.0 1.0e-12))
      (t/is (m/delta-eq (* s (sut/spherical-bessel-1-i n x)) (sut/spherical-bessel-1-i n (- x)) 0.0 1.0e-12))
      (t/is (m/delta-eq (* -1.0 s (sut/spherical-bessel-2-i n x)) (sut/spherical-bessel-2-i n (- x)) 0.0 1.0e-12)))
    (t/is (m/nan? (sut/spherical-bessel-j 1.5 -1.0)))
    (t/is (m/nan? (sut/spherical-bessel-k 1.0 -1.0)))
    (t/is (m/nan? (sut/spherical-bessel-k 4.0 -1.0)))))

(t/deftest spherical-hankel
  (doseq [n [0.0 1.0 2.0 4.0] x [0.5 3.0 12.0]
          :let [h1 (sut/spherical-hankel-1 n x)
                h2 (sut/spherical-hankel-2 n x)]]
    (t/is (m/delta-eq (sut/spherical-bessel-j n x) (.x h1) 1.0e-15))
    (t/is (m/delta-eq (sut/spherical-bessel-y n x) (.y h1) 1.0e-15))
    (t/is (m/delta-eq (.x h1) (.x h2) 1.0e-15))
    (t/is (m/delta-eq (- (.y h1)) (.y h2) 1.0e-15))))

(t/deftest bessel-k-half-odd
  (t/testing "matches bessel-K, K is even in the order"
    (doseq [n [1 3 5 7 9 -3 -7] x [0.1 1.0 5.0]]
      (t/is (m/delta-eq (sut/bessel-K (/ (m/abs n) 2.0) x) (sut/bessel-K-half-odd n x) 0.0 1.0e-13))
      (t/is (m/delta-eq (* (m/exp x) (sut/bessel-K (/ (m/abs n) 2.0) x)) (sut/bessel-K-half-odd-scaled n x) 0.0 1.0e-13))))
  (t/testing "even order is a caller error"
    (doseq [n [0 2 4 -2]]
      (t/is (thrown? IllegalArgumentException (sut/bessel-K-half-odd n 1.0)))
      (t/is (thrown? IllegalArgumentException (sut/bessel-K-half-odd-scaled n 1.0)))
      (t/is (thrown? IllegalArgumentException (sut/bessel-K-half-odd-complex n (cplx/complex 1.0 0.5))))
      (t/is (thrown? IllegalArgumentException (sut/bessel-K-half-odd-scaled-complex n (cplx/complex 1.0 0.5))))))
  (t/testing "complex forms agree with the real ones on the real axis"
    (doseq [n [1 3 5 9] x [0.3 2.0]]
      (t/is (m/delta-eq (sut/bessel-K-half-odd n x) (.x (sut/bessel-K-half-odd-complex n (cplx/complex x 0.0))) 0.0 1.0e-13))
      (t/is (m/delta-eq (sut/bessel-K-half-odd-scaled n x) (.x (sut/bessel-K-half-odd-scaled-complex n (cplx/complex x 0.0))) 0.0 1.0e-13)))))

(t/deftest jinc-definition
  (t/is (m/delta-eq 1.0 (sut/jinc 0.0) 0.0))
  (t/is (m/zero? (sut/jinc ##Inf)))
  (t/is (m/delta-eq 0.18119175498741528 (sut/jinc 1.0) 1.0e-14)) ;; 2 J1(pi)/pi
  (t/is (m/delta-eq (sut/jinc 0.7) (sut/jinc -0.7) 0.0))
  (t/is (m/delta-eq 0.0 (sut/jinc 1.2196698912665045) 1.0e-12))) ;; first zero: j_{1,1}/pi

(t/deftest spherical-bessel-fixed-orders
  (t/testing "the closed forms of orders 0, 1, 2 agree with the Bessel-function definition"
    (doseq [[n j y i1 i2 k] [[0 sut/spherical-bessel-j0 sut/spherical-bessel-y0 sut/spherical-bessel-1-i0 sut/spherical-bessel-2-i0 sut/spherical-bessel-k0]
                             [1 sut/spherical-bessel-j1 sut/spherical-bessel-y1 sut/spherical-bessel-1-i1 sut/spherical-bessel-2-i1 sut/spherical-bessel-k1]
                             [2 sut/spherical-bessel-j2 sut/spherical-bessel-y2 sut/spherical-bessel-1-i2 sut/spherical-bessel-2-i2 sut/spherical-bessel-k2]]
            x [0.001 0.3 0.7 1.5 4.0 12.0]
            :let [f (m/sqrt (/ m/HALF_PI x))
                  close (fn [a b] (m/delta-eq a b 0.0 1.0e-11))]]
      (t/is (close (* f (sut/bessel-J (+ n 0.5) x)) (j x)) (str "j" n " " x))
      (t/is (close (* f (sut/bessel-Y (+ n 0.5) x)) (y x)) (str "y" n " " x))
      (t/is (close (* f (sut/bessel-I (+ n 0.5) x)) (i1 x)) (str "1-i" n " " x))
      (t/is (close (* f (sut/bessel-I (- (+ n 0.5)) x)) (i2 x)) (str "2-i" n " " x))
      (t/is (close (* f (sut/bessel-K (+ n 0.5) x)) (k x)) (str "k" n " " x)))))

;; Airy

(t/deftest airy-edges
  (t/is (m/zero? (sut/airy-Ai ##Inf)))
  (t/is (m/pos-inf? (sut/airy-Bi ##Inf)))
  (t/is (m/pos-inf? (sut/airy-Bi 800.0)))
  (t/is (m/nan? (sut/airy-Ai ##NaN)))
  (t/testing "x <= -1e8 is documented as NaN (loss of the phase)"
    (t/is (m/nan? (sut/airy-Ai -1.0e8)))
    (t/is (m/nan? (sut/airy-Bi ##-Inf))))
  (t/testing "Wronskian Ai Bi' - Ai' Bi = 1/pi"
    (doseq [x [-20.0 -3.0 0.0 2.0 10.0]]
      (t/is (m/delta-eq m/INV_PI (- (* (sut/airy-Ai x) (sut/airy-Bi' x)) (* (sut/airy-Ai' x) (sut/airy-Bi x))) 1.0e-12 1.0e-12)))))

;; elliptic

(t/deftest elliptic-D-limits
  (t/is (m/delta-eq m/QUARTER_PI (sut/elliptic-D 0.0) 1.0e-15))
  (t/is (m/pos-inf? (sut/elliptic-D 1.0)))
  (t/is (m/nan? (sut/elliptic-D 1.5)))
  (t/testing "m -> 0: D(phi,0) = phi/2 - sin(2 phi)/4"
    (doseq [phi [0.3 1.0 1.5 -0.7 2.0 5.0]
            :let [ref (- (/ phi 2.0) (/ (Math/sin (* 2.0 phi)) 4.0))]]
      (t/is (m/delta-eq ref (sut/elliptic-D phi 0.0) 1.0e-14))))
  (t/testing "series and quotient agree at the switch |m| = 0.5"
    (doseq [phi [0.5 1.2] m [0.4999999 -0.4999999]]
      (t/is (m/delta-eq (sut/elliptic-D phi m) (sut/elliptic-D phi (+ m 2.0e-7)) 1.0e-6))))
  (t/testing "quasi-periodic: D(phi + k pi) = D(phi) + 2 k D(pi/2)"
    (t/is (m/delta-eq (+ (sut/elliptic-D 0.4 0.3) (* 4.0 (sut/elliptic-D 0.3))) (sut/elliptic-D (+ 0.4 (* 2.0 m/PI)) 0.3) 1.0e-13))))

(t/deftest elliptic-PI-edges
  (t/is (m/nan? (sut/elliptic-PI 1.0 0.5))) ;; pole n = 1
  (t/is (m/nan? (sut/elliptic-PI 1.0 0.3 0.5)))
  (t/is (m/nan? (sut/elliptic-PI 2.0 0.5))) ;; n > 1
  (t/is (m/nan? (sut/elliptic-PI 0.3 1.5)))
  (t/is (m/pos-inf? (sut/elliptic-PI 0.3 1.0)))
  (t/is (m/delta-eq (sut/elliptic-K 0.5) (sut/elliptic-PI 0.0 0.5) 1.0e-14))
  (t/testing "incomplete with m = 1 is finite below pi/2 (was Inf)"
    (t/is (m/delta-eq -1.3616887091204068 (sut/elliptic-PI 0.3 -1.0 1.0) 1.0e-13)))
  (t/testing "quasi-periodicity in phi"
    (doseq [phi [0.7 -1.1 3.0] k [-2 1 3]]
      (t/is (m/delta-eq (+ (sut/elliptic-PI 0.3 phi 0.4) (* 2.0 k (sut/elliptic-PI 0.3 0.4)))
                        (sut/elliptic-PI 0.3 (+ phi (* k m/PI)) 0.4) 1.0e-12)))))

(t/deftest elliptic-edges
  (t/is (m/pos-inf? (sut/elliptic-K 1.0)))
  (t/is (m/nan? (sut/elliptic-K 1.5)))
  (t/is (m/delta-eq 1.0 (sut/elliptic-E 1.0) 1.0e-15))
  (t/is (m/nan? (sut/elliptic-E 2.0)))
  (t/testing "Rc"
    (t/is (m/pos-inf? (sut/elliptic-Rc 1.0 0.0))) ;; was an exception
    (t/is (m/pos-inf? (sut/elliptic-Rc 0.0 0.0)))
    (t/is (m/nan? (sut/elliptic-Rc -1.0 1.0)))
    (t/is (m/delta-eq (m// 1.0 (m/sqrt 3.0)) (sut/elliptic-Rc 3.0 3.0) 1.0e-15)))
  (t/testing "Rf, Rd, Rj, Rg"
    (t/is (m/pos-inf? (sut/elliptic-Rf 0.0 0.0 1.0)))
    (t/is (m/nan? (sut/elliptic-Rf -1.0 1.0 1.0)))))

(t/deftest jacobi-inverse-roundtrip
  (let [m 0.5
        pairs [[sut/jacobi-sn sut/jacobi-asn] [sut/jacobi-cn sut/jacobi-acn] [sut/jacobi-dn sut/jacobi-adn]
               [sut/jacobi-sc sut/jacobi-asc] [sut/jacobi-sd sut/jacobi-asd] [sut/jacobi-cs sut/jacobi-acs]
               [sut/jacobi-cd sut/jacobi-acd] [sut/jacobi-ds sut/jacobi-ads] [sut/jacobi-dc sut/jacobi-adc]
               [sut/jacobi-ns sut/jacobi-ans] [sut/jacobi-nc sut/jacobi-anc] [sut/jacobi-nd sut/jacobi-and]]]
    (doseq [[f inv] pairs u [0.2 0.5 0.9 1.4]]
      (t/is (m/delta-eq u (inv (f u m) m) 1.0e-9 1.0e-9) (str inv " " u)))))

(t/deftest jacobi-identities
  (doseq [u [0.3 1.1 2.7 -0.8] m [-1.0 0.2 0.8 1.5]
          :let [sn (sut/jacobi-sn u m) cn (sut/jacobi-cn u m) dn (sut/jacobi-dn u m)]]
    (t/is (m/delta-eq 1.0 (+ (* sn sn) (* cn cn)) 1.0e-13))
    (t/is (m/delta-eq (- 1.0 (* m sn sn)) (* dn dn) 1.0e-13))))

;; Lambert, Minkowski, Owen

(t/deftest lambert-w-edges
  (t/is (m/pos-inf? (sut/lambert-W ##Inf)))
  (t/is (m/nan? (sut/lambert-W ##NaN)))
  (t/is (m/nan? (sut/lambert-W -0.5)))
  (t/is (m/delta-eq -1.0 (sut/lambert-W (- (m/exp -1.0))) 0.0))
  (t/is (m/zero? (sut/lambert-W 0.0)))
  (t/is (m/neg-inf? (sut/lambert-W-1 0.0)))
  (t/is (m/neg-inf? (sut/lambert-W-1 -0.0)))
  (t/is (m/nan? (sut/lambert-W-1 0.5)))
  (t/is (m/nan? (sut/lambert-W-1 -0.5)))
  (t/testing "W(x e^x) = x"
    (doseq [x [-0.9 -0.1 0.5 1.0 5.0 40.0]]
      (t/is (m/delta-eq x (sut/lambert-W (* x (m/exp x))) 1.0e-12 1.0e-12) (str "W0 " x)))
    (doseq [x [-1.1 -2.0 -5.0 -20.0]]
      (t/is (m/delta-eq x (sut/lambert-W-1 (* x (m/exp x))) 1.0e-12 1.0e-12) (str "W-1 " x)))))

(t/deftest minkowski-function
  (t/testing "known values"
    (t/is (m/zero? (sut/minkowski 0.0)))
    (t/is (m/delta-eq 1.0 (sut/minkowski 1.0) 1.0e-15))
    (t/is (m/delta-eq 0.5 (sut/minkowski 0.5) 1.0e-15))
    (t/is (m/delta-eq 0.25 (sut/minkowski (/ 1.0 3.0)) 1.0e-12))
    (t/is (m/delta-eq 0.75 (sut/minkowski (/ 2.0 3.0)) 1.0e-12)))
  (t/testing "extended to all reals by ?(x+1) = ?(x) + 1"
    (t/is (m/delta-eq 1.5 (sut/minkowski 1.5) 1.0e-15))
    (t/is (m/delta-eq -0.21875 (sut/minkowski -0.3) 1.0e-12))
    (t/is (m/delta-eq 10.125 (sut/minkowski 10.25) 1.0e-15))
    (t/is (m/delta-eq (+ (sut/minkowski 0.3) 7.0) (sut/minkowski 7.3) 1.0e-12))
    (t/is (m/delta-eq (- (sut/minkowski 0.3) 2.0) (sut/minkowski -1.7) 1.0e-12)))
  (t/testing "special values and monotonicity"
    (t/is (m/pos-inf? (sut/minkowski ##Inf)))
    (t/is (m/neg-inf? (sut/minkowski ##-Inf)))
    (t/is (m/nan? (sut/minkowski ##NaN)))
    (let [xs (range -3.0 3.0 0.0371)
          ys (map sut/minkowski xs)]
      (t/is (every? true? (map <= ys (rest ys)))))))

(t/deftest owens-t-properties
  (t/is (m/zero? (sut/owens-t 1.0 0.0)))
  (t/is (m/nan? (sut/owens-t ##NaN 1.0)))
  (t/testing "T(0, a) = atan(a) / (2 pi)"
    (doseq [a [0.1 1.0 3.0 -2.0]]
      (t/is (m/delta-eq (/ (Math/atan a) m/TWO_PI) (sut/owens-t 0.0 a) 1.0e-15))))
  (t/testing "T(h, +-Inf) = +-erfc(|h|/sqrt 2)/4, T(+-Inf, a) = 0"
    (doseq [h [0.0 0.5 3.0 -2.0]]
      (t/is (m/delta-eq (* 0.25 (sut/erfc (/ (m/abs h) m/SQRT2))) (sut/owens-t h ##Inf) 1.0e-14))
      (t/is (m/delta-eq (* -0.25 (sut/erfc (/ (m/abs h) m/SQRT2))) (sut/owens-t h ##-Inf) 1.0e-14)))
    (t/is (m/zero? (sut/owens-t ##Inf 2.0)))
    (t/is (m/zero? (sut/owens-t ##-Inf 0.5))))
  (t/testing "even in h, odd in a"
    (t/is (m/delta-eq (sut/owens-t 1.3 0.7) (sut/owens-t -1.3 0.7) 0.0))
    (t/is (m/delta-eq (sut/owens-t 1.3 0.7) (- (sut/owens-t 1.3 -0.7)) 0.0)))
  (t/testing "far tail keeps the relative accuracy (the lower-tail formula lost it): T(10,2) = 3.8099e-24"
    (t/is (m/delta-eq 3.809926512080263E-24 (sut/owens-t 10.0 2.0) 0.0 1.0e-12))
    (t/is (pos? (sut/owens-t 10.0 2.0)))
    (t/is (m/delta-eq 1.4332578593959696E-7 (sut/owens-t 5.0 2.0) 0.0 1.0e-12)))
  (t/testing "T(h,a) + T(ah,1/a) = (Phi(h) + Phi(ah))/2 - Phi(h) Phi(ah) for h >= 0"
    (doseq [h [0.2 1.0 2.5] a [0.4 2.0 5.0]
            :let [ph (* 0.5 (sut/erfc (/ (- h) m/SQRT2)))
                  pah (* 0.5 (sut/erfc (/ (- (* a h)) m/SQRT2)))]]
      (t/is (m/delta-eq (- (+ (* 0.5 ph) (* 0.5 pah)) (* ph pah))
                        (+ (sut/owens-t h a) (sut/owens-t (* a h) (/ 1.0 a))) 1.0e-14)))))

;; hypergeometric 2F1

(t/deftest hypergeometric-2F1-regressions
  (t/testing "c-a-b < 0 with c-b a non-positive integer: the Euler transformation terminates (was NaN)"
    (t/is (m/delta-eq 483.3197655095485 (sut/hypergeometric-2F1 0.1 20.0 2.0 0.5) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq 483.3197655095485 (sut/hypergeometric-2F1 20.0 0.1 2.0 0.5) 1.0e-12 1.0e-12))
        (t/is (m/delta-eq 4045969.145305963 (sut/hypergeometric-2F1 30.0 0.5 2.0 0.5) 1.0e-12 1.0e-12)))
  (t/testing "|x| <= 0.72 is the Maclaurin series"
    (t/is (m/delta-eq 71.57554993740976 (sut/hypergeometric-2F1 10.0 20.0 30.0 0.5) 1.0e-11 1.0e-11)))
  (t/testing "integer c-a-b uses the accurate digamma (Commons Math digamma(1) is off by 3e-9)"
    (t/is (m/delta-eq 1.145289847807781 (sut/hypergeometric-2F1 2.5 1.5 15.0 0.5) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq 3.3525532382524896 (sut/hypergeometric-2F1 5.0 6.0 15.0 0.5) 1.0e-9 1.0e-9))
    (t/is (m/delta-eq 3.722811731536524 (sut/hypergeometric-2F1 8.0 4.0 15.0 0.5) 1.0e-8 1.0e-8)))
  (t/testing "negative non-integer a: Gamma(a+m+eps) at a pole means 1/Gamma = 0, not NaN"
    (t/is (m/delta-eq 1.0606601717798212 (sut/hypergeometric-2F1 -2.5 0.5 -0.5 0.5) 1.0e-12 1.0e-12))
    (t/is (m/delta-eq 0.37565942021996457 (sut/hypergeometric-2F1 -2.5 0.5 -0.5 0.8) 1.0e-9 1.0e-9))
    (t/is (m/delta-eq 0.1535963434938927 (sut/hypergeometric-2F1 -2.5 -1.5 -3.5 0.9) 1.0e-9 1.0e-9))))

(t/deftest gamma-ratio-poles
  (t/testing "a pole gives the plain Inf quotient (no log-space fallback, whose sign of a pole is undefined)"
    (t/is (m/pos-inf? (gam/gamma-ratio [0.0 5.0] [1.0 1.0])))
    (t/is (m/zero? (gam/gamma-ratio [5.0 1.0] [-2.0 1.0])))
    (t/is (m/nan? (gam/gamma-ratio [-1.0] [-2.0])))))

(t/deftest beta-tiny-arguments
  (t/testing "Commons Math logBeta overflows for tiny arguments"
    (t/is (m/delta-eq 2.0e200 (sut/beta 1.0e-200 1.0e-200) 0.0 1.0e-12))
        (t/is (m/delta-eq 461.210165779369 (sut/log-beta 1.0e-200 1.0e-200) 0.0 1.0e-13))))

(t/deftest hurwitz-zeta-large-integer-offset
  (t/testing "the exact primitive sum is used for an offset up to 1000"
    (t/is (m/delta-eq (- (sut/zeta -3.0) (reduce + (map #(Math/pow % 3.0) (range 1 1000)))) (sut/zeta -3.0 1000.0) 0.0 1.0e-10))
    (t/is (m/valid-double? (sut/zeta -3.0 5000.0)))))

(t/deftest hypergeometric-2F1-terminating-condition-number
  (t/testing "x <= 0.5 and c' > -n: the direct alternating sum lost 1.8e-10 (checked vs mpmath)"
    (t/is (m/delta-eq 2.4322779220811614 (sut/hypergeometric-2F1 -40.0 3.0 -2.5 0.2) 0.0 1.0e-12))))

(t/deftest hypergeometric-2F1-terminating-high-degree
  (t/testing "terminating series 2F1(1/2, -n; 3/2; x): the direct alternating sum loses everything (x = 0.5: 4e-3 at n = 89); fixed by the transformation to 1-x"
    (t/is (m/delta-eq 0.13229452383074777 (sut/hypergeometric-2F1 0.5 -89.0 1.5 0.5) 0.0 1.0e-12))
    (t/is (m/delta-eq 0.09307607390171233 (sut/hypergeometric-2F1 0.5 -90.0 1.5 0.999) 0.0 1.0e-12))
    (t/is (m/delta-eq 0.0977672141388646 (sut/hypergeometric-2F1 0.5 -81.5 1.5 0.999) 0.0 1.0e-12))
    (t/is (m/delta-eq (sut/hypergeometric-2F1 -89.0 0.5 1.5 0.7) (sut/hypergeometric-2F1 0.5 -89.0 1.5 0.7) 0.0 1.0e-15)))
  (t/testing "boundaries of the transformation: x = 0, x = 1, n = 0, negative x, c a pole of the transformed series"
    (t/is (m/delta-eq 1.0 (sut/hypergeometric-2F1 0.5 -7.0 1.5 0.0) 1.0e-15))
    (t/is (m/delta-eq 1.0 (sut/hypergeometric-2F1 -7.0 2.0 3.0 0.0) 1.0e-15))
    (t/is (m/delta-eq 1.0 (sut/hypergeometric-2F1 0.0 2.0 3.0 0.5) 1.0e-15))
    (t/is (m/delta-eq (* (/ (* (- 3.0 0.5) (- 3.0 0.5 -1.0)) (* 3.0 4.0))) (sut/hypergeometric-2F1 -2.0 0.5 3.0 1.0) 1.0e-14)) ;; Chu-Vandermonde: (c-b)_n / (c)_n
    (t/is (m/delta-eq 32.0 (sut/hypergeometric-2F1 -3.0 1.5 0.5 -1.0) 0.0 1.0e-13)))) ;; negative x: direct sum 1 + 9 + 15 + 7

(t/deftest pfq-variants
  (t/testing "ratio (exact rational arithmetic) and the Maclaurin series agree with mpmath"
    (t/is (m/delta-eq 1.0503759676008382 (sut/hypergeometric-pFq-ratio [1/2 1/3] [3/4] 1/5) 1.0e-14))
    (t/is (m/delta-eq 1.0503759676008382 (sut/hypergeometric-pFq [0.5 (/ 1.0 3.0)] [0.75] 0.2) 1.0e-13)))
  (t/testing "with-reason variants return the value and why the iteration stopped"
    (let [r (hg/hypergeometric-pFq-weniger-with-reason [0.5 0.3] [0.75] 0.2)]
      (t/is (= :converged (:reason r)))
      (t/is (m/delta-eq 1.0451929019609443 (:value r) 1.0e-13)))
    (let [r (hg/hypergeometric-pFq-weniger-complex-with-reason [(cplx/complex 0.5 0.0) (cplx/complex 0.3 0.0)] [(cplx/complex 0.75 0.0)] (cplx/complex 0.2 0.0))]
      (t/is (= :converged (:reason r)))
      (t/is (m/delta-eq 1.0451929019609443 (.x ^fastmath.vector.Vec2 (:value r)) 1.0e-13))
      (t/is (m/delta-eq 0.0 (.y ^fastmath.vector.Vec2 (:value r)) 1.0e-13))))
  (t/testing "Maclaurin on complex numbers"
    (let [v (hg/hypergeometric-pFq-maclaurin-complex [(cplx/complex 0.5 0.0) (cplx/complex 0.3 0.0)] [(cplx/complex 0.75 0.0)] (cplx/complex 0.2 0.0))]
      (t/is (m/delta-eq 1.0451929019609443 (.x ^fastmath.vector.Vec2 v) 1.0e-13)))))

;; native gamma namespace

(t/deftest gamma-helpers
  (t/testing "sin-pi / cos-pi are exact at integers / half-integers"
    (doseq [k (range -5 6)]
      (t/is (m/zero? (gam/sin-pi (double k))))
      (t/is (m/zero? (gam/cos-pi (+ k 0.5)))))
    (t/is (m/delta-eq 1.0 (gam/sin-pi 0.5) 1.0e-15))
    (t/is (m/delta-eq -1.0 (gam/sin-pi 1.5) 1.0e-15))
    (t/is (m/delta-eq -1.0 (gam/cos-pi 1.0) 1.0e-15))
    (t/is (m/neg? (m/copy-sign 1.0 (gam/sin-pi -0.0)))))
  (t/testing "unsafe-gamma is +Inf at the poles (1/Gamma = 0)"
    (t/is (m/pos-inf? (gam/unsafe-gamma 0.0)))
    (t/is (m/pos-inf? (gam/unsafe-gamma -4.0)))
    (t/is (m/delta-eq 2.0 (gam/unsafe-gamma 3.0) 1.0e-14)))
  (t/testing "gamma-ratio falls back to log space when the products overflow"
    (t/is (m/delta-eq 150.0 (gam/gamma-ratio [151.0 200.0] [150.0 200.0]) 0.0 1.0e-12))
    (t/is (m/delta-eq (/ 1.0 150.0) (gam/gamma-ratio [150.0 200.0] [151.0 200.0]) 0.0 1.0e-12))))
