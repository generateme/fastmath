(ns fastmath.kernel.vector-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.kernel :as k]
            [fastmath.kernel.vector :as vk]
            [fastmath.random :as r]
            [fastmath.vector :as v])
  (:import [org.apache.commons.math3.linear Array2DRowRealMatrix EigenDecomposition]))

;; Reference data: `utils/fastmath/dev/generate_kernel_vector_reference.py` (scikit-learn 1.9, scipy 1.18,
;; mpmath; file `kernel_vector_reference.edn`) and `generate_kernel_vector_kernlab_reference.R` (kernlab 0.9.33;
;; file `kernel_vector_kernlab_reference.edn`). The reference cases are built through the registry
;; (`fastmath.kernel/kernel`) because they use the presets (`:matern-32`) too.

(def ^:private reference (delay (edn/read-string (slurp (io/resource "kernel/kernel_vector_reference.edn")))))
(def ^:private kernlab (delay (edn/read-string (slurp (io/resource "kernel/kernel_vector_kernlab_reference.edn")))))

(defn- close?
  "Absolute error for `b` below 1, relative otherwise; false for NaN."
  [a b tol]
  (<= (m/abs (- (double a) (double b))) (* (double tol) (m/max 1.0 (m/abs (double b))))))

(defn- case-kernel
  [{:keys [kernel params]}]
  (k/kernel kernel (cond-> params (:distance params) (update :distance v/distances))))

(defn- case-error
  "Maximum error (as in `close?`) of the case, NaN if any value is NaN."
  ^double [{:keys [x y expected] :as c}]
  (let [kf (case-kernel c)]
    (reduce (fn [^double acc [xi yi e]]
              (let [a (double (kf xi yi))
                    err (/ (m/abs (- a (double e))) (m/max 1.0 (m/abs (double e))))]
                (if (or (m/nan? err) (> err acc)) err acc)))
            0.0 (map vector x y expected))))

;; Measured maximum errors (scikit-learn / scipy / mpmath / kernlab, all cases below): everything at or below 2.2e-15
;; except :b-spline 6.5e-13 (alternating sum of a degree 11 spline), :bessel 1.0e-14 (mpmath besselj vs the
;; implementation of the Bessel function), :geometric 5.0e-15, :dirichlet 1.8e-15.
(def ^:private tolerance {:b-spline 2e-12 :bessel 1e-13 :geometric 1e-13})

(deftest reference-test
  (testing "scikit-learn, scipy and mpmath"
    (let [cases (:cases @reference)]
      (is (> (count cases) 100))
      (is (= 33 (count (distinct (map :kernel cases)))))
      (doseq [c cases
              :let [err (case-error c)
                    tol (get tolerance (:kernel c) 1e-14)]]
        (is (<= err tol) (str (:kernel c) " " (:params c) " " (:oracle c) ": " err)))))
  (testing "kernlab"
    (let [cases (:cases @kernlab)]
      (is (= 17 (count cases)))
      (doseq [c cases
              :let [err (case-error c)
                    ;; kernlab laplacedot computes the distance from expanded dot products: error 1.1e-10 at identical points
                    tol (if (= :laplacian (:kernel c)) 1e-9 1e-14)]]
        (is (<= err tol) (str (:kernel c) " " (:params c) " " (:oracle c) ": " err))))))

(deftest constructors-test
  (testing "no argument, nil and an empty map give the same kernel"
    (doseq [[nm f] {:linear vk/linear :polynomial vk/polynomial :gaussian vk/gaussian :exponential vk/exponential
                    :laplacian vk/laplacian :anova vk/anova :hyperbolic-tangent vk/hyperbolic-tangent
                    :rational-quadratic vk/rational-quadratic :multiquadratic vk/multiquadratic
                    :inverse-multiquadratic vk/inverse-multiquadratic :wave vk/wave :periodic vk/periodic
                    :power vk/power :log vk/log :spline vk/spline :b-spline vk/b-spline :bessel vk/bessel
                    :bessel2 vk/bessel2 :cauchy vk/cauchy :chi-square vk/chi-square :chi-square2 vk/chi-square2
                    :histogram vk/histogram :generalized-histogram vk/generalized-histogram
                    :generalized-t-student vk/generalized-t-student :dirichlet vk/dirichlet
                    :hellinger vk/hellinger :pearson vk/pearson :hyperbolic-secant vk/hyperbolic-secant
                    :matern vk/matern}
            :let [x [0.2 0.7] y [0.5 0.1]
                  a ((f) x y) b ((f nil) x y) c ((f {}) x y)]]
      (is (= a b c) (str nm " " a " " b " " c)))))

(deftest pearson-test
  (testing "0.5 at half of sigma for every omega (sigma is the full width at half maximum)"
    (is (close? ((vk/pearson {:sigma 1.0 :omega 2.0}) [0.0] [0.5]) 0.5 1e-15))
    (doseq [sigma [0.3 1.0 4.0] omega [0.5 1.0 2.0 5.0 20.0]]
      (is (close? ((vk/pearson {:sigma sigma :omega omega}) [0.0] [(* 0.5 sigma)]) 0.5 2e-15) (str sigma " " omega))))
  (testing "the form of Ustun (2006) and Weka Puk, measured error 2.2e-16"
    (doseq [sigma [0.5 1.0 2.5] omega [0.5 1.0 2.0 3.0 7.0] d [0.0 0.1 0.5 1.0 3.0]
            :let [t (/ (* 2.0 d (Math/sqrt (- (Math/pow 2.0 (/ 1.0 omega)) 1.0))) sigma)
                  expected (/ 1.0 (Math/pow (+ 1.0 (* t t)) omega))]]
      (is (close? ((vk/pearson {:sigma sigma :omega omega}) [0.0] [d]) expected 1e-15) (str sigma " " omega " " d))))
  (testing "omega 1 is the Cauchy kernel with sigma / 2"
    (doseq [sigma [0.5 2.0 3.0] d [0.0 0.4 2.0]]
      (is (close? ((vk/pearson {:sigma sigma}) [0.0] [d]) ((vk/cauchy {:sigma (* 0.5 sigma)}) [0.0] [d]) 1e-15))))
  (testing "defaults and K(x, x)"
    (is (= 1.0 ((vk/pearson) [1.0 2.0] [1.0 2.0])))))

(deftest geometric-test
  (testing "K(0) = 1 for every dimension and radius"
    (doseq [n (range 1 13) rad [1.0 0.3 5.0]]
      (is (close? ((vk/geometric {:n n :r rad}) [0.0] [0.0]) 1.0 1e-15) (str n " " rad))))
  (testing "values of the intersection volume of two balls (scipy betainc, mpmath, 1e-16)"
    ;; the old values: n 4 at 0.5: 0.6369 (ratio 1.08), n 4 at 1.9: 0.0325 (31x), n 6 at 1.9: 0.155 (1774x)
    (doseq [[n d expected] [[4 0.5 0.5887244480896832] [4 1.9 0.0010450705922466436] [4 1.99 3.3862153069294155e-06]
                            [6 0.5 0.5164895523012262] [6 1.9 8.762523965086255e-05] [6 1.99 2.8961364057779423e-08]]
            :let [a ((vk/geometric {:n n}) [0.0] [d])]]
      (is (< (m/abs (- a expected)) 1e-14) (str n " " d " " a))
      (is (< (m/abs (- a expected)) (* 1e-8 (m/abs expected)))  (str "relative " n " " d " " a))))
  (testing "triangular, circular and spherical closed forms"
    (doseq [s [0.0 0.1 0.45 0.8 0.999] :let [d (* 2.0 s)]]
      (is (close? ((vk/geometric {:n 1}) [0.0] [d]) (- 1.0 s) 1e-15))
      (is (close? ((vk/geometric {:n 2}) [0.0] [d]) (* (/ 2.0 m/PI) (- (m/acos s) (* s (m/sqrt (- 1.0 (* s s)))))) 1e-15))
      (is (close? ((vk/geometric {:n 3}) [0.0] [d]) (+ (- 1.0 (* 1.5 s)) (* 0.5 s s s)) 1e-15))))
  (testing "n 2 is normalised by 2/pi: not pi/2 at zero"
    (is (close? ((vk/geometric {:n 2}) [0.0] [0.0]) 1.0 1e-15))
    (is (close? ((vk/geometric {:n 2}) [0.0] [0.5]) 0.6850 1e-4)))
  (testing "support radius is 2 r, exactly zero at and beyond"
    (doseq [n [1 2 3 4 5 8] rad [1.0 0.7 3.0]
            :let [kf (vk/geometric {:n n :r rad})]]
      (is (pos? (kf [0.0] [(* 1.99 rad)])) (str n " " rad))
      (is (== 0.0 (kf [0.0] [(* 2.0 rad)])))
      (is (== 0.0 (kf [0.0] [(* 2.5 rad)])))
      (is (== 0.0 (kf [0.0] [1e300])))))
  (testing "non-negative, non-increasing in the distance"
    (doseq [n [1 2 3 4 5 6 7 9 12 20]
            :let [kf (vk/geometric {:n n})
                  vals (map #(kf [0.0] [(* 0.01 %)]) (range 0 220))]]
      (is (every? #(>= (double %) 0.0) vals) (str n))
      (is (every? (fn [[a b]] (<= (double b) (+ (double a) 1e-15))) (partition 2 1 vals)) (str n))))
  (testing "numbers, a custom distance and symmetry"
    (is (= ((vk/geometric {:n 3}) 0.0 0.5) ((vk/geometric {:n 3}) [0.0] [0.5])))
    (is (= ((vk/geometric {:n 4}) 0.3 0.1) ((vk/geometric {:n 4}) 0.1 0.3)))
    (is (close? ((vk/geometric {:n 1 :distance v/dist-abs}) [0.0 0.0] [0.5 0.5]) 0.5 1e-15))
    (is (close? ((vk/geometric {:n 1 :distance v/dist-abs}) [0.0 0.0] [0.25 0.25]) 0.75 1e-15)))
  (testing "NaN distance gives NaN, empty vectors have distance zero"
    (is (m/nan? ((vk/geometric {:n 4}) [##NaN] [0.0])))
    (is (== 1.0 ((vk/geometric {:n 4}) [] []))))
  (testing "the dimension is required and must be an integer >= 1"
    (doseq [params [nil {} {:r 1.0} {:n nil}]]
      (let [e (try (vk/geometric params) nil (catch clojure.lang.ExceptionInfo e e))]
        (is (some? e) (str params))
        (is (= :n (:option (ex-data e))))
        (is (= :invalid-input (:type (ex-data e))))))
    (doseq [n [0 -1 2.5 "3" :a ##NaN]]
      (is (= :n (:option (ex-data (try (vk/geometric {:n n}) (catch clojure.lang.ExceptionInfo e e))))) (str n)))
    (is (= ((vk/geometric {:n 3.0}) [0.0] [0.5]) ((vk/geometric {:n 3}) [0.0] [0.5])))
    (is (fn? (vk/geometric {:n 1})))
    (is (= :n (:option (ex-data (try (k/kernel :geometric) (catch clojure.lang.ExceptionInfo e e))))))
    (doseq [rad [0.0 -1.0 ##NaN]]
      (is (= :r (:option (ex-data (try (vk/geometric {:n 2 :r rad}) (catch clojure.lang.ExceptionInfo e e))))) (str rad)))))

(deftest log-and-t-student-test
  (testing "defaults: euclidean distance, p 2.0 (NullPointerException before)"
    (let [l (vk/log) g (vk/generalized-t-student)
          y [(m/sqrt 0.2)]]
      (is (close? (l [0.0] y) (- (Math/log 1.2)) 1e-14))
      (is (close? (g [0.0] y) (/ 1.0 1.2) 1e-14))
      (is (close? (l [0.0] y) -0.1823 1e-4))
      (is (close? (g [0.0] y) 0.8333 1e-4))
      (is (close? ((vk/log {}) [0.0] y) (l [0.0] y) 0.0))
      (is (close? ((vk/generalized-t-student nil) [0.0] y) (g [0.0] y) 0.0))))
  (testing "formulas for several p and distances"
    (doseq [p [0.5 1.0 2.0 3.5] d [0.0 0.01 1.0 12.0]]
      (is (close? ((vk/log {:p p}) [0.0] [d]) (- (Math/log1p (Math/pow d p))) 1e-15))
      (is (close? ((vk/generalized-t-student {:p p}) [0.0] [d]) (/ 1.0 (+ 1.0 (Math/pow d p))) 1e-15))
      (is (close? ((vk/power {:p p}) [0.0] [d]) (- (Math/pow d p)) 1e-15))))
  (testing "a custom distance and numbers"
    (is (close? ((vk/log {:p 1.0 :distance v/dist-abs}) [0.0 0.0] [1.0 1.0]) (- (Math/log 3.0)) 1e-15))
    (is (close? ((vk/generalized-t-student {:distance v/dist-abs}) [0.0 0.0] [1.0 1.0]) 0.2 1e-15))
    (is (= ((vk/log) 0.0 2.0) ((vk/log) [0.0] [2.0]))))
  (testing "identical points"
    (is (== 0.0 ((vk/log) [1.0 2.0] [1.0 2.0])))
    (is (== 1.0 ((vk/generalized-t-student) [1.0 2.0] [1.0 2.0])))))

(deftest chi-square-test
  (testing "empty bins (x + y = 0) contribute 0, not NaN (NaN before)"
    (is (== 1.0 ((vk/chi-square) [0 1] [0 1])))
    (is (== 1.0 ((vk/chi-square2) [0 1] [0 1])))
    (is (== 0.0 ((vk/chi-square) [0 0 0] [0 0 0])))
    (is (== 1.0 ((vk/chi-square2) [0 0 0] [0 0 0])))
    (is (== 0.0 ((vk/chi-square) 0 0)))
    (is (== 1.0 ((vk/chi-square2) 0 0)))
    (is (== 0.0 ((vk/chi-square) [1.0 0.0] [-1.0 0.0])) "x + y = 0 with x not 0"))
  (testing "values: sklearn additive_chi2_kernel mapping (reference-test) and hand values"
    (is (close? ((vk/chi-square) [1.0 3.0] [3.0 1.0]) 3.0 1e-15))
    (is (close? ((vk/chi-square2) [1.0 3.0] [3.0 1.0]) -3.0 1e-15))
    (is (close? ((vk/chi-square) 2.0 6.0) 3.0 1e-15)))
  (testing "empty vectors and NaN"
    (is (== 0.0 ((vk/chi-square) [] [])))
    (is (== 1.0 ((vk/chi-square2) [] [])))
    (is (m/nan? ((vk/chi-square) [##NaN 1.0] [1.0 1.0])))
    (is (m/nan? ((vk/chi-square2) [##NaN 1.0] [1.0 1.0])))))

(deftest b-spline-test
  (testing "exactly zero at and beyond the edge of the support, n + 1 (cancellation noise before)"
    (doseq [n (range 0 7) dist [(inc n) (+ n 1.5) (+ n 5.0) 30.0 100.0 1e6 1e300]]
      (is (== 0.0 ((vk/b-spline {:n n}) [0.0] [dist])) (str n " " dist))
      (is (== 0.0 ((vk/b-spline {:n n}) [dist] [0.0])) (str n " " dist))
      (is (== 0.0 ((vk/b-spline {:n n}) [0.0 0.0] [0.1 dist])) (str n " " dist))))
  (testing "central values: hat, cubic (2/3), degree 5 (11/20), degree 7 (151/315)"
    (is (close? ((vk/b-spline {:n 0}) [1.0] [1.0]) 1.0 1e-15))
    (is (close? ((vk/b-spline {:n 1}) [1.0] [1.0]) (/ 2.0 3.0) 1e-15))
    (is (close? ((vk/b-spline {:n 2}) [1.0] [1.0]) 0.55 1e-15))
    (is (close? ((vk/b-spline {:n 3}) [1.0] [1.0]) (/ 151.0 315.0) 1e-15))
    (is (close? ((vk/b-spline {:n 0}) [0.0] [0.25]) 0.75 1e-15))
    (is (close? ((vk/b-spline) [1.0] [1.0]) 0.55 1e-15) "default n is 2"))
  (testing "partition of unity of the cubic spline (integer shifts sum to 1)"
    (let [kf (vk/b-spline {:n 1})]
      (doseq [t [0.0 0.13 0.5 0.77]]
        (is (close? (reduce + (map #(kf [0.0] [(- t %)]) (range -4 5))) 1.0 1e-14)))))
  (testing "product over components, numbers, symmetry"
    (let [kf (vk/b-spline {:n 1})]
      (is (close? (kf [0.0 0.0] [0.3 0.9]) (* (kf [0.0] [0.3]) (kf [0.0] [0.9])) 1e-15))
      (is (= (kf 0.0 0.7) (kf [0.0] [0.7])))
      (is (= (kf [0.2 0.1] [0.7 0.5]) (kf [0.7 0.5] [0.2 0.1])))
      (is (== 1.0 (kf [] [])))))
  (testing "NaN and invalid degree"
    (is (m/nan? ((vk/b-spline {:n 1}) [##NaN] [0.0])))
    (is (= :n (:option (ex-data (try (vk/b-spline {:n -1}) (catch clojure.lang.ExceptionInfo e e))))))))

(deftest matern-test
  (testing "closed forms with theta = sqrt(2) ell: orders 1, 3 and 5 (measured error 3e-16)"
    (doseq [theta [0.5 1.0 3.3] d [0.0 0.01 0.3 1.0 2.7 12.0]
            :let [ell (/ theta (m/sqrt 2.0))
                  t (/ d ell)]]
      (is (close? ((vk/matern {:order 1 :theta theta}) [0.0] [d]) (Math/exp (- t)) 1e-14))
      (is (close? ((vk/matern {:order 3 :theta theta}) [0.0] [d]) (* (+ 1.0 (* (m/sqrt 3.0) t)) (Math/exp (- (* (m/sqrt 3.0) t)))) 1e-14))
      (is (close? ((vk/matern {:order 5 :theta theta}) [0.0] [d])
                  (* (+ 1.0 (* (m/sqrt 5.0) t) (/ (* 5.0 t t) 3.0)) (Math/exp (- (* (m/sqrt 5.0) t)))) 1e-14))))
  (testing "matern-12 is exp(-sqrt(2) d / theta), the exponential kernel"
    (is (close? ((vk/matern) [0.0] [1.0]) (Math/exp (- (m/sqrt 2.0))) 1e-15))
    (is (close? ((vk/matern) [0.0] [1.0]) 0.2431167344342142 1e-15)))
  (testing "only odd positive integer orders are accepted (even orders were silently wrong)"
    (doseq [order [1 3 5 7 9 3.0 11]]
      (is (fn? (vk/matern {:order order})) (str order)))
    (is (= ((vk/matern {:order 3}) [0.0] [1.3]) ((vk/matern {:order 3.0}) [0.0] [1.3])))
    (doseq [order [0 2 4 6 -1 -3 -2 0.0 2.5 3.5 1.5 "3" :a nil ##NaN ##Inf]
            :let [e (try (vk/matern {:order order}) nil (catch clojure.lang.ExceptionInfo e e))]]
      (is (some? e) (str order))
      (is (= :order (:option (ex-data e))) (str order))
      (is (= :invalid-input (:type (ex-data e))))))
  (testing "limits: zero and huge distances, NaN"
    (doseq [order [1 3 5 7 9 11]
            :let [kf (vk/matern {:order order})]]
      (is (== 1.0 (kf [0.0] [0.0])))
      (is (close? (kf [0.0] [1e-12]) 1.0 1e-10))
      (is (== 0.0 (kf [0.0] [1e150])) (str order))
      (is (== 0.0 (kf [0.0] [##Inf])) (str order))
      (is (m/nan? (kf [##NaN] [0.0])))))
  (testing "1-D numbers and vectors"
    (is (= ((vk/matern {:order 5}) 0.2 0.9) ((vk/matern {:order 5}) [0.2] [0.9])))))

(deftest units-and-mappings-test
  (testing "exponential(sigma) is laplacian(2 sigma^2); Genton's exp(-d / theta) is laplacian(theta)"
    (doseq [sigma [0.4 1.0 1.7] d [0.0 0.3 2.0]]
      (is (close? ((vk/exponential {:sigma sigma}) [0.0] [d]) ((vk/laplacian {:sigma (* 2.0 sigma sigma)}) [0.0] [d]) 1e-15))
      (is (close? ((vk/laplacian {:sigma sigma}) [0.0] [d]) (Math/exp (- (/ d sigma))) 1e-15))
      (is (close? ((vk/exponential {:sigma sigma}) [0.0] [d]) (Math/exp (- (/ d (* 2.0 sigma sigma)))) 1e-15))))
  (testing "gaussian is exp(-d^2 / (2 sigma^2))"
    (is (close? ((vk/gaussian {:sigma 2.0}) [0.0] [3.0]) (Math/exp (- (/ 9.0 8.0))) 1e-15)))
  (testing "cauchy(sigma) is rational-quadratic(c = sigma^2); with p 2 the generalized t-student"
    (doseq [sigma [0.5 1.0 3.0] d [0.0 0.7 5.0]]
      (is (close? ((vk/cauchy {:sigma sigma}) [0.0] [d]) ((vk/rational-quadratic {:c (* sigma sigma)}) [0.0] [d]) 1e-15)))
    (is (close? ((vk/cauchy) [0.0] [2.0]) ((vk/generalized-t-student) [0.0] [2.0]) 1e-15)))
  (testing "rational-quadratic c is a squared distance: c / (d^2 + c); scaling points by a and c by a^2 keeps the value"
    (is (close? ((vk/rational-quadratic {:c 2.0}) [0.0] [1.0]) (/ 2.0 3.0) 1e-15))
    (doseq [a [0.1 3.0 40.0]]
      (is (close? ((vk/rational-quadratic {:c (* a a 1.5)}) [0.0] [(* a 0.8)]) ((vk/rational-quadratic {:c 1.5}) [0.0] [0.8]) 1e-14))))
  (testing "rational-quadratic at huge distances keeps the relative accuracy (1 - d/(d + c) gave 0 and NaN)"
    (is (close? (/ ((vk/rational-quadratic) [0.0] [1e8]) 1e-16) 1.0 1e-12))
    (is (== 0.0 ((vk/rational-quadratic) [0.0] [##Inf])))
    (is (== 1.0 ((vk/rational-quadratic) [0.0] [0.0]))))
  (testing "multiquadratic and inverse multiquadratic: c is a distance, K(x, x) = abs(c)"
    (is (close? ((vk/multiquadratic {:c 2.0}) [1.0] [1.0]) 2.0 1e-15))
    (is (close? ((vk/inverse-multiquadratic {:c 2.0}) [1.0] [1.0]) 0.5 1e-15))
    (is (close? ((vk/multiquadratic {:c -2.0}) [1.0] [1.0]) 2.0 1e-15))
    (doseq [a [0.1 3.0 40.0]]
      (is (close? ((vk/multiquadratic {:c (* a 1.5)}) [0.0] [(* a 0.8)]) (* a ((vk/multiquadratic {:c 1.5}) [0.0] [0.8])) 1e-14))
      (is (close? ((vk/inverse-multiquadratic {:c (* a 1.5)}) [0.0] [(* a 0.8)]) (/ ((vk/inverse-multiquadratic {:c 1.5}) [0.0] [0.8]) a) 1e-14))))
  (testing "laplacian with the manhattan distance equals the scikit-learn laplacian_kernel (reference-test)"
    (is (close? ((vk/laplacian {:sigma 2.0 :distance v/dist-abs}) [0.0 0.0] [1.0 1.0]) (Math/exp -1.0) 1e-15))
    (is (close? ((vk/laplacian {:sigma 2.0}) [0.0 0.0] [1.0 1.0]) (Math/exp (- (/ (m/sqrt 2.0) 2.0))) 1e-15))))

(deftest bessel-test
  (testing "bessel is zero at identical points for v > -1 (crsouza form), one for the default v -1"
    (doseq [v [-0.5 0.0 0.5 1.0 2.0]]
      (is (== 0.0 ((vk/bessel {:v v}) [0.3 0.4] [0.3 0.4])) (str v)))
    (is (== 1.0 ((vk/bessel) [0.3 0.4] [0.3 0.4])))
    (is (== 1.0 ((vk/bessel {:v -1.0}) [0.3 0.4] [0.3 0.4]))))
  (testing "bessel2: K(x, x) = 1, default order 0 (kernlab 1), limit below the switch 1e-5"
    (is (== 1.0 ((vk/bessel2) [0.3 0.4] [0.3 0.4])))
    (is (== 1.0 ((vk/bessel2 {:order 1.0}) [0.3 0.4] [0.3 0.4])))
    (is (== 1.0 ((vk/bessel2) [0.0] [1e-6])))
    (is (close? ((vk/bessel2) [0.0] [1.1e-5]) 1.0 1e-9))
    (is (close? ((vk/bessel2 {:order 1.0}) [0.0] [1.1e-5]) 1.0 1e-9))
    (is (not (close? ((vk/bessel2) [0.0] [2.0]) ((vk/bessel2 {:order 1.0}) [0.0] [2.0]) 1e-3)))
    (is (close? ((vk/bessel2) [0.0] [2.0]) 0.22389077914123567 1e-14) "J0(2), mpmath")
    (is (close? ((vk/bessel2 {:order 1.0}) [0.0] [2.0]) (* 2.0 0.5767248077568734 0.5) 1e-14) "2 J1(2) / 2, mpmath")))

(deftest various-kernels-test
  (testing "linear and polynomial"
    (is (== 11.0 ((vk/linear) [1.0 2.0] [3.0 4.0])))
    (is (== 11.0 ((vk/linear {}) [1.0 2.0] [3.0 4.0])))
    (is (== 121.0 ((vk/polynomial) [1.0 2.0] [3.0 4.0])))
    (is (close? ((vk/polynomial {:alpha 0.5 :c 1.0 :p 3.0}) [1.0 2.0] [3.0 4.0]) (Math/pow 6.5 3.0) 1e-15))
    (is (== 0.0 ((vk/linear) [] []))))
  (testing "wave: 1 at zero, sin(d) / d, first zero at pi sigma"
    (is (== 1.0 ((vk/wave) [1.0] [1.0])))
    (is (close? ((vk/wave {:sigma 2.0}) [0.0] [3.0]) (* (/ 2.0 3.0) (Math/sin 1.5)) 1e-15))
    (is (close? ((vk/wave {:sigma 2.0}) [0.0] [(* m/PI 2.0)]) 0.0 1e-15)))
  (testing "periodic: 1 at multiples of the period, periodic in each argument"
    (doseq [p [0.5 1.0 2.5] kk [0 1 2 3]]
      (is (close? ((vk/periodic {:periodicity p}) [0.0] [(* kk p)]) 1.0 1e-14)))
    (is (close? ((vk/periodic {:periodicity 2.0 :sigma 0.7}) [0.3] [1.1]) ((vk/periodic {:periodicity 2.0 :sigma 0.7}) [0.3] [3.1]) 1e-14)))
  (testing "histogram, generalized histogram and hellinger: empty vectors, numbers, negative values"
    (is (== 0.0 ((vk/histogram) [] [])))
    (is (== 0.0 ((vk/generalized-histogram) [] [])))
    (is (== 0.0 ((vk/hellinger) [] [])))
    (is (== 0.5 ((vk/histogram) 0.5 2.0)))
    (is (== 2.0 ((vk/hellinger) 1.0 4.0)))
    (is (== 0.0 ((vk/hellinger) [1.0] [-1.0])) "a negative product contributes 0")
    (is (== 1.0 ((vk/generalized-histogram {:p 1.0}) [-1.0 3.0] [0.5 -0.5])) "abs values 0.5 + 0.5"))
  (testing "spline: empty product is 1, symmetric"
    (is (== 1.0 ((vk/spline) [] [])))
    (is (= ((vk/spline) [0.2 0.5] [0.7 0.1]) ((vk/spline) [0.7 0.1] [0.2 0.5]))))
  (testing "dirichlet: n + 1/2 per component at identical points"
    (is (== 1.5 ((vk/dirichlet) [0.3] [0.3])))
    (is (== 2.25 ((vk/dirichlet) [0.3 0.1] [0.3 0.1])))
    (is (close? ((vk/dirichlet {:n 2.0}) 0.0 0.6) (/ (Math/sin 1.5) (* 2.0 (Math/sin 0.3))) 1e-15)))
  (testing "hyperbolic tangent, secant"
    (is (close? ((vk/hyperbolic-tangent {:alpha 0.5 :c 0.1}) [1.0 2.0] [3.0 4.0]) (Math/tanh 5.6) 1e-15))
    (is (close? ((vk/hyperbolic-secant {:a 2.0}) [0.0] [0.5]) (/ 1.0 (Math/cosh 1.0)) 1e-15))))

(deftest anova-reading-test
  (testing "K(x, x) is the dimension"
    (is (close? ((vk/anova) [0.3 0.4 0.5] [0.3 0.4 0.5]) 3.0 1e-15)))
  (testing "d acts on every term, not on the sum"
    (let [x [0.0 0.0] y [1.0 2.0]
          per-term (+ (Math/exp -2.0) (Math/exp -8.0))
          on-sum (Math/pow (+ (Math/exp -1.0) (Math/exp -4.0)) 2.0)]
      (is (close? ((vk/anova {:d 2.0}) x y) per-term 1e-15))
      (is (> (m/abs (- per-term on-sum)) 1e-2))))
  (testing "k is a power of the components before the difference"
    (is (close? ((vk/anova {:k 2.0}) [0.0] [2.0]) (Math/exp -16.0) 1e-15))
    (is (close? ((vk/anova {:sigma 0.5}) [0.0] [2.0]) (Math/exp -2.0) 1e-15))))

;; Gram matrices: the claims of the docstrings about positive definiteness.

(defn- gram
  "Symmetrised Gram matrix as `double[][]`."
  ^"[[D" [kf pts]
  (let [n (count pts)
        g (into-array (map double-array (for [x pts] (for [y pts] (double (kf x y))))))]
    (dotimes [i n]
      (dotimes [j i]
        (let [a (* 0.5 (+ (aget ^doubles (aget g i) j) (aget ^doubles (aget g j) i)))]
          (aset ^doubles (aget g i) j a)
          (aset ^doubles (aget g j) i a))))
    g))

(defn- min-eigenvalue
  ^double [^"[[D" g]
  (double (reduce min (.getRealEigenvalues (EigenDecomposition. (Array2DRowRealMatrix. g))))))

(defn- center
  "Double centering `P G P`, the Gram matrix on vectors with a zero sum."
  ^"[[D" [^"[[D" g]
  (let [n (alength g)
        rows (vec (for [i (range n)] (/ (reduce + (aget g i)) n)))
        total (/ (reduce + rows) n)]
    (into-array (for [i (range n)]
                  (double-array (for [j (range n)]
                                  (+ (- (aget ^doubles (aget g i) j) (rows i) (rows j)) total)))))))

(defn- points
  [seed n d lo hi]
  (let [rng (r/rng :isaac seed)]
    (vec (repeatedly n (fn [] (vec (repeatedly d #(r/drandom rng lo hi))))))))

(deftest positive-definiteness-test
  ;; 30 seeded random points; measured minimum eigenvalues of the PSD kernels: at least -7.6e-14
  (testing "positive semi-definite kernels"
    (doseq [[kw params dims] [[:linear {} [1 3]]
                              [:polynomial {:p 2.0 :c 1.0} [1 3]]
                              [:polynomial {:p 3.0 :c 0.5 :alpha 0.3} [3]]
                              [:gaussian {:sigma 0.8} [1 3 5]]
                              [:exponential {} [1 3 5]]
                              [:laplacian {} [1 3 5]]
                              [:rational-quadratic {} [1 3 5]]
                              [:inverse-multiquadratic {} [1 3 5]]
                              [:cauchy {} [1 3 5]]
                              [:matern-12 {} [1 3 5]]
                              [:matern-32 {} [1 3 5]]
                              [:matern-52 {} [1 3 5]]
                              [:matern {:order 7 :theta 2.0} [1 3]]
                              [:pearson {:omega 2.0} [1 3 5]]
                              [:generalized-t-student {} [1 3]]
                              [:hyperbolic-secant {} [1 3]]
                              [:wave {} [1 2 3]]
                              [:periodic {} [1]]
                              [:b-spline {:n 1} [1 3]]
                              [:b-spline {:n 2} [1 3]]
                              [:bessel2 {} [1 2]]
                              [:bessel2 {:order 1.0} [1 2 3]]
                              [:dirichlet {:n 2.0} [1 3]]
                              [:anova {} [1 3]]
                              [:triangular {} [1]]
                              [:circular {} [1 2]]
                              [:spherical {} [1 2 3]]
                              [:geometric {:n 4} [1 2 3 4]]
                              [:geometric {:n 5 :r 1.5} [5]]]
            d dims
            :let [pts (points 20260514 30 d -1.5 1.5)
                  eig (min-eigenvalue (gram (k/kernel kw params) pts))]]
      (is (>= eig -1e-10) (str kw " " params " dim " d ": " eig))))
  (testing "positive semi-definite kernels on non-negative input"
    (doseq [kw [:spline :histogram :chi-square :hellinger :generalized-histogram]
            d [1 3 5]
            :let [pts (points 20260514 30 d 0.0 2.0)
                  eig (min-eigenvalue (gram (k/kernel kw) pts))]]
      (is (>= eig -1e-10) (str kw " dim " d ": " eig))))
  (testing "not positive semi-definite (the documented counterexamples, margin 1e-3)"
    (doseq [[kw params d] [[:hyperbolic-tangent {} 3]
                           [:wave {} 5]
                           [:periodic {} 3]
                           [:bessel {} 3]
                           [:bessel {} 5]
                           [:bessel2 {} 3]
                           [:bessel2 {} 5]
                           [:bessel2 {:order 1.0} 5]
                           [:chi-square2 {} 3]
                           [:power {} 3]
                           [:log {} 3]
                           [:multiquadratic {} 3]]
            :let [lo (if (= kw :chi-square2) 0.0 -1.5)
                  hi (if (= kw :chi-square2) 2.0 1.5)
                  pts (points 20260514 30 d lo hi)
                  eig (min-eigenvalue (gram (k/kernel kw params) pts))]]
      (is (< eig -1e-3) (str kw " dim " d ": " eig))))
  (testing "conditionally positive definite: power, log, chi-square2 (and minus multiquadratic) on zero-sum vectors"
    (doseq [[kw params d sign] [[:power {} 3 1.0] [:power {:p 1.0} 3 1.0] [:log {} 3 1.0] [:log {:p 1.5} 5 1.0]
                                [:chi-square2 {} 3 1.0] [:multiquadratic {} 3 -1.0]]
            :let [pts (points 20260514 30 d (if (= kw :chi-square2) 0.0 -1.5) (if (= kw :chi-square2) 2.0 1.5))
                  kf (k/kernel kw params)
                  eig (min-eigenvalue (center (gram (fn [x y] (* sign (double (kf x y)))) pts)))]]
      ;; magnitudes are up to 1e2: the rounding level is 1e-13
      (is (>= eig -1e-9) (str kw " dim " d ": " eig)))))

(deftest wrong-arguments-test
  (testing "points which are not numbers or vectors throw"
    (doseq [kf [(vk/gaussian) (vk/laplacian) (vk/pearson) (vk/geometric {:n 3}) (vk/matern) (vk/linear)]]
      (is (thrown? Exception (kf "a" 1.0)))
      (is (thrown? Exception (kf :a 1.0)))))
  (testing "NaN in an argument gives NaN for the distance kernels"
    (doseq [kf [(vk/gaussian) (vk/laplacian) (vk/pearson) (vk/cauchy) (vk/wave) (vk/hyperbolic-secant)
                (vk/exponential) (vk/rational-quadratic) (vk/periodic)]]
      (is (m/nan? (kf [##NaN 1.0] [0.0 1.0]))))))
