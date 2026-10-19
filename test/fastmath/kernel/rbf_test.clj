(ns fastmath.kernel.rbf-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.kernel.rbf :as rbf]
            [fastmath.random :as r])
  (:import [org.apache.commons.math3.linear Array2DRowRealMatrix QRDecomposition RealMatrix]))

;; Reference data: `utils/fastmath/dev/generate_rbf_kernel_reference.py` (sympy exact rationals for wu, wendland and
;; gneiting; mpmath for shifted surface splines, poisson, matern, whittaker and gaussians-laguerre; the R package
;; fields for wendland through `generate_rbf_kernel_reference.R`), file `rbf_kernel_reference.edn`.
;; The CPD sign table is checked on seeded random points (see `cpd-sign`).

(def ^:private reference (delay (edn/read-string (slurp (io/resource "kernel/rbf_kernel_reference.edn")))))

(defn- max-err
  "Maximum absolute error of `f` on a reference case (parallel vectors `:r` and `:v`); NaN if `f` gives NaN."
  ^double [f c]
  (reduce (fn [^double acc [rr v]]
            (let [e (Math/abs (- (double (f rr)) (double v)))]
              (if (or (Double/isNaN e) (> e acc)) e acc)))
          0.0 (map vector (:r c) (:v c))))

(defn- max-rel-err
  "Maximum error scaled by `max(1, abs(v))`."
  ^double [f c]
  (reduce (fn [^double acc [rr v]]
            (let [e (/ (Math/abs (- (double (f rr)) (double v))) (Math/max 1.0 (Math/abs (double v))))]
              (if (or (Double/isNaN e) (> e acc)) e acc)))
          0.0 (map vector (:r c) (:v c))))

(defn- max-true-rel-err
  "Maximum relative error (the reference values are positive)."
  ^double [f c]
  (reduce (fn [^double acc [rr v]]
            (let [e (/ (Math/abs (- (double (f rr)) (double v))) (Math/abs (double v)))]
              (if (or (Double/isNaN e) (> e acc)) e acc)))
          0.0 (map vector (:r c) (:v c))))

(defn- thrown-data
  "Data of the `ex-info` thrown by `f` (nil when nothing or another exception is thrown)."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e (assoc (ex-data e) ::message (.getMessage e)))
       (catch Throwable _ nil)))

(defn- invalid-input?
  "True when `f` throws `ex-info` with `:type :invalid-input` and the given data."
  [f expected]
  (let [d (thrown-data f)]
    (and (= :invalid-input (:type d))
         (= expected (select-keys d (keys expected))))))

(defn- bits=
  "Bitwise equality of two doubles (NaN equals NaN)."
  [a b]
  (= (Double/doubleToLongBits (double a)) (Double/doubleToLongBits (double b))))

(def ^:private rs [0.0 1.0e-300 1.0e-9 1.0e-3 0.1 0.25 0.5 0.9 0.999999 1.0 1.0001 1.5 3.0 10.0 1.0e5 1.0e100])

(def ^:private full-params
  "A valid parameter map for every family (also the registry of the families)."
  {:linear [rbf/linear {}]
   :gaussian [rbf/gaussian {}]
   :truncated-power [rbf/truncated-power {:k 2.5}]
   :gaussians-laguerre [rbf/gaussians-laguerre {:dimension 3.0 :degree 2}]
   :poisson [rbf/poisson {:d 5.0}]
   :matern [rbf/matern {:order 5}]
   :generalized-multiquadratic [rbf/generalized-multiquadratic {:beta 2.5}]
   :radial-powers [rbf/radial-powers {:beta 2.5}]
   :thin-plate-splines [rbf/thin-plate-splines {:beta 2.0}]
   :whittaker [rbf/whittaker {:alpha 0 :k 2 :beta 1.0}]
   :shifted-surface-splines [rbf/shifted-surface-splines {:s 3 :beta 2.5}]
   :wendland [rbf/wendland {:s 3 :k 2}]
   :gneiting [rbf/gneiting {:s 3 :l 5.0}]
   :wu [rbf/wu {:l 3.0 :k 2.0}]})

;; ------------------------------------------------------------------------------------------------ wu

(deftest wu-test
  (testing "every closed form (l 0..3, k 0..l) equals the exact rational value of sympy, long and double arguments"
    (let [cases (filter #(<= (:l %) 3) (filter #(integer? (:l %)) (:wu @reference)))]
      (is (= 10 (count cases)))
      (doseq [c cases
              :let [p (select-keys c [:l :k])
                    pd {:l (double (:l c)) :k (double (:k c))}
                    f (rbf/wu p)
                    fd (rbf/wu pd)]]
        ;; measured 3.3e-16 over all 10 closed forms and 12 distances
        (is (< (max-err f c) 2.0e-14) (str p))
        (is (< (max-err fd c) 2.0e-14) (str pd))
        (is (< (max-err (rbf/wu {:l (double (:l c)) :k (:k c)}) c) 2.0e-14) "double l, long k")
        (doseq [x rs]
          (is (bits= (f x) (fd x)) (str p " long and double agree at " x))))))
  (testing "wu-33 at r = 1e-3: 0.9978 (the numeric path gave -20958)"
    (doseq [p [{:l 3.0 :k 3.0} {:l 3 :k 3} {:l 3.0 :k 3}]]
      (is (m/delta-eq 0.9978125021874987 ((rbf/wu p) 1.0e-3) 1.0e-13))
      (is (m/delta-eq 0.9978 ((rbf/wu p) 1.0e-3) 1.0e-4))))
  (testing "numeric path for the non-closed-form cases"
    ;; measured 8.1e-10 (l 4, k 1) and 1.6e-10 (l 1.5, k 0); k 0 is accurate, k >= 2 is unstable and is not tested
    (doseq [c (:wu @reference)
            :when (or (= 4 (:l c)) (= 1.5 (:l c)))
            :when (<= (:k c) 1)]
      (is (< (max-err (rbf/wu (select-keys c [:l :k])) c) 1.0e-8) (str (select-keys c [:l :k])))))
  (testing "boundaries"
    (doseq [p [{:l 1 :k 0} {:l 2.0 :k 1.0} {:l 3 :k 3} {:l 3.0 :k 1}]
            :let [f (rbf/wu p)]]
      (is (= 1.0 (f 0.0)) (str p))
      (is (= 0.0 (f 1.0)) "exactly at the support radius")
      (is (= 0.0 (f 1.0001)))
      (is (= 0.0 (f 1.0e100)))
      (is (= 0.0 (f 1.0e300)) "no 0 * Inf = NaN beyond the support")
      (is (= 0.0 (f ##Inf)))
      (is (bits= (f 0.37) (f -0.37)) "symmetric")
      (is (< 0.0 (f 0.999999) 1.0e-5))
      (is (m/delta-eq 1.0 (f 1.0e-300) 1.0e-15))
      (is (m/delta-eq 1.0 (f 1.0e-9) 1.0e-8))))
  (testing "scale is the support radius"
    (doseq [sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/wu {:l 2 :k 1})
                  fs (rbf/wu {:l 2 :k 1 :scale sc})]
            x [0.0 0.3 0.999 1.0 1.2]]
      (is (bits= (f1 (/ x sc)) (fs x)) (str sc " " x)))
    (is (= 0.0 ((rbf/wu {:l 2 :k 1 :scale 1.0e-300}) 1.0e-300))))
  (testing "invalid input"
    (is (invalid-input? #(rbf/wu) {:family :wu :property :l}))
    (is (invalid-input? #(rbf/wu {:l 2.0}) {:family :wu :property :k}))
    (is (invalid-input? #(rbf/wu {:k 2.0}) {:family :wu :property :l}))
    (is (invalid-input? #(rbf/wu {:l 2.0 :k 1.5}) {:family :wu :property :k}))
    (is (invalid-input? #(rbf/wu {:l 2.0 :k -1}) {:family :wu :property :k}))
    (is (invalid-input? #(rbf/wu {:l -1.0 :k 0}) {:family :wu :property :l}))))

;; ------------------------------------------------------------------------------------------------ whittaker

(deftest whittaker-test
  (testing "closed forms (alpha 0 and 1, k 2 and 3, any numeric type) vs mpmath quadrature, r up to 1e5"
    (let [cases (filter #(integer? (:alpha %)) (:whittaker @reference))]
      (is (= 12 (count cases)))
      (doseq [c cases
              :let [p (select-keys c [:alpha :k :beta])
                    fd (rbf/whittaker {:alpha (double (:alpha c)) :k (double (:k c)) :beta (double (:beta c))})
                    f (rbf/whittaker p)]]
        ;; measured 4.3e-16 relative (the series for r > 0.67 beta removes the cancellation)
        (is (< (max-true-rel-err f c) 1.0e-14) (str p))
        (is (< (max-true-rel-err fd c) 1.0e-14) (str p " doubles"))
        (doseq [x rs]
          (is (bits= (f x) (fd x)) (str p " long and double agree at " x))))))
  (testing "numerical integration path (alpha 0.5)"
    ;; measured 3.8e-5 relative to max(1, v) for beta 2.5 (quadrature tolerance 1e-6 relative to its own scale), 3.4e-6 for beta 0.5
    (doseq [c (remove #(integer? (:alpha %)) (:whittaker @reference))]
      (is (< (max-rel-err (rbf/whittaker (select-keys c [:alpha :k :beta])) c) 1.0e-4))))
  (testing "boundaries"
    (doseq [p [{:alpha 0 :k 2 :beta 1.0} {:alpha 1 :k 2 :beta 0.5} {:alpha 0 :k 3 :beta 2.0} {:alpha 1 :k 3 :beta 1.0}]
            :let [f (rbf/whittaker p)]]
      (is (= 1.0 (f 0.0)) (str p))
      (is (m/delta-eq 1.0 (f 1.0e-300) 1.0e-15))
      (is (bits= (f 0.37) (f -0.37)))
      (is (<= 0.0 (f 1.0e150) 1.0e-100) "no Inf - Inf = NaN for a huge r")
      (is (<= 0.0 (f 1.0e300) 1.0e-250))
      (is (<= 0.0 (f ##Inf) 1.0e-300))
      (is (apply > (map f [0.0 0.5 1.0 5.0 50.0])) "decreasing")))
  (testing "scale"
    (doseq [sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/whittaker {:alpha 1 :k 3 :beta 1.0})
                  fs (rbf/whittaker {:alpha 1 :k 3 :beta 1.0 :scale sc})]
            x [0.0 0.3 1.0 7.0]]
      (is (bits= (f1 (/ x sc)) (fs x)) (str sc " " x))))
  (testing "invalid input"
    (is (invalid-input? #(rbf/whittaker) {:family :whittaker :property :alpha}))
    (is (invalid-input? #(rbf/whittaker {:alpha 0 :beta 1.0}) {:family :whittaker :property :k}))
    (is (invalid-input? #(rbf/whittaker {:alpha 0 :k 2}) {:family :whittaker :property :beta}))
    (is (invalid-input? #(rbf/whittaker {:alpha 0 :k 2 :beta 0.0}) {:family :whittaker :property :beta}))
    (is (invalid-input? #(rbf/whittaker {:alpha -1.0 :k 2 :beta 1.0}) {:family :whittaker :property :alpha}))
    (is (invalid-input? #(rbf/whittaker {:alpha 0 :k 0 :beta 1.0}) {:family :whittaker :property :k}))))

;; ------------------------------------------------------------------------------------------------ shifted surface splines

(defn- sss-closed
  "Luh equations (1) and (2): the oracle in closed form for an integer s and beta."
  ^double [^long s ^double beta ^double x]
  (let [t (- beta (* 0.5 s))
        u (+ 1.0 (* x x))]
    (if (odd? s)
      (* (if (even? (long (Math/ceil t))) 1.0 -1.0) (Math/pow u t))
      (* (if (even? (long (+ 1 t))) 1.0 -1.0) (Math/pow u t) (Math/log (Math/sqrt u))))))

(deftest shifted-surface-splines-test
  (testing "vs mpmath for s 1..5 and several beta"
    (let [cases (:shifted-surface-splines @reference)]
      (is (= 18 (count cases)))
      (doseq [c cases
              :let [p (select-keys c [:s :beta])
                    f (rbf/shifted-surface-splines p)]]
        ;; measured 2.2e-16 relative to max(1, v)
        (is (< (max-rel-err f c) 1.0e-14) (str p))
        (is (< (max-rel-err (rbf/shifted-surface-splines {:s (double (:s c)) :beta (:beta c)}) c) 1.0e-14) "s as a double"))))
  (testing "closed form oracle"
    (doseq [[s betas] {1 [1.0 1.5 2.0 4.0] 2 [2.0 3.0 5.0] 3 [2.0 2.5 3.5] 4 [3.0 4.0] 5 [3.0 4.5 5.5]}
            beta betas
            x [0.0 0.01 0.3 1.0 2.5 7.0]
            :let [f (rbf/shifted-surface-splines {:s s :beta beta})
                  expected (sss-closed s beta x)]]
      (is (< (Math/abs (- expected (f x))) (* 1.0e-14 (Math/max 1.0 (Math/abs expected)))) (str s " " beta " " x))))
  (testing "s 3, beta 2.5, r 1 is the sign times 2 (the old code gave 2.828)"
    (is (= -2.0 ((rbf/shifted-surface-splines {:s 3 :beta 2.5}) 1.0)))
    (is (= -2.0 ((rbf/shifted-surface-splines {:s 3 :beta 2.5}) -1.0))))
  (testing "s 2, beta 2, r 1: log sqrt(2) * 2 (the old code gave 1.665)"
    (is (m/delta-eq 0.6931471805599453 ((rbf/shifted-surface-splines {:s 2 :beta 2}) 1.0) 1.0e-15)))
  (testing "boundaries"
    (is (= -1.0 ((rbf/shifted-surface-splines {:s 3 :beta 2.5}) 0.0)) "odd s: the sign")
    (is (= 1.0 ((rbf/shifted-surface-splines {:s 1 :beta 2.0}) 0.0)))
    (is (= 0.0 ((rbf/shifted-surface-splines {:s 2 :beta 2}) 0.0)) "even s: 0 at 0")
    ;; the logarithm is of 1 + r^2: smooth, 1/2 r^2 near 0 (not abs(r) as the old code)
    (let [f (rbf/shifted-surface-splines {:s 2 :beta 2})]
      (is (m/delta-eq 1.0 (/ (f 1.0e-5) (* 0.5 1.0e-10)) 1.0e-9))
      (is (m/delta-eq 1.0 (/ (f 1.0e-9) (* 0.5 1.0e-18)) 1.0e-9))
      (is (bits= (f 0.2) (f -0.2))))
    (is (= 0.0 ((rbf/shifted-surface-splines {:s 2 :beta 2}) 1.0e-300)))
    (is (Double/isInfinite ((rbf/shifted-surface-splines {:s 3 :beta 2.5}) 1.0e200)) "grows without bound")
    (is (m/delta-eq 1.0 (/ (- ((rbf/shifted-surface-splines {:s 3 :beta 2.5}) 1.0e100)) 1.0e200) 1.0e-12) "sign * r^2t"))
  (testing "scale"
    (doseq [s [3 4]
            sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/shifted-surface-splines {:s s :beta 3.0})
                  fs (rbf/shifted-surface-splines {:s s :beta 3.0 :scale sc})]
            x [0.0 0.3 1.0 7.0]]
      (is (bits= (f1 (/ x sc)) (fs x)) (str s " " sc " " x))))
  (testing "invalid input"
    (is (invalid-input? #(rbf/shifted-surface-splines) {:family :shifted-surface-splines :property :s}))
    (is (invalid-input? #(rbf/shifted-surface-splines {:s 3}) {:family :shifted-surface-splines :property :beta}))
    (is (invalid-input? #(rbf/shifted-surface-splines {:beta 3.0}) {:family :shifted-surface-splines :property :s}))
    (is (invalid-input? #(rbf/shifted-surface-splines {:s 0 :beta 3.0}) {:property :s}))
    (is (invalid-input? #(rbf/shifted-surface-splines {:s 2.5 :beta 3.0}) {:property :s}))
    (is (invalid-input? #(rbf/shifted-surface-splines {:s 3 :beta 1.5}) {:property :beta}) "beta must be above s / 2")
    (is (invalid-input? #(rbf/shifted-surface-splines {:s 4 :beta 2}) {:property :beta}))
    (is (invalid-input? #(rbf/shifted-surface-splines {:s 4 :beta 3.5}) {:property :beta}) "even s needs an integer beta")
    (is (some? (rbf/shifted-surface-splines {:s 3 :beta 1.5000001})))))

;; ------------------------------------------------------------------------------------------------ wendland

(deftest wendland-test
  (testing "vs sympy I^k for s 1..8, k 0..4"
    (is (= 40 (count (:wendland @reference))))
    (doseq [c (:wendland @reference)
            :let [p (select-keys c [:s :k])]]
      ;; measured 5.6e-16
      (is (< (max-err (rbf/wendland p) c) 1.0e-14) (str p))
      (is (< (max-err (rbf/wendland {:s (double (:s c)) :k (double (:k c))}) c) 1.0e-14) "doubles")))
  (testing "vs R fields::Wendland(dimension s, k) for k 1..4"
    (is (= 32 (count (:fields-wendland @reference))))
    (doseq [c (:fields-wendland @reference)]
      ;; measured 2.2e-16
      (is (< (max-err (rbf/wendland (select-keys c [:s :k])) c) 1.0e-13) (str (select-keys c [:s :k])))))
  (testing "phi(3, 1)(0.5) = 0.1875 (the old exponent gave 0.1436)"
    (is (m/delta-eq 0.1875 ((rbf/wendland {:s 3 :k 1}) 0.5) 1.0e-15))
    (is (m/delta-eq 0.1875 ((rbf/wendland {:s 3.0 :k 1.0}) -0.5) 1.0e-15)))
  (testing "only floor(s / 2) matters"
    (doseq [x [0.0 0.2 0.7 0.99]
            k (range 5)]
      (is (bits= ((rbf/wendland {:s 2 :k k}) x) ((rbf/wendland {:s 3 :k k}) x)))
      (is (bits= ((rbf/wendland {:s 2.0 :k k}) x) ((rbf/wendland {:s 2.9 :k k}) x)))
      (is (not (bits= ((rbf/wendland {:s 1 :k k}) 0.5) ((rbf/wendland {:s 2 :k k}) 0.5))))))
  (testing "boundaries"
    (doseq [s [1 2 3 7] k (range 5)
            :let [f (rbf/wendland {:s s :k k})]]
      (is (= 1.0 (f 0.0)) (str s " " k))
      (is (m/delta-eq 1.0 (f 1.0e-300) 1.0e-15))
      (is (= 0.0 (f 1.0)) "exactly at the support radius")
      (is (= 0.0 (f 1.0000001)))
      (is (= 0.0 (f 1.0e150)))
      (is (= 0.0 (f ##Inf)))
      (is (< 0.0 (f 0.9999999) 1.0e-6))
      (is (bits= (f 0.37) (f -0.37)))))
  (testing "scale is the support radius"
    (doseq [sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/wendland {:s 3 :k 3})
                  fs (rbf/wendland {:s 3 :k 3 :scale sc})]
            x [0.0 0.3 0.999 1.0 1.2]]
      (is (bits= (f1 (/ x sc)) (fs x)) (str sc " " x))))
  (testing "invalid input"
    (is (invalid-input? #(rbf/wendland) {:family :wendland :property :s}))
    (is (invalid-input? #(rbf/wendland {:k 1}) {:family :wendland :property :s}))
    (is (invalid-input? #(rbf/wendland {:s 3}) {:family :wendland :property :k}))
    (is (invalid-input? #(rbf/wendland {:s 3 :k 5}) {:family :wendland :property :k}) "k >= 5")
    (is (invalid-input? #(rbf/wendland {:s 3 :k 6}) {:property :k}))
    (is (invalid-input? #(rbf/wendland {:s 3 :k -1}) {:property :k}))
    (is (invalid-input? #(rbf/wendland {:s 3 :k 1.5}) {:property :k}))
    (is (invalid-input? #(rbf/wendland {:s 0 :k 1}) {:property :s}))
    (is (some? (rbf/wendland {:s 3 :k 4})))))

;; ------------------------------------------------------------------------------------------------ gneiting

(deftest gneiting-test
  (testing "vs the turning bands operator applied to the Wendland function (sympy), s 1..4, l 3.5 .. 7"
    (is (= 16 (count (:gneiting @reference))))
    (doseq [c (:gneiting @reference)]
      ;; measured 2.2e-16
      (is (< (max-err (rbf/gneiting (select-keys c [:s :l])) c) 1.0e-14) (str (select-keys c [:s :l])))))
  (testing "the return type is hinted (a primitive function)"
    (is (instance? clojure.lang.IFn$DD (rbf/gneiting {:s 3 :l 5.0}))))
  (testing "boundaries"
    (doseq [p [{:s 1 :l 4.0} {:s 3.0 :l 5.5} {:s 2 :l 7}]
            :let [f (rbf/gneiting p)]]
      (is (= 1.0 (f 0.0)) (str p))
      (is (= 0.0 (f 1.0)))
      (is (= 0.0 (f 1.5)))
      (is (= 0.0 (f 1.0e300)))
      (is (bits= (f 0.37) (f -0.37)))
      (is (some neg? (map f (range 0.0 1.0 0.01))) "oscillates: negative values")
      (is (m/delta-eq 1.0 (f 1.0e-300) 1.0e-15))))
  (testing "scale"
    (doseq [sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/gneiting {:s 3 :l 5.0})
                  fs (rbf/gneiting {:s 3 :l 5.0 :scale sc})]
            x [0.0 0.3 0.999 1.0 1.2]]
      (is (bits= (f1 (/ x sc)) (fs x)))))
  (testing "invalid input: s = 0 is guarded (the formula divides by s)"
    (is (invalid-input? #(rbf/gneiting {:s 0 :l 5.0}) {:family :gneiting :property :s}))
    (is (invalid-input? #(rbf/gneiting {:s -1.0 :l 5.0}) {:property :s}))
    (is (invalid-input? #(rbf/gneiting) {:family :gneiting :property :s}))
    (is (invalid-input? #(rbf/gneiting {:s 3}) {:family :gneiting :property :l}))))

;; ------------------------------------------------------------------------------------------------ poisson

(deftest poisson-test
  (testing "vs mpmath, normalised Gamma(p + 1) 2^p J_p(r) / r^p"
    (is (= 8 (count (:poisson @reference))))
    (doseq [c (:poisson @reference)]
      ;; measured 4.0e-15
      (is (< (max-err (rbf/poisson (select-keys c [:d])) c) 1.0e-14) (str (select-keys c [:d])))))
  (testing "normalised to 1 at r = 0 for every d (the old code gave 0.798, 0.5, 0.266, 0.125 for d 3..6)"
    (doseq [d [1.0 1.5 2.0 2.5 3.0 4.0 5.0 6.0 7.0 12.0 40.0 300.0]]
      (is (= 1.0 ((rbf/poisson {:d d}) 0.0)) (str d))
      (is (m/delta-eq 1.0 ((rbf/poisson {:d d}) 1.0e-9) 1.0e-12) (str d))
      (is (m/delta-eq 1.0 ((rbf/poisson {:d d}) 1.0e-300) 1.0e-12) (str d))))
  (testing "the series near 0 is continuous with the Bessel function"
    (doseq [d [2.5 5.0 6.0 7.0]
            :let [f (rbf/poisson {:d d})]]
      (is (m/delta-eq (f 0.99999e-5) (f 1.00001e-5) 1.0e-14))))
  (testing "closed forms"
    (is (m/delta-eq (Math/cos 2.0) ((rbf/poisson {:d 1}) 2.0) 1.0e-15))
    (is (m/delta-eq (/ (Math/sin 2.0) 2.0) ((rbf/poisson {:d 3}) 2.0) 1.0e-15))
    (is (m/delta-eq (/ (Math/sin 2.0) 2.0) ((rbf/poisson {:d 3.0}) -2.0) 1.0e-15)))
  (testing "boundaries"
    (doseq [d [2.0 3.0 4.0 5.0 6.5]
            :let [f (rbf/poisson {:d d})]]
      (is (bits= (f 0.37) (f -0.37)))
      (is (<= (Math/abs (double (f 1.0e5))) 1.0e-2))
      (is (not (Double/isNaN (f 1.0e100))))
      (is (<= (Math/abs (double (f 1.0e10))) 1.0e-2))))
  (testing "scale"
    (doseq [sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/poisson {:d 5.0})
                  fs (rbf/poisson {:d 5.0 :scale sc})]
            x [0.0 0.3 1.0 7.0]]
      (is (bits= (f1 (/ x sc)) (fs x)))))
  (testing "invalid input"
    (is (invalid-input? #(rbf/poisson) {:family :poisson :property :d}))
    (is (invalid-input? #(rbf/poisson {:d 0.5}) {:family :poisson :property :d}))))

;; ------------------------------------------------------------------------------------------------ matern

(deftest matern-test
  (testing "vs mpmath for orders 1, 3, 5, 7, 9 up to r = 800"
    (is (= 5 (count (:matern @reference))))
    (doseq [c (:matern @reference)]
      ;; measured 5.6e-16 absolute
      (is (< (max-err (rbf/matern (select-keys c [:order])) c) 1.0e-14) (str (select-keys c [:order])))
      (is (< (max-err (rbf/matern {:order (double (:order c))}) c) 1.0e-14) "order as a double")))
  (testing "closed forms of the half-odd orders"
    (doseq [x [1.0e-3 0.5 1.0 4.0 20.0]]
      (is (m/delta-eq (Math/exp (- x)) ((rbf/matern {:order 1}) x) 1.0e-15 1.0e-14))
      (is (m/delta-eq (* (+ 1.0 x) (Math/exp (- x))) ((rbf/matern {:order 3}) x) 1.0e-15 1.0e-14))
      (is (m/delta-eq (* (+ 1.0 x (/ (* x x) 3.0)) (Math/exp (- x))) ((rbf/matern {:order 5}) x) 1.0e-15 1.0e-14))))
  (testing "boundaries"
    (doseq [order [1 3 5 9]
            :let [f (rbf/matern {:order order})]]
      (is (= 1.0 (f 0.0)))
      (is (m/delta-eq 1.0 (f 1.0e-300) 1.0e-15))
      (is (m/delta-eq 1.0 (f 1.0e-9) 1.0e-8))
      (is (bits= (f 0.37) (f -0.37)))
      (is (= 0.0 (f 1.0e200)) "no 0 * Inf = NaN for a huge r (the old code gave NaN for order 5)")
      (is (= 0.0 (f 1.0e300)))
      (is (= 0.0 (f ##Inf)))
      (is (apply > (map f [0.0 0.1 1.0 5.0 50.0])) "decreasing")))
  (testing "scale is the length scale (no sqrt(order) factor)"
    (doseq [sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/matern {:order 3})
                  fs (rbf/matern {:order 3 :scale sc})]
            x [0.0 0.3 1.0 7.0]]
      (is (bits= (f1 (/ x sc)) (fs x)))))
  (testing "even, non-positive and non-integer order throw (even order gave garbage)"
    (doseq [o [0 2 4 6 -1 -3 2.5]]
      (is (invalid-input? #(rbf/matern {:order o}) {:family :matern :property :order}) (str o)))
    (is (invalid-input? #(rbf/matern) {:family :matern :property :order}))
    (is (some? (rbf/matern {:order 1})))))

;; ------------------------------------------------------------------------------------------------ gaussians-laguerre

(deftest gaussians-laguerre-test
  (testing "vs mpmath exp(-r^2) L_n^(d/2)(r^2), degree up to 20"
    (is (= 24 (count (:gaussians-laguerre @reference))))
    (doseq [c (:gaussians-laguerre @reference)
            :let [p {:dimension (:dimension c) :degree (:degree c)}]]
      ;; measured 2.4e-14 relative to max(1, v) (degree 20)
      (is (< (max-rel-err (rbf/gaussians-laguerre p) c) 1.0e-13) (str p))
      (is (< (max-rel-err (rbf/gaussians-laguerre {:dimension (double (:dimension c)) :degree (double (:degree c))}) c) 1.0e-13) "doubles")))
  (testing "unnormalised: L_n^(d/2)(0) = C(n + d/2, n) at r = 0"
    (is (= 1.5 ((rbf/gaussians-laguerre {:dimension 1.0 :degree 1}) 0.0)))
    (is (= 2.0 ((rbf/gaussians-laguerre {:dimension 2.0 :degree 1}) 0.0)))
    (is (= 3.0 ((rbf/gaussians-laguerre {:dimension 4.0 :degree 1}) 0.0)))
    (is (m/delta-eq 4.375 ((rbf/gaussians-laguerre {:dimension 3.0 :degree 2}) 0.0) 1.0e-15))
    (is (= 1.0 ((rbf/gaussians-laguerre {:dimension 5.0 :degree 0}) 0.0))))
  (testing "boundaries"
    (let [f (rbf/gaussians-laguerre {:dimension 3.0 :degree 4})]
      (is (bits= (f 0.37) (f -0.37)))
      (is (= 0.0 (f 1.0e100)))
      (is (= 0.0 (f 30.0)))
      (is (some neg? (map f (range 0.0 5.0 0.05))) "oscillates for a positive degree")))
  (testing "scale"
    (doseq [sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/gaussians-laguerre {:dimension 2.0 :degree 2})
                  fs (rbf/gaussians-laguerre {:dimension 2.0 :degree 2 :scale sc})]
            x [0.0 0.3 1.0 7.0]]
      (is (bits= (f1 (/ x sc)) (fs x)))))
  (testing "invalid input"
    (is (invalid-input? #(rbf/gaussians-laguerre) {:family :gaussians-laguerre :property :dimension}))
    (is (invalid-input? #(rbf/gaussians-laguerre {:dimension 2.0}) {:family :gaussians-laguerre :property :degree}))
    (is (invalid-input? #(rbf/gaussians-laguerre {:degree 2}) {:property :dimension}))
    (is (invalid-input? #(rbf/gaussians-laguerre {:dimension 2.0 :degree -1}) {:property :degree}))
    (is (invalid-input? #(rbf/gaussians-laguerre {:dimension 2.0 :degree 1.5}) {:property :degree}))
    (is (invalid-input? #(rbf/gaussians-laguerre {:dimension 0.0 :degree 1}) {:property :dimension}))))

;; ------------------------------------------------------------------------------------------------ truncated power

(deftest truncated-power-test
  (testing "closed forms and the general path agree with (1 - r)^k"
    (doseq [k [m/THIRD 0.5 1.0 2.0 3.0 2.5 4.0 0.7]
            :let [f (rbf/truncated-power {:k k})]
            x [0.0 1.0e-9 0.3 0.9 0.999999]]
      (is (m/delta-eq (Math/pow (- 1.0 x) k) (f x) 1.0e-15 1.0e-14) (str k " " x))
      (is (bits= (f x) (f (- x))))))
  (testing "boundaries"
    (doseq [k [m/THIRD 0.5 1.0 2.0 3.0 2.5]
            :let [f (rbf/truncated-power {:k k})]]
      (is (= 1.0 (f 0.0)) (str k))
      (is (= 0.0 (f 1.0)) "exactly at the support radius")
      (is (= 0.0 (f 1.0000001)))
      (is (= 0.0 (f 1.0e300)))
      (is (= 0.0 (f ##Inf)))))
  (testing "k 0 is the indicator of the closed ball"
    (let [f (rbf/truncated-power {:k 0.0})]
      (is (= 1.0 (f 0.0)))
      (is (= 1.0 (f 0.999)))
      (is (= 1.0 (f 1.0)))
      (is (= 0.0 (f 1.0001)))))
  (testing "scale is the support radius"
    (doseq [sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/truncated-power {:k 2.0})
                  fs (rbf/truncated-power {:k 2.0 :scale sc})]
            x [0.0 0.3 1.0 1.2]]
      (is (bits= (f1 (/ x sc)) (fs x)))))
  (testing "invalid input"
    (is (invalid-input? #(rbf/truncated-power) {:family :truncated-power :property :k}))
    (is (invalid-input? #(rbf/truncated-power {:scale 2.0}) {:property :k}))))

;; ------------------------------------------------------------------------------------------------ linear, gaussian

(deftest linear-gaussian-test
  (testing "linear: abs(r / scale)"
    (let [f (rbf/linear)]
      (is (= 0.0 (f 0.0)))
      (is (= 2.5 (f -2.5)))
      (is (= 1.0e300 (f 1.0e300)))
      (is (= 1.0e-300 (f -1.0e-300))))
    (is (= 2.5 ((rbf/linear {:scale 2.0}) 5.0)))
    (is (= 2.5 ((rbf/linear {:scale 2.0 :unknown 1}) 5.0)))
    (is (= 2.5 ((rbf/linear {:scale nil}) 2.5))))
  (testing "gaussian: exp(-(r / scale)^2)"
    (let [f (rbf/gaussian)]
      (is (= 1.0 (f 0.0)))
      (is (m/delta-eq (Math/exp -1.0) (f 1.0) 1.0e-15))
      (is (m/delta-eq (Math/exp -1.0) (f -1.0) 1.0e-15))
      (is (= 0.0 (f 1.0e100)))
      (is (= 0.0 (f ##Inf)))
      (is (= 1.0 (f 1.0e-300))))
    (is (m/delta-eq (Math/exp -1.0) ((rbf/gaussian {:scale 3.0}) 3.0) 1.0e-15))
    (is (= 1.0 ((rbf/gaussian {:scale 1.0e300}) 1.0e10)))
    (is (= 0.0 ((rbf/gaussian {:scale 1.0e-300}) 1.0e-100)))))

;; ------------------------------------------------------------------------------------------------ signs

(defn- random-points
  [seed n d]
  (let [g (r/rng :mersenne seed)]
    (vec (repeatedly n (fn [] (vec (repeatedly d #(- (* 2.0 ^double (r/drandom g)) 1.0))))))))

(defn- exponents
  "Exponent vectors of the monomials of total degree <= deg in d variables."
  [d deg]
  (if (zero? d)
    [[]]
    (vec (for [e (range (inc deg))
               tail (exponents (dec d) (- deg e))]
           (into [e] tail)))))

(defn- jacobi-eigenvalues
  "Eigenvalues of the symmetric matrix `a` (a `RealMatrix`) by cyclic Jacobi rotations (robust for the
  clusters of zero eigenvalues of a semi-definite matrix, where the QL algorithm of Commons Math fails)."
  [^RealMatrix a]
  (let [n (.getRowDimension a)
        a (.copy a)
        off (fn [] (reduce + (for [i (range n) j (range n) :when (not= i j)] (let [x (.getEntry a i j)] (* x x)))))
        scale (.getFrobeniusNorm a)]
    (loop [sweep 0]
      (when (and (< sweep 100) (> (Math/sqrt (double (off))) (* 1.0e-15 scale)))
        (doseq [p (range (dec n)) q (range (inc p) n)
                :let [apq (.getEntry a p q)]
                :when (not (zero? apq))]
          (let [app (.getEntry a p p)
                aqq (.getEntry a q q)
                theta (/ (- aqq app) (* 2.0 apq))
                t (/ (Math/signum (+ theta 0.0)) (+ (Math/abs theta) (Math/sqrt (+ 1.0 (* theta theta)))))
                t (if (zero? theta) 1.0 t)
                c (/ 1.0 (Math/sqrt (+ 1.0 (* t t))))
                s (* t c)]
            (doseq [k (range n) :when (and (not= k p) (not= k q))]
              (let [akp (.getEntry a k p)
                    akq (.getEntry a k q)
                    np (- (* c akp) (* s akq))
                    nq (+ (* s akp) (* c akq))]
                (.setEntry a k p np) (.setEntry a p k np)
                (.setEntry a k q nq) (.setEntry a q k nq)))
            (.setEntry a p p (- app (* t apq)))
            (.setEntry a q q (+ aqq (* t apq)))
            (.setEntry a p q 0.0)
            (.setEntry a q p 0.0)))
        (recur (inc sweep))))
    (mapv #(.getEntry a % %) (range n))))

(defn- cpd-sign
  "Sign of the Gram matrix of the radial function `f` at the points `pts` restricted to the coefficient vectors
  which annihilate the polynomials of degree `deg` (all vectors for nil): 1 when every eigenvalue of the projected
  matrix is non-negative (and the largest is positive), -1 when every one is non-positive, 0 when mixed.
  Eigenvalues below 1e-9 of the largest in magnitude count as zero."
  [f pts deg]
  (let [n (count pts)
        d (count (first pts))
        dist (fn [a b] (Math/sqrt (reduce + (map (fn [x y] (let [t (- (double x) (double y))] (* t t))) a b))))
        gram (Array2DRowRealMatrix. (into-array (map double-array
                                                            (for [a pts] (for [b pts] (double (f (dist a b))))))))
        z (if (nil? deg)
            (Array2DRowRealMatrix. (into-array (map double-array (for [i (range n)] (for [j (range n)] (if (== i j) 1.0 0.0))))))
            (let [ex (exponents d deg)
                  q (count ex)
                  p (Array2DRowRealMatrix. (into-array (map double-array
                                                                   (for [pt pts]
                                                                     (for [e ex] (reduce * (map (fn [x k] (Math/pow (double x) (double k))) pt e)))))))
                  qq (.getQ (QRDecomposition. p))]
              (.getSubMatrix ^RealMatrix qq 0 (dec n) q (dec n))))
        a (.multiply (.transpose ^RealMatrix z) (.multiply gram ^RealMatrix z))
        amax (.getNorm ^RealMatrix a)
        sym (.scalarMultiply (.add ^RealMatrix a (.transpose ^RealMatrix a)) (/ 0.5 amax))
        ev (jacobi-eigenvalues sym)
        top (reduce max (map #(Math/abs (double %)) ev))
        tol (* 1.0e-9 top)]
    (cond
      (and (every? #(>= (double %) (- tol)) ev) (some #(> (double %) tol) ev)) 1
      (and (every? #(<= (double %) tol) ev) (some #(< (double %) (- tol)) ev)) -1
      :else 0)))

(defn- sign-in-dims
  "Signs of the Gram matrix on 40 seeded random points in R^1, R^2 and R^3."
  [f deg]
  (mapv (fn [d] (cpd-sign f (random-points (+ 100 d) 40 d) deg)) [1 2 3]))

(def ^:private gmq-auto-signs
  "Sign of `generalized-multiquadratic` with `:negate? :auto`, beta -> sign of (1 + r^2)^beta, ceil(beta) odd means -1."
  {-4.0 1 -3.0 1 -2.0 1 -1.0 1 -0.5 1 0.5 -1 1.0 -1 1.5 1 2.0 1 2.5 -1 3.0 -1 3.5 1 4.0 1 4.5 -1})

(deftest generalized-multiquadratic-sign-test
  (testing "the :auto sign table is (-1)^ceil(beta) for beta > 0 and +1 otherwise"
    (doseq [[beta expected] gmq-auto-signs]
      (is (= (double expected) ((rbf/generalized-multiquadratic {:beta beta}) 0.0)) (str beta))
      (is (= (double expected) (Math/signum ((rbf/generalized-multiquadratic {:beta beta}) 0.7))) (str beta))))
  (testing "the sign makes the kernel conditionally positive definite (seeded numerical check, R^1..R^3, 40 points)"
    ;; the Gram matrix is projected on the coefficient vectors annihilating the polynomials of degree
    ;; ceil(beta) - 1 (order ceil(beta)), no projection for beta < 0 (positive definite); beta up to 3.5 is checked,
    ;; for the highest beta of the table the check is not made (conditioning)
    (doseq [beta (filter #(<= % 3.5) (keys gmq-auto-signs))
            :let [deg (when (pos? beta) (dec (long (Math/ceil beta))))
                  auto (rbf/generalized-multiquadratic {:beta beta})
                  neg (rbf/generalized-multiquadratic {:beta beta :negate? (= 1 (gmq-auto-signs beta))})]]
      (is (= [1 1 1] (sign-in-dims auto deg)) (str "auto beta " beta))
      (is (= [-1 -1 -1] (sign-in-dims neg deg)) (str "the opposite sign is not conditionally positive definite, beta " beta))))
  (testing "negate? true, false, nil"
    (is (= -1.5 ((rbf/generalized-multiquadratic {:beta 1.0 :negate? true}) (Math/sqrt 0.5))))
    (is (= 1.5 ((rbf/generalized-multiquadratic {:beta 1.0 :negate? false}) (Math/sqrt 0.5))))
    (is (= -1.0 ((rbf/generalized-multiquadratic {:beta 1.0 :negate? true}) 0.0)))
    (is (= 1.0 ((rbf/generalized-multiquadratic {:beta 2.0 :negate? false}) 0.0)))
    (is (= -1.0 ((rbf/generalized-multiquadratic {:beta 0.5 :negate? nil}) 0.0)) "nil is :auto")
    (is (= 1.0 ((rbf/generalized-multiquadratic {:beta 0.5 :negate? false}) 0.0)))))

(deftest generalized-multiquadratic-test
  (testing "closed forms and the general path"
    (doseq [beta [m/THIRD (- m/THIRD) 0.5 -0.5 1.0 -1.0 2.0 -2.0 3.0 -3.0 1.7 -1.7 0.0 6.0]
            :let [f (rbf/generalized-multiquadratic {:beta beta :negate? false})]
            x [0.0 1.0e-9 0.3 1.0 5.0 -5.0]]
      (is (m/delta-eq (Math/pow (+ 1.0 (* x x)) beta) (f x) 1.0e-14 1.0e-14) (str beta " " x))))
  (testing "boundaries"
    (is (= 1.0 ((rbf/generalized-multiquadratic {:beta -0.5}) 0.0)))
    (is (= 0.0 ((rbf/generalized-multiquadratic {:beta -0.5}) 1.0e300)) "inverse multiquadric decays")
    (is (= 1.0 ((rbf/generalized-multiquadratic {:beta -0.5}) 1.0e-300)))
    (is (Double/isInfinite ((rbf/generalized-multiquadratic {:beta 0.5 :negate? false}) 1.0e200)) "grows, overflows")
    (is (m/delta-eq 1.0 (/ ((rbf/generalized-multiquadratic {:beta 0.5 :negate? false}) 1.0e100) 1.0e100) 1.0e-12)))
  (testing "scale (the shape parameter is 1 / scale)"
    (doseq [sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/generalized-multiquadratic {:beta 1.5})
                  fs (rbf/generalized-multiquadratic {:beta 1.5 :scale sc})]
            x [0.0 0.3 1.0 7.0]]
      (is (bits= (f1 (/ x sc)) (fs x)))))
  (testing "invalid input"
    (is (invalid-input? #(rbf/generalized-multiquadratic) {:family :generalized-multiquadratic :property :beta}))
    (is (invalid-input? #(rbf/generalized-multiquadratic {:negate? true}) {:property :beta}))))

(deftest radial-powers-test
  (testing "the sign is (-1)^ceil(beta / 2), unchanged"
    (doseq [[beta expected] {-1.0 1.0 0.5 -1.0 1.0 -1.0 2.0 -1.0 2.5 1.0 3.0 1.0 4.0 1.0 4.5 -1.0 5.0 -1.0 6.0 -1.0}]
      (is (= expected (Math/signum ((rbf/radial-powers {:beta beta}) 2.0))) (str beta))))
  (testing "values"
    (doseq [beta [1.0 3.0 2.5 0.5 -1.0 -2.5]
            x [0.0 0.3 1.0 5.0 -5.0]
            :when (not (and (zero? x) (neg? beta)))]
      (is (m/delta-eq (Math/pow (Math/abs x) beta) ((rbf/radial-powers {:beta beta :negate? false}) x) 1.0e-14 1.0e-14))
      (is (m/delta-eq (- (Math/pow (Math/abs x) beta)) ((rbf/radial-powers {:beta beta :negate? true}) x) 1.0e-14 1.0e-14))))
  (testing "the auto sign gives conditionally positive definite kernels in R^1..R^3 (order ceil(beta / 2))"
    (doseq [beta [1.0 1.5 2.5 3.0 3.5 5.0]
            :let [deg (dec (long (Math/ceil (/ beta 2.0))))]]
      (is (= [1 1 1] (sign-in-dims (rbf/radial-powers {:beta beta}) deg)) (str beta))))
  (testing "boundaries"
    (is (= 0.0 (Math/abs ^double ((rbf/radial-powers {:beta 2.5}) 0.0))))
    (is (Double/isInfinite ((rbf/radial-powers {:beta -1.0}) 0.0)))
    (is (= 1.0e-300 ((rbf/radial-powers {:beta -1.0}) 1.0e300)))
    (is (= 1.0 ((rbf/radial-powers {:beta 0.0 :negate? false}) 0.0)))
    (is (= -1.0e-300 ((rbf/radial-powers {:beta 1.0}) 1.0e-300))))
  (testing "scale"
    (doseq [sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/radial-powers {:beta 2.5})
                  fs (rbf/radial-powers {:beta 2.5 :scale sc})]
            x [0.0 0.3 1.0 7.0]]
      (is (bits= (f1 (/ x sc)) (fs x)))))
  (testing "invalid input"
    (is (invalid-input? #(rbf/radial-powers) {:family :radial-powers :property :beta}))))

(deftest thin-plate-splines-test
  (testing "the :auto sign is (-1)^(ceil(beta) + 1) for integer beta: + - + -"
    (doseq [[beta expected] {1.0 1.0 2.0 -1.0 3.0 1.0 4.0 -1.0}]
      ;; at r = e the logarithm is 1
      (is (= expected (Math/signum ((rbf/thin-plate-splines {:beta beta}) Math/E))) (str beta))))
  (testing "non-integer beta: the rule (-1)^(ceil(beta) + 1) is applied but is uncertain"
    (doseq [[beta expected] {0.5 1.0 1.5 -1.0 2.5 1.0 3.5 -1.0 4.5 1.0}]
      (is (= expected (Math/signum ((rbf/thin-plate-splines {:beta beta}) Math/E))) (str beta))))
  (testing "the sign makes the kernel conditionally positive definite for beta 1, 2, 3 (seeded, R^1..R^3, order beta + 1)"
    (doseq [beta [1.0 2.0 3.0]]
      (is (= [1 1 1] (sign-in-dims (rbf/thin-plate-splines {:beta beta}) (long beta))) (str beta))
      (is (= [-1 -1 -1] (sign-in-dims (rbf/thin-plate-splines {:beta beta :negate? (odd? (long beta))}) (long beta)))
          (str "the opposite sign is not conditionally positive definite, beta " beta))))
  (testing "values and boundaries"
    (doseq [beta [1.0 2.0 1.5]
            x [1.0e-9 0.3 1.0 5.0 -5.0]]
      (is (m/delta-eq (* (Math/pow (Math/abs x) (* 2.0 beta)) (Math/log (Math/abs x)))
                      ((rbf/thin-plate-splines {:beta beta :negate? false}) x) 1.0e-14 1.0e-13)))
    (is (= 0.0 ((rbf/thin-plate-splines {:beta 1.0}) 0.0)))
    (is (= 0.0 ((rbf/thin-plate-splines {:beta 1.0}) 1.0)))
    (is (= 0.0 ((rbf/thin-plate-splines {:beta 1.0}) -1.0)))
    (is (neg? ((rbf/thin-plate-splines {:beta 1.0}) 0.5)) "r^2 log r < 0 for r < 1, positive sign")
    (is (m/delta-eq 0.0 ((rbf/thin-plate-splines {:beta 1.0}) 1.0e-100) 1.0e-190))
    (is (Double/isInfinite ((rbf/thin-plate-splines {:beta 1.0}) 1.0e200))))
  (testing "scale"
    (doseq [sc [1.0e-10 3.7 1.0e10]
            :let [f1 (rbf/thin-plate-splines {:beta 2.0})
                  fs (rbf/thin-plate-splines {:beta 2.0 :scale sc})]
            x [0.0 0.3 1.0 7.0]]
      (is (bits= (f1 (/ x sc)) (fs x)))))
  (testing "invalid input"
    (is (invalid-input? #(rbf/thin-plate-splines) {:family :thin-plate-splines :property :beta}))
    (is (invalid-input? #(rbf/thin-plate-splines {:beta 0.0}) {:property :beta}))
    (is (invalid-input? #(rbf/thin-plate-splines {:beta -1.0}) {:property :beta}))))

;; ------------------------------------------------------------------------------------------------ cross-cutting

(deftest required-parameters-test
  (let [required {:truncated-power [rbf/truncated-power {:k 1.0} [:k]]
                  :gaussians-laguerre [rbf/gaussians-laguerre {:dimension 1.0 :degree 1} [:dimension :degree]]
                  :poisson [rbf/poisson {:d 3.0} [:d]]
                  :matern [rbf/matern {:order 3} [:order]]
                  :generalized-multiquadratic [rbf/generalized-multiquadratic {:beta 0.5} [:beta]]
                  :radial-powers [rbf/radial-powers {:beta 1.0} [:beta]]
                  :thin-plate-splines [rbf/thin-plate-splines {:beta 1.0} [:beta]]
                  :whittaker [rbf/whittaker {:alpha 0 :k 2 :beta 1.0} [:alpha :k :beta]]
                  :shifted-surface-splines [rbf/shifted-surface-splines {:s 3 :beta 2.5} [:s :beta]]
                  :wendland [rbf/wendland {:s 3 :k 1} [:s :k]]
                  :gneiting [rbf/gneiting {:s 3 :l 5.0} [:s :l]]
                  :wu [rbf/wu {:l 2.0 :k 1.0} [:l :k]]}]
    (is (= 12 (count required)))
    (doseq [[family [ctor valid ks]] required]
      (is (fn? (ctor valid)) (str family " constructs with its minimal map"))
      (is (invalid-input? ctor {:family family :property (first ks)}) (str family " without arguments"))
      (is (invalid-input? #(ctor nil) {:family family :property (first ks)}))
      (is (invalid-input? #(ctor {}) {:family family :property (first ks)}))
      (is (invalid-input? #(ctor {:scale 2.0}) {:family family :property (first ks)}))
      (doseq [k ks]
        (is (invalid-input? #(ctor (dissoc valid k)) {:family family :property k}) (str family " " k))
        (is (invalid-input? #(ctor (assoc valid k nil)) {:family family :property k}) (str family " " k " nil"))
        (is (invalid-input? #(ctor (assoc valid k "x")) {:family family :property k}) (str family " " k " string"))
        (is (invalid-input? #(ctor (assoc valid k ##NaN)) {:family family :property k}) (str family " " k " NaN"))
        (is (invalid-input? #(ctor (assoc valid k ##Inf)) {:family family :property k}) (str family " " k " Inf"))
        (is (str/includes? (str (::message (thrown-data #(ctor (dissoc valid k))))) (name k))
            "the message names the key")))
    (testing "linear and gaussian have no required parameters"
      (is (fn? (rbf/linear)))
      (is (fn? (rbf/gaussian)))
      (is (fn? (rbf/linear nil)))
      (is (fn? (rbf/gaussian {}))))))

(deftest all-families-test
  (doseq [[family [ctor params]] full-params
          :let [f (ctor params)]]
    (testing (str family " returns a primitive function of a double")
      (is (instance? clojure.lang.IFn$DD f)))
    (testing (str family " is symmetric and never NaN for finite r")
      (doseq [x rs]
        (is (bits= (f x) (f (- x))) (str family " " x))
        (is (not (Double/isNaN (f x))) (str family " " x))))
    (testing (str family " scale divides the distance")
      (doseq [sc [1.0e-10 3.7 1.0e10]
              :let [fs (ctor (assoc params :scale sc))]
              x [0.0 0.3 1.0 7.0 1.0e-3]]
        (is (bits= (f (/ x sc)) (fs x)) (str family " " sc " " x))))
    (testing (str family " nil scale is the default")
      (is (bits= (f 0.37) ((ctor (assoc params :scale nil)) 0.37)))))
  (testing "constructors accept no arguments only when they have no required parameter"
    (is (some? (rbf/linear)))
    (is (some? (rbf/gaussian)))))

(deftest docstrings-test
  (let [publics (ns-publics 'fastmath.kernel.rbf)]
    (is (= 14 (count publics)))
    (doseq [[sym v] publics
            :let [doc (:doc (meta v))]]
      (is (string? doc) (str sym))
      (is (str/includes? doc "Parameters:") (str sym " has Parameters"))
      (is (str/includes? doc "Returns") (str sym " has Returns"))
      (is (str/includes? doc "See also [[") (str sym " has See also"))
      (is (str/includes? doc "`:scale`") (str sym " documents scale"))
      (is (re-find #"r = 0" doc) (str sym " documents the value at r = 0"))
      (is (not (str/includes? doc "\"")) (str sym " no double quotes"))
      (is (not (str/includes? doc "\\")) (str sym " no backslashes"))
      (is (not (str/includes? doc "```")) (str sym " no code blocks"))
      (is (not (re-find #"Wittaker|Shifter" doc)) (str sym " typos")))
    (is (str/includes? (:doc (meta (publics 'matern))) "modified Bessel function of the second kind"))
    (is (str/includes? (:doc (meta (publics 'whittaker))) "Whittaker"))
    (is (str/includes? (:doc (meta (publics 'shifted-surface-splines))) "Shifted surface splines"))))
