(ns fastmath.special.en-test
  "Tests of the generalized exponential integral `fastmath.special/En`.

  Reference values come from `mpmath.expint(n, x)` (50 digits, cross-checked at
  70 digits) and are stored in `test/resources/special/en_reference.edn`; see
  `utils/fastmath/dev/generate_en_reference.py`.  Arguments in the file are the
  exact doubles used by the oracle.

  Measured relative error of the current implementation over the whole file
  (18239 grid points and 2059 further points): at most 1.6e-14 for x <= 1.5,
  at most 5.3e-13 overall (worst at n = 3.45, x = 2.96, where the power series
  of the `x < 3` branch loses digits to alternating terms; the continued
  fraction for x >= 3 reaches 2.5e-13, worst at n = 2, x = 50)."
  (:require [fastmath.special :as sut]
            [fastmath.core :as m]
            [clojure.test :as t]
            [clojure.java.io :as io]
            [clojure.edn :as edn]))

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "special/en_reference.edn")))))

(defn- tol
  "Relative tolerance by region: tight where the series is well conditioned
  (observed 1.6e-14), looser close to the series / continued fraction switch
  at x = 3 (observed 5.3e-13)."
  ^double [^double x]
  (if (m/<= x 1.5) 1.0e-13 2.0e-12))

(defn- block-rows
  "Rows `[order arg ref]` of a parallel-vector block."
  [k]
  (let [{:keys [order arg ref]} (get @reference k)]
    (mapv vector order arg ref)))

(defn- grid-rows
  "Rows `[order arg ref]` of the row-major order x argument grid."
  []
  (let [{:keys [orders xs ref]} (:grid @reference)
        nx (count xs)]
    (vec (for [i (range (count orders))
               j (range nx)]
           [(orders i) (xs j) (ref (+ (* i nx) j))]))))

(defn- bad-points
  "Rows where `En` differs from the reference by more than the relative `(tol-fn x)`
  (and, when given, more than the absolute `abs-tol`, for results at a real zero)."
  ([rows tol-fn] (bad-points rows tol-fn 0.0))
  ([rows tol-fn ^double abs-tol]
   (into []
         (comp (map (fn [[n x r]] {:order n :arg x :ref r :got (sut/En n x)}))
               (remove (fn [{:keys [^double got ^double ref ^double arg]}]
                         (m/delta-eq got ref abs-tol (tol-fn arg)))))
         rows)))

(defn- rel-err
  ^double [^double a ^double b]
  (if (m/== a b) 0.0 (m// (m/abs (m/- a b)) (m/abs b))))

(t/deftest en-edge-cases
  (t/testing "NaN"
    (t/is (m/nan? (sut/En ##NaN 1.0)))
    (t/is (m/nan? (sut/En 1.5 ##NaN)))
    (t/is (m/nan? (sut/En ##NaN ##NaN)))
    (t/is (m/nan? (sut/En ##NaN 0.0))))
  (t/testing "x = 0: infinite for n < 1 (n = 0 and n = 1 included), 1/(n-1) for n > 1"
    (t/is (m/pos-inf? (sut/En 0.0 0.0)))
    (t/is (m/pos-inf? (sut/En 1.0 0.0)))
    (t/is (m/pos-inf? (sut/En 0.5 0.0)))
    (t/is (m/pos-inf? (sut/En 0.999999 0.0)))
    (t/is (m/pos-inf? (sut/En -0.5 0.0)))
    (t/is (m/pos-inf? (sut/En -3.0 0.0)))
    (t/is (m/delta-eq 2.0 (sut/En 1.5 0.0) 1.0e-15))
    (t/is (m/delta-eq 1.0 (sut/En 2.0 0.0) 1.0e-15))
    (t/is (m/delta-eq 0.5 (sut/En 3.0 0.0) 1.0e-15))
    (t/is (m/delta-eq 0.2 (sut/En 6.0 0.0) 1.0e-15))
    (t/is (m/delta-eq 1.0e6 (sut/En 1.000001 0.0) 1.0e-3 1.0e-9))
    (t/is (m/delta-eq 1.0 (sut/En 2.0 -0.0) 1.0e-15)))
  (t/testing "x tiny: finite limit 1/(n-1) for n > 1, divergence for n < 1"
    (t/is (m/delta-eq 2.0 (sut/En 1.5 1.0e-300) 1.0e-15))
    (t/is (m/delta-eq 1.0 (sut/En 2.0 5.0e-324) 1.0e-15))
    (t/is (m/delta-eq 0.4 (sut/En 3.5 1.0e-300) 1.0e-15))
    (t/is (m/delta-eq 1.7724538509055159E150 (sut/En 0.5 1.0e-300) 0.0 1.0e-13))
    (t/is (m/pos-inf? (sut/En -1.5 1.0e-300))))
  (t/testing "large x: underflows to 0"
    (t/is (m/zero? (sut/En 1.5 746.0)))
    (t/is (m/zero? (sut/En 2.0 1000.0)))
    (t/is (m/zero? (sut/En -0.5 1.0e5)))
    (t/is (m/zero? (sut/En 1.5 ##Inf))))
  (t/testing "negative x: real only for a non-positive integer order, NaN otherwise"
    (t/is (m/nan? (sut/En 1.0 -1.0)))
    (t/is (m/nan? (sut/En 2.0 -1.0)))
    (t/is (m/nan? (sut/En 5.0 -0.001)))
    (t/is (m/nan? (sut/En 0.5 -1.0)))
    (t/is (m/nan? (sut/En 1.5 -0.001)))
    (t/is (m/nan? (sut/En -0.5 -1.0)))
    (t/is (m/nan? (sut/En -1.5 -2.0)))
    (t/is (m/nan? (sut/En 1.0e-9 -1.0)))
    (t/is (m/valid-double? (sut/En 0.0 -1.0)))
    (t/is (m/valid-double? (sut/En -1.0 -1.0)))
    (t/is (m/valid-double? (sut/En -2.0 -5.0)))))

(t/deftest en-dispatch-to-e0-and-e1
  (doseq [[n x _] (block-rows :dispatch)]
    (if (m/zero? n)
      (t/is (== (sut/E0 x) (sut/En n x)) (str "n = 0, x = " x))
      (t/is (== (sut/E1 x) (sut/En n x)) (str "n = 1, x = " x)))))

(t/deftest en-spike-points
  ;; Isolated points where the former near-integer Taylor branch fired by accident
  ;; for orders 0.45 away from an integer (relative errors 5.28, 2.11, 0.81, 0.62).
  ;; References: mpmath.expint at 50 digits.
  (doseq [[n x r] [[1.55 2.88 0.01349696706404477]
                   [1.6 2.66 0.0176199827368727]
                   [1.65 2.46 0.02248144474201001]
                   [1.65 2.47 0.022195482367631924]]]
    (t/is (< (rel-err (sut/En n x) r) 1.0e-12) (str "n = " n ", x = " x))))

(t/deftest en-vs-mpmath-grid
  ;; 61 orders (0.5 .. 3.5 step 0.05, integers go through the positive-integer
  ;; series) x 299 arguments (0.01 .. 2.99): the former code had 87 points above
  ;; 1e-8, all in the near-integer Taylor branch.
  (let [rows (grid-rows)]
    (t/is (== 18239 (count rows)))
    (t/is (empty? (bad-points rows tol)))
    (t/testing "no point above 1e-10 anywhere in the grid"
      (t/is (empty? (bad-points rows (constantly 1.0e-10)))))))

(t/deftest en-vs-mpmath-near-integer
  ;; orders k +- d for d from 0.4 down to 1e-12, exactly k, and around 0
  ;; (k = 1, 2, 3, 6, 11); x up to 2.999999. The delicate region is |d| < 1e-2,
  ;; where the series term with the vanishing denominator cancels the Gamma term.
  (t/is (empty? (bad-points (block-rows :near-integer) tol)))
  (t/testing "offsets 1e-3, 1e-6, 1e-9, 1e-12 on both sides of k, individually, all x < 3"
    (let [lookup (into {} (map (fn [[n x r]] [[n x] r])) (block-rows :near-integer))]
      (doseq [k [1.0 2.0 3.0 6.0 11.0]
              d [1.0e-3 1.0e-6 1.0e-9 1.0e-12]
              sgn [-1.0 1.0]
              :let [n (m/+ k (m/* sgn d))]
              x [0.01 0.5 2.0 2.9]
              :let [r (get lookup [n x])]]
        (t/is (some? r) (str "reference present for n = " n ", x = " x))
        (t/is (< (rel-err (sut/En n x) r) 2.0e-12) (str "n = " n ", x = " x))))))

(t/deftest en-vs-mpmath-regions
  (t/testing "x tiny (1e-300 .. 1e-5)"
    (t/is (empty? (bad-points (block-rows :tiny-x) (constantly 1.0e-13)))))
  (t/testing "around the series / continued fraction switch at x = 3"
    (t/is (empty? (bad-points (block-rows :branch-switch) tol))))
  (t/testing "x from 3 to 700 (continued fraction)"
    (t/is (empty? (bad-points (block-rows :large-x) (constantly 1.0e-12)))))
  (t/testing "large fractional orders (4.5 .. 1200.5)"
    (t/is (empty? (bad-points (block-rows :large-order) tol))))
  (t/testing "negative orders, fractional and integer, and negative x for integer orders"
    ;; absolute 1e-15: E-1(-1) is exactly 0 (observed 1.1e-16)
    (t/is (empty? (bad-points (block-rows :negative-order) (constantly 1.0e-12) 1.0e-15))))
  (t/testing "positive integer orders"
    (t/is (empty? (bad-points (block-rows :posint) tol)))))

(t/deftest en-continuity
  (t/testing "across the order n = k (series with vanishing denominator vs integer series), x < 3"
    (doseq [k [1.0 2.0 3.0 4.0 7.0]
            x [0.01 0.5 1.5 2.9]
            d [1.0e-9 -1.0e-9]]
      ;; |dEn/dn| is of the order of En itself, so the relative step is about 1e-9
      (t/is (< (rel-err (sut/En (m/+ k d) x) (sut/En k x)) 1.0e-8) (str "k = " k ", x = " x))))
  (t/testing "across the switch between the n = 0.5, 1.5, 2.5 rounding cells"
    (doseq [n [0.5 1.5 2.5 3.5]
            x [0.1 1.0 2.9]]
      (t/is (< (rel-err (sut/En (m/- n 1.0e-9) x) (sut/En (m/+ n 1.0e-9) x)) 1.0e-8)
            (str "n = " n ", x = " x))))
  (t/testing "across the series / continued fraction switch at x = 3"
    (doseq [n [0.5 1.5 2.0 2.5 3.7 -0.5 -2.0]]
      (t/is (< (rel-err (sut/En n (m/- 3.0 1.0e-9)) (sut/En n (m/+ 3.0 1.0e-9))) 1.0e-8)
            (str "n = " n)))))

(t/deftest en-structure
  (t/testing "En decreases in x for fixed n on the grid"
    (doseq [n [0.5 1.55 2.0 2.45 3.5]
            :let [ys (mapv #(sut/En n %) (map #(m/* 0.01 %) (range 1 300)))]]
      (t/is (every? true? (map > ys (rest ys))) (str "n = " n))))
  (t/testing "En decreases in n for fixed x > 0"
    (doseq [x [0.01 0.5 1.0 2.88 3.0 5.0]
            :let [ys (mapv #(sut/En (m/* 0.05 %) x) (range 10 70))]]
      (t/is (every? true? (map > ys (rest ys))) (str "x = " x))))
  (t/testing "recurrence n * En+1(x) = exp(-x) - x * En(x), including the former spike orders"
    (doseq [n [0.5 1.55 1.6 1.65 2.3 2.7 3.45 -0.5 -1.5]
            x [0.05 0.5 1.0 2.46 2.66 2.88 2.99]
            :let [lhs (m/* n (sut/En (m/inc n) x))
                  rhs (m/- (m/exp (m/- x)) (m/* x (sut/En n x)))]]
      ;; cancellation in rhs: absolute 1e-13 or relative 1e-11
      (t/is (m/delta-eq lhs rhs 1.0e-13 1.0e-11) (str "n = " n ", x = " x))))
  (t/testing "derivative d/dx En(x) = -E(n-1)(x), central difference"
    (doseq [n [0.5 1.55 1.6 2.3]
            x [0.5 1.0 2.66 2.88]
            :let [h 1.0e-5
                  dd (m// (m/- (sut/En n (m/+ x h)) (sut/En n (m/- x h))) (m/* 2.0 h))]]
      ;; central difference truncation error ~ h^2 En'''/6 ~ 1e-11
      (t/is (m/delta-eq dd (m/- (sut/En (m/dec n) x)) 1.0e-8 1.0e-8) (str "n = " n ", x = " x)))))
