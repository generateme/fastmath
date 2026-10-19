(ns fastmath.interpolation.monotone-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.interpolation.monotone :as sut]))

(set! *warn-on-reflection* true)

;; Oracles (utils/fastmath/dev/generate_interp1d_reference.py):
;; - :pchip     scipy PchipInterpolator
;; - :monotone  an independent Fritsch-Carlson clone (arithmetic mean tangents, sign test, all intervals
;;              limited to the circle alpha^2 + beta^2 <= 9)
;; - :random-monotone (200) and :random-general (100) datasets for the shape properties

(def ^:private ref-data
  (delay (edn/read-string (slurp (io/resource "interpolation/interp1d_reference.edn")))))

(defn- property
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (= :invalid-input (:type d)) (:property d))))))

(defn- ymax ^double [ys] (double (reduce max 0.0 (map #(Math/abs (double %)) ys))))

(defn- grid
  "10000 points covering the domain of `xs` (the knots are not necessarily on the grid)."
  [xs]
  (let [lo (double (first xs))
        hi (double (last xs))]
    (map #(+ lo (* (- hi lo) (/ (double %) 9999.0))) (range 10000))))

(defn- violations
  "Number of intervals where `f` leaves [min(y_i, y_i+1), max(y_i, y_i+1)] by more than `tol` on the 10000 point grid,
  and the largest overshoot."
  [f xs ys tol]
  (let [xs (vec xs)
        ys (vec ys)
        n (count xs)
        pts (grid xs)
        interval (fn [x] (loop [i 0] (if (and (< i (- n 2)) (>= x (xs (inc i)))) (recur (inc i)) i)))
        res (map (fn [x]
                   (let [i (interval x)
                         lo (min (ys i) (ys (inc i)))
                         hi (max (ys i) (ys (inc i)))
                         v (f x)]
                     [i (max 0.0 (- lo v) (- v hi))]))
                 pts)]
    {:violating-intervals (count (distinct (map first (filter #(> (double (second %)) tol) res))))
     :overshoot (reduce max 0.0 (map second res))}))

(defn- decreases
  "Largest decrease of a non-decreasing (increase of a non-increasing) function between consecutive grid points."
  [f xs increasing?]
  (let [vs (map f (grid xs))
        sgn (if increasing? 1.0 -1.0)]
    (reduce max 0.0 (map (fn [a b] (* sgn (- (double a) (double b)))) vs (rest vs)))))

;;

(t/deftest monotone-edge-cases
  (t/testing "two points: the straight line, also outside"
    (doseq [ctor [sut/monotone sut/pchip]]
      (let [f (ctor [0 2] [1 5])]
        (t/is (m/delta-eq 1.0 (f 0) 1e-15))
        (t/is (m/delta-eq 3.0 (f 1) 1e-15))
        (t/is (m/delta-eq 5.0 (f 2) 1e-15))
        (t/is (m/delta-eq -1.0 (f -1) 1e-14))
        (t/is (m/delta-eq 9.0 (f 4) 1e-14)))))
  (t/testing "knots return the data exactly"
    (doseq [ctor [sut/monotone sut/pchip]]
      (let [xs [-3.0 -1.0 0.5 2.0 10.0]
            ys [4.0 -2.0 7.5 0.25 -8.0]
            f (ctor xs ys)]
        (t/is (= ys (mapv f xs))))))
  (t/testing "constant data stay constant"
    (doseq [ctor [sut/monotone sut/pchip]]
      (let [f (ctor [0 1 2 3] [5 5 5 5])]
        (t/is (every? #(m/delta-eq 5.0 (f %) 1e-15) [-1.0 0.3 1.5 2.9 4.0])))))
  (t/testing "plateaus stay flat, extrema are at the data points"
    (doseq [ctor [sut/monotone sut/pchip]]
      (let [f (ctor [0 1 2 3 4 5] [0 1 1 1 3 3])]
        (t/is (every? #(m/delta-eq 1.0 (f %) 1e-15) [1.0 1.25 1.5 1.75 2.0 2.5 3.0]))
        (t/is (every? #(m/delta-eq 3.0 (f %) 1e-15) [4.0 4.5 5.0])))))
  (t/testing "NaN query gives NaN; NaN in ys is allowed"
    (doseq [ctor [sut/monotone sut/pchip]]
      (t/is (m/nan? ((ctor [0 1 2] [0 1 4]) ##NaN)))
      (t/is (m/nan? ((ctor [0 1 2 3] [0 1 ##NaN 4]) 1.5)))))
  (t/testing "huge and tiny ys scale the interpolant (scale invariance)"
    (doseq [ctor [sut/monotone sut/pchip]
            c [1e150 1e-150]]
      (let [xs [0 1 2 3 4]
            ys [0.0 1.0 1.0 4.0 4.5]
            f (ctor xs ys)
            g (ctor xs (map #(* c %) ys))]
        (doseq [x [0.3 1.5 2.2 3.7 -0.5 5.0]]
          (t/is (m/delta-eq (* c (f x)) (g x) 0.0 1e-14)))))))

(t/deftest monotone-validation
  (doseq [ctor [sut/monotone sut/pchip]]
    (t/is (= :length-mismatch (property #(ctor [0 1 2] [0 1]))))
    (t/is (= :too-few-points (property #(ctor [] []))))
    (t/testing "one point threw NullPointerException"
      (t/is (= :too-few-points (property #(ctor [1] [1])))))
    (t/is (= :not-finite (property #(ctor [0 ##NaN 2] [0 1 2]))))
    (t/is (= :not-finite (property #(ctor [0 1 ##Inf] [0 1 2]))))
    (t/testing "unsorted xs gave a silently wrong value, duplicates NaN"
      (t/is (= :not-increasing (property #(ctor [0 2 1 3] [0 1 2 3]))))
      (t/is (= :not-increasing (property #(ctor [0 1 1 3] [0 1 2 3]))))
      (t/is (= :not-increasing (property #(ctor [3 2 1 0] [0 1 2 3])))))
    (t/is (nil? (property #(ctor [0 1] [0 1]))))))

(t/deftest monotone-shape-preservation
  ;; 200 seeded monotone datasets with plateaus and steep jumps, 10000 point grid:
  ;; before the repair 28 (no limiting of the last interval) and 85 (no sign test) of the datasets and
  ;; 50 of 200 in the audit of the old implementation violated; now 0.
  (t/testing "random monotone datasets: no decrease, no overshoot"
    (doseq [[k d] (map-indexed vector (:random-monotone @ref-data))]
      (let [xs (:xs d)
            ys (:ys d)
            up? (<= (first ys) (last ys))]
        ;; measured: overshoot 8.5e-14 and decrease 1.4e-13 (data range up to 400), absolute
        (doseq [ctor [sut/monotone sut/pchip]]
          (let [f (ctor xs ys)
                v (violations f xs ys 1e-12)]
            (t/is (zero? (:violating-intervals v)) (str "dataset " k " overshoot " (:overshoot v)))
            (t/is (<= (:overshoot v) 1e-12))
            (t/is (<= (decreases f xs up?) 1e-12) (str "dataset " k " not monotone"))))))))

(t/deftest monotone-no-overshoot-general
  ;; random data with sign changes and plateaus: every interval stays between its two data values
  ;; measured overshoot 7.1e-15 (absolute)
  (doseq [[k d] (map-indexed vector (:random-general @ref-data))]
    (let [xs (:xs d)
          ys (:ys d)]
      (doseq [ctor [sut/monotone sut/pchip]]
        (let [v (violations (ctor xs ys) xs ys 1e-12)]
          (t/is (zero? (:violating-intervals v)) (str "dataset " k " overshoot " (:overshoot v))))))))

(t/deftest monotone-known-failures
  ;; the secants 1 1 0.01 (monotone data): the last interval was not limited, maximum 2.07672 (limit 2.01)
  (t/testing "[0 1 2 2.01] does not exceed 2.01"
    (doseq [ctor [sut/monotone sut/pchip]]
      (let [f (ctor [0 1 2 3] [0 1 2 2.01])
            mx (reduce max (map f (grid [0 3])))]
        (t/is (<= mx (+ 2.01 1e-12)) (str "max " mx)))))
  ;; the secants 3 -1 -2: a local maximum at x = 1 needs the zero tangent (sign test), maximum 3.07655 (limit 3)
  (t/testing "[0 3 2 0] does not exceed 3"
    (doseq [ctor [sut/monotone sut/pchip]]
      (let [f (ctor [0 1 2 3] [0 3 2 0])
            mx (reduce max (map f (grid [0 3])))]
        (t/is (<= mx (+ 3.0 1e-12)) (str "max " mx)))))
  (t/testing "the zero tangent at the local extremum: derivative at the node is zero"
    (let [f (sut/monotone [0 1 2 3] [0 3 2 0])
          h 1e-7]
      ;; measured 2e-7; the arithmetic mean tangent without the sign test would give 1.0
      (t/is (m/delta-eq 0.0 (/ (- (f (+ 1.0 h)) (f (- 1.0 h))) (* 2.0 h)) 1e-5))))
  (t/testing "the demonstration data discriminate: the clone without each repair violates"
    (let [needed (:fc-needed @ref-data)]
      (t/is (zero? (:full needed)))
      (t/is (pos? (:no-last-interval needed)))
      (t/is (pos? (:no-sign-test needed)))
      (t/is (pos? (:neither needed))))))

(t/deftest monotone-values
  (t/testing "equals the Fritsch-Carlson clone (closed form up to rounding)"
    ;; measured: 5.3e-16 * max(ys) inside, 5.1e-15 relative outside
    (doseq [d (:monotone @ref-data)]
      (let [f (sut/monotone (:xs d) (:ys d))
            sc (ymax (:ys d))]
        (doseq [[q e] (map vector (:qs d) (:expected d))]
          (t/is (m/delta-eq e (f q) (* 1e-12 sc) 1e-12) (str (:name d) " at " q))))))
  (t/testing "hand-computed: uniform data 0 1 3 with both secants positive (xs 0 1 2)"
    ;; secants 1 and 2, m = [1 1.5 2], alpha/beta inside the circle -> no limiting;
    ;; Hermite on [0,1]: p(0.5) = (y0+y1)/2 + h (m0 - m1)/8 = 0.5 + (1 - 1.5)/8
    (let [f (sut/monotone [0 1 2] [0 1 3])]
      (t/is (m/delta-eq (+ 0.5 (/ (- 1.0 1.5) 8.0)) (f 0.5) 1e-15))
      (t/is (m/delta-eq (+ 2.0 (/ (- 1.5 2.0) 8.0)) (f 1.5) 1e-15)))))

(t/deftest pchip-values
  (t/testing "equals scipy PchipInterpolator"
    ;; measured: 5.3e-16 * max(ys) inside, 3.4e-11 relative outside (many200, up to two widths away: a cubic
    ;; amplifies the rounding of the tangents)
    (doseq [d (:pchip @ref-data)]
      (let [f (sut/pchip (:xs d) (:ys d))
            lo (first (:xs d))
            hi (last (:xs d))
            sc (ymax (:ys d))]
        (doseq [[q e] (map vector (:qs d) (:expected d))]
          (if (<= lo q hi)
            (t/is (m/delta-eq e (f q) (* 1e-13 sc) 1e-13) (str (:name d) " at " q " inside"))
            (t/is (m/delta-eq e (f q) 0.0 1e-9) (str (:name d) " at " q " outside"))))))))

(t/deftest pchip-hand-computed
  ;; uniform spacing 1: harmonic mean of secants 1 and 2 is 2/(1/1+1/2) = 4/3 at node 1;
  ;; end tangents of the three point formula: d0 = (3*1 - 2)/2 = 0.5, d2 = (3*2 - 1)/2 = 2.5 -> the
  ;; sign test of the second kind: sign(m0) == sign(m1) and abs(d) > 3 abs(m) is not reached
  (let [f (sut/pchip [0 1 2] [0 1 3])]
    (t/is (m/delta-eq (+ 0.5 (/ (- 0.5 (/ 4.0 3.0)) 8.0)) (f 0.5) 1e-15))
    (t/is (m/delta-eq (+ 2.0 (/ (- (/ 4.0 3.0) 2.5) 8.0)) (f 1.5) 1e-15))))

(t/deftest monotone-docstrings
  (t/is (re-find #"Fritsch-Carlson" (:doc (meta #'sut/monotone))))
  (t/is (re-find #"arithmetic mean" (:doc (meta #'sut/monotone))))
  (t/is (re-find #"harmonic mean" (:doc (meta #'sut/pchip))))
  (t/is (re-find #"PchipInterpolator" (:doc (meta #'sut/pchip)))))
