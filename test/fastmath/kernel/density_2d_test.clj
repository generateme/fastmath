(ns fastmath.kernel.density-2d-test
  "Tests of the 2d kernel density estimate (`fastmath.kernel/kernel-density-2d`, `fastmath.kernel.density/kernel-density-2d+`).
  Reference: `MASS::kde2d` of R with its default bandwidth (utils/fastmath/dev/kde2d_reference.R)."
  (:require [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.kernel :as k]
            [fastmath.kernel.density :as d]))

(def ^:private xs [1.2 2.5 3.1 0.4 2.2 4.0 3.3 1.9 2.8 0.9])
(def ^:private ys [5.0 3.2 4.4 6.1 2.7 3.9 5.5 4.8 3.0 5.9])

;; kde2d(x, y, n = 3, lims = c(1, 3, 3, 5)): z at x in (1 2 3) varying fastest, y in (3 4 5)
(def ^:private kde2d-z [0.0146642877356993 0.0634118509782272 0.0713087853329909
                        0.0275996421766086 0.0599303312367149 0.0725221667279065
                        0.0596464768923657 0.0606284207894995 0.0573334883620581])

(defn- cloud
  [^long n seed]
  (let [rng (java.util.Random. seed)
        xs (vec (repeatedly n #(.nextGaussian rng)))]
    [xs (mapv #(+ (* 0.5 (double %)) (.nextGaussian rng)) xs)]))

(defn- integral
  "Midpoint rule of `f` over its support, `n` x `n` cells."
  [{:keys [kde mn mx]} ^long n]
  (let [[x0 y0] mn
        [x1 y1] mx
        dx (/ (- (double x1) (double x0)) n)
        dy (/ (- (double y1) (double y0)) n)]
    (* dx dy (reduce + (for [i (range n) j (range n)]
                         (kde (+ (double x0) (* dx (+ i 0.5))) (+ (double y0) (* dy (+ j 0.5)))))))))

(deftest mass-kde2d-test
  (testing "default :nrd per axis equals MASS::kde2d"
    (let [{:keys [kde h]} (k/kernel-density-2d :gaussian xs ys nil true)]
      (is (m/delta-eq 0.763250229428961 (first h)))
      (is (m/delta-eq 0.815631139896643 (second h)))
      (is (every? true? (map (fn [[x y] z] (m/delta-eq z (kde x y)))
                             (for [y [3.0 4.0 5.0] x [1.0 2.0 3.0]] [x y])
                             kde2d-z))))))

(deftest product-kernel-test
  (testing "one point: the product of the axis kernels"
    (let [f (k/kernel-density-2d :gaussian [1.0] [2.0] [0.5 2.0])
          normal (fn [^double v ^double s] (/ (m/exp (* -0.5 (m/sq (/ v s)))) (* s (m/sqrt m/TWO_PI))))]
      (is (m/delta-eq (* (normal 0.3 0.5) (normal 0.5 2.0)) (f 1.3 2.5)))
      (is (m/delta-eq (* (normal 0.0 0.5) (normal 0.0 2.0)) (f 1.0 2.0)))))
  (testing "integrates to 1, every kernel"
    (let [[cx cy] (cloud 200 3)]
      (doseq [kernel [:gaussian :epanechnikov :uniform :quartic :laplace]]
        (is (< (m/abs (- 1.0 (double (integral (k/kernel-density-2d kernel cx cy nil true) 300)))) 1.0e-3) (str kernel)))))
  (testing "exact sum with the cutoff equals a naive sum"
    (let [[cx cy] (cloud 300 7)
          {:keys [kde h]} (k/kernel-density-2d :epanechnikov cx cy :nrd true)
          [hx hy] h
          naive (fn [x y] (/ (reduce + (map (fn [xi yi] (* (d/epanechnikov (/ (- x xi) hx)) (d/epanechnikov (/ (- y yi) hy)))) cx cy))
                             (* (count cx) hx hy)))]
      (doseq [[x y] [[0.0 0.0] [1.5 -0.3] [-2.0 1.0] [5.0 5.0]]]
        (is (< (m/abs (- (double (kde x y)) (double (naive x y)))) 1.0e-12) (str [x y]))))))

(deftest bandwidth-2d-test
  (let [[cx cy] (cloud 100 11)
        h (fn [b] (:h (k/kernel-density-2d :gaussian cx cy b true)))]
    (testing "a number, a pair, a selector per axis"
      (is (= [0.3 0.3] (h 0.3)))
      (is (= [0.5 2.0] (h [0.5 2.0])))
      (is (= [(k/bandwidth :gaussian cx :nrd0) (k/bandwidth :gaussian cy :sj)] (h [:nrd0 :sj])))
      (is (= [(k/bandwidth :gaussian cx :nrd) (k/bandwidth :gaussian cy :nrd)] (h :nrd) (h {:bandwidth nil}))))
    (testing ":scott - sd * n^(-1/6)"
      (let [[hx] (h :scott)
            mean (/ (reduce + cx) (count cx))
            sd (m/sqrt (/ (reduce + (map #(m/sq (- (double %) mean)) cx)) (dec (count cx))))]
        (is (m/delta-eq (* sd (m/pow 100 (/ -1.0 6.0))) hx))))
    (testing "support: the ranges extended by 1.01 h extent"
      (let [{:keys [mn mx] [hx hy] :h} (k/kernel-density-2d :epanechnikov cx cy :nrd true)]
        (is (m/delta-eq (- (double (reduce min cx)) (* 1.01 hx)) (first mn)))
        (is (m/delta-eq (+ (double (reduce max cy)) (* 1.01 hy)) (second mx)))))))

(deftest binned-2d-test
  (testing "binned within 1 % of exact, faster for many points"
    (let [[cx cy] (cloud 10000 13)
          exact (k/kernel-density-2d :gaussian cx cy)
          binned (k/kernel-density-2d :gaussian cx cy {:binned? true})]
      (doseq [[x y] [[0.0 0.0] [1.0 0.5] [-1.5 -1.0]]]
        (is (< (m/abs (- 1.0 (/ (double (binned x y)) (double (exact x y))))) 0.01) (str [x y])))))
  (testing "a divisor"
    (let [[cx cy] (cloud 500 17)]
      (is (< (m/abs (- 1.0 (/ (double ((k/kernel-density-2d :gaussian cx cy {:binned? 10}) 0.2 0.1))
                              (double ((k/kernel-density-2d :gaussian cx cy) 0.2 0.1)))))
             0.01)))))

(deftest validation-2d-test
  (let [err (fn [f] (try (f) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))]
    (testing "invalid input throws :invalid-input"
      (is (= :invalid-input (:type (err #(k/kernel-density-2d :gaussian [] [])))))
      (is (= :ys (:property (err #(k/kernel-density-2d :gaussian [1 2 3] [1 2])))))
      (is (= :xs (:property (err #(k/kernel-density-2d :gaussian [1 ##NaN] [1 2])))))
      (is (= :bandwidth (:property (err #(k/kernel-density-2d :gaussian [1 2] [1 2] -1.0)))))
      (is (= :bandwidth (:property (err #(k/kernel-density-2d :gaussian [1 2] [1 2] [0.1 0.2 0.3])))))
      (is (= :bandwidth (:property (err #(k/kernel-density-2d :gaussian [1 2] [1 2] :bogus)))))
      (is (= :data (:property (err #(k/kernel-density-2d :gaussian [1] [1])))) "a selector on a single point")
      (is (= :kernel (:property (err #(k/kernel-density-2d :bogus [1 2] [1 2])))))
      (is (= :binned? (:property (err #(k/kernel-density-2d :gaussian [1 2] [1 2] {:binned? -1}))))))
    (testing "one point with an explicit bandwidth; a constant axis falls back as the 1d rule"
      (is (pos? (double ((k/kernel-density-2d :gaussian [1] [1] 0.5) 1 1))))
      (is (pos? (double (first (:h (k/kernel-density-2d :gaussian [2 2 2] [1 2 3] nil true)))))))
    (testing "NaN query"
      (is (m/nan? ((k/kernel-density-2d :gaussian [1 2] [1 2] 0.5) ##NaN 1.0))))
    (testing "the inputs are not modified"
      (let [a (double-array [3.0 1.0 2.0]) b (double-array [1.0 2.0 3.0])]
        (k/kernel-density-2d :gaussian a b 0.5)
        (is (= [3.0 1.0 2.0] (vec a)))))))
