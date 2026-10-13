(ns fastmath.calculus.vegas-test
  (:require [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.special :as special]
            [fastmath.calculus.vegas :as sut]))

;; Every call below is seeded, so the tests are deterministic. The statistical tests do not compare with a tolerance of
;; one run, they check that the reported standard deviation describes the errors of 30 seeded runs.

;; https://github.com/JohannesBuchner/cuba/blob/master/demo/demo-c.c
;; + wolfram alpha results

(def fun1-res 0.664669679781377268445214)
(defn fun1 [[^double x ^double y ^double z]] (m/* (m/sin x) (m/cos y) (m/exp z)))

(def fun1+-res 3.019450507398802024611853)
(defn fun1+ [[^double x ^double y ^double z]] (m/+ (m/sin x) (m/cos y) (m/exp z)))

(def fun2-res 5.26852)
(defn fun2 [[^double x ^double y ^double z]] (m/* (m// (m/+ (m/sq (m/+ x y)) 0.003))
                                               (m/cos y) (m/exp z)))

(def fun3-res 0.307807)
(defn fun3 [[^double x ^double y ^double z]] (m// (m/- 3.75 (m/cos (m/* m/PI x))
                                                    (m/cos (m/* m/PI y))
                                                    (m/cos (m/* m/PI z)))))

(def fun4-res 0.877314)
(defn fun4 [v] (m/abs (m/- (v/magsq v) 0.125)))

(def fun5-res 0.416538385886638169609660243601001780010941391101)
(defn fun5 [v] (m/exp (m/- (v/magsq v))))

(def fun6-res 1.202056902995100878845769442835673986997509412384)
(defn fun6 [[^double x ^double y ^double z]] (m// 1.0 (m/- 1.0 (m/* x y z) -1.0e-10)))

(def fun7-res 0.709615885636934)
(defn fun7 [[^double x ^double y ^double z]] (m/sqrt (m/abs (m/- x y z))))

(def fun8-res 0.891213)
(defn fun8 [[^double x ^double y ^double z]] (m/exp (m/- (m/* x y z))))

(def fun9-res 0.0801856)
(defn fun9 [[^double x ^double y ^double z]] (m// (m/sq x)
                                               (m/+ 5.0 (m/cos (m/+ 1.0 x y z)))))

(def fun11-res 0.5235987755982988)
(defn fun11 [v] (if (m/< (v/magsq v) 1.0) 1.0 0.0))

(def normal-pdf (partial r/pdf (r/distribution :multi-normal {:means [0 0 0]
                                                            :covariances [[1 0 0]
                                                                          [0 0.5 0]
                                                                          [0 0 2]]})))

(defn gaussian [[^double _phi ^double r]] (m/* r (m/exp (m/* -0.5 r r))))

;; https://tutorial.math.lamar.edu/classes/calciii/TISphericalCoords.aspx
(defn polar-3d
  [[^double ro ^double _theta ^double phi]]
  (m/* 16.0 ro ro ro (m/sin phi) (m/cos phi)))

(defn- seeded
  "Options with a generator seeded with `seed`."
  ([seed] (seeded seed nil))
  ([seed options] (assoc options :rng (r/rng :jdk seed))))

(defn- ex-data-of [thunk]
  (try (thunk) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))

(t/deftest multivariate-integrals
  (t/testing "integrals on the unit cube"
    (t/are [f res] (m/delta-eq res (sut/vegas f [0 0 0] [1 1 1]
                                              (seeded 1 {:max-iters 8 :warmup 3 :rel 1.0e-5 :abs 1.0e-5
                                                         :nstrats 12 :nevals 8000})) 1.0e-2)
      fun1+ fun1+-res
      fun1 fun1-res
      fun2 fun2-res
      fun3 fun3-res
      fun4 fun4-res
      fun5 fun5-res
      fun6 fun6-res
      fun7 fun7-res
      fun8 fun8-res
      fun9 fun9-res
      fun11 fun11-res))
  (t/testing "improper integrals"
    (t/are [lower upper res] (m/delta-eq res (sut/vegas normal-pdf lower upper
                                                        (seeded 1 {:max-iters 8 :warmup 3 :rel 1.0e-5 :abs 1.0e-5
                                                                   :nstrats 12 :nevals 8000})) 1.0e-2)
      [##-Inf ##-Inf ##-Inf] [##Inf ##Inf ##Inf] 1.0
      [##-Inf ##-Inf ##-Inf] [0 0 0] 0.125
      [0 0 0] [##Inf ##Inf ##Inf] 0.125
      [##-Inf 0 ##-Inf] [0 ##Inf ##Inf] 0.25)
    (t/is (m/delta-eq m/TWO_PI (sut/vegas gaussian [0.0 0.0] [m/TWO_PI ##Inf]
                                          (seeded 1 {:max-iters 20 :warmup 5 :random-sequence :sobol})) 1.0e-2))
    (t/is (m/delta-eq (m/* 4.0 m/PI)
                      (sut/vegas polar-3d [0 0 0] [1.0 m/TWO_PI m/HALF_PI]
                                 (seeded 1 {:max-iters 8 :warmup 3 :nstrats 12 :nevals 8000})) 1.0e-2))))

;; calibration: for the exact integral `I`, z = (result - I) / sd is a standard normal variable

(def ^:private exact-peak4 (m/pow (m/* (m// (m/sqrt m/PI) 10.0) (special/erf 5.0)) 4.0))

(def ^:private calibration-cases
  {:peak4 [(fn [^doubles x] (m/exp (m/- (m/* 100.0 (m/+ (m/sq (m/- (aget x 0) 0.5)) (m/sq (m/- (aget x 1) 0.5))
                                                        (m/sq (m/- (aget x 2) 0.5)) (m/sq (m/- (aget x 3) 0.5)))))))
           (repeat 4 0.0) (repeat 4 1.0) exact-peak4 {:nevals 4000 :max-iters 6}]
   :prod3 [(fn [^doubles x] (m/* 8.0 (aget x 0) (aget x 1) (aget x 2)))
           (repeat 3 0.0) (repeat 3 1.0) 1.0 {:nevals 4000 :max-iters 6}]
   :prod3-warmup [(fn [^doubles x] (m/* 8.0 (aget x 0) (aget x 1) (aget x 2)))
                  (repeat 3 0.0) (repeat 3 1.0) 1.0 {:nevals 4000 :max-iters 4 :warmup 3}]
   :prod3-classic [(fn [^doubles x] (m/* 8.0 (aget x 0) (aget x 1) (aget x 2)))
                   (repeat 3 0.0) (repeat 3 1.0) 1.0 {:nevals 4000 :max-iters 4 :warmup 3 :nstrats 1}]
   :ball3 [(fn [^doubles x] (if (m/< (m/+ (m/sq (aget x 0)) (m/sq (aget x 1)) (m/sq (aget x 2))) 1.0) 1.0 0.0))
           (repeat 3 0.0) (repeat 3 1.0) (m// m/PI 6.0) {:nevals 4000 :max-iters 6}]
   :square1d [(fn [^doubles x] (m/sq (aget x 0))) [0.0] [1.0] m/THIRD {:nevals 2000 :max-iters 5}]
   :gauss-plane [(fn [x] (m/exp (m/- (m/+ (m/sq (x 0)) (m/sq (x 1)))))) [##-Inf ##-Inf] [##Inf ##Inf] m/PI
                 {:nevals 4000 :max-iters 6}]
   :reversed [(fn [^doubles x] (m/* 8.0 (aget x 0) (aget x 1) (aget x 2)))
              [1.0 0.0 1.0] [0.0 1.0 0.0] 1.0 {:nevals 4000 :max-iters 6}]})

(defn- mean [xs] (m// (reduce + xs) (count xs)))
(defn- sd [xs] (let [mu (mean xs)] (m/sqrt (m// (reduce + (map #(m/sq (m/- % mu)) xs)) (m/dec (count xs))))))

(t/deftest calibration
  ;; measured over seeds 1..30: sdz 0.81-1.25, mean z -0.54..0.16, coverage of |z|<2 0.80-1.0, mean Q 0.42-0.54,
  ;; and over 150 seeds with nevals 20000 sdz 0.92-1.02 on peak4, prod3, ball3 (Lepage's vegas: 0.98-1.08)
  (doseq [[id [f lower upper exact options]] calibration-cases]
    (let [runs (mapv #(sut/vegas f lower upper (seeded % (merge {:rel 0.0 :abs 0.0 :info? true} options))) (range 1 31))
          zs (mapv #(m// (m/- (:result %) exact) (:sd %)) runs)]
      (t/testing (name id)
        (t/is (m/< 0.6 (sd zs) 1.6) "the standard deviation of z is about 1")
        (t/is (m/< (m/abs (mean zs)) 0.9) "no bias")
        (t/is (m/>= (m// (count (filter #(m/< (m/abs %) 2.0) zs)) 30.0) 0.75) "about 95% of the results within 2 sd")
        (t/is (m/< 0.25 (mean (map :Q runs)) 0.75) "Q is spread over (0, 1)")
        (t/is (every? #(m/== (:dof %) (m/dec (:iterations %))) runs))))))

(t/deftest scale-invariance
  ;; the standard deviation of the result scales with the integrand (an absolute floor of the variance did not allow it)
  (let [f (fn [^doubles x] (m/* 8.0 (aget x 0) (aget x 1) (aget x 2)))
        run (fn [scale] (sut/vegas (fn [x] (m/* scale (f x))) [0 0 0] [1 1 1] (seeded 5 {:nevals 4000 :max-iters 4 :rel 0.0 :abs 0.0 :info? true})))
        base (run 1.0)]
    (doseq [scale [1.0e-9 1.0e-6 1.0e6 1.0e9]]
      (let [res (run scale)]
        (t/is (m/delta-eq (m/* scale (:result base)) (:result res) 0.0 1.0e-12) (str scale))
        (t/is (m/delta-eq (m/* scale (:sd base)) (:sd res) 0.0 1.0e-9) (str scale))))))

(t/deftest pinned-results
  ;; one seeded run per sequence, to notice any change of the algorithm (the tolerance covers 1-4 ulp differences of Math functions)
  (doseq [[sequence expected] {:uniform {:result 0.666566619009684 :sd 9.372031390058366E-4}
                               :r2 {:result 0.6654201197982242 :sd 9.129793848006404E-4}
                               :sobol {:result 0.6645637960388204 :sd 9.12797234352088E-4}
                               :halton {:result 0.6651706487303184 :sd 8.536889357120502E-4}}]
    (let [res (sut/vegas fun1 [0 0 0] [1 1 1] (seeded 1 {:max-iters 3 :nevals 2000 :info? true :random-sequence sequence}))]
      (t/testing (str sequence)
        (t/is (m/delta-eq (:result expected) (:result res) 0.0 1.0e-9))
        (t/is (m/delta-eq (:sd expected) (:sd res) 0.0 1.0e-9))
        (t/is (m/delta-eq fun1-res (:result res) (m/* 5.0 (:sd res))) "the exact value is within 5 sd")))))

(t/deftest rng
  (let [run (fn [rng-fn extra] (sut/vegas fun1 [0 0 0] [1 1 1]
                                          (merge {:max-iters 3 :nevals 400 :info? true :rng (rng-fn)} extra)))
        results (fn [rng-fn extra] (select-keys (run rng-fn extra) [:result :sd :evaluations]))]
    (doseq [extra [{} {:random-sequence :r2} {:random-sequence :sobol} {:random-sequence :halton :jitter 0.5}]]
      (t/testing (str extra)
        (t/is (= (results #(r/rng :jdk 1) extra) (results #(r/rng :jdk 1) extra)) "equal seeds give equal results")
        (t/is (not= (results #(r/rng :jdk 1) extra) (results #(r/rng :jdk 2) extra)) "different seeds differ")))
    (t/testing "a missing and a nil rng create a new generator"
      (t/is (number? (:result (sut/vegas fun1 [0 0 0] [1 1 1] {:max-iters 2 :nevals 200 :info? true}))))
      (t/is (number? (:result (sut/vegas fun1 [0 0 0] [1 1 1] {:max-iters 2 :nevals 200 :info? true :rng nil
                                                                 :random-sequence :sobol})))))
    (t/testing "a value which is not a generator throws ex-info"
      (doseq [bad [5 :jdk "x" (r/distribution :normal)]
              extra [{} {:random-sequence :r2}]]
        (t/is (= {:rng bad} (ex-data-of #(sut/vegas fun1 [0 0 0] [1 1 1] (merge {:max-iters 1 :nevals 200 :rng bad} extra)))))))
    (t/testing "the shared generator is not used"
      (r/set-seed! 1)
      (let [expected (r/drand)]
        (r/set-seed! 1)
        (sut/vegas fun1 [0 0 0] [1 1 1] {:max-iters 2 :nevals 200 :rng (r/rng :jdk 3)})
        (sut/vegas fun1 [0 0 0] [1 1 1] {:max-iters 2 :nevals 200 :random-sequence :sobol :rng (r/rng :jdk 3)})
        (sut/vegas fun1 [0 0 0] [1 1 1] {:max-iters 2 :nevals 200})
        (t/is (= expected (r/drand)))))))

(t/deftest degenerate-integrands
  (let [opts (fn [& kv] (merge (seeded 42 {:info? true :max-iters 3 :nevals 1000}) (apply hash-map kv)))]
    (t/testing "a zero integrand"
      (let [res (sut/vegas (fn [_] 0.0) [0 0] [1 1] (opts))]
        (t/is (= 0.0 (:result res)))
        (t/is (m/< (:sd res) 1.0e-100))))
    (t/testing "constant integrands"
      (doseq [c [1.0 -3.5 1.0e6 1.0e-9]]
        (let [res (sut/vegas (fn [_] c) [0 0] [2 3] (opts))]
          (t/is (m/delta-eq (m/* 6.0 c) (:result res) 0.0 1.0e-12) (str c))
          (t/is (m/< (:sd res) (m/* 1.0e-6 (m/abs (:result res)))) (str c)))))
    (t/testing "an indicator which is zero on a half of the domain"
      (let [res (sut/vegas (fn [^doubles x] (if (m/< (aget x 0) 0.5) 1.0 0.0)) [0 0] [1 1] (opts :max-iters 5))]
        (t/is (m/delta-eq 0.5 (:result res) (m/* 5.0 (:sd res))))))
    (t/testing "an indicator of a tiny region (most of the hypercubes see zeros only)"
      (let [res (sut/vegas (fn [^doubles x] (if (and (m/< (aget x 0) 0.01) (m/< (aget x 1) 0.01)) 1.0 0.0))
                           [0 0] [1 1] (opts :max-iters 5 :nevals 5000))]
        (t/is (m/valid-double? (:result res)))
        (t/is (m/valid-double? (:sd res)))))
    (t/testing "an integrand which is nonzero in a single bin"
      (let [res (sut/vegas (fn [^doubles x] (if (m/< (aget x 0) 0.001) 1.0 0.0)) [0] [1] (opts :max-iters 4 :nevals 5000))]
        (t/is (m/valid-double? (:result res)))))
    (t/testing "NaN and infinite values throw ex-info with the point and the value"
      (doseq [bad [##NaN ##Inf ##-Inf]]
        (let [data (ex-data-of #(sut/vegas (fn [_] bad) [0 0] [1 1] (opts)))]
          (t/is (= 2 (count (:x data))))
          (t/is (= (str bad) (str (:value data))))))
      (let [data (ex-data-of #(sut/vegas (fn [^doubles x] (if (m/< (aget x 0) 0.5) 1.0 ##Inf)) [0 0] [1 1] (opts)))]
        (t/is (m/>= ^double (first (:x data)) 0.5) "the point is the one where the integrand failed"))
      (let [data (ex-data-of #(sut/vegas (fn [_] ##NaN) [##-Inf 0] [##Inf 1] (opts)))]
        (t/is (every? m/valid-double? (:x data)) "the point is given in the variables of the integrand")))))

(t/deftest bounds
  (let [f (fn [^doubles x] (m/* 8.0 (aget x 0) (aget x 1) (aget x 2)))
        opts (fn [] (seeded 7 {:nevals 4000 :max-iters 4 :info? true}))
        straight (sut/vegas f [0 0 0] [1 1 1] (opts))]
    (t/testing "orientation"
      (let [one (sut/vegas f [1 0 0] [0 1 1] (opts))
            two (sut/vegas f [1 1 0] [0 0 1] (opts))
            three (sut/vegas f [1 1 1] [0 0 0] (opts))]
        (t/is (= (m/- (:result straight)) (:result one)))
        (t/is (= (:result straight) (:result two)))
        (t/is (= (m/- (:result straight)) (:result three)))
        (t/is (= (:sd straight) (:sd one) (:sd two) (:sd three)) "the standard deviation is positive")
        (t/is (= (m/- (:result straight)) (sut/vegas f [1 1 1] [0 0 0] (dissoc (opts) :info?)))))
      (t/testing "also with infinite bounds"
        (let [g (fn [x] (m/exp (m/- (m/+ (m/sq (x 0)) (m/sq (x 1))))))
              base (sut/vegas g [##-Inf ##-Inf] [##Inf ##Inf] (opts))
              reversed (sut/vegas g [##Inf ##-Inf] [##-Inf ##Inf] (opts))]
          (t/is (m/delta-eq (m/- (:result base)) (:result reversed) 0.0 1.0e-12)))))
    (t/testing "equal bounds in any dimension give zero"
      (doseq [[lower upper] [[[1 1] [1 1]] [[0 1] [1 1]] [[0 0 0] [1 0 1]] [[2] [2]] [[##Inf 0] [##Inf 1]]]]
        (t/is (= 0.0 (sut/vegas f lower upper (dissoc (opts) :info?))))
        (t/is (= {:result 0.0 :sd 0.0 :iterations 0 :evaluations 0}
                 (select-keys (sut/vegas f lower upper (opts)) [:result :sd :iterations :evaluations])))))
    (t/testing "empty and mismatched bounds throw ex-info"
      (doseq [[lower upper] [[[] []] [[0 0] [1]] [[0] [1 1]] [nil nil]]]
        (t/is (some? (ex-data-of #(sut/vegas f lower upper (opts)))))))))

(t/deftest dimensions
  (t/testing "one dimension"
    (doseq [rs [:uniform :r2 :sobol :halton]]
      (let [res (sut/vegas (fn [^doubles x] (m/sq (aget x 0))) [0.0] [1.0] (seeded 3 {:nevals 2000 :max-iters 4 :info? true :random-sequence rs}))]
        (t/is (m/delta-eq m/THIRD (:result res) (m/* 5.0 (:sd res))) (str rs))))
    (t/is (m/delta-eq m/SQRTPI (sut/vegas (fn [x] (m/exp (m/- (m/sq (x 0))))) [##-Inf] [##Inf] (seeded 3 {:nevals 4000 :max-iters 5})) 1.0e-2))
    (t/is (m/delta-eq 8.0 (sut/vegas (fn [^doubles x] (aget x 0)) [0.0] [4.0] (seeded 3 {:nevals 2000 :max-iters 4})) 1.0e-2))))

(t/deftest arguments-of-the-integrand
  (let [seen (atom #{})
        remember (fn [f] (fn [x] (swap! seen conj (type x)) (f x)))]
    (sut/vegas (remember (fn [_] 1.0)) [0 0] [1 1] (seeded 1 {:nevals 100 :max-iters 1}))
    (t/is (= #{(type (double-array 0))} @seen) "a double array for finite bounds")
    (reset! seen #{})
    (sut/vegas (remember (fn [_] 1.0)) [0 ##-Inf] [1 1] (seeded 1 {:nevals 100 :max-iters 1}))
    (t/is (= #{clojure.lang.PersistentVector} @seen) "a vector for infinite bounds")))

(t/deftest stopping-and-evaluations
  (let [f (fn [^doubles x] (m/+ 1.0e6 (m/* 100.0 (aget x 0))))
        base {:nevals 1000 :info? true}]
    (t/testing "evaluations are the calls of the integrand, also with warmup"
      (doseq [extra [{:max-iters 3} {:max-iters 3 :warmup 2} {:max-iters 2 :random-sequence :sobol :nstrats 1}]]
        (let [calls (atom 0)
              res (sut/vegas (fn [x] (swap! calls inc) (f x)) [0 0] [1 1] (seeded 1 (merge base extra {:rel 0.0 :abs 0.0})))]
          (t/is (= @calls (:evaluations res)) (str extra))
          (t/is (= (:max-iters extra) (:iterations res)) (str extra)))))
    (t/testing "the iterations stop when the sd is below the larger of the relative and the absolute error"
      (let [res (sut/vegas (fn [_] 1.0e6) [0 0] [1 1] (seeded 1 (assoc base :rel 1.0e-3 :abs 1.0e-3 :max-iters 10)))]
        (t/is (= 1 (:iterations res)) "the relative error is met, the absolute one is not"))
      (let [res (sut/vegas f [0 0] [1 1] (seeded 1 (assoc base :rel 0.0 :abs 1.0e9 :max-iters 10)))]
        (t/is (= 1 (:iterations res)) "the absolute error is met, the relative one is not"))
      (let [res (sut/vegas f [0 0] [1 1] (seeded 1 (assoc base :rel 0.0 :abs 0.0 :max-iters 4)))]
        (t/is (= 4 (:iterations res)) "neither is met")))
    (t/testing "after a single iteration there are no degrees of freedom"
      (let [res (sut/vegas f [0 0] [1 1] (seeded 1 (assoc base :max-iters 1)))]
        (t/is (zero? (:dof res)))
        (t/is (m/nan? (:chi2-avg res)))
        (t/is (m/nan? (:Q res)))))
    (t/testing "recorded data"
      (let [res (sut/vegas f [0 0] [1 1] (seeded 1 (assoc base :max-iters 2 :record-data? true :rel 0.0 :abs 0.0)))]
        (t/is (= 2 (count (:data res))))
        (t/is (contains? res :hcubes))))))

(t/deftest option-validation
  (let [f (fn [_] 1.0)
        data (fn [options] (ex-data-of #(sut/vegas f [0 0] [1 1] (seeded 1 options))))]
    (doseq [[option bad] {:nevals [0 1 -5 2.5 "a" ##NaN]
                          :max-iters [0 -1 1.5 "a"]
                          :warmup [-1 1.5 "a"]
                          :nintervals [0 -1 1.5]
                          :rel [-1 -1.0e-9 "a" ##NaN]
                          :abs [-1 "a" ##NaN]
                          :alpha [-0.5 "a"]
                          :beta [-0.5 "a"]
                          :random-sequence [:bogus "sobol" 1]}
            value bad]
      (t/testing (str option " " (pr-str value))
        (let [d (data {option value})]
          (t/is (= option (:option d)))
          (t/is (= (str value) (str (:value d))))
          (t/is (string? (:expected d))))))
    (t/testing "valid values and unknown options are accepted"
      (doseq [options [{:nevals 2 :max-iters 1} {:warmup 0} {:rel 0 :abs 0} {:alpha 0} {:beta 0}
                       {:niters 20 :res 1.0e-5 :bogus-option true} {:nevals nil :alpha nil}]]
        (t/is (nil? (data (merge {:nevals 100 :max-iters 1} options))) (pr-str options))))))

(t/deftest default-alpha
  ;; measured (60 seeds, warmup 5): with one stratum per dimension alpha 1.5 has 1.5-3x smaller errors than 0.5 on
  ;; peak4, prod3 and the Gaussian; with stratified sampling the two are close (0.5 better on 3 of 5 integrands)
  (let [f (fn [^doubles x] (m/* 8.0 (aget x 0) (aget x 1) (aget x 2)))
        result (fn [options] (:result (sut/vegas f [0 0 0] [1 1 1] (seeded 11 (merge {:nevals 1000 :max-iters 2 :rel 0.0 :abs 0.0 :info? true} options)))))]
    (t/testing "one stratum per dimension: 1.5"
      (doseq [nstrats [1 [1 1 1] [1]]]
        (t/is (= (result {:nstrats 1 :alpha 1.5}) (result {:nstrats nstrats})) (pr-str nstrats))
        (t/is (not= (result {:nstrats 1 :alpha 0.5}) (result {:nstrats nstrats})) (pr-str nstrats))))
    (t/testing "otherwise 0.5"
      (doseq [nstrats [[1 2 1] 2 [3 1] nil]]
        (t/is (= (result {:nstrats nstrats :alpha 0.5}) (result {:nstrats nstrats})) (pr-str nstrats))
        (t/is (not= (result {:nstrats nstrats :alpha 1.5}) (result {:nstrats nstrats})) (pr-str nstrats))))
    (t/testing "an explicit alpha wins"
      (t/is (= (result {:nstrats 1 :alpha 0.5}) (result {:nstrats [1 1 1] :alpha 0.5}))))))
