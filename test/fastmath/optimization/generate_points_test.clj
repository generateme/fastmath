(ns fastmath.optimization.generate-points-test
  "Scan of the search domain (`generate-points`) with functions which are not finite."
  (:require [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.optimization :as o]
            [fastmath.random :as r]
            [fastmath.vector :as v]))


(def ^:private generate-points @#'o/generate-points)

(def ^:private timeout-ms 10000)

(defn- psum [p] (reduce + 0.0 p))

(defn- call-with-timeout
  "Result of `thunk` (a thrown exception is returned) or `::timeout` when it runs longer than `timeout-ms`. A hang of the code under test fails the test instead of the whole suite."
  [thunk]
  (let [fut (future (try (thunk) (catch Throwable e e)))
        res (deref fut timeout-ms ::timeout)]
    (when (= res ::timeout) (future-cancel fut))
    res))

(defn- invalid-input?
  [res attempts N]
  (and (instance? clojure.lang.ExceptionInfo res)
       (let [d (ex-data res)]
         (and (= :invalid-input (:type d))
              (= attempts (:attempts d))
              (= N (:N d))))))

(defn- close?
  [^double a ^double b ^double tol]
  (m/<= (m/abs (m/- a b)) tol))

(defn- old-generate-points
  "Version of `generate-points` before the attempt bound."
  [f bounds goal N jitter rng]
  (let [dim (count bounds)
        lo (map first bounds)
        hi (map second bounds)
        N (long (m/max (m/+ 4.5 (m/* dim (m/log2 dim))) (long N)))
        gen (r/jittered-sequence-generator (if (m/< dim 15) :r2 :sobol) dim jitter rng)]
    (->> (if (m/one? dim) (map vector gen) gen)
         (map (fn [v] (let [p (v/einterpolate lo hi v)] [(f p) p])))
         (filter (comp m/valid-double? first))
         (take N)
         (sort-by first (if (= goal :minimize) m/< m/>))
         (map second))))

(deftest all-not-finite-throws
  (testing "generate-points: NaN, +Inf and -Inf everywhere, 1000 N attempts"
    (doseq [[label v] [["NaN" ##NaN] ["+Inf" ##Inf] ["-Inf" ##-Inf]]]
      (let [calls (atom 0)
            f (fn [_] (swap! calls inc) v)
            res (call-with-timeout #(generate-points f [[-1.0 1.0] [-1.0 1.0]] :minimize 100 0.25 (r/rng :isaac 1)))]
        (is (invalid-input? res 100000 100) label)
        (is (= 100000 @calls) label))))
  (testing "N is raised to 4.5 + d log2 d before the bound is computed"
    (is (invalid-input? (call-with-timeout #(generate-points (fn [_] ##NaN) [[0.0 1.0]] :minimize 1 0.25 (r/rng :isaac 1)))
                        4000 4))
    (is (invalid-input? (call-with-timeout #(generate-points (fn [_] ##NaN) (repeat 5 [0.0 1.0]) :maximize 1 0.25 (r/rng :isaac 1)))
                        16000 16)))
  (testing "ex-info has the message about the scanned points"
    (let [res (call-with-timeout #(generate-points (fn [_] ##NaN) [[0.0 1.0]] :minimize 10 0.25 (r/rng :isaac 1)))]
      (is (instance? clojure.lang.ExceptionInfo res))
      (is (re-find #"no finite value" (ex-message res))))))

(deftest all-not-finite-scan-and-optimize
  (testing "scan-and-minimize / maximize / optimize throw, they do not hang"
    (is (invalid-input? (call-with-timeout #(o/scan-and-minimize :lbfgsb (fn [_] ##NaN)
                                                                 {:bounds [[-1 1] [-1 1]] :rng (r/rng :isaac 1)}))
                        100000 100))
    (is (invalid-input? (call-with-timeout #(o/scan-and-maximize :brent (fn [^double _] ##Inf)
                                                                 {:bounds [-1 1] :N 10 :rng (r/rng :isaac 1)}))
                        10000 10))
    (is (invalid-input? (call-with-timeout #(o/scan-and-optimize :cmaes (fn [_] ##-Inf)
                                                                 {:bounds [[-1 1] [-1 1]] :N 300 :goal :minimize
                                                                  :vector-arg? true :rng (r/rng :isaac 1)}))
                        300000 300))))

(deftest all-not-finite-bayesian-optimization
  (testing "a NaN objective gives NaN Gaussian process and utility function, the scan throws"
    (let [res (call-with-timeout #(doall (take 2 (o/bayesian-optimization (fn [_] ##NaN)
                                                                          {:bounds [[-1 1] [-1 1]] :warm-up 50
                                                                           :rng (r/rng :isaac 1)}))))]
      (is (invalid-input? res 50000 50)))
    (let [res (call-with-timeout #(doall (take 2 (o/bayesian-optimization (fn [^double _] ##NaN)
                                                                          {:bounds [[-1 1]] :warm-up 20
                                                                           :vector-arg? false
                                                                           :rng (r/rng :isaac 1)}))))]
      (is (invalid-input? res 20000 20)))))

(deftest part-of-the-domain-not-finite
  (testing "NaN on half of the domain: N finite points, none is NaN"
    (let [f (fn [[^double x ^double y]] (if (m/< x 0.0) ##NaN (m/+ (m/sq (m/- x 1.0)) (m/sq (m/- y 0.5)))))
          pts (generate-points f [[-3.0 3.0] [-3.0 3.0]] :minimize 100 0.25 (r/rng :isaac 2))]
      (is (= 100 (count pts)))
      (is (every? (fn [[^double x _]] (m/>= x 0.0)) pts))
      (is (apply <= (map f pts)))))
  (testing "scan-and-minimize finds the minimum"
    (let [f (fn [[^double x ^double y]] (if (m/< x 0.0) ##NaN (m/+ (m/sq (m/- x 1.0)) (m/sq (m/- y 0.5)))))
          [[^double x ^double y] ^double v] (o/scan-and-minimize :lbfgsb f {:bounds [[-3 3] [-3 3]] :N 50 :rng (r/rng :isaac 3)})]
      ;; the optimum is interior, the optimizer converges to the gradient tolerance (observed error below 1e-6)
      (is (close? x 1.0 1.0e-5))
      (is (close? y 0.5 1.0e-5))
      (is (close? v 0.0 1.0e-9))))
  (testing "the best points come first for both goals"
    (let [f (fn [^double x] (m/sin (m/* 3.0 x)))
          pmin (map first (generate-points (fn [[x]] (f x)) [[0.0 3.0]] :minimize 20 0.25 (r/rng :isaac 4)))
          pmax (map first (generate-points (fn [[x]] (f x)) [[0.0 3.0]] :maximize 20 0.25 (r/rng :isaac 4)))]
      (is (= 20 (count pmin) (count pmax)))
      (is (apply <= (map f pmin)))
      (is (apply >= (map f pmax))))))

(deftest exactly-n-and-fewer-valid-points
  (testing "exactly N finite values: N points, the stream is not consumed further"
    (let [calls (atom 0)
          f (fn [_] (if (m/<= (long (swap! calls inc)) 100) 1.0 ##NaN))
          pts (generate-points f [[-1.0 1.0] [-1.0 1.0]] :minimize 100 0.25 (r/rng :isaac 5))]
      (is (= 100 (count pts)))
      (is (= 100 @calls))))
  (testing "N - 1 finite values: the found ones are returned after all 1000 N attempts"
    (let [calls (atom 0)
          f (fn [_] (if (m/<= (long (swap! calls inc)) 99) 1.0 ##NaN))
          pts (generate-points f [[-1.0 1.0] [-1.0 1.0]] :minimize 100 0.25 (r/rng :isaac 5))]
      (is (= 99 (count pts)))
      (is (= 100000 @calls))))
  (testing "one finite value only"
    (let [calls (atom 0)
          f (fn [_] (if (== 1 (long (swap! calls inc))) 1.0 ##Inf))
          pts (generate-points f [[-1.0 1.0]] :maximize 1 0.25 (r/rng :isaac 5))]
      (is (= 1 (count pts)))
      (is (= 4000 @calls))))
  (testing "a small finite region: fewer than N points, scan-and-minimize still runs"
    (let [f (fn [[^double x ^double y]] (if (m/> x 2.999) (m/+ (m/sq (m/- x 3.0)) (m/sq y)) ##NaN))
          pts (generate-points f [[-3.0 3.0] [-3.0 3.0]] :minimize 100 0.25 (r/rng :isaac 6))
          n (count pts)]
      (is (< 0 n 100))
      (is (every? (fn [[^double x _]] (m/> x 2.999)) pts))
      (let [res (o/scan-and-minimize :lbfgsb f {:bounds [[-3 3] [-3 3]] :N 100 :n 3 :rng (r/rng :isaac 6)})]
        (is (some? res))))))

(deftest dimensions-and-small-n
  (testing "N = 1 is raised to 4.5 + d log2 d, 1-D and 5-D"
    (doseq [[bounds expected] [[[[0.0 1.0]] 4]
                               [[[0.0 1.0] [0.0 2.0]] 6]
                               [(vec (repeat 5 [0.0 1.0])) 16]]]
      (let [pts (generate-points (fn [p] (psum p)) bounds :minimize 1 0.25 (r/rng :isaac 7))]
        (is (= expected (count pts)) (str (count bounds) "-D"))
        (is (every? #(= (count bounds) (count %)) pts)))))
  (testing "5-D with half of the domain not finite"
    (let [f (fn [p] (if (m/< (double (first p)) 0.5) ##NaN (psum p)))
          pts (generate-points f (vec (repeat 5 [0.0 1.0])) :maximize 1 0.25 (r/rng :isaac 8))]
      (is (= 16 (count pts)))
      (is (every? #(m/>= (double (first %)) 0.5) pts))))
  (testing "N = 1 in scan-and-minimize, 1-D"
    (let [[^double x ^double v] (o/scan-and-minimize :brent (fn [^double x] (m/sq (m/- x 0.3)))
                                                     {:bounds [-1 1] :N 1 :rng (r/rng :isaac 9)})]
      (is (close? x 0.3 1.0e-6))
      (is (close? v 0.0 1.0e-10)))))

(deftest same-as-before-for-finite-functions
  (testing "the same points as the version without the bound and the same state of the generator afterwards"
    (let [f (fn [[^double x ^double y]] (m/+ (m/sq (m/- x 1.0)) (m/* 2.0 (m/sq (m/- y 0.5))) (m/* 0.3 (m/sin (m/* 5.0 x)))))
          half (fn [[^double x ^double y]] (if (m/< x 0.0) ##NaN (f [x y])))]
      (doseq [[fun bounds N goal] [[f [[-3.0 3.0] [-3.0 3.0]] 50 :minimize]
                                   [f [[-3.0 3.0] [-3.0 3.0]] 30 :maximize]
                                   [half [[-3.0 3.0] [-3.0 3.0]] 40 :minimize]
                                   [(fn [p] (f (concat p p))) [[-3.0 3.0]] 7 :minimize]
                                   [(fn [p] (psum p)) (vec (repeat 16 [0.0 1.0])) 20 :maximize]]]
        (let [rng-old (r/rng :isaac 42)
              rng-new (r/rng :isaac 42)
              old (doall (old-generate-points fun bounds goal N 0.25 rng-old))
              new (doall (generate-points fun bounds goal N 0.25 rng-new))]
          (is (= old new))
          (is (= (r/irandom rng-old 1000000000) (r/irandom rng-new 1000000000)))))))
  (testing "scan-and-minimize with seeded generators, reference values taken before the change"
    (let [f (fn [[^double x ^double y]] (m/+ (m/sq (m/- x 1.0)) (m/* 2.0 (m/sq (m/- y 0.5))) (m/* 0.3 (m/sin (m/* 5.0 x)))))
          [[^double x ^double y] ^double v] (o/scan-and-minimize :lbfgsb f {:bounds [[-3 3] [-3 3]] :N 50 :rng (r/rng :isaac 42)})
          rng (r/rng :isaac 9)
          _ (o/scan-and-minimize :lbfgsb f {:bounds [[-3 3] [-3 3]] :N 50 :rng rng})]
      (is (close? x 0.954593582936739 1.0e-9))
      (is (close? y 0.5 1.0e-8))
      (is (close? v -0.2973879545203336 1.0e-12))
      ;; the generator state after the call is the same as before the change
      (is (= 756402 (r/irandom rng 1000000))))))

