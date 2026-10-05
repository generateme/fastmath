(ns fastmath.optimization-test
  (:require [fastmath.optimization :as sut]
            [fastmath.optimization.problems :as p]
            [clojure.test :as t]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.core :as m])
  (:import [org.apache.commons.math3.exception TooManyEvaluationsException TooManyIterationsException]
           [org.apache.commons.math3.optim.linear NoFeasibleSolutionException UnboundedSolutionException]))

(t/deftest linear-programming
  (t/are [res target constrains opts] (let [[p v] res
                                            [rp rv] (sut/linear-optimization target constrains opts)]
                                        (and (v/delta-eq (vec rp) p)
                                             (m/delta-eq v rv)))
    ;; scipy
    [[10 -3] -22] [-1 4 0] [[-3 1] :<= 6
                            [-1 -2] :>= -4
                            [0 1] :>= -3] nil
    
    ;; http://people.brunel.ac.uk/~mastjjb/jeb/or/morelp.html
    [[45.0 6.25] 1.25] [1 1 -50] [[50 24] :leq 2400
                                  [30 33] :leq 2100
                                  [1 0] :geq 45
                                  [0 1] :geq 5] {:goal :maximize}
    
    ;; transportation problem https://dewwool.com/linear-programming-examples/
    [[0.0 0.0 10.0 30.0 20.0 20.0] 290.0] [5 4 6 3 2 5 0] '[[1 1 1 0 0 0] <= 50
                                                            [0 0 0 1 1 1] <= 70
                                                            [1 0 0 1 0 0] >= 30
                                                            [0 1 0 0 1 0] >= 20
                                                            [0 0 1 0 0 1] >= 30] {:non-negative? true}))

;; sudoku solver

(defn sudoku-idx ^long [^long v ^long r ^long c] (+ (* 81 v) (* 9 r) c))

(defn sudoku-solver
  [input]
  (let [zeros (vec (repeat (* 9 9 9) 0))
        c1 (for [r (range 9)
                 c (range 9)]
             (reduce (fn [z v]
                       (assoc z (sudoku-idx v r c) 1)) zeros (range 9)))
        c2 (for [v (range 9)
                 c (range 9)]
             (reduce (fn [z r]
                       (assoc z (sudoku-idx v r c) 1)) zeros (range 9)))
        c3 (for [v (range 9)
                 r (range 9)]
             (reduce (fn [z c]
                       (assoc z (sudoku-idx v r c) 1)) zeros (range 9)))
        c4 (for [v (range 9)
                 p (range 3)
                 q (range 3)
                 :let [rc (for [r (range (* 3 p) (* 3 (inc p)))
                                c (range (* 3 q) (* 3 (inc q)))]
                            [r c])]]
             (reduce (fn [z [r c]]
                       (assoc z (sudoku-idx v r c) 1)) zeros rc))
        c5 (for [[v r c] input]
             (assoc zeros (sudoku-idx (dec v) r c) 1))
        all (reduce (fn [buff c]
                      (conj buff c :eq 1)) [] (mapcat identity [c1 c2 c3 c4 c5]))
        res (vec (first (sut/linear-optimization (conj zeros 0) all {:non-negative? true})))
        solution (vec (repeat 81 0))]
    (->> (for [r (range 9)
               c (range 9)
               v (range 9)
               :let [id (sudoku-idx v r c)]
               :when (m/one? (m/round (res id)))]
           [v r c])
         (reduce (fn [s [v r c]]
                   (assoc s (+ c (* r 9)) (inc v))) solution)
         (partition 9))))

(t/deftest sudoku
  ;; v - value 1-9
  ;; r - row 0-8
  ;; c - column 0-8
  ;;
  ;;                        v r c
  (t/is (= (sudoku-solver [[8 0 2]
                           [4 0 4]
                           [2 0 5]
                           [3 0 8]

                           [6 1 1]
                           [1 1 7]
                           
                           [3 2 1]
                           [7 2 5]

                           [5 3 5]
                           [3 3 6]

                           [4 4 1]
                           [9 4 3]
                           [1 4 5]
                           [7 4 7]

                           [5 5 2]
                           [2 5 3]

                           [3 6 3]
                           [4 6 7]

                           [5 7 1]
                           [6 7 7]

                           [7 8 0]
                           [6 8 3]
                           [1 8 4]
                           [9 8 6]])
           [[9 7 8 1 4 2 6 5 3]
            [2 6 4 5 9 3 7 1 8]
            [5 3 1 8 6 7 4 2 9]
            [1 2 7 4 8 5 3 9 6]
            [8 4 6 9 3 1 5 7 2]
            [3 9 5 2 7 6 1 8 4]
            [6 1 9 3 5 8 2 4 7]
            [4 5 3 7 2 9 8 6 1]
            [7 8 2 6 1 4 9 3 5]])))




;; general optimizers
;;
;; Reference: Himmelblau's function has four minima with value 0:
;; (3, 2), (-2.805118, 3.131313), (-3.779310, -3.283186), (3.584428, -1.848126) (Wikipedia).
;; problem02 (sin x + sin(10x/3) on [2.7, 7.5]) has its global minimum -1.899599 at 5.145735
;; (https://infinity77.net/global_optimization/test_functions_1d.html).
;; The optimizers are checked for sanity only: tolerance 1e-3 on values and points close to a known minimum.

(def ^:private val-tol 1.0e-3)
(def ^:private himmelblau-minima [[3.0 2.0] [-2.805118 3.131313] [-3.779310 -3.283186] [3.584428 -1.848126]])
(def ^:private hb (p/himmelblau-bounds))
(def ^:private multivariate-methods [:lbfgsb :bobyqa :cmaes :sceua :nelder-mead :multidirectional-simplex :powell :gradient :non-linear-gradient])
;; :sceua is stochastic. With its defaults it stops after 7 loops without an improvement, which happened in about 10%
;; of 300 seeded runs on himmelblau (27 runs with values up to 0.63), with :stop-loops 25 in 1 of 300. The tests of the
;; method matrices test the wiring, so they use a seeded generator and :stop-loops 25 for it.
(defn- seeded [method opts] (cond-> opts (= :sceua method) (assoc :rng (r/rng :jdk 1) :stop-loops 25)))
(def ^:private all-methods (conj multivariate-methods :brent))

(defn- near-minimum? [pt] (boolean (some #(v/delta-eq (vec pt) % 1.0e-3) himmelblau-minima)))
(defn- neg-himmelblau [x] (- (p/himmelblau x)))
(defn- value-of [res] (if (map? res) (:value res) (second res)))
(defn- point-of [res] (if (map? res) (:point res) (first res)))

(defn- himmelblau-gradient [[^double x ^double y]]
  [(m/+ (m/* 4.0 x (m/+ (m/* x x) y -11.0)) (m/* 2.0 (m/+ x (m/* y y) -7.0)))
   (m/+ (m/* 2.0 (m/+ (m/* x x) y -11.0)) (m/* 4.0 y (m/+ x (m/* y y) -7.0)))])

(defn- ex-data-of [thunk]
  (try (thunk) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))

(t/deftest unknown-method
  (doseq [bad [:foo :bfgs nil "lbfgsb" :LBFGSB]
          call [#(sut/minimize bad p/himmelblau {:bounds hb})
                #(sut/maximize bad p/himmelblau {:bounds hb})
                #(sut/optimize bad p/himmelblau {:bounds hb})
                #(sut/minimizer bad p/himmelblau {:bounds hb})
                #(sut/maximizer bad p/himmelblau {:bounds hb})
                #(sut/scan-and-minimize bad p/himmelblau {:bounds hb})
                #(sut/scan-and-maximize bad p/himmelblau {:bounds hb})
                #(sut/scan-and-optimize bad p/himmelblau {:bounds hb})]]
    (let [d (ex-data-of call)]
      (t/is (= bad (:method d)) (pr-str bad))
      (t/is (= (set all-methods) (:allowed d)) (pr-str bad)))))

(t/deftest minimize-maximize-matrix
  (doseq [method multivariate-methods
          vector-arg? [true false nil]
          :let [f (if (false? vector-arg?) (fn [x y] (p/himmelblau [x y])) p/himmelblau)
                nf (if (false? vector-arg?) (fn [x y] (- (p/himmelblau [x y]))) neg-himmelblau)
                opts (seeded method {:bounds hb :vector-arg? vector-arg?})
                label (str method " " vector-arg?)]]
    (let [res (sut/minimize method f opts)]
      (t/is (vector? res) label)
      (t/is (vector? (first res)) label)
      (t/is (near-minimum? (first res)) label)
      (t/is (< (second res) val-tol) label))
    (let [[pt val] (sut/maximize method nf opts)]
      (t/is (near-minimum? pt) label)
      (t/is (<= val 0.0) label)
      (t/is (< (- val) val-tol) label))
    ;; optimize takes the goal from the options, minimize is the default
    (t/is (near-minimum? (point-of (sut/optimize method f opts))) label)
    (t/is (near-minimum? (point-of (sut/optimize method nf (assoc opts :goal :maximize)))) label)
    ;; minimize and maximize override the goal
    (t/is (< (value-of (sut/minimize method f (assoc opts :goal :maximize))) val-tol) label)
    (t/is (< (- (value-of (sut/maximize method nf (assoc opts :goal :minimize)))) val-tol) label)
    ;; stats
    (let [s (sut/minimize method f (assoc opts :stats? true))]
      (t/is (= (case method
                 :lbfgsb #{:point :value :iterations :gradient :status}
                 :sceua #{:point :value :evaluations :iterations :status :complexes}
                 #{:point :value :evaluations :iterations})
               (set (keys s))) label)
      (t/is (near-minimum? (:point s)) label)
      (t/is (integer? (:iterations s)) label))))

(t/deftest invalid-goal
  (doseq [method all-methods
          bad [:min "maximize"]]
    (t/is (= {:goal bad :allowed #{:minimize :maximize}}
             (ex-data-of #(sut/optimize method p/himmelblau {:bounds [[0 1]] :goal bad}))) (str method bad))))

(t/deftest bounds-guards
  (doseq [method all-methods
          [bounds initial part] [[[[5.0 -5.0]] nil "greater"]
                                 [[[0.0 ##NaN]] nil "NaN"]
                                 [[[0.0 1.0] [0.0 1.0]] [0.0] "differs"]]]
    ;; brent has one dimension only and rejects a second pair earlier
    (when-not (and (= :brent method) (= 2 (count bounds)))
      (let [d (ex-data-of #(sut/minimize method p/himmelblau {:bounds bounds :initial initial}))]
        ;; :gradient is an alias of :non-linear-gradient
        (t/is (contains? (set [method :non-linear-gradient]) (:method d)) (str method (pr-str bounds)))
        (t/is (re-find (re-pattern part) (str (:reason d))) (str method (pr-str bounds)))))))

(t/deftest one-dimension
  (let [f (fn [[x]] (p/problem02 x))
        bounds [[2.7 7.5]]]
    (doseq [method (remove #{:bobyqa} multivariate-methods)
            goal [:minimize :maximize]
            :let [f (if (= goal :minimize) f (fn [v] (- (f v))))]]
      (let [[pt val] (sut/optimize method f (seeded method {:bounds bounds :goal goal}))]
        (t/is (vector? pt) (str method))
        (t/is (= 1 (count pt)) (str method))
        (t/is (m/delta-eq 5.1457 (first pt) 1.0e-1) (str method))
        (t/is (< (Math/abs (+ 1.8996 (if (= goal :minimize) val (- val)))) 1.0e-1) (str method))))
    ;; brent: a number, the point is a number
    (let [[pt val] (sut/minimize :brent p/problem02 {:bounds bounds})]
      (t/is (m/delta-eq 5.145735 pt 1.0e-4))
      (t/is (m/delta-eq -1.899599 val 1.0e-5)))
    (t/is (= (sut/minimize :brent p/problem02 {:bounds bounds})
             (sut/minimize :brent p/problem02 {:bounds [2.7 7.5]})))
    (let [[pt] (sut/minimize :brent f {:bounds bounds :vector-arg? true})]
      (t/is (vector? pt)))
    (t/is (= 1 (:n (ex-data-of #(sut/minimize :bobyqa f {:bounds bounds})))))))

(t/deftest limits
  (t/is (thrown? TooManyEvaluationsException (sut/minimize :nelder-mead p/himmelblau {:bounds hb :max-evals 5})))
  (t/is (thrown? TooManyIterationsException (sut/minimize :powell p/himmelblau {:bounds hb :max-iters 1})))
  ;; lbfgsb and cmaes stop at the iteration limit and return the best point
  (t/is (= :max-iterations (:status (sut/minimize :lbfgsb p/rosenbrock {:bounds (p/rosenbrock-bounds 2) :initial [-1.2 1.0] :max-iters 2 :stats? true}))))
  (t/is (vector? (sut/minimize :cmaes p/himmelblau {:bounds hb :max-iters 2}))))

(t/deftest sceua-method
  (t/testing "the evaluation limit throws ex-info, the iteration limit stops the run"
    (t/is (= {:max-evals 5 :evaluations 5} (ex-data-of #(sut/minimize :sceua p/himmelblau {:bounds hb :max-evals 5}))))
    (let [s (sut/minimize :sceua p/himmelblau {:bounds hb :max-iters 1 :stats? true :rng (r/rng :jdk 1)})]
      (t/is (= :max-iterations (:status s)))
      (t/is (= 1 (:iterations s)))))
  (t/testing "bounds are required, finite and of a positive range"
    (t/is (= :sceua (:method (ex-data-of #(sut/minimize :sceua p/himmelblau {})))))
    (doseq [bounds [[[##-Inf 5] [-5 5]] [[-5 ##Inf] [-5 5]] [[1 1] [-5 5]] [[5 -5] [-5 5]]]]
      (t/is (= :sceua (:method (ex-data-of #(sut/minimize :sceua p/himmelblau {:bounds bounds})))) (pr-str bounds))))
  (t/testing "minimizer: the initial point joins the population, nil and a point outside the bounds"
    (let [mz (sut/minimizer :sceua p/himmelblau {:bounds hb :rng (r/rng :jdk 2) :stop-loops 25})]
      (t/is (near-minimum? (first (mz [3.0 2.0]))))
      (t/is (<= (second (mz [3.0 2.0])) 1.0e-12) "the best value is not worse than the value at the initial point")
      (t/is (near-minimum? (first (mz nil))))
      (t/is (= :initial (:option (ex-data-of #(mz [10.0 0.0])))))
      ;; the number of bounds differs from the length of the point
      (t/is (= :sceua (:method (ex-data-of #(mz [1.0])))))))
  (t/testing "one dimension and a flat pair of bounds"
    (let [[pt val] (sut/minimize :sceua p/problem02 {:bounds [2.7 7.5] :vector-arg? false :rng (r/rng :jdk 3) :stop-loops 25})]
      (t/is (m/delta-eq 5.145735 (first pt) 1.0e-2))
      (t/is (m/delta-eq -1.899599 val 1.0e-4))))
  (t/testing "scan-and-minimize starts :sceua from scanned points"
    (let [[pt val] (sut/scan-and-minimize :sceua p/himmelblau {:bounds hb :N 20 :n 2 :rng (r/rng :jdk 4) :stop-loops 25})]
      (t/is (near-minimum? pt))
      (t/is (< val val-tol)))))

;; minimizer and maximizer

(t/deftest minimizer-maximizer
  (let [calls (atom 0)
        counting (fn [x] (swap! calls inc) (p/himmelblau x))
        mz (sut/minimizer :lbfgsb counting {:bounds hb})]
    (t/is (zero? @calls) "creation does not run the function")
    (t/is (fn? mz))
    (t/is (near-minimum? (first (mz [1 1]))))
    (t/is (pos? @calls))
    ;; nil is the default initial point: the middle of the bounds
    (t/is (= (sut/minimize :lbfgsb p/himmelblau {:bounds hb}) (mz nil)))
    (t/is (= (sut/minimize :lbfgsb p/himmelblau {:bounds hb}) (mz)))
    ;; the initial point of the call replaces the :initial option, and gives the same result as minimize
    (let [mz2 (sut/minimizer :lbfgsb p/himmelblau {:bounds hb :initial [-3 -3]})]
      (t/is (= (sut/minimize :lbfgsb p/himmelblau {:bounds hb :initial [2 2]}) (mz2 [2 2])))
      (t/is (= (sut/minimize :lbfgsb p/himmelblau {:bounds hb}) (mz2 nil)))
      (t/is (= (sut/minimize :lbfgsb p/himmelblau {:bounds hb}) (mz2)))
      (t/is (= (mz2) (mz2 nil)))))
  ;; different minima from different initial points
  (let [mz (sut/minimizer :nelder-mead p/himmelblau {:bounds hb})]
    (t/is (v/delta-eq [3.0 2.0] (first (mz [2.5 1.5])) 1.0e-3))
    (t/is (v/delta-eq [-3.77931 -3.283186] (first (mz [-4.0 -3.0])) 1.0e-3)))
  ;; goal is overridden
  (let [mx (sut/maximizer :powell neg-himmelblau {:bounds hb :goal :minimize})]
    (t/is (near-minimum? (first (mx [1 1]))))
    (t/is (< (- (second (mx nil))) val-tol)))
  (let [mz (sut/minimizer :powell p/himmelblau {:bounds hb :goal :maximize})]
    (t/is (near-minimum? (first (mz [1 1])))))
  ;; options are validated when the function is called
  (t/is (near-minimum? (first ((sut/minimizer :lbfgsb p/himmelblau {:bounds hb :goal :foo}) [1 1]))) "an invalid goal is overridden")
  (t/is (= :lbfgsb (:method (ex-data-of #((sut/minimizer :lbfgsb p/himmelblau {:bounds [[5 -5] [0 1]]}) [1 1])))))
  ;; every method, with the initial point and nil
  (doseq [method multivariate-methods]
    (let [mz (sut/minimizer method p/himmelblau (seeded method {:bounds hb}))]
      (t/is (near-minimum? (first (mz [1 1]))) (str method))
      (t/is (near-minimum? (first (mz nil))) (str method))))
  ;; brent: number, one number sequence or nil
  (let [mz (sut/minimizer :brent p/problem02 {:bounds [[2.7 7.5]]})]
    (doseq [initial [5.0 [5.0] nil]]
      (t/is (m/delta-eq 5.145735 (first (mz initial)) 1.0e-4) (pr-str initial)))
    (t/is (m/delta-eq 4.1966 (first ((sut/maximizer :brent p/problem02 {:bounds [[2.7 7.5]]}) 4.0)) 1.0e-3))))

;; gradient

(t/deftest gradient-option
  (let [calls (fn [method goal f g]
                (let [c (atom 0)
                      res (sut/optimize method (fn [x] (swap! c inc) (f x))
                                        (cond-> {:bounds hb :goal goal} g (assoc :gradient g)))]
                  [@c res]))]
    (doseq [method [:lbfgsb :gradient :non-linear-gradient]]
      (let [[c-num _res-num] (calls method :minimize p/himmelblau nil)
            [c-grad res-grad] (calls method :minimize p/himmelblau himmelblau-gradient)]
        (t/is (near-minimum? (first res-grad)) (str method))
        (t/is (< (second res-grad) val-tol) (str method))
        (t/is (< c-grad c-num) (str method " analytic gradient needs fewer evaluations"))
        ;; maximization: the gradient of f itself
        (let [[_ res-max] (calls method :maximize neg-himmelblau (fn [x] (mapv - (himmelblau-gradient x))))]
          (t/is (near-minimum? (first res-max)) (str method))
          (t/is (< (- (second res-max)) val-tol) (str method)))))
    ;; through minimizer
    (t/is (near-minimum? (first ((sut/minimizer :lbfgsb p/himmelblau {:bounds hb :gradient himmelblau-gradient}) [1 1]))))
    ;; separate arguments of the function, one sequence for the gradient
    (doseq [method [:lbfgsb :gradient]]
      (t/is (near-minimum? (first (sut/minimize method (fn [x y] (p/himmelblau [x y]))
                                                {:bounds hb :vector-arg? false :gradient himmelblau-gradient}))) (str method)))
    ;; methods which do not use the gradient ignore it
    (doseq [method [:bobyqa :cmaes :nelder-mead :multidirectional-simplex :powell]]
      (t/is (near-minimum? (first (sut/minimize method p/himmelblau {:bounds hb :gradient (fn [_] (throw (RuntimeException. "not used")))}))) (str method)))
    ;; gradient with scan
    (t/is (near-minimum? (first (sut/scan-and-minimize :lbfgsb p/himmelblau {:bounds hb :gradient himmelblau-gradient}))))))

;; scan

(t/deftest scan-and-minimize-matrix
  (doseq [method multivariate-methods
          vector-arg? [true false nil]
          :let [f (if (false? vector-arg?) (fn [x y] (p/himmelblau [x y])) p/himmelblau)
                nf (if (false? vector-arg?) (fn [x y] (- (p/himmelblau [x y]))) neg-himmelblau)
                opts (seeded method {:bounds hb :vector-arg? vector-arg?})
                label (str method " " vector-arg?)]]
    (let [res (sut/scan-and-minimize method f opts)]
      (t/is (vector? res) label)
      (t/is (near-minimum? (first res)) label)
      (t/is (< (second res) val-tol) label))
    (let [[pt val] (sut/scan-and-maximize method nf opts)]
      (t/is (near-minimum? pt) label)
      (t/is (< (- val) val-tol) label))
    (let [s (sut/scan-and-minimize method f (assoc opts :stats? true))]
      (t/is (map? s) label)
      (t/is (< (:value s) val-tol) label)
      (t/is (near-minimum? (:point s)) label))))

(t/deftest scan-and-optimize-goal
  (t/is (= :minimize-like (if (< (second (sut/scan-and-optimize :lbfgsb p/himmelblau {:bounds hb})) val-tol) :minimize-like :other)))
  (let [[pt val] (sut/scan-and-optimize :lbfgsb p/himmelblau {:bounds hb :goal :maximize})]
    (t/is (= [5.0 5.0] pt))
    (t/is (m/delta-eq 890.0 val 1.0e-6)))
  (t/is (= {:goal :foo :allowed #{:minimize :maximize}} (ex-data-of #(sut/scan-and-optimize :lbfgsb p/himmelblau {:bounds hb :goal :foo})))))

(t/deftest scan-simplex-in-parallel
  ;; regression: one mutable simplex was shared by the threads and the scan failed
  (doseq [method [:nelder-mead :multidirectional-simplex]
          parallel? [true false]
          i (range 8)]
    (let [res (sut/scan-and-minimize method p/himmelblau {:bounds hb :N 40 :n 10 :parallel? parallel?})]
      (t/is (some? res) (str method parallel? i))
      (t/is (< (second res) val-tol) (str method parallel? i)))))

(t/deftest scan-brent
  (let [bounds (p/problem02-bounds)
        [pt val] (sut/scan-and-minimize :brent p/problem02 {:bounds bounds})]
    (t/is (m/delta-eq 5.145735 pt 1.0e-4))
    (t/is (m/delta-eq -1.899599 val 1.0e-5)))
  ;; the flat form of bounds
  (let [[pt val] (sut/scan-and-minimize :brent p/problem02 {:bounds [2.7 7.5]})]
    (t/is (double? pt))
    (t/is (m/delta-eq -1.899599 val 1.0e-5)))
  ;; problem10 bounds: a global minimum of -x sin x on [0, 10] is at 7.978666, value -7.916727
  (let [[pt val] (sut/scan-and-minimize :brent p/problem10 {:bounds (p/problem10-bounds) :N 200 :n 10})]
    (t/is (m/delta-eq 7.978666 pt 1.0e-3))
    (t/is (m/delta-eq -7.916727 val 1.0e-3)))
  ;; vector-arg?: points are vectors
  (let [[pt val] (sut/scan-and-minimize :brent (fn [[x]] (p/problem02 x)) {:bounds [2.7 7.5] :vector-arg? true})]
    (t/is (vector? pt))
    (t/is (m/delta-eq -1.899599 val 1.0e-5)))
  ;; stats, regression: brent did not return them and the scan failed on them
  (let [s (sut/scan-and-minimize :brent p/problem02 {:bounds [2.7 7.5] :stats? true})]
    (t/is (= #{:point :value :evaluations :iterations :lo :hi} (set (keys s))))
    (t/is (m/delta-eq -1.899599 (:value s) 1.0e-5)))
  ;; the global maximum on the interval (a dense grid search gives 6.217309, 0.888315),
  ;; while the local maximum at 4.1966 is found by a plain minimize :brent from its neighbourhood
  (let [[pt val] (sut/scan-and-maximize :brent p/problem02 {:bounds [2.7 7.5]})]
    (t/is (m/delta-eq 6.217309 pt 1.0e-4))
    (t/is (m/delta-eq 0.888315 val 1.0e-5))))

(t/deftest scan-number-of-results
  (let [run (fn [opts] (sut/scan-and-minimize :lbfgsb p/himmelblau (merge {:bounds hb :N 20 :take-last-n 100} opts)))]
    ;; n >= 1 is the number of runs, a fraction is a part of N, at least one run
    (t/is (= 3 (count (run {:n 3}))))
    (t/is (= 3 (count (run {:n 3.0}))))
    (t/is (= 20 (count (run {:n 1.0}))))
    (t/is (= 10 (count (run {:n 0.5}))))
    (t/is (= 1 (count (run {:n 0.01}))))
    (t/is (= 1 (count (run {:n 0.0}))))
    ;; more runs than points: not more than the points
    (t/is (<= (count (run {:n 1000})) 20)))
  (let [res (sut/scan-and-minimize :lbfgsb p/himmelblau {:bounds hb :N 20 :n 1.0 :take-last-n 5})]
    (t/is (= 5 (count res)))
    (t/is (every? vector? res))
    ;; sorted from the best
    (t/is (apply <= (map second res))))
  ;; ... also for maximization and stats
  (let [res (sut/scan-and-maximize :lbfgsb p/himmelblau {:bounds hb :N 20 :n 1.0 :take-last-n 5 :stats? true})]
    (t/is (= 5 (count res)))
    (t/is (every? map? res))
    (t/is (apply >= (map :value res))))
  ;; 0 and 1 mean the single best result
  (doseq [take-n [0 1 nil]]
    (t/is (vector? (sut/scan-and-minimize :lbfgsb p/himmelblau {:bounds hb :take-last-n take-n :N 20}))))
  ;; N smaller than the minimum is raised to it
  (t/is (some? (sut/scan-and-minimize :lbfgsb p/himmelblau {:bounds hb :N 1 :n 1.0})))
  (t/is (some? (sut/scan-and-minimize :lbfgsb p/himmelblau {:bounds hb :N 0})))
  (t/is (some? (sut/scan-and-minimize :lbfgsb p/himmelblau {:bounds hb :N 20 :jitter 0.0}))))

(t/deftest scan-failed-runs
  ;; failed runs are skipped: nothing left gives nil or an empty sequence
  (t/is (nil? (sut/scan-and-minimize :nelder-mead p/himmelblau {:bounds hb :max-evals 3})))
  (t/is (empty? (sut/scan-and-minimize :nelder-mead p/himmelblau {:bounds hb :max-evals 3 :take-last-n 3})))
  (t/is (nil? (sut/scan-and-minimize :nelder-mead p/himmelblau {:bounds hb :max-evals 3 :parallel? false})))
  ;; some runs fail, the rest is returned
  ;; the first 40 calls scan the domain (exceptions there are not swallowed), later calls are optimization runs
  (let [calls (atom 0)
        f (fn [x] (let [c (swap! calls inc)]
                    (when (and (> c 40) (zero? (rem c 200))) (throw (RuntimeException. "sometimes")))
                    (p/himmelblau x)))
        res (sut/scan-and-minimize :nelder-mead f {:bounds hb :N 40 :n 10 :parallel? false})]
    (t/is (> @calls 150) "many runs, some of them failed")
    (t/is (some? res))
    (t/is (< (second res) val-tol)))
  ;; exceptions of the function at the scanned points are thrown
  (t/is (thrown? RuntimeException (sut/scan-and-minimize :lbfgsb (fn [_] (throw (RuntimeException. "boom"))) {:bounds hb})))
  (t/is (thrown? clojure.lang.ArityException (sut/scan-and-minimize :lbfgsb (fn [_x _y] 0.0) {:bounds hb}))))

(t/deftest scan-bounds
  (t/is (= {:method :lbfgsb :bounds nil :initial nil}
           (select-keys (ex-data-of #(sut/scan-and-minimize :lbfgsb p/himmelblau {})) [:method :bounds :initial])))
  (doseq [method [:powell :gradient :nelder-mead]]
    (t/is (some? (ex-data-of #(sut/scan-and-minimize method p/himmelblau {}))) (str method)))
  (doseq [method all-methods]
    (t/is (= method (:method (ex-data-of #(sut/scan-and-minimize method p/himmelblau {:bounds [[5 -5] [0 1]]})))) (str method))
    (t/is (= method (:method (ex-data-of #(sut/scan-and-minimize method p/himmelblau {:bounds [[##-Inf ##Inf] [0 1]]})))) (str method))))

;; bayesian optimization

;; rng option

(defn- scan [method rng-fn opts]
  (sut/scan-and-minimize method p/himmelblau (merge {:bounds hb :N 40 :n 6 :stats? true :rng (rng-fn)} opts)))

(t/deftest scan-rng
  (doseq [method [:cmaes :lbfgsb :nelder-mead :bobyqa]
          rng-fn [#(r/rng :jdk 1) #(r/rng :mersenne 1) #(r/synced-rng :isaac 1)]]
    (let [label (str method " " (class (rng-fn)))
          seq-1 (scan method rng-fn {:parallel? false})]
      (t/is (some? seq-1) label)
      (t/is (= seq-1 (scan method rng-fn {:parallel? false})) (str label ": equal seeds, sequential"))
      (t/is (= seq-1 (scan method rng-fn {:parallel? true})) (str label ": parallel equals sequential"))))
  (t/testing "a not thread safe generator gives the same result in repeated parallel runs"
    (let [first-run (scan :cmaes #(r/rng :mersenne 5) {:parallel? true :n 12})]
      (dotimes [_ 5]
        (t/is (= first-run (scan :cmaes #(r/rng :mersenne 5) {:parallel? true :n 12}))))))
  (t/testing "the stochastic method depends on the seed"
    (t/is (not= (scan :cmaes #(r/rng :jdk 1) {}) (scan :cmaes #(r/rng :jdk 2) {}))))
  (t/testing "results with take-last-n are reproducible"
    (let [run #(sut/scan-and-minimize :cmaes p/himmelblau {:bounds hb :N 30 :n 5 :take-last-n 3 :rng (r/rng :jdk 3)})]
      (t/is (= (run) (run)))))
  (t/testing "a missing and a nil rng create a new generator"
    (t/is (some? (sut/scan-and-minimize :cmaes p/himmelblau {:bounds hb :N 20})))
    (t/is (some? (sut/scan-and-minimize :cmaes p/himmelblau {:bounds hb :N 20 :rng nil}))))
  (t/testing "a value which is not a generator throws ex-info"
    (doseq [bad [5 :jdk "x" (r/distribution :normal)]
            f [sut/scan-and-minimize sut/scan-and-maximize sut/scan-and-optimize]]
      (t/is (= {:rng bad} (ex-data-of #(f :cmaes p/himmelblau {:bounds hb :rng bad}))))))
  (t/testing "the shared generator is not used"
    (r/set-seed! 1)
    (let [expected (r/drand)]
      (r/set-seed! 1)
      (sut/scan-and-minimize :cmaes p/himmelblau {:bounds hb :N 20 :rng (r/rng :jdk 1)})
      (sut/scan-and-minimize :cmaes p/himmelblau {:bounds hb :N 20})
      (t/is (= expected (r/drand))))))

(defn- bayesian-steps [f opts n]
  (mapv #(select-keys % [:xs :ys :x :y]) (take n (sut/bayesian-optimization f (merge {:warm-up 50} opts)))))

(t/deftest bayesian-optimization-rng
  (let [f2 (fn [v] (- (p/himmelblau v)))
        f1 (fn [[x]] (- (p/problem02 x)))
        cases [[f2 {:bounds hb :optimizer :lbfgsb}] [f1 {:bounds [[2.7 7.5]] :optimizer :cmaes}] [f2 {:bounds hb}]]]
    (doseq [[f opts] cases]
      (let [run (fn [seed extra] (bayesian-steps f (merge opts {:rng (r/rng :jdk seed)} extra) 3))]
        (t/is (= (run 1 {}) (run 1 {})) (str opts ": equal seeds"))
        (t/is (not= (run 1 {}) (run 2 {})) (str opts ": different seeds"))
        (t/is (= (run 1 {}) (run 1 {:optimizer-params {:rng (r/rng :jdk 99)}}))
              (str opts ": rng of optimizer-params is overridden")))))
  (t/testing "initial points given as a sequence leave the steps seeded"
    (let [run #(bayesian-steps (fn [v] (- (p/himmelblau v))) {:bounds hb :init-points [[0.0 0.0] [1.0 1.0]] :rng (r/rng :jdk 4)} 2)]
      (t/is (= (run) (run)))))
  (t/testing "a missing and a nil rng create a new generator"
    (t/is (= 2 (count (bayesian-steps (fn [v] (- (p/himmelblau v))) {:bounds hb} 2))))
    (t/is (= 2 (count (bayesian-steps (fn [v] (- (p/himmelblau v))) {:bounds hb :rng nil} 2)))))
  (t/testing "a value which is not a generator throws ex-info"
    (doseq [bad [5 :jdk "x" (r/distribution :normal)]]
      (t/is (= {:rng bad} (ex-data-of #(sut/bayesian-optimization (fn [v] (- (p/himmelblau v))) {:bounds hb :rng bad}))))))
  (t/testing "the shared generator is not used"
    (r/set-seed! 1)
    (let [expected (r/drand)]
      (r/set-seed! 1)
      (bayesian-steps (fn [v] (- (p/himmelblau v))) {:bounds hb :rng (r/rng :jdk 1)} 2)
      (bayesian-steps (fn [v] (- (p/himmelblau v))) {:bounds hb} 2)
      (t/is (= expected (r/drand))))))

(t/deftest bayesian-optimization-steps
  (let [f (fn [v] (- (p/himmelblau v)))
        steps (sut/bayesian-optimization f {:bounds hb :warm-up 50 :init-points 3})
        taken (vec (take 4 steps))]
    (t/is (= 4 (count taken)))
    (doseq [[i s] (map-indexed vector taken)]
      (t/is (= #{:x :y :xs :ys :gp :util-fn :util-best} (set (keys s))))
      (t/is (= (+ 3 (inc i)) (count (:xs s))))
      (t/is (= (count (:xs s)) (count (:ys s))))
      (t/is (every? (fn [[x y]] (and (<= -5.0 x 5.0) (<= -5.0 y 5.0))) (:xs s)) "points are in the bounds")
      (t/is (= (:y s) (f (:x s))) "y is the value at x")
      (t/is (= (:y s) (apply max (:ys s))) "y is the best value")
      (t/is (= (f (:util-best s)) (first (:ys s))) "the maximum of the utility function is evaluated (the newest value is the first one)"))
    ;; the best value does not decrease
    (t/is (apply <= (map :y taken)))))

(t/deftest bayesian-optimization-options
  (let [f (fn [v] (- (p/himmelblau v)))
        step (fn [opts] (nth (sut/bayesian-optimization f (merge {:bounds hb :warm-up 50} opts)) 1))]
    (doseq [opts [{:utility-function-type :ucb} {:utility-function-type :ei} {:utility-function-type :poi}
                  {:utility-function-type :ucb :utility-param 1.0} {:utility-function-type :ei :utility-param 0.1}
                  {:kernel :gaussian} {:kscale 2.0} {:jitter 0.1} {:noise 1.0e-4} {:normalize? false}
                  {:init-points 1} {:init-points 5} {:init-points [[0 0] [1 1] [2 2]]} {:optimizer-params {:max-iters 100}}]]
      (let [s (step opts)]
        (t/is (map? s) (pr-str opts))
        (t/is (<= (:y s) 0.0) (pr-str opts))))
    (t/is (= 7 (count (:xs (step {:init-points [[0 0] [1 1] [2 2] [3 2] [-3 2]]})))) "5 given points and 2 steps"))
  ;; every optimizer of the utility function, also these without constraints and one dimensional brent
  (doseq [optimizer (remove #{:bobyqa} multivariate-methods)
          :let [s (nth (sut/bayesian-optimization (fn [v] (- (p/himmelblau v))) {:bounds hb :warm-up 50 :optimizer optimizer}) 2)]]
    (t/is (every? (fn [[x y]] (and (<= -5.0 x 5.0) (<= -5.0 y 5.0))) (:xs s)) (str optimizer " points are in the bounds")))
  (t/is (map? (nth (sut/bayesian-optimization (fn [v] (- (p/himmelblau v))) {:bounds hb :warm-up 50 :optimizer :bobyqa}) 1)))
  (doseq [optimizer [nil :brent :powell :nelder-mead :cmaes :lbfgsb]
          :let [s (nth (sut/bayesian-optimization (fn [[x]] (- (p/problem02 x))) {:bounds [[2.7 7.5]] :warm-up 50 :optimizer optimizer}) 2)]]
    (t/is (every? (fn [[x]] (<= 2.7 x 7.5)) (:xs s)) (str optimizer " points are in the bounds"))
    (t/is (vector? (:x s)) (str optimizer))))

(t/deftest bayesian-optimization-vector-arg
  (let [f (fn [x y] (- (p/himmelblau [x y])))
        s (nth (sut/bayesian-optimization f {:bounds hb :warm-up 50 :vector-arg? false}) 2)]
    (t/is (= (:y s) (apply f (:x s))))
    (t/is (= 2 (count (:x s)))))
  ;; the default is one sequence
  (let [s (nth (sut/bayesian-optimization (fn [v] (- (p/himmelblau v))) {:bounds hb :warm-up 50 :vector-arg? nil}) 1)]
    (t/is (map? s)))
  (t/is (thrown? clojure.lang.ArityException (doall (take 1 (sut/bayesian-optimization (fn [_x _y] 0.0) {:bounds hb :warm-up 50}))))))

(t/deftest bayesian-optimization-guards
  ;; unknown values throw immediately, before any step is taken
  (t/is (= {:utility-function-type :foo :allowed #{:ucb :ei :poi}}
           (ex-data-of #(sut/bayesian-optimization (fn [_] 0.0) {:bounds hb :utility-function-type :foo}))))
  (t/is (= :bfgs (:method (ex-data-of #(sut/bayesian-optimization (fn [_] 0.0) {:bounds hb :optimizer :bfgs})))))
  (doseq [bounds [nil [[5 -5] [0 1]] [[##-Inf ##Inf] [0 1]] [[0 1] [0 ##NaN]] []]]
    (t/is (some? (ex-data-of #(sut/bayesian-optimization (fn [_] 0.0) {:bounds bounds}))) (pr-str bounds)))
  ;; the utility function can not be maximized when every run of the optimizer fails
  (let [d (ex-data-of #(first (sut/bayesian-optimization (fn [v] (- (p/himmelblau v)))
                                                         {:bounds hb :warm-up 50 :optimizer :nelder-mead :optimizer-params {:max-evals 1}})))]
    (t/is (= :nelder-mead (:optimizer d)))))

;; linear optimization

(t/deftest linear-optimization-relations
  ;; every token of a relation gives the same result
  (doseq [[tokens relation] [[[:<= '<= :leq] :leq] [[:>= '>= :geq] :geq] [[:= '= :eq] :eq]]
          token tokens]
    (let [constraints (case relation
                        :leq [[1 1] token 4 [1 0] :<= 3]
                        :geq [[1 1] token 4 [1 0] :<= 3]
                        :eq [[1 1] token 4 [1 0] :<= 3])
          [pt val] (sut/linear-optimization [1 2 0] constraints {:non-negative? true
                                                                  :goal (if (= relation :leq) :maximize :minimize)})]
      (t/is (vector? pt) (pr-str token))
      (t/is (every? double? pt) (pr-str token))
      (t/is (double? val) (pr-str token))
      (t/is (case relation
              :leq (m/delta-eq 8.0 val) ;; maximize x + 2y, x + y <= 4, x <= 3: (0, 4)
              :geq (m/delta-eq 5.0 val) ;; minimize x + 2y, x + y >= 4, x <= 3: (3, 1)
              :eq (m/delta-eq 5.0 val))
            (pr-str token val)))))

(t/deftest linear-optimization-relation-guards
  (doseq [bad [:< :> :== :neq '< nil "=" "<=" 3 [:<=]]]
    (let [d (ex-data-of #(sut/linear-optimization [1 1 0] [[1 1] bad 4] {:non-negative? true}))]
      (t/is (= bad (:relation d)) (pr-str bad))
      (t/is (contains? (set (:allowed d)) :<=) (pr-str bad))))
  ;; triplets only
  (doseq [n [1 2 4 5 7]]
    (t/is (= {:count n} (ex-data-of #(sut/linear-optimization [1 1 0] (take n (cycle [[1 1] :<= 4])) {:non-negative? true}))) (str n)))
  ;; the unknown relation is found anywhere in the constraints
  (t/is (= :bad (:relation (ex-data-of #(sut/linear-optimization [1 1 0] [[1 1] :<= 4 [1 0] :bad 2] {:non-negative? true}))))))

(t/deftest linear-optimization-options
  (let [target [-1 4 0]
        constraints [[-3 1] :<= 6 [-1 -2] :>= -4 [0 1] :>= -3]
        [pt val] (sut/linear-optimization target constraints)]
    (t/is (v/delta-eq [10.0 -3.0] pt 1.0e-6))
    (t/is (m/delta-eq -22.0 val 1.0e-6))
    ;; the form with the constant term on both sides: x1 + 2 >= 3 - x2 + ... as left [a... c] R right [b... c]
    (let [[pt2 val2] (sut/linear-optimization [1 1 0] [[1 0 1] :<= [0 -1 5] [0 1] :<= 2 [1 0] :>= 0] {:goal :maximize :non-negative? true})]
      ;; x + 1 <= -y + 5 (x + y <= 4), y <= 2: maximize x + y = 4
      (t/is (m/delta-eq 4.0 val2 1.0e-6))
      (t/is (m/delta-eq 4.0 (v/sum pt2) 1.0e-6)))
    ;; rule
    (doseq [rule [:dantzig :bland]]
      (t/is (m/delta-eq -22.0 (second (sut/linear-optimization target constraints {:rule rule})) 1.0e-6) (str rule)))
    (t/is (= {:rule :foo :allowed #{:dantzig :bland}} (ex-data-of #(sut/linear-optimization target constraints {:rule :foo}))))
    (t/is (= {:rule nil :allowed #{:dantzig :bland}} (ex-data-of #(sut/linear-optimization target constraints {:rule nil}))))
    ;; goal
    (t/is (m/delta-eq 22.0 (second (sut/linear-optimization [1 -4 0] constraints {:goal :maximize})) 1.0e-6))
    (t/is (= (sut/linear-optimization target constraints) (sut/linear-optimization target constraints {:goal nil})))
    (t/is (= {:goal :foo :allowed #{:minimize :maximize}} (ex-data-of #(sut/linear-optimization target constraints {:goal :foo}))))
    ;; tolerances
    (t/is (m/delta-eq -22.0 (second (sut/linear-optimization target constraints {:epsilon 1.0e-8 :max-ulps 5 :cut-off 1.0e-9})) 1.0e-6))
    ;; the limit is :max-iters
    (t/is (thrown? TooManyIterationsException (sut/linear-optimization target constraints {:max-iters 1})))
    (t/is (vector? (sut/linear-optimization target constraints {:max-iters 100})))
    ;; stats
    (let [s (sut/linear-optimization target constraints {:stats? true})]
      (t/is (= #{:point :value :evaluations :iterations} (set (keys s))))
      (t/is (vector? (:point s)))
      (t/is (= 0 (:evaluations s)))
      (t/is (pos? (:iterations s))))
    (t/is (vector? (sut/linear-optimization target constraints {:stats? nil}))))
  ;; non-negative variables change the solution
  (let [constraints [[1 1] :<= 4]]
    (t/is (m/delta-eq 4.0 (second (sut/linear-optimization [1 1 0] constraints {:goal :maximize :non-negative? true})) 1.0e-6))
    (t/is (thrown? UnboundedSolutionException (sut/linear-optimization [1 1 0] constraints {:goal :minimize :non-negative? false})))
    (t/is (m/delta-eq 0.0 (second (sut/linear-optimization [1 1 0] constraints {:goal :minimize :non-negative? true})) 1.0e-6)))
  ;; infeasible and unbounded problems
  (t/is (thrown? NoFeasibleSolutionException (sut/linear-optimization [1 1 0] [[1 1] :>= 4 [1 1] :<= 2] {:non-negative? true})))
  (t/is (thrown? UnboundedSolutionException (sut/linear-optimization [1 1 0] [[1 1] :>= 4] {:goal :maximize :non-negative? true})))
  ;; a constant of the objective
  (t/is (m/delta-eq 6.0 (second (sut/linear-optimization [1 1 5] [[1 1] :>= 1] {:non-negative? true})) 1.0e-6)))
