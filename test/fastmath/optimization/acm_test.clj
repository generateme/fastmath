(ns fastmath.optimization.acm-test
  (:require [fastmath.optimization.acm :as sut]
            [fastmath.optimization.problems :as p]
            [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.test :as t])
  (:import [org.apache.commons.math3.exception TooManyEvaluationsException TooManyIterationsException OutOfRangeException]))

;; Reference: Himmelblau's function has four minima with value 0:
;; (3, 2), (-2.805118, 3.131313), (-3.779310, -3.283186), (3.584428, -1.848126) (Wikipedia).
;; problem02 (sin x + sin(10x/3) on [2.7, 7.5]) has its global minimum -1.899599 at 5.145735 and a local
;; maximum 0.119 at 4.1966 (https://infinity77.net/global_optimization/test_functions_1d.html).
;; Tolerances: the optimizers are checked only for sanity, 1e-3 on values (cmaes and the fast-stopping
;; settings 1e-2), 1e-3 on points close to a known minimum.

(def ^:private val-tol 1.0e-3)
(def ^:private himmelblau-minima [[3.0 2.0] [-2.805118 3.131313] [-3.779310 -3.283186] [3.584428 -1.848126]])
(def ^:private hb (p/himmelblau-bounds))

(defn- near-minimum? [pt] (boolean (some #(v/delta-eq (vec pt) % 1.0e-3) himmelblau-minima)))

(defn- himmelblau-gradient [[^double x ^double y]]
  [(m/+ (m/* 4.0 x (m/+ (m/* x x) y -11.0)) (m/* 2.0 (m/+ x (m/* y y) -7.0)))
   (m/+ (m/* 2.0 (m/+ (m/* x x) y -11.0)) (m/* 4.0 y (m/+ x (m/* y y) -7.0)))])

(defn- ex-data-of [thunk]
  (try (thunk) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))

(def ^:private multivariate
  {:bobyqa sut/bobyqa
   :cmaes sut/cmaes
   :nelder-mead sut/nelder-mead
   :multidirectional-simplex sut/multidirectional-simplex
   :powell sut/powell
   :non-linear-gradient sut/non-linear-gradient})

(defn- himmelblau-min [opt opts] (opt p/himmelblau (merge {:bounds hb} opts)))

(t/deftest public-api
  (t/is (= #{'brent 'bobyqa 'cmaes 'nelder-mead 'multidirectional-simplex 'powell 'non-linear-gradient}
           (set (keys (ns-publics 'fastmath.optimization.acm))))))

;; shapes

(t/deftest result-shape
  (doseq [[nm opt] multivariate]
    (let [res (himmelblau-min opt {})
          [pt val] res]
      (t/is (vector? res) (str nm))
      (t/is (= 2 (count res)) (str nm))
      (t/is (vector? pt) (str nm))
      (t/is (every? double? pt) (str nm))
      (t/is (double? val) (str nm))
      (t/is (near-minimum? pt) (str nm))
      (t/is (< val val-tol) (str nm)))
    (let [s (himmelblau-min opt {:stats? true})]
      (t/is (= #{:point :value :evaluations :iterations} (set (keys s))) (str nm))
      (t/is (vector? (:point s)) (str nm))
      (t/is (pos? (:evaluations s)) (str nm))
      (t/is (integer? (:iterations s)) (str nm)))
    ;; falsy :stats?
    (t/is (vector? (himmelblau-min opt {:stats? nil})) (str nm))
    (t/is (vector? (himmelblau-min opt {:stats? false})) (str nm))))

(t/deftest goal
  (doseq [[nm opt] multivariate
          :let [neg-himmelblau (fn [x] (- (p/himmelblau x)))]]
    ;; maximization: the value is the value of f
    (let [[pt val] (opt neg-himmelblau {:bounds hb :goal :maximize})]
      (t/is (near-minimum? pt) (str nm))
      (t/is (<= val 0.0) (str nm))
      (t/is (< (- val) val-tol) (str nm)))
    ;; nil is minimization
    (t/is (< (second (opt p/himmelblau {:bounds hb :goal nil})) val-tol) (str nm))
    (doseq [bad [:min :MAXIMIZE "minimize" 0]]
      (t/is (= {:goal bad :allowed #{:minimize :maximize}}
               (ex-data-of #(opt p/himmelblau {:bounds hb :goal bad}))) (str nm " " (pr-str bad))))))

(t/deftest powell-maximize-regression
  ;; Apache Commons Math stops Powell after one iteration when asked to maximize (value -1.76, -5.48, ... for -f):
  ;; the optimizer minimizes the negated function instead
  (doseq [initial [[1.0 1.0] [0.0 0.0] [-1.0 1.0] [3.3 1.7]]
          stats? [false true]]
    (let [res (sut/powell (fn [x] (- (p/himmelblau x))) {:bounds hb :initial initial :goal :maximize :stats? stats?})
          val (if stats? (:value res) (second res))
          pt (if stats? (:point res) (first res))]
      (t/is (m/delta-eq 0.0 val 1.0e-6) (pr-str initial stats?))
      (t/is (near-minimum? pt) (pr-str initial stats?))))
  ;; the same with a function with a positive maximum
  (let [[pt val] (sut/powell (fn [[x y]] (- 5.0 (m/sq (- x 1.0)) (m/sq (+ y 2.0)))) {:initial [0 0] :goal :maximize})]
    (t/is (v/delta-eq [1.0 -2.0] pt 1.0e-4))
    (t/is (m/delta-eq 5.0 val 1.0e-6))))

(t/deftest vector-arg
  (let [f2 (fn [x y] (p/himmelblau [x y]))]
    (doseq [[nm opt] multivariate]
      ;; separate arguments, also for the finite difference gradient and the hessian preconditioner
      (let [[pt val] (opt f2 {:bounds hb :vector-arg? false})]
        (t/is (near-minimum? pt) (str nm))
        (t/is (< val val-tol) (str nm)))
      ;; nil is the default: one sequence
      (t/is (near-minimum? (first (himmelblau-min opt {:vector-arg? nil}))) (str nm))
      (t/is (near-minimum? (first (himmelblau-min opt {:vector-arg? true}))) (str nm))
      ;; the function receives exactly the declared shape
      (t/is (thrown? clojure.lang.ArityException (opt f2 {:bounds hb})) (str nm))))
  (t/is (near-minimum? (first (sut/non-linear-gradient (fn [x y] (p/himmelblau [x y]))
                                                       {:bounds hb :vector-arg? false :preconditioner :hessian}))))
  (t/is (near-minimum? (first (sut/non-linear-gradient (fn [x y] (p/himmelblau [x y]))
                                                       {:bounds hb :vector-arg? false :gradient-acc 4})))))

;; initial point

(t/deftest initial-point
  (doseq [[nm opt] (dissoc multivariate :bobyqa)] ;; BOBYQA shifts the first evaluation to the bounds
    (let [calls (atom [])
          f (fn [x] (swap! calls conj (vec x)) (p/himmelblau x))
          first-call (fn [opts] (reset! calls []) (try (opt f (merge {:bounds hb :max-evals 3} opts)) (catch Exception _ nil)) (first @calls))
          ;; the gradient method evaluates around the initial point
          close? (fn [expected actual] (v/delta-eq expected actual 1.0e-5))]
      (t/is (close? [0.0 0.0] (first-call {})) (str nm " middle of the bounds"))
      (t/is (close? [1.0 2.0] (first-call {:initial [1 2]})) (str nm))
      (t/is (close? [1.0 2.0] (first-call {:initial (double-array [1 2])})) (str nm))
      (t/is (close? [1.0 2.0] (first-call {:initial '(1 2)})) (str nm))))
  ;; no bounds, a number as the initial point of a one dimensional problem
  (let [f (fn [[x]] (p/problem02 x))]
    (doseq [nm [:nelder-mead :multidirectional-simplex :powell :non-linear-gradient]]
      (t/is (v/delta-eq [5.1457] (first ((multivariate nm) f {:initial 5.0})) 1.0e-3) (str nm))
      (t/is (v/delta-eq [5.1457] (first ((multivariate nm) f {:initial [5.0]})) 1.0e-3) (str nm)))))

(t/deftest bounds-and-initial-required
  ;; unconstrained methods need the initial point or the bounds, constrained ones the bounds
  (doseq [nm [:nelder-mead :multidirectional-simplex :powell :non-linear-gradient]
          :let [opt (multivariate nm)]]
    (t/is (map? (ex-data-of #(opt p/himmelblau {}))) (str nm))
    (t/is (map? (ex-data-of #(opt p/himmelblau {:bounds nil :initial nil}))) (str nm))
    (t/is (near-minimum? (first (opt p/himmelblau {:initial [1 1]}))) (str nm)))
  (doseq [nm [:bobyqa :cmaes]
          :let [opt (multivariate nm)]]
    (let [d (ex-data-of #(opt p/himmelblau {:initial [1 1]}))]
      (t/is (= nm (:method d)) (str nm))
      (t/is (re-find #"required" (str (:reason d))) (str nm)))))

(t/deftest bounds-guards
  ;; the guards of fastmath.optimization.common are active in every optimizer
  (doseq [[nm opt] multivariate]
    (doseq [[bounds initial part] [[[[5.0 -5.0] [0.0 1.0]] nil "greater"]
                                   [[[0.0 ##NaN] [0.0 1.0]] nil "NaN"]
                                   [[[nil 1.0] [0.0 1.0]] nil "pair of numbers"]
                                   [[[0.0 1.0] [0.0 1.0]] [0.0 0.0 0.0] "differs"]
                                   [[[0.0 1.0] [0.0 1.0]] [0.0] "differs"]]]
      (let [d (ex-data-of #(opt p/himmelblau {:bounds bounds :initial initial}))]
        (t/is (= nm (:method d)) (str nm " " (pr-str bounds initial)))
        (t/is (re-find (re-pattern part) (str (:reason d))) (str nm " " (pr-str bounds initial))))))
  ;; the objective is not called for invalid input
  (let [calls (atom 0)]
    (doseq [[_ opt] multivariate]
      (ex-data-of #(opt (fn [_] (swap! calls inc) 0.0) {:bounds [[1 0] [0 1]]})))
    (t/is (zero? @calls))))

(t/deftest bounds-finite
  ;; infinite bounds are rejected where the method needs the range
  (doseq [nm [:bobyqa :cmaes :nelder-mead :multidirectional-simplex]
          initial [nil [0.0 0.0]]]
    (let [d (ex-data-of #((multivariate nm) p/himmelblau {:bounds [[##-Inf ##Inf] [-5.0 5.0]] :initial initial}))]
      (t/is (some? d) (str nm initial))
      (t/is (= nm (:method d)))))
  ;; ... and accepted where the bounds are not used, with the initial point
  (doseq [nm [:powell :non-linear-gradient]]
    (t/is (near-minimum? (first ((multivariate nm) p/himmelblau {:bounds [[##-Inf ##Inf] [##-Inf ##Inf]] :initial [1 1]}))) (str nm)))
  ;; a flat range is legal where it works: bobyqa, cmaes; it is not for the simplex methods
  (doseq [nm [:bobyqa :cmaes]]
    (let [[pt] ((multivariate nm) (fn [[x y]] (p/himmelblau [x y])) {:bounds [[1.0 1.0] [-5.0 5.0]] :initial [1.0 0.0]})]
      (t/is (= 1.0 (first pt)) (str nm))))
  (doseq [nm [:nelder-mead :multidirectional-simplex]]
    (t/is (some? (ex-data-of #((multivariate nm) p/himmelblau {:bounds [[1.0 1.0] [-5.0 5.0]] :initial [1.0 0.0]}))) (str nm))))

;; limits

(t/deftest limits
  (doseq [nm [:bobyqa :nelder-mead :multidirectional-simplex :powell :non-linear-gradient]]
    (t/is (thrown? TooManyEvaluationsException (himmelblau-min (multivariate nm) {:max-evals 8})) (str nm)))
  (doseq [nm [:nelder-mead :multidirectional-simplex :powell :non-linear-gradient]]
    (t/is (thrown? TooManyIterationsException (himmelblau-min (multivariate nm) {:max-iters 1})) (str nm)))
  ;; BOBYQA does not count iterations, CMA-ES stops at the limit and returns the best point
  (t/is (= 0 (:iterations (himmelblau-min sut/bobyqa {:max-iters 1 :stats? true}))))
  (let [s (himmelblau-min sut/cmaes {:max-iters 3 :stats? true})]
    (t/is (<= (:iterations s) 3))
    (t/is (vector? (:point s))))
  ;; generous limits are fine
  (doseq [[nm opt] multivariate]
    (t/is (< (second (himmelblau-min opt {:max-evals 100000 :max-iters 100000})) val-tol) (str nm))))

;; bobyqa

(t/deftest bobyqa-options
  (doseq [opts [{:number-of-points 4} {:number-of-points 5} {:number-of-points 6}
                {:initial-radius :default} {:initial-radius :inferred} {:initial-radius 1.0} {:initial-radius 2}
                {:stopping-radius 1.0e-6} {:initial [1 1]} {:max-iters 5}]]
    (let [[pt val] (himmelblau-min sut/bobyqa opts)]
      (t/is (near-minimum? pt) (pr-str opts))
      (t/is (< val val-tol) (pr-str opts))))
  ;; a coarse stopping radius stops early
  (t/is (> (:evaluations (himmelblau-min sut/bobyqa {:stats? true :stopping-radius 1.0e-8}))
           (:evaluations (himmelblau-min sut/bobyqa {:stats? true :stopping-radius 1.0e-1}))))
  ;; the number of interpolation points is limited by Apache Commons Math: [n+2, (n+1)(n+2)/2] = [4, 6]
  (doseq [bad [3 7 100]]
    (t/is (thrown? OutOfRangeException (himmelblau-min sut/bobyqa {:number-of-points bad})) (str bad)))
  ;; at least two dimensions
  (t/is (= {:n 1 :bounds [[0.0 5.0]]} (ex-data-of #(sut/bobyqa (fn [[x]] x) {:bounds [0 5]}))))
  ;; three dimensions, default number of points 2n+1 = 7
  (let [[pt] (sut/bobyqa p/sphere {:bounds (p/sphere-bounds 3)})]
    (t/is (v/delta-eq [0.0 0.0 0.0] pt 1.0e-4))))

;; cmaes

(t/deftest cmaes-options
  (doseq [opts [{:sigma 0.5} {:sigma 0.05} {:population-size 6} {:population-size 40} {:active-cma? false} {:diagonal-only 5}
                {:check-feasible-count 3} {:rel 1.0e-6 :abs 1.0e-6} {:initial [1 1]}]]
    (let [[pt val] (himmelblau-min sut/cmaes opts)]
      (t/is (near-minimum? pt) (pr-str opts))
      (t/is (< val 1.0e-2) (pr-str opts))))
  ;; stop-fitness: stops as soon as the value is below it
  (let [[_ val] (himmelblau-min sut/cmaes {:stop-fitness 1.0})]
    (t/is (< val 1.0)))
  ;; the generator makes the run reproducible: same seed, same result, other seed, other result
  (let [run (fn [seed] (himmelblau-min sut/cmaes {:rng (r/rng :jdk seed) :stats? true}))]
    (t/is (= (run 42) (run 42)))
    (t/is (not= (:point (run 42)) (:point (run 43)))))
  (let [run #(himmelblau-min sut/cmaes {:rng (r/rng :isaac 1) :stats? true})]
    (t/is (= (run) (run)))))

;; simplex

(t/deftest simplex-options
  (doseq [[nm opt] (select-keys multivariate [:nelder-mead :multidirectional-simplex])]
    (doseq [opts [{:length 0.1} {:length 0.9} {:rel 1.0e-6 :abs 1.0e-6} {:initial [1 1]} {:khi 1.5} {:gamma 0.3}]]
      (let [[pt val] (himmelblau-min opt opts)]
        (t/is (near-minimum? pt) (str nm (pr-str opts)))
        (t/is (< val val-tol) (str nm (pr-str opts)))))
    ;; without bounds the length is absolute
    (let [[pt] (opt p/himmelblau {:initial [1 1] :length 0.5})]
      (t/is (near-minimum? pt) (str nm)))
    ;; the initial simplex depends on the length
    (t/is (not= (himmelblau-min opt {:length 0.1 :initial [1 1] :stats? true})
                (himmelblau-min opt {:length 0.8 :initial [1 1] :stats? true})) (str nm))
    ;; no state is shared between calls
    (t/is (= (himmelblau-min opt {:stats? true}) (himmelblau-min opt {:stats? true})) (str nm)))
  (doseq [opts [{:rho 2.0} {:khi 3.0} {:gamma 0.4} {:sigma 0.3}]]
    (t/is (near-minimum? (first (himmelblau-min sut/nelder-mead opts))) (pr-str opts))))

;; powell

(t/deftest powell-options
  (doseq [opts [{:line-rel 1.0e-4 :line-abs 1.0e-4} {:rel 1.0e-6 :abs 1.0e-6} {:initial [1 1]} {:bounds nil :initial [1 1]}]]
    (let [[pt val] (sut/powell p/himmelblau (merge {:bounds hb} opts))]
      (t/is (near-minimum? pt) (pr-str opts))
      (t/is (< val val-tol) (pr-str opts)))))

;; non-linear gradient

(t/deftest gradient-options
  (doseq [opts [{:formula :polak-ribiere} {:formula :fletcher-reeves} {:preconditioner :identity} {:preconditioner :hessian}
                {:preconditioner :hessian :hessian-h 1.0e-2} {:gradient-h 1.0e-4} {:gradient-acc 2} {:gradient-acc 4}
                {:bracketing-range 1.0e-3} {:line-rel 1.0e-6 :line-abs 1.0e-6} {:rel 1.0e-8 :abs 1.0e-8} {:initial [1 1]}
                {:bounds nil :initial [1 1]}]]
    (let [[pt val] (himmelblau-min sut/non-linear-gradient opts)]
      (t/is (near-minimum? pt) (pr-str opts))
      (t/is (< val val-tol) (pr-str opts)))))

(t/deftest gradient-validation
  (doseq [[opts k] [[{:formula :foo} :formula] [{:formula nil} :formula] [{:preconditioner :foo} :preconditioner]
                    [{:gradient-acc 3} :gradient-acc] [{:gradient-acc 0} :gradient-acc] [{:gradient-acc nil} :gradient-acc]
                    [{:gradient-h 0.0} :gradient-h] [{:gradient-h -1.0e-6} :gradient-h]]]
    (t/is (contains? (ex-data-of #(himmelblau-min sut/non-linear-gradient opts)) k) (pr-str opts)))
  (t/is (= #{2 4} (:allowed (ex-data-of #(himmelblau-min sut/non-linear-gradient {:gradient-acc 3}))))))

(t/deftest user-gradient
  (let [counter (atom 0)
        args (atom [])
        g (fn [x] (swap! counter inc) (swap! args conj x) (himmelblau-gradient x))
        [pt val] (himmelblau-min sut/non-linear-gradient {:gradient g :stats? false})]
    (t/is (near-minimum? pt))
    (t/is (< val val-tol))
    (t/is (pos? @counter) "gradient is used")
    (t/is (every? #(= 2 (count %)) @args)))
  ;; the same gradient as the numerical one
  (t/is (v/delta-eq (first (himmelblau-min sut/non-linear-gradient {}))
                    (first (himmelblau-min sut/non-linear-gradient {:gradient himmelblau-gradient})) 1.0e-4))
  ;; with fewer calls of f than the numerical gradient (which calls f 2n times per gradient)
  (let [calls (fn [opts] (let [c (atom 0)]
                           (sut/non-linear-gradient (fn [x] (swap! c inc) (p/himmelblau x)) (merge {:bounds hb} opts))
                           @c))]
    (t/is (< (calls {:gradient himmelblau-gradient}) (calls {}))))
  ;; it is the gradient of f also when maximizing
  (let [[pt val] (sut/non-linear-gradient (fn [x] (- (p/himmelblau x)))
                                          {:bounds hb :goal :maximize :gradient (fn [x] (mapv - (himmelblau-gradient x)))})]
    (t/is (near-minimum? pt))
    (t/is (< (- val) val-tol)))
  ;; any sequence of numbers, separate arguments of f
  (doseq [conv [vec seq double-array #(apply list %)]]
    (let [[pt] (sut/non-linear-gradient (fn [x y] (p/himmelblau [x y]))
                                        {:bounds hb :vector-arg? false :gradient (fn [x] (conv (himmelblau-gradient x)))})]
      (t/is (near-minimum? pt))))
  ;; wrong length
  (doseq [bad [[1.0] [1.0 2.0 3.0] [] nil]]
    (t/is (= {:expected 2 :actual (count bad)}
             (ex-data-of #(himmelblau-min sut/non-linear-gradient {:gradient (fn [_] bad)}))) (pr-str bad)))
  ;; the gradient is ignored by the methods which do not use it
  (doseq [nm [:bobyqa :cmaes :nelder-mead :multidirectional-simplex :powell]]
    (t/is (near-minimum? (first (himmelblau-min (multivariate nm) {:gradient (fn [_] (throw (RuntimeException. "not used")))}))) (str nm))))

;; one dimension

(t/deftest one-dimension
  (let [f (fn [[x]] (p/problem02 x))
        bounds (p/problem02-bounds)]
    (doseq [nm [:cmaes :nelder-mead :multidirectional-simplex :powell :non-linear-gradient]]
      (let [[pt val] ((multivariate nm) f {:bounds bounds})]
        (t/is (vector? pt) (str nm))
        (t/is (= 1 (count pt)) (str nm))
        (t/is (m/delta-eq 5.1457 (first pt) 1.0e-1) (str nm)) ;; cmaes is stochastic and approximate
        (t/is (< val -1.88) (str nm))))
    ;; the flat form of bounds
    (t/is (< (second (sut/powell f {:bounds [2.7 7.5]})) -1.89))
    ;; bobyqa needs two dimensions
    (t/is (= 1 (:n (ex-data-of #(sut/bobyqa f {:bounds bounds})))))))

;; brent

(t/deftest brent-shape
  (let [f p/problem02
        bounds (p/problem02-bounds)
        [pt val :as res] (sut/brent f {:bounds bounds})]
    (t/is (vector? res))
    (t/is (double? pt) "number for the default vector-arg? false")
    (t/is (m/delta-eq 5.145735 pt 1.0e-4))
    (t/is (m/delta-eq -1.899599 val 1.0e-5))
    ;; flat bounds, nil vector-arg?
    (t/is (= res (sut/brent f {:bounds [2.7 7.5]})))
    (t/is (= res (sut/brent f {:bounds bounds :vector-arg? nil})))
    ;; vector-arg?: the function and the result use a sequence with one number
    (let [[vpt vval] (sut/brent (fn [x] (f (first x))) {:bounds bounds :vector-arg? true})]
      (t/is (vector? vpt))
      (t/is (= [pt] vpt))
      (t/is (= val vval)))
    (t/is (thrown? Exception (sut/brent (fn [x] (f (first x))) {:bounds bounds :vector-arg? false})))
    ;; stats
    (let [s (sut/brent f {:bounds bounds :stats? true})
          sv (sut/brent (fn [x] (f (first x))) {:bounds bounds :stats? true :vector-arg? true})]
      (t/is (= #{:point :value :evaluations :iterations :lo :hi} (set (keys s))))
      (t/is (= pt (:point s)))
      (t/is (= val (:value s)))
      (t/is (pos? (:evaluations s)))
      (t/is (pos? (:iterations s)))
      (t/is (= [pt] (:point sv))))
    (t/is (vector? (sut/brent f {:bounds bounds :stats? nil})))))

(t/deftest brent-options
  (let [f p/problem02
        bounds (p/problem02-bounds)]
    ;; initial point: a number or a sequence with one number
    (doseq [initial [5.0 [5.0] (double-array [5.0]) 3 '(5.0)]]
      (t/is (m/delta-eq 5.145735 (first (sut/brent f {:bounds bounds :initial initial})) 1.0e-4) (pr-str initial)))
    ;; goal
    (let [[pt val] (sut/brent f {:bounds bounds :goal :maximize})]
      (t/is (m/delta-eq 4.1966 pt 1.0e-3))
      (t/is (m/delta-eq 0.119 val 1.0e-3)))
    (t/is (= (sut/brent f {:bounds bounds}) (sut/brent f {:bounds bounds :goal nil})))
    (t/is (= {:goal :foo :allowed #{:minimize :maximize}} (ex-data-of #(sut/brent f {:bounds bounds :goal :foo}))))
    ;; tolerances: looser tolerances need fewer evaluations
    (t/is (> (:evaluations (sut/brent f {:bounds bounds :stats? true}))
             (:evaluations (sut/brent f {:bounds bounds :stats? true :rel 1.0e-2 :abs 1.0e-2}))))
    ;; bracket: true or a map; a local minimum of an interval is found
    (doseq [fb [true {} {:grow-limit 50.0 :max-evals 100}]]
      (t/is (m/delta-eq 3.387252 (first (sut/brent f {:bounds [[3 4]] :find-bracket fb})) 1.0e-4) (pr-str fb)))
    ;; limits are errors
    (t/is (thrown? TooManyEvaluationsException (sut/brent f {:bounds bounds :max-evals 3})))
    (t/is (thrown? TooManyIterationsException (sut/brent f {:bounds bounds :max-iters 2})))
    (t/is (thrown? OutOfRangeException (sut/brent f {:bounds bounds :initial 10.0})))
    ;; an interval with a boundary minimum
    (t/is (m/delta-eq 3.0 (first (sut/brent (fn [x] x) {:bounds [3 7]})) 1.0e-3))))

(t/deftest brent-bounds-guards
  (doseq [[bounds initial part] [[nil nil "required"]
                                 [[] nil "non-empty"]
                                 [[[2.7 7.5] [1.0 2.0]] nil "exactly one"]
                                 [[[7.5 2.7]] nil "greater"]
                                 [[[2.7 2.7]] nil "less than"]
                                 [[[2.7 ##Inf]] nil "finite"]
                                 [[[##-Inf ##Inf]] 1.0 "finite"]
                                 [[[2.7 ##NaN]] nil "NaN"]
                                 [[[2.7 7.5]] [4.0 5.0] "differs"]]]
    (let [d (ex-data-of #(sut/brent p/problem02 {:bounds bounds :initial initial}))]
      (t/is (= :brent (:method d)) (pr-str bounds initial))
      (t/is (re-find (re-pattern part) (str (:reason d))) (pr-str bounds initial)))))

;; concurrency

(t/deftest parallel-runs
  ;; every call creates its own Apache Commons Math objects: the mutable simplex used to be shared between
  ;; the threads of a scan and made it fail
  (doseq [[nm opt] multivariate]
    (let [results (doall (pmap (fn [i] (opt p/himmelblau {:bounds hb :initial [(- (mod i 7) 3) (- (mod i 5) 2)] :stats? true})) (range 60)))]
      (t/is (= 60 (count results)) (str nm))
      (t/is (every? #(< (:value %) val-tol) results) (str nm))))
  (let [results (doall (pmap (fn [i] (sut/brent p/problem02 {:bounds (p/problem02-bounds) :initial (+ 3.0 (* 0.05 i))})) (range 60)))]
    (t/is (= 60 (count results)))))
