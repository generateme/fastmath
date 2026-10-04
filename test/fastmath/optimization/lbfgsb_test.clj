(ns fastmath.optimization.lbfgsb-test
  (:require [fastmath.optimization.lbfgsb :as sut]
            [fastmath.optimization.problems :as p]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.test :as t])
  (:import [org.generateme.lbfgsb Parameters Parameters$LINESEARCH LBFGSBException]))

;; Reference: closed-form minimum of the convex quadratic below is (1, -2) with value 0,
;; Himmelblau's function has minima with value 0 (one of them at (3, 2)).
;; Tolerances: L-BFGS-B stops when the projected gradient norm is below 1e-8, the observed error
;; of the point is about 1e-9; the tests use 1e-5 on points and 1e-3 on values (sanity, not accuracy).

(def ^:private pt-tol 1.0e-5)
(def ^:private val-tol 1.0e-3)

(defn- quad ^double [[^double x ^double y]] (m/+ (m/sq (m/- x 1.0)) (m/sq (m/+ y 2.0))))
(defn- quad-gradient [[^double x ^double y]] [(m/* 2.0 (m/- x 1.0)) (m/* 2.0 (m/+ y 2.0))])
(def ^:private quad-bounds [[-10.0 10.0] [-10.0 10.0]])

(defn- ex-data-of [thunk]
  (try (thunk) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))

(defn- point-close? [a b] (v/delta-eq (vec a) (vec b) pt-tol))

;; parameters

(t/deftest parameters-defaults
  (let [^Parameters p (sut/parameters {})]
    (t/is (= 6 (.-m p)))
    (t/is (= 1.0e-8 (.-epsilon p)))
    (t/is (= 1.0e-8 (.-epsilon_rel p)))
    (t/is (= 3 (.-past p)))
    (t/is (= 1.0e-10 (.-delta p)))
    (t/is (= 1000 (.-max_iterations p)))
    (t/is (= 10 (.-max_submin p)))
    (t/is (= 20 (.-max_linesearch p)))
    (t/is (= Parameters$LINESEARCH/MORETHUENTE_ORIG (.-linesearch p)))
    (t/is (= 1.0e-8 (.-xtol p)))
    (t/is (= 1.0e-20 (.min_step p)))
    (t/is (= 1.0e20 (.max_step p)))
    (t/is (= 1.0e-4 (.ftol p)))
    (t/is (= 0.9 (.wolfe p)))
    (t/is (true? (.weak_wolfe p)))))

(t/deftest parameters-values
  (let [^Parameters p (sut/parameters {:m 3 :abs 1.0e-6 :rel 1.0e-5 :past 0 :delta 0.0 :max-iters 0 :max-submin 0
                                       :max-linesearch 7 :linesearch :lewis-overton :xtol 1.0e-6 :min-step 1.0e-10
                                       :max-step 1.0e-10 :ftol 0.25 :wolfe 0.5 :weak-wolfe? false})]
    (t/is (= 3 (.-m p)))
    (t/is (= 1.0e-6 (.-epsilon p)))
    (t/is (= 1.0e-5 (.-epsilon_rel p)))
    (t/is (= 0 (.-past p)))
    (t/is (= 0.0 (.-delta p)))
    (t/is (= 0 (.-max_iterations p)))
    (t/is (= 0 (.-max_submin p)))
    (t/is (= 7 (.-max_linesearch p)))
    (t/is (= Parameters$LINESEARCH/LEWISOVERTON (.-linesearch p)))
    (t/is (= 1.0e-6 (.-xtol p)))
    (t/is (= 1.0e-10 (.min_step p)))
    (t/is (= 1.0e-10 (.max_step p)))
    (t/is (= 0.25 (.ftol p)))
    (t/is (= 0.5 (.wolfe p)))
    (t/is (false? (.weak_wolfe p)))))

(t/deftest parameters-validation
  ;; invalid: just outside of the allowed range
  (doseq [bad [{:m 0} {:m -1} {:abs 0.0} {:abs -1.0e-9} {:rel 0.0} {:rel -1.0} {:max-linesearch 0} {:max-linesearch -3}
               {:min-step 0.0} {:min-step -1.0} {:past -1} {:delta -1.0e-12} {:max-iters -1} {:max-submin -1}
               {:max-step 1.0e-30} {:ftol 0.0} {:ftol 0.5} {:ftol -0.1} {:wolfe 1.0} {:wolfe 1.0e-4} {:wolfe 2.0}
               {:ftol 0.3 :wolfe 0.3}]]
    (t/is (thrown? clojure.lang.ExceptionInfo (sut/parameters bad)) (pr-str bad)))
  ;; valid: exactly on the allowed edge
  (doseq [ok [{:m 1} {:past 0} {:delta 0.0} {:max-iters 0} {:max-submin 0} {:max-step 1.0e-20} {:max-linesearch 1}
              {:ftol 0.4999 :wolfe 0.9999} {:abs Double/MIN_VALUE}]]
    (t/is (instance? Parameters (sut/parameters ok)) (pr-str ok))))

(t/deftest linesearch-names
  (doseq [[nm expected] {:more-thuente Parameters$LINESEARCH/MORETHUENTE_ORIG
                         :orig Parameters$LINESEARCH/MORETHUENTE_ORIG
                         :more-thuente-lbfgspp Parameters$LINESEARCH/MORETHUENTE_LBFGSPP
                         :lbfgsb Parameters$LINESEARCH/MORETHUENTE_LBFGSPP
                         :lewis-overton Parameters$LINESEARCH/LEWISOVERTON}]
    (t/is (= expected (.-linesearch (sut/parameters {:linesearch nm}))) (str nm))
    ;; ... and every line search finds the minimum
    (let [[pt val] (sut/lbfgsb quad {:bounds quad-bounds :linesearch nm})]
      (t/is (point-close? [1.0 -2.0] pt) (str nm))
      (t/is (< val val-tol) (str nm))))
  ;; the old misspelled name, nil, strings are rejected with the allowed values listed
  (doseq [bad [:levis-overton nil "lewis-overton" :foo]]
    (let [d (ex-data-of #(sut/parameters {:linesearch bad}))]
      (t/is (= bad (:linesearch d)) (pr-str bad))
      (t/is (contains? (set (:allowed d)) :lewis-overton) (pr-str bad)))))

;; basic runs

(t/deftest minimize-basic
  (let [res (sut/lbfgsb quad {:bounds quad-bounds})
        [pt val] res]
    (t/is (vector? res))
    (t/is (= 2 (count res)))
    (t/is (vector? pt))
    (t/is (every? double? pt))
    (t/is (double? val))
    (t/is (point-close? [1.0 -2.0] pt))
    (t/is (m/delta-eq 0.0 val 1.0e-9))))

(t/deftest goal
  (let [neg-quad (fn [v] (- (quad v)))
        [pt val] (sut/lbfgsb neg-quad {:bounds quad-bounds :goal :maximize})]
    (t/is (point-close? [1.0 -2.0] pt))
    (t/is (m/delta-eq 0.0 val 1.0e-9))
    ;; maximization returns the value of f, not of -f
    (t/is (m/delta-eq (quad [-10.0 10.0])
                      (second (sut/lbfgsb quad {:bounds quad-bounds :goal :maximize}))
                      1.0e-6) "maximum of the convex quadratic is in a corner"))
  (t/is (= (sut/lbfgsb quad {:bounds quad-bounds})
           (sut/lbfgsb quad {:bounds quad-bounds :goal nil})
           (sut/lbfgsb quad {:bounds quad-bounds :goal :minimize})))
  (doseq [bad [:min :MAXIMIZE "maximize" 1]]
    (t/is (= {:goal bad :allowed #{:minimize :maximize}}
             (ex-data-of #(sut/lbfgsb quad {:bounds quad-bounds :goal bad}))) (pr-str bad))))

(t/deftest vector-arg
  (let [f2 (fn [x y] (quad [x y]))
        expected (sut/lbfgsb quad {:bounds quad-bounds})]
    (t/is (= expected (sut/lbfgsb f2 {:bounds quad-bounds :vector-arg? false})))
    ;; nil is the default, true
    (t/is (= expected (sut/lbfgsb quad {:bounds quad-bounds :vector-arg? nil})))
    (t/is (= expected (sut/lbfgsb quad {:bounds quad-bounds :vector-arg? true})))
    ;; 1d, one argument
    (t/is (point-close? [2.0] (first (sut/lbfgsb (fn [x] (m/sq (m/- x 2.0))) {:bounds [[0.0 5.0]] :vector-arg? false}))))
    ;; the function receives exactly the declared shape
    (t/is (thrown? clojure.lang.ArityException (sut/lbfgsb f2 {:bounds quad-bounds})))))

(t/deftest dimensions
  ;; 1d with flat bounds, many dimensions
  (t/is (point-close? [2.0] (first (sut/lbfgsb (fn [[x]] (m/sq (m/- x 2.0))) {:bounds [0.0 5.0]}))))
  (let [n 50
        [pt val] (sut/lbfgsb p/sphere {:bounds (p/sphere-bounds n)})]
    (t/is (= n (count pt)))
    (t/is (every? #(m/delta-eq 0.0 % 1.0e-6) pt))
    (t/is (< val 1.0e-9)))
  ;; himmelblau: one of four minima
  (let [[_ val] (sut/lbfgsb p/himmelblau {:bounds (p/himmelblau-bounds)})]
    (t/is (< val val-tol))))

;; initial point

(t/deftest initial-point
  (let [calls (atom [])
        f (fn [v] (swap! calls conj (vec v)) (quad v))
        ;; with the gradient given the first call of the objective is at the initial point
        ;; (finite differences would call it first at a neighbouring point)
        first-call (fn [opts] (reset! calls []) (sut/lbfgsb f (merge {:bounds [[-4.0 6.0] [0.0 2.0]] :gradient quad-gradient} opts)) (first @calls))]
    (t/is (= [1.0 1.0] (first-call {})) "middle of the bounds")
    (t/is (= [3.0 0.5] (first-call {:initial [3 0.5]})))
    (t/is (= [3.0 0.5] (first-call {:initial (double-array [3 0.5])})))
    (t/is (= [3.0 0.5] (first-call {:initial '(3 0.5)})))
    ;; outside the bounds: moved to the nearest bound, objective never sees the outside
    (t/is (= [6.0 0.0] (first-call {:initial [100 -100]})))
    (t/is (= [-4.0 2.0] (first-call {:initial [-100 100]})))
    (t/is (= [6.0 0.0] (first-call {:initial [6.0 0.0]})) "on the bounds"))
  ;; middle of huge bounds does not overflow
  (let [calls (atom [])]
    (sut/lbfgsb (fn [[x]] (swap! calls conj x) (m/sq x)) {:bounds [[(- Double/MAX_VALUE) Double/MAX_VALUE]] :gradient (fn [[x]] [(m/* 2.0 x)]) :max-iters 1})
    (t/is (= 0.0 (first @calls)))))

;; bounds

(t/deftest bounds-guards
  (doseq [[bounds initial reason-part]
          [[nil nil "required"]
           [[] nil "non-empty"]
           [[[5.0 -5.0] [0.0 1.0]] nil "greater"]
           [[[0.0 ##NaN] [0.0 1.0]] nil "NaN"]
           [[[nil 1.0] [0.0 1.0]] nil "pair of numbers"]
           [[[0.0 1.0]] [0.0 0.0] "differs"]
           [[[0.0 1.0] [0.0 1.0]] [0.5] "differs"]
           [[[##-Inf ##Inf] [0.0 1.0]] nil "initial"]
           [[[0.0 ##Inf]] nil "initial"]]]
    (let [d (ex-data-of #(sut/lbfgsb quad {:bounds bounds :initial initial}))]
      (t/is (some? d) (pr-str bounds initial))
      (t/is (= :lbfgsb (:method d)) (pr-str bounds initial))
      (t/is (re-find (re-pattern reason-part) (str (:reason d))) (pr-str bounds initial))))
  ;; the objective is never called when the input is invalid
  (let [calls (atom 0)]
    (ex-data-of #(sut/lbfgsb (fn [_] (swap! calls inc) 0.0) {:bounds [[1 0]]}))
    (t/is (zero? @calls))))

(t/deftest infinite-bounds
  ;; unconstrained, constrained from one side, with the initial point
  (let [[pt] (sut/lbfgsb quad {:bounds [[##-Inf ##Inf] [##-Inf ##Inf]] :initial [5 5]})]
    (t/is (point-close? [1.0 -2.0] pt)))
  (let [[pt] (sut/lbfgsb quad {:bounds [[2.0 ##Inf] [##-Inf -3.0]] :initial [5 5]})]
    (t/is (point-close? [2.0 -3.0] pt) "optimum on the finite bounds")))

(t/deftest bounds-are-respected
  ;; the objective is never evaluated outside the box, also by the finite differences, also at the bounds
  (let [calls (atom [])
        bounds [[0.0 1.0] [0.0 1.0]]
        f (fn [v] (swap! calls conj (vec v)) (m/+ (v/sum v)))
        [pt val] (sut/lbfgsb f {:bounds bounds :initial [0.5 0.5]})]
    (t/is (point-close? [0.0 0.0] pt))
    (t/is (m/delta-eq 0.0 val 1.0e-6))
    (t/is (every? (fn [[x y]] (and (<= 0.0 x 1.0) (<= 0.0 y 1.0))) @calls))
    (t/is (> (count @calls) 5)))
  ;; the same on the upper bound when maximizing
  (let [calls (atom [])
        [pt] (sut/lbfgsb (fn [[x]] (swap! calls conj x) x) {:bounds [[0.0 1.0]] :goal :maximize :initial [0.5]})]
    (t/is (point-close? [1.0] pt))
    (t/is (every? #(<= 0.0 % 1.0) @calls))))

(t/deftest domain-restricted-objective
  ;; regression: the finite difference gradient used to step outside of the bounds, the objective returned
  ;; NaN there and the whole optimization ended with NaN after 1000 iterations
  (let [f (fn [[x]] (if (neg? x) ##NaN (m/+ (m/pow x 1.5) (m/* 0.5 x))))
        [pt val] (sut/lbfgsb f {:bounds [[0.0 4.0]] :initial [2.0]})]
    (t/is (point-close? [0.0] pt))
    (t/is (m/delta-eq 0.0 val 1.0e-9)))
  (let [f (fn [[x]] (if (or (neg? x) (> x 4.0)) ##NaN (m/- (m/sqrt x) x)))
        [pt] (sut/lbfgsb f {:bounds [[0.0 4.0]] :goal :maximize :initial [4.0]})]
    (t/is (point-close? [0.25] pt) "maximum of sqrt(x) - x is at 1/4, inside, from the upper bound")))

(t/deftest narrow-and-degenerate-bounds
  ;; box narrower than the finite difference step (1e-6)
  (let [calls (atom [])
        f (fn [[x y]] (swap! calls conj [x y]) (m/+ (m/* 3.0 x) y))
        [pt] (sut/lbfgsb f {:bounds [[0.0 1.0e-8] [0.0 1.0]] :initial [5.0e-9 0.5]})]
    (t/is (point-close? [0.0 0.0] pt))
    (t/is (every? (fn [[x y]] (and (<= 0.0 x 1.0e-8) (<= 0.0 y 1.0))) @calls)))
  ;; zero width: the variable is fixed
  (let [[pt val] (sut/lbfgsb (fn [[x y]] (m/+ (m/sq (m/- x 4.0)) (m/sq (m/- y 1.0)))) {:bounds [[2.0 2.0] [0.0 3.0]]})]
    (t/is (point-close? [2.0 1.0] pt))
    (t/is (m/delta-eq 4.0 val 1.0e-9))))

;; gradient

(t/deftest user-gradient
  (let [counter (atom 0)
        args (atom [])
        g (fn [v] (swap! counter inc) (swap! args conj v) (quad-gradient v))
        [pt val] (sut/lbfgsb quad {:bounds quad-bounds :gradient g})]
    (t/is (point-close? [1.0 -2.0] pt))
    (t/is (m/delta-eq 0.0 val 1.0e-9))
    (t/is (pos? @counter) "gradient is used")
    ;; one sequence of the point
    (t/is (every? #(= 2 (count %)) @args))
    (t/is (every? #(every? number? %) @args)))
  ;; any sequence of numbers is accepted as the gradient
  (doseq [conv [vec seq double-array (fn [v] (map double v)) (fn [v] (list* v)) (fn [v] (into-array Double v))]]
    (let [[pt] (sut/lbfgsb quad {:bounds quad-bounds :gradient (fn [v] (conv (quad-gradient v)))})]
      (t/is (point-close? [1.0 -2.0] pt))))
  ;; integers as well
  (t/is (point-close? [3.0] (first (sut/lbfgsb (fn [[x]] (m/sq (m/- x 3.0))) {:bounds [[-10 10]] :gradient (fn [[x]] [(* 2 (- x 3))])}))))
  ;; it is the gradient of f also for maximization
  (let [[pt val] (sut/lbfgsb (fn [v] (- (quad v))) {:bounds quad-bounds :goal :maximize
                                                    :gradient (fn [v] (mapv - (quad-gradient v)))})]
    (t/is (point-close? [1.0 -2.0] pt))
    (t/is (m/delta-eq 0.0 val 1.0e-9)))
  ;; gradient follows the single sequence convention also when f takes separate arguments
  (let [[pt] (sut/lbfgsb (fn [x y] (quad [x y])) {:bounds quad-bounds :vector-arg? false :gradient quad-gradient})]
    (t/is (point-close? [1.0 -2.0] pt)))
  ;; analytic and numerical gradients give the same optimum
  (t/is (point-close? (first (sut/lbfgsb quad {:bounds quad-bounds}))
                      (first (sut/lbfgsb quad {:bounds quad-bounds :gradient quad-gradient})))))

(t/deftest user-gradient-guards
  ;; wrong length: too short, too long, empty, nil
  (doseq [bad [[1.0] [1.0 2.0 3.0] [] nil]]
    (t/is (= {:expected 2 :actual (count bad)}
             (ex-data-of #(sut/lbfgsb quad {:bounds quad-bounds :gradient (fn [_] bad)}))) (pr-str bad)))
  ;; non-finite gradient stops the optimization
  (t/is (thrown? LBFGSBException (sut/lbfgsb quad {:bounds quad-bounds :gradient (fn [_] [##NaN 0.0])})))
  (t/is (thrown? LBFGSBException (sut/lbfgsb quad {:bounds quad-bounds :gradient (fn [_] [##Inf 0.0])})))
  ;; exceptions of the user gradient are not swallowed
  (t/is (thrown-with-msg? RuntimeException #"boom" (sut/lbfgsb quad {:bounds quad-bounds :gradient (fn [_] (throw (RuntimeException. "boom")))}))))

(t/deftest gradient-h
  (let [ref (first (sut/lbfgsb quad {:bounds quad-bounds}))]
    (t/is (point-close? ref (first (sut/lbfgsb quad {:bounds quad-bounds :gradient-h nil}))))
    (t/is (v/delta-eq (vec ref) (vec (first (sut/lbfgsb quad {:bounds quad-bounds :gradient-h 1.0e-3}))) 1.0e-3)
          "a coarse step still finds the optimum (central differences are exact for a quadratic)")
    (t/is (v/delta-eq (vec ref) (vec (first (sut/lbfgsb quad {:bounds quad-bounds :gradient-h 1})))  1.0e-6)))
  (doseq [bad [0 0.0 -1.0e-6 ##NaN]]
    (t/is (contains? (ex-data-of #(sut/lbfgsb quad {:bounds quad-bounds :gradient-h bad})) :gradient-h)
          (pr-str bad))))

;; stats

(t/deftest stats
  (let [s (sut/lbfgsb quad {:bounds quad-bounds :stats? true})]
    (t/is (= #{:point :value :iterations :gradient :status} (set (keys s))))
    (t/is (point-close? [1.0 -2.0] (:point s)))
    (t/is (vector? (:point s)))
    (t/is (integer? (:iterations s)))
    (t/is (pos? (:iterations s)))
    (t/is (= :converged (:status s)))
    (t/is (= 2 (count (:gradient s))))
    (t/is (every? #(m/delta-eq 0.0 % 1.0e-6) (:gradient s)) "gradient at the optimum")
    (t/is (= [(:point s) (:value s)] (sut/lbfgsb quad {:bounds quad-bounds}))))
  ;; falsy values
  (t/is (vector? (sut/lbfgsb quad {:bounds quad-bounds :stats? nil}))))

(t/deftest stats-gradient-sign
  ;; :gradient is always the gradient of f, whatever the goal; on a bound it is not zero
  (let [bounds [[0.0 1.0] [0.0 1.0]]
        mn (sut/lbfgsb (fn [[x y]] (m/+ x y)) {:bounds bounds :stats? true})
        mx (sut/lbfgsb (fn [[x y]] (m/- (m/+ x y))) {:bounds bounds :goal :maximize :stats? true})]
    (t/is (point-close? [0.0 0.0] (:point mn)))
    (t/is (v/delta-eq [1.0 1.0] (:gradient mn) 1.0e-5))
    (t/is (point-close? [0.0 0.0] (:point mx)))
    (t/is (v/delta-eq [-1.0 -1.0] (:gradient mx) 1.0e-5))
    (t/is (m/delta-eq 0.0 (:value mx) 1.0e-9))))

(t/deftest stats-status
  (let [rosen (fn [opts] (sut/lbfgsb p/rosenbrock (merge {:bounds (p/rosenbrock-bounds 2) :initial [-1.2 1.0] :stats? true} opts)))]
    (t/is (= :converged (:status (rosen {}))))
    (t/is (= :max-iterations (:status (rosen {:max-iters 2}))))
    (t/is (= 2 (:iterations (rosen {:max-iters 2}))))
    (t/is (= :max-iterations (:status (rosen {:max-iters 1}))))
    (t/is (= :stalled (:status (rosen {:delta 0.5 :past 1 :abs 1.0e-30 :rel 1.0e-30}))))
    ;; 0 means unlimited
    (t/is (= :converged (:status (rosen {:max-iters 0})))))
  ;; starting at the optimum: no iterations
  (let [s (sut/lbfgsb quad {:bounds quad-bounds :initial [1.0 -2.0] :stats? true})]
    (t/is (= 0 (:iterations s)))
    (t/is (= :converged (:status s)))))

;; failures

(t/deftest non-finite-objective
  (doseq [bad [##NaN ##Inf ##-Inf]
          opts [{} {:max-iters 0} {:linesearch :lbfgsb} {:linesearch :lewis-overton}]]
    (t/is (thrown-with-msg? LBFGSBException #"non-finite"
                            (sut/lbfgsb (fn [_] bad) (merge {:bounds quad-bounds} opts)))
          (pr-str bad opts)))
  ;; not finite only far away from the visited points is not a problem
  (t/is (point-close? [1.0 -2.0] (first (sut/lbfgsb (fn [v] (if (> (v/mag v) 1000.0) ##NaN (quad v))) {:bounds quad-bounds}))))
  ;; exceptions of the objective are not swallowed
  (t/is (thrown-with-msg? RuntimeException #"boom" (sut/lbfgsb (fn [_] (throw (RuntimeException. "boom"))) {:bounds quad-bounds}))))

(t/deftest invalid-options-do-not-start-optimization
  (let [calls (atom 0)
        f (fn [_] (swap! calls inc) 0.0)]
    (doseq [opts [{:m 0} {:linesearch :foo} {:goal :foo} {:gradient-h 0}]]
      (t/is (thrown? clojure.lang.ExceptionInfo (sut/lbfgsb f (merge {:bounds quad-bounds} opts))) (pr-str opts)))
    (t/is (zero? @calls))))

;; options

(t/deftest accepted-options
  (doseq [opts [{:m 1} {:m 20} {:abs 1.0e-12 :rel 1.0e-12} {:past 0} {:past 10} {:delta 0.0} {:max-submin 1} {:max-linesearch 50}
                {:xtol 1.0e-12} {:ftol 1.0e-3 :wolfe 0.8} {:weak-wolfe? false} {:max-step 1.0e5} {:min-step 1.0e-12}
                {:max-iters 0} {:max-iters 10000}]]
    (let [[pt val] (sut/lbfgsb quad (merge {:bounds quad-bounds} opts))]
      (t/is (point-close? [1.0 -2.0] pt) (pr-str opts))
      (t/is (< val val-tol) (pr-str opts)))))

(t/deftest debug-flag
  ;; the Java flag is global and set on every call
  (t/is (not (empty? (with-out-str (sut/lbfgsb quad {:bounds quad-bounds :debug? true})))))
  (t/is (empty? (with-out-str (sut/lbfgsb quad {:bounds quad-bounds :debug? false}))))
  (t/is (empty? (with-out-str (sut/lbfgsb quad {:bounds quad-bounds})))))

;; concurrency

(t/deftest parallel-runs
  ;; a single call shares nothing with the other calls
  (let [results (doall (pmap (fn [i] (sut/lbfgsb quad {:bounds quad-bounds :initial [(- (mod i 7) 3) (- (mod i 5) 2)]})) (range 200)))]
    (t/is (= 200 (count results)))
    (t/is (every? (fn [[pt val]] (and (point-close? [1.0 -2.0] pt) (< val val-tol))) results))))
