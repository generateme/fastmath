(ns fastmath.optimization.sceua-test
  (:require [fastmath.optimization.sceua :as sut]
            [fastmath.optimization.problems :as p]
            [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.stats :as stats]
            [fastmath.matrix :as mat]
            [clojure.test :as t])
  (:import [java.util Arrays]
           [java.util.concurrent.atomic AtomicLong]))

;; References: Himmelblau's function has four minima with value 0: (3, 2), (-2.805118, 3.131313),
;; (-3.779310, -3.283186), (3.584428, -1.848126) (Wikipedia). Sphere, Rosenbrock (minimum 0 at (1, ..., 1)) and Ackley
;; (minimum 0 at the origin) are described at https://www.sfu.ca/~ssurjano/optimization.html. problem02 (sin x + sin(10x/3)
;; on [2.7, 7.5]) has its global minimum -1.899599 at 5.145735 (https://infinity77.net/global_optimization/test_functions_1d.html).
;;
;; The optimizer is stochastic, so all runs use a seeded generator. Tolerances of the values below are the worst error
;; over 30 seeds with the default options, times at least 10: sphere 2D 8e-8 and 5D 7.5e-7 (tolerance 1e-5), Himmelblau
;; 3.7e-6, Rosenbrock 2D 3.6e-7 and 3D 9.7e-7 (tolerance 1e-4), problem02 1.5e-10 for the value and 1.6e-6 for the point
;; (tolerances 1e-8 and 1e-4). The worst of Ackley (2D), 1.5e-3, is limited by :stop-range (1e-3 of the range 65.5 is
;; 0.065), so its tolerance is 1e-2.

(def ^:private himmelblau-minima [[3.0 2.0] [-2.805118 3.131313] [-3.779310 -3.283186] [3.584428 -1.848126]])
(def ^:private hb (p/himmelblau-bounds))
(def ^:private ackley (p/->ackley))

(defn- ex-data-of [thunk]
  (try (thunk) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))

(defn- near-minimum? [pt] (boolean (some #(v/delta-eq (vec pt) % 1.0e-3) himmelblau-minima)))

(defn- run
  "Runs the optimizer with a seeded generator."
  ([f opts] (run f 1 opts))
  ([f seed opts] (sut/sceua f (merge {:rng (r/rng :jdk seed)} opts))))

(defn- never-stop
  "Options which end a run only by :max-iters."
  [max-iters]
  {:stop-range 0.0 :stop-loops 100000 :max-iters max-iters})

;; API

(t/deftest public-api
  (t/is (= #{'sceua} (set (keys (ns-publics 'fastmath.optimization.sceua))))))

;; convergence

(t/deftest convergence
  (doseq [[nm f opts tol] [[:sphere-2d p/sphere {:bounds (p/sphere-bounds 2)} 1.0e-5]
                           [:sphere-5d p/sphere {:bounds (p/sphere-bounds 5)} 1.0e-5]
                           [:himmelblau p/himmelblau {:bounds hb} 1.0e-4]
                           [:rosenbrock-2d p/rosenbrock {:bounds (p/rosenbrock-bounds 2)} 1.0e-4]
                           [:rosenbrock-3d p/rosenbrock {:bounds (p/rosenbrock-bounds 3)} 1.0e-4]
                           [:ackley-2d ackley {:bounds (p/ackley-bounds 2)} 1.0e-2]]
          seed (range 1 6)
          :let [[pt val] (run f seed opts)
                label (str nm " seed " seed)]]
    (t/is (vector? pt) label)
    (t/is (every? double? pt) label)
    (t/is (< (Math/abs ^double val) tol) label)
    (t/is (= val (f pt)) (str label ": the value is the value at the point")))
  (t/testing "the minimum of Himmelblau is one of four"
    (doseq [seed (range 1 6)]
      (t/is (near-minimum? (first (run p/himmelblau seed {:bounds hb}))) (str "seed " seed))))
  (t/testing "known points of the other problems"
    (t/is (v/delta-eq (vec (repeat 3 1.0)) (first (run p/rosenbrock {:bounds (p/rosenbrock-bounds 3)})) 1.0e-2))
    (t/is (v/delta-eq [0.0 0.0] (first (run ackley {:bounds (p/ackley-bounds 2)})) 5.0e-2)))
  (t/testing "one dimension, separate arguments, a flat pair of bounds"
    (doseq [seed (range 1 6)
            :let [[pt val] (run p/problem02 seed {:bounds [2.7 7.5] :vector-arg? false})]]
      (t/is (m/delta-eq 5.145735 (first pt) 1.0e-4) (str "seed " seed))
      (t/is (m/delta-eq -1.899599349 val 1.0e-8) (str "seed " seed)))))

(t/deftest maximize
  (let [[pt val] (run (fn [x] (- (p/sphere x))) {:bounds (p/sphere-bounds 2) :goal :maximize})]
    (t/is (< (- val) 1.0e-5) "the value keeps the sign of the function")
    (t/is (v/delta-eq [0.0 0.0] pt 1.0e-2)))
  (t/testing "the maximum in a corner of the bounds is reached, points are clipped"
    (let [[pt val] (run p/sphere {:bounds (p/sphere-bounds 2) :goal :maximize})]
      (t/is (m/delta-eq (* 2 5.12 5.12) val 1.0e-3))
      (t/is (every? #(m/delta-eq 5.12 (Math/abs ^double %) 1.0e-3) pt))))
  (t/testing "a minimum and a maximum of the same function"
    (let [[_ vmin] (run p/problem02 {:bounds [[2.7 7.5]] :vector-arg? false})
          [_ vmax] (run p/problem02 {:bounds [[2.7 7.5]] :vector-arg? false :goal :maximize})]
      (t/is (< vmin 0.0 vmax))))
  (t/testing "goal nil is minimize, an unknown goal throws"
    (t/is (= (run p/sphere {:bounds (p/sphere-bounds 2)}) (run p/sphere {:bounds (p/sphere-bounds 2) :goal nil})))
    (t/is (= {:goal :foo :allowed #{:minimize :maximize}} (ex-data-of #(run p/sphere {:bounds (p/sphere-bounds 2) :goal :foo}))))))

;; result

(t/deftest result-shapes
  (let [opts {:bounds hb}
        res (run p/himmelblau opts)
        s (run p/himmelblau (assoc opts :stats? true))]
    (t/is (vector? res))
    (t/is (= 2 (count res)))
    (t/is (vector? (first res)))
    (t/is (double? (second res)))
    (t/is (= #{:point :value :evaluations :iterations :status :complexes} (set (keys s))))
    (t/is (= res [(:point s) (:value s)]) "the same run gives the same result in both forms")
    (t/is (vector? (:point s)))
    (t/is (int? (:evaluations s)))
    (t/is (pos? (:iterations s)))
    (t/is (#{:converged-range :converged-improvement :max-iterations} (:status s)))
    (t/is (= 5 (:complexes s)))
    (t/is (= res (run p/himmelblau (assoc opts :stats? nil))) ":stats? nil is false")
    (t/is (= res (run p/himmelblau (assoc opts :stats? false)))))
  (t/testing "one dimension gives a vector point"
    (t/is (vector? (first (run (fn [[x]] (* x x)) {:bounds [[-1 1]]}))))
    (t/is (= 1 (count (first (run (fn [[x]] (* x x)) {:bounds [[-1 1]]})))))))

;; initial point and population

(t/deftest initial-point
  (let [sphere-opts {:bounds (p/sphere-bounds 2)}]
    (t/testing "the result is never worse than the initial point"
      (doseq [seed (range 1 6)
              initial [[0.0 0.0] [1.0 1.0] [-5.12 5.12] [5.12 5.12] [4.0 -3.0]]
              :let [[_ val] (run p/sphere seed (assoc sphere-opts :initial initial :max-iters 1))]]
        (t/is (<= val (p/sphere initial)) (str seed " " initial))))
    (t/testing "the optimum as the initial point is kept"
      (t/is (= 0.0 (second (run p/sphere (assoc sphere-opts :initial [0.0 0.0] :max-iters 1))))))
    (t/testing "the number of one dimension, any sequence of numbers"
      (doseq [initial [2.0 [2.0] '(2.0) (double-array [2.0])]]
        (t/is (<= (second (run (fn [[x]] (* x x)) {:bounds [[-5 5]] :initial initial :max-iters 1})) 4.0) (pr-str initial))))
    (t/testing "invalid initial points"
      (doseq [bad [[6.0 0.0] [0.0 -6.0] [##NaN 0.0] [##Inf 0.0] ["a" 0.0] [nil 0.0] [:a :b]]]
        (t/is (= :initial (:option (ex-data-of #(run p/sphere (assoc sphere-opts :initial bad))))) (pr-str bad)))
      (doseq [bad [[1.0] [1.0 2.0 3.0] []]]
        (t/is (= :sceua (:method (ex-data-of #(run p/sphere (assoc sphere-opts :initial bad))))) (pr-str bad))))
    (t/testing "a point on the bound is inside the bounds"
      (t/is (some? (run p/sphere (assoc sphere-opts :initial [5.12 -5.12])))))))

(t/deftest initial-population
  (let [ev (#'sut/evaluator p/sphere 1.0 100000)
        lo (double-array [-1.0 -2.0])
        hi (double-array [1.0 2.0])
        build (fn [initial size jitter seed]
                (#'sut/initial-population {:lo lo :hi hi :initial initial :evaluate ev :jitter jitter
                                           :rng (r/rng :jdk seed) :complexes 1 :complex-size size}))]
    (doseq [size [1 2 25 100]
            initial [nil (double-array [0.5 0.5])]
            jitter [0.0 0.25 1.0]
            :let [pop (build initial size jitter 1)
                  label (str size " " (some? initial) " " jitter)]]
      (t/is (vector? pop) label)
      (t/is (= size (count pop)) label)
      (t/is (every? #(= 2 (count (:x %))) pop) label)
      (t/is (apply <= (map :value pop)) (str label ": sorted by value"))
      (t/is (every? (fn [{[x y] :x}] (and (<= -1.0 x 1.0) (<= -2.0 y 2.0))) pop) (str label ": in the bounds"))
      (t/is (every? #(= (:value %) (p/sphere (:x %))) pop) (str label ": values"))
      (when initial
        (t/is (some #(= [0.5 0.5] (vec (:x %))) pop) (str label ": the initial point is there"))))
    (t/testing "the generator decides the points"
      (t/is (= (map (comp vec :x) (build nil 10 0.25 1)) (map (comp vec :x) (build nil 10 0.25 1))))
      (t/is (not= (map (comp vec :x) (build nil 10 0.25 1)) (map (comp vec :x) (build nil 10 0.25 2)))))
    (t/testing "more than 14 dimensions use a Sobol sequence"
      (doseq [n [14 15 16]
              :let [pop (#'sut/initial-population {:lo (double-array (repeat n 0.0)) :hi (double-array (repeat n 1.0))
                                                   :evaluate (#'sut/evaluator (fn [_] 0.0) 1.0 1000) :jitter 0.25
                                                   :rng (r/rng :jdk 1) :complexes 1 :complex-size 10})]]
        (t/is (= 10 (count pop)) (str n))
        (t/is (every? #(= n (count (:x %))) pop) (str n))))))

(t/deftest evaluator
  (let [x (double-array [1.0 2.0])]
    (t/testing "a point keeps the array and the value, calls are counted"
      (let [ev (#'sut/evaluator p/sphere 1.0 10)
            point (ev x)]
        (t/is (zero? ((#'sut/evaluator p/sphere 1.0 10))) "no calls yet")
        (t/is (identical? x (:x point)))
        (t/is (= 5.0 (:value point)))
        (t/is (= 1 (ev)))
        (ev x)
        (t/is (= 2 (ev)))))
    (t/testing "the sign turns a maximum into a minimum"
      (t/is (= -5.0 (:value ((#'sut/evaluator p/sphere -1.0 10) x)))))
    (t/testing "NaN is the worst value, infinities are kept"
      (t/is (= ##Inf (:value ((#'sut/evaluator (fn [_] ##NaN) 1.0 10) x))))
      (t/is (= ##Inf (:value ((#'sut/evaluator (fn [_] ##NaN) -1.0 10) x))))
      (t/is (= ##-Inf (:value ((#'sut/evaluator (fn [_] ##-Inf) 1.0 10) x))))
      (t/is (= 3.0 (:value ((#'sut/evaluator (fn [_] 3) 1.0 10) x))) "an integer becomes a double"))
    (t/testing "the call number max-evals plus one throws"
      (doseq [limit [1 2 10]
              :let [ev (#'sut/evaluator p/sphere 1.0 limit)]]
        (dotimes [_ limit] (ev x))
        (t/is (= limit (ev)))
        (t/is (= {:max-evals limit :evaluations limit} (ex-data-of #(ev x))) (str limit))))
    (t/testing "the counter is exact in many threads"
      (let [ev (#'sut/evaluator p/sphere 1.0 1000000)]
        (doall (pmap (fn [_] (dotimes [_ 1000] (ev x))) (range 8)))
        (t/is (= 8000 (ev)))))))

;; limits

(t/deftest evaluation-limit
  (let [opts {:bounds hb}
        s (run p/himmelblau (assoc opts :stats? true))
        e (:evaluations s)]
    (t/testing "exactly the needed number of evaluations is enough, one less throws"
      (t/is (= (:point s) (:point (run p/himmelblau (assoc opts :stats? true :max-evals e)))))
      (t/is (= {:max-evals (dec e) :evaluations (dec e)} (ex-data-of #(run p/himmelblau (assoc opts :max-evals (dec e)))))))
    (t/testing "limits of the first population (25 points) and below"
      (doseq [limit [1 2 24 25 26]]
        (t/is (= {:max-evals limit :evaluations limit} (ex-data-of #(run p/himmelblau (assoc opts :max-evals limit)))) (str limit))))
    (t/testing "a huge limit"
      (t/is (some? (run p/himmelblau (assoc opts :max-evals Long/MAX_VALUE)))))
    (t/testing "the exception is ex-info in parallel mode as well"
      (let [e (try (run p/himmelblau (assoc opts :parallel? true :max-evals 300)) (catch Throwable t t))]
        (t/is (instance? clojure.lang.ExceptionInfo e))
        (t/is (= {:max-evals 300 :evaluations 300} (ex-data e)))))))

(t/deftest iteration-limit
  (doseq [k [1 2 5 20]
          :let [s (run p/himmelblau (merge {:bounds hb :stats? true} (never-stop k)))]]
    (t/is (= :max-iterations (:status s)) (str k))
    (t/is (= k (:iterations s)) (str k)))
  (t/testing "the best point of the last loop is returned, not an error"
    (t/is (vector? (run p/himmelblau (merge {:bounds hb} (never-stop 1)))))))

(t/deftest stop-criteria
  (t/testing "a range of the population within the bounds stops after the first loop"
    (let [s (run p/himmelblau {:bounds hb :stats? true :stop-range 1.0})]
      (t/is (= :converged-range (:status s)))
      (t/is (= 1 (:iterations s)))))
  (t/testing "a huge improvement tolerance stops after :stop-loops loops"
    (doseq [loops [1 2 5]
            :let [s (run p/himmelblau {:bounds hb :stats? true :stop-range 0.0 :stop-loops loops :stop-improvement 1.0e300})]]
      (t/is (= :converged-improvement (:status s)) (str loops))
      (t/is (= loops (:iterations s)) (str loops))))
  (t/testing "the range test goes first when both hold"
    (t/is (= :converged-range (:status (run p/himmelblau {:bounds hb :stats? true :stop-range 1.0 :stop-improvement 1.0e300 :stop-loops 1})))))
  (t/testing "zero tolerances: a constant function never improves, so it stops by :stop-loops"
    (let [s (run (fn [_] 1) {:bounds hb :stats? true :stop-range 0.0 :stop-improvement 0.0 :stop-loops 3})]
      (t/is (= :converged-improvement (:status s)))
      (t/is (= 3 (:iterations s)))
      (t/is (= 1.0 (:value s)) "an integer value becomes a double")))
  (t/testing "the default patience is 40 loops, and it is enough for Himmelblau (7 loops were not: about 10% of runs stopped early)"
    (let [s (run (fn [_] 1) {:bounds hb :stats? true :stop-range 0.0})]
      (t/is (= :converged-improvement (:status s)))
      (t/is (= 40 (:iterations s))))
    (doseq [seed (range 1 41)]
      (t/is (< (second (run p/himmelblau seed {:bounds hb})) 1.0e-3) (str "seed " seed))))
  (t/testing "at least one loop is always made"
    (t/is (= 1 (:iterations (run p/himmelblau (merge {:bounds hb :stats? true} (never-stop 1)))))))
  (t/testing "a tighter range gives a more accurate result"
    (let [val-of (fn [range-tol] (second (run p/himmelblau {:bounds hb :stop-range range-tol :stop-loops 100})))]
      (t/is (< (val-of 1.0e-6) 1.0e-8))
      (t/is (< (val-of 1.0e-6) (val-of 1.0e-1))))))

;; complexes

(t/deftest complex-reduction
  (doseq [[complexes minimum k expected] [[4 2 5 2] [4 4 5 4] [4 3 2 3] [3 1 10 1] [1 1 3 1] [2 1 1 2] [6 2 3 4] [4 nil 5 4]]
          :let [s (run p/himmelblau (merge {:bounds hb :stats? true :complexes complexes :min-complexes minimum} (never-stop k)))]]
    (t/is (= expected (:complexes s)) (str complexes " " minimum " " k)))
  (t/testing "the private reduction drops the worst complex"
    (let [pop (vec (range 15))
          reduced (fn [complex-size minimum] (#'sut/reduce-complexes pop {:complex-size complex-size :min-complexes minimum}))]
      (t/is (= (vec (range 10)) (reduced 5 2)) "the worst points are the last ones")
      (t/is (= pop (reduced 5 3)) "at the minimum nothing is removed")
      (t/is (= (vec (range 10)) (reduced 5 1)) "one complex is removed")
      (t/is (= pop (reduced 15 1)) "the only complex stays")))
  (t/testing "a run with reduction is valid and cheaper per loop"
    (let [[_ val] (run p/himmelblau {:bounds hb :complexes 6 :min-complexes 2})]
      (t/is (double? val)))))

(t/deftest minimal-sizes
  (doseq [opts [{:complexes 1 :complex-size 2 :subcomplex-size 2 :evolution-steps 1}
                {:complexes 1 :complex-size 3}
                {:complexes 1}
                {:complex-size 2}
                {:complex-size 2 :subcomplex-size 2}
                {:evolution-steps 1}
                {:subcomplex-size 2}
                {:subcomplex-size 5 :complex-size 5}
                {:complexes 20 :complex-size 3}]
          bounds [[[2.7 7.5]] hb (p/sphere-bounds 4)]
          :let [f (if (= 1 (count bounds)) (fn [[x]] (p/problem02 x)) (fn [x] (reduce + (map #(* % %) x))))
                s (run f (merge {:bounds bounds :stats? true :max-iters 30} opts))
                label (str opts " " (count bounds))]]
    (t/is (every? (fn [[x [lo hi]]] (<= lo x hi)) (map vector (:point s) bounds)) (str label ": in the bounds"))
    (t/is (= (:value s) (f (:point s))) label)
    (t/is (pos? (:evaluations s)) label)))

;; determinism, random generators

(t/deftest generators
  (let [f p/himmelblau
        opts {:bounds hb :stats? true}]
    (t/testing "equal seeds give equal results, different seeds differ"
      (doseq [parallel? [false true]
              pca? [false true]]
        (let [with-seed (fn [seed] (run f seed (assoc opts :parallel? parallel? :pca-recovery? pca?)))]
          (t/is (= (with-seed 1) (with-seed 1)) (str parallel? pca?))
          (t/is (not= (with-seed 1) (with-seed 2)) (str parallel? pca?)))))
    (t/testing "other generators, a synchronized one too"
      (doseq [make [#(r/rng :mersenne 3) #(r/rng :isaac 3) #(r/synced-rng :well512a 3)]
              parallel? [false true]]
        (let [res (fn [] (sut/sceua f (assoc opts :rng (make) :parallel? parallel?)))]
          (t/is (= (res) (res)) (str (class (make)) parallel?)))))
    (t/testing "repeated parallel runs with a generator which is not thread safe give the same result"
      (let [res (fn [] (sut/sceua f (assoc opts :rng (r/rng :mersenne 9) :parallel? true :complexes 8)))
            first-res (res)]
        (dotimes [_ 5] (t/is (= first-res (res))))))
    (t/testing "a missing and a nil generator create a new one"
      (t/is (map? (sut/sceua f opts)))
      (t/is (map? (sut/sceua f (assoc opts :rng nil)))))
    (t/testing "a value which is not a generator throws"
      (doseq [bad [5 :jdk "x" [1] (r/distribution :normal)]]
        (t/is (= {:rng bad} (ex-data-of #(sut/sceua f (assoc opts :rng bad)))) (pr-str bad))))
    (t/testing "the shared generator is not used"
      (doseq [extra [{} {:parallel? true} {:pca-recovery? true}]]
        (r/set-seed! 1)
        (let [expected (r/drand)]
          (r/set-seed! 1)
          (sut/sceua f (merge opts extra {:rng (r/rng :jdk 1)}))
          (sut/sceua f (merge opts extra))
          (t/is (= expected (r/drand)) (str extra)))))))

(t/deftest child-generators
  ;; parallel mode asks fastmath.random/child-rngs for one generator per complex and loop, sequential mode does not
  (let [calls (atom [])
        child-rngs r/child-rngs
        ;; the primitive hint is needed: sceua calls child-rngs with a primitive long
        record-calls (fn [rng ^long n] (swap! calls conj n) (child-rngs rng n))]
    (with-redefs [r/child-rngs record-calls]
      (run p/himmelblau (merge {:bounds hb :parallel? true :complexes 4} (never-stop 3)))
      (t/is (= [4 4 4] @calls))
      (reset! calls [])
      (run p/himmelblau (merge {:bounds hb :parallel? true :complexes 4 :min-complexes 2} (never-stop 4)))
      (t/is (= [4 3 2 2] @calls) "the number of complexes shrinks")
      (reset! calls [])
      (run p/himmelblau (merge {:bounds hb :parallel? true :complexes 1} (never-stop 2)))
      (t/is (= [1 1] @calls) "one complex")
      (reset! calls [])
      (run p/himmelblau (merge {:bounds hb :complexes 4} (never-stop 3)))
      (t/is (empty? @calls) "sequential mode uses the generator itself"))))

;; evaluations

(t/deftest evaluation-count
  (doseq [extra [{} {:parallel? true} {:pca-recovery? true} {:goal :maximize} {:complexes 4 :min-complexes 2} {:initial [1.0 1.0]}
                 {:vector-arg? false} {:max-iters 1}]
          :let [calls (AtomicLong. 0)
                f (if (false? (:vector-arg? extra))
                    (fn [x y] (.incrementAndGet calls) (p/himmelblau [x y]))
                    (fn [pt] (.incrementAndGet calls) (p/himmelblau pt)))
                s (run f (merge {:bounds hb :stats? true} extra))]]
    (t/is (= (.get calls) (:evaluations s)) (str extra))
    (t/is (<= 25 (:evaluations s)) (str extra ": at least the first population"))))

(t/deftest objective-contract
  (t/testing "the function gets a new double array of the coordinates for every call"
    (let [seen (atom [])]
      (run (fn [x] (swap! seen conj x) (p/himmelblau x)) {:bounds hb :max-iters 2})
      (t/is (every? #(instance? (Class/forName "[D") %) @seen))
      (t/is (every? #(= 2 (count %)) @seen))
      (t/is (= (count @seen) (count (set (map #(System/identityHashCode %) @seen)))))))
  (t/testing "the sequence functions and destructuring work on the argument"
    (t/is (< (second (run (fn [[x y]] (+ (* x x) (* y y))) {:bounds [[-1 1] [-1 1]]})) 1.0e-5))
    (t/is (< (second (run (fn [x] (reduce + (map #(* % %) x))) {:bounds [[-1 1] [-1 1] [-1 1]]})) 1.0e-5))
    (t/is (< (second (run (fn [x] (Math/abs ^double (nth x 0))) {:bounds [[-1 1]]})) 1.0e-3)))
  (t/testing "NaN and infinite values on a part of the domain: the optimum of the rest is found"
    (doseq [bad [##NaN ##Inf]
            :let [[pt val] (run (fn [[x y]] (if (< x 0.0) bad (p/himmelblau [x y]))) {:bounds hb})]]
      (t/is (<= 0.0 (first pt)) (str bad))
      (t/is (< val 1.0e-4) (str bad))))
  (t/testing "a function which is NaN or infinite everywhere never converges, the evaluation limit ends the run"
    (doseq [bad [##NaN ##Inf]]
      (t/is (= {:max-evals 600 :evaluations 600} (ex-data-of #(run (fn [_] bad) {:bounds hb :max-evals 600}))) (str bad))))
  (t/testing "the objective may return any number"
    (t/is (= 1.0 (second (run (fn [_] 1) {:bounds hb :max-iters 2}))))
    (t/is (= 2.5 (second (run (fn [_] 5/2) {:bounds hb :max-iters 2})))))
  (t/testing "an exception of the function goes out as it is, also from the threads of the parallel mode"
    (t/is (thrown? ArithmeticException (run (fn [_] (/ 1 0)) {:bounds hb})))
    (doseq [parallel? [false true]
            :let [calls (atom 0)
                  f (fn [x] (if (> (swap! calls inc) 40) (throw (ArithmeticException. "boom")) (p/himmelblau x)))]]
      (t/is (thrown-with-msg? ArithmeticException #"boom" (run f {:bounds hb :parallel? parallel?})) (str parallel?))))
  (t/testing "wrong arity of the function"
    (t/is (thrown? clojure.lang.ArityException (run (fn [_x _y _z] 0.0) {:bounds hb :vector-arg? false})))))

(t/deftest best-value-never-increases
  (doseq [extra [{} {:pca-recovery? true} {:complexes 4 :min-complexes 2} {:parallel? false :complexes 1}]
          f [p/himmelblau p/sphere]
          :let [vals (mapv #(second (run f 3 (merge {:bounds hb} extra (never-stop %)))) (range 1 9))]]
    (t/is (apply >= vals) (str extra " " vals))))

;; validation

(t/deftest option-validation
  (doseq [option [:complexes :complex-size :subcomplex-size :evolution-steps :min-complexes :stop-loops :max-evals :max-iters]
          bad [0 -1 1.5 2.0 "3" ##NaN true :a [1] 1/2]]
    (let [d (ex-data-of #(run p/himmelblau {:bounds hb option bad}))]
      (t/is (= option (:option d)) (str option " " (pr-str bad)))
      (t/is (= (pr-str bad) (pr-str (:value d))) (str option " " (pr-str bad)))
      (t/is (string? (:reason d)) (str option " " (pr-str bad)))))
  (doseq [option [:stop-improvement :stop-range :jitter]
          bad [-1 -1.0e-300 ##NaN ##Inf ##-Inf "a" :a [1] true]]
    (t/is (= option (:option (ex-data-of #(run p/himmelblau {:bounds hb option bad})))) (str option " " (pr-str bad))))
  (t/testing "relations between the options"
    (t/is (= :complex-size (:option (ex-data-of #(run p/himmelblau {:bounds hb :complex-size 1})))))
    (t/is (= :subcomplex-size (:option (ex-data-of #(run p/himmelblau {:bounds hb :subcomplex-size 1})))))
    (t/is (= :subcomplex-size (:option (ex-data-of #(run p/himmelblau {:bounds hb :subcomplex-size 6 :complex-size 5})))))
    (t/is (= :min-complexes (:option (ex-data-of #(run p/himmelblau {:bounds hb :min-complexes 6})))))
    (t/is (= :min-complexes (:option (ex-data-of #(run p/himmelblau {:bounds hb :complexes 2 :min-complexes 3})))))
    (t/is (= :jitter (:option (ex-data-of #(run p/himmelblau {:bounds hb :jitter 1.0000001}))))))
  (t/testing "limit values are valid"
    (doseq [extra [{:complexes 1} {:complex-size 2} {:subcomplex-size 2} {:subcomplex-size 5 :complex-size 5}
                   {:min-complexes 5} {:jitter 0.0} {:jitter 1.0} {:stop-range 0.0} {:stop-improvement 0.0}
                   {:stop-loops 1} {:max-iters 1} {:evolution-steps 1} {:stop-range 1.0e300} {:jitter 1} {:stop-range 1}
                   {:complexes (int 2)} {:complexes (short 2)} {:max-evals Long/MAX_VALUE}]]
      (t/is (some? (run p/himmelblau (merge {:bounds hb :max-iters 3} extra))) (pr-str extra))))
  (t/testing "an explicit nil is the default"
    (doseq [option [:complexes :complex-size :subcomplex-size :evolution-steps :min-complexes :stop-loops
                    :stop-improvement :stop-range :max-evals :max-iters :jitter :parallel? :pca-recovery? :stats? :vector-arg? :goal]]
      (t/is (= (run p/himmelblau {:bounds hb})
               (run p/himmelblau {:bounds hb option nil})) (str option)))))

(t/deftest bounds-validation
  (t/is (= :sceua (:method (ex-data-of #(sut/sceua p/himmelblau nil)))))
  (t/is (= :sceua (:method (ex-data-of #(sut/sceua p/himmelblau {})))))
  (doseq [bounds [[[##-Inf 5] [-5 5]] [[-5 ##Inf] [-5 5]] [[##-Inf ##Inf]] [[1 1] [-5 5]] [[-5 5] [2 2]] [[5 -5] [-5 5]]
                  [[##NaN 5]] [[nil 5]] [] [[]] [[1]] [[1 2 3]] 5 "ab" [["a" "b"]]]]
    (t/is (= :sceua (:method (ex-data-of #(sut/sceua p/himmelblau {:bounds bounds})))) (pr-str bounds)))
  (t/testing "more than 1000 dimensions"
    (let [d (ex-data-of #(sut/sceua (fn [_] 0.0) {:bounds (repeat 1001 [0 1])}))]
      (t/is (= :bounds (:option d)))))
  (t/testing "a tiny range of a dimension is valid"
    (t/is (some? (run (fn [[x _]] x) {:bounds [[0 Double/MIN_VALUE] [0 1]] :max-iters 2}))))
  (t/testing "bounds of integers and ratios, lazy sequences"
    (t/is (some? (run p/sphere {:bounds '((-5 5) (-1/2 1/2)) :max-iters 2})))
    (t/is (some? (run p/sphere {:bounds (map vector [-5 -5] [5 5]) :max-iters 2})))))

(t/deftest dimensions
  (t/testing "from one dimension up to the switch of the sequence (14 and 15 dimensions) and above"
    (doseq [n [1 2 3 14 15 16]
            :let [s (run p/sphere {:bounds (p/sphere-bounds n) :stats? true :max-iters 2})]]
      (t/is (= n (count (:point s))) (str n))
      (t/is (every? #(<= -5.12 % 5.12) (:point s)) (str n)))))

;; competitive complex evolution

(defn- points-of
  "A population of points of the function `f` at the coordinates, sorted by value, and the evaluator which made it."
  [f coords]
  (let [ev (#'sut/evaluator f 1.0 1000000)]
    [(#'sut/sorted-by-value (map #(ev (double-array %)) coords)) ev]))

(defn- values-of [points] (map :value points))

(defn- select-parents-reference
  "The previous implementation of parent selection: indices drawn into a sorted set until there are `q` different ones."
  [m q rng]
  (vec (loop [chosen (sorted-set)]
         (if (== (count chosen) q)
           chosen
           (recur (conj chosen (#'sut/triangular-index m (r/drandom rng))))))))

(defn- triangular-index-reference
  "The previous implementation of the triangular index, with `m/floor`."
  [^long m ^double u]
  (let [mh (+ m 0.5)
        root (m/safe-sqrt (- (* mh mh) (* m (inc m) u)))]
    (m/constrain (long (m/floor (- mh root))) 0 (dec m))))

(defn- centroid-reference
  "The previous implementation of the centroid of the parents except the last one."
  [complex parents]
  (v/average-vectors (map #(:x (complex %)) (butlast parents))))

(defn- random-point-reference
  "The previous implementation of the random point of a box."
  [lo hi rng]
  (v/einterpolate lo hi (double-array (repeatedly (count lo) #(r/drandom rng)))))

(t/deftest parent-selection
  (t/testing "the triangular index: the best index is the most probable one"
    (doseq [m [1 2 7 100]]
      (t/is (= 0 (#'sut/triangular-index m 0.0)) (str m))
      (t/is (= (dec m) (#'sut/triangular-index m 0.9999999999999999)) (str m))
      (t/is (every? #(<= 0 (#'sut/triangular-index m %) (dec m)) (range 0.0 1.0 0.001)) (str m))
      (t/is (apply <= (map #(#'sut/triangular-index m %) (range 0.0 1.0 0.001))) (str m ": monotone"))))
  (t/testing "frequencies of the indices follow the triangular distribution"
    (let [m 7
          n 20000
          rng (r/rng :jdk 11)
          counts (frequencies (repeatedly n #(first (#'sut/select-parents m 1 rng))))]
      (doseq [i (range m)]
        (t/is (m/delta-eq (/ (* 2.0 (- m i)) (* m (inc m))) (/ (get counts i 0) (double n)) 0.02) (str i)))))
  (t/testing "different indices in the ascending order, for any size of the sub-complex"
    (let [rng (r/rng :jdk 5)]
      (doseq [m [2 3 7 20]
              q (range 1 (inc m))
              _ (range 20)
              :let [idx (#'sut/select-parents m q rng)]]
        (t/is (instance? (Class/forName "[J") idx) (str m " " q))
        (t/is (= q (count idx)) (str m " " q))
        (t/is (apply < idx) (str m " " q " " idx))
        (t/is (every? #(<= 0 % (dec m)) idx) (str m " " q " " idx)))))
  (t/testing "the whole complex is selected when the sub-complex is the complex"
    (t/is (= [0 1 2 3 4] (vec (#'sut/select-parents 5 5 (r/rng :jdk 1)))))
    (t/is (= [0] (vec (#'sut/select-parents 1 1 (r/rng :jdk 1))))))
  (t/testing "the same indices and the same draws as the sorted set of the previous implementation"
    (doseq [m [2 3 5 11 41]
            q (range 2 (inc m))
            seed (range 1 21)
            :let [rng-a (r/rng :jdk seed)
                  rng-b (r/rng :jdk seed)
                  idx (vec (#'sut/select-parents m q rng-a))]]
      (t/is (= (select-parents-reference m q rng-b) idx) (str m " " q " " seed))
      (t/is (= (r/drandom rng-a) (r/drandom rng-b)) (str m " " q " " seed ": the same number of draws")))))

(t/deftest complexes-of-population
  (t/testing "a complex takes every ngs-th point starting at k"
    (let [pop (vec (range 15))]
      (t/is (= [0 3 6 9 12] (#'sut/complex-of pop 0 3)))
      (t/is (= [1 4 7 10 13] (#'sut/complex-of pop 1 3)))
      (t/is (= [2 5 8 11 14] (#'sut/complex-of pop 2 3)))
      (t/is (= pop (#'sut/complex-of pop 0 1)) "one complex is the population")
      (t/is (vector? (#'sut/complex-of pop 0 3)))))
  (t/testing "the complexes are evolved and merged into a sorted population of the same size"
    (let [[pop ev] (points-of p/sphere (for [x (range -3.0 3.1 1.0) y (range -3.0 3.1 1.0)] [x y]))
          settings {:complex-size 7 :subcomplex-size 3 :evolution-steps 5 :evaluate ev
                    :lo (double-array [-4.0 -4.0]) :hi (double-array [4.0 4.0])}]
      (t/is (= 49 (count pop)))
      (doseq [parallel? [false true]
              :let [evolved (#'sut/evolve-population pop (assoc settings :parallel? parallel?) (r/rng :jdk 1))]]
        (t/is (vector? evolved) (str parallel?))
        (t/is (= 49 (count evolved)) (str parallel?))
        (t/is (apply <= (values-of evolved)) (str parallel? ": sorted"))
        (t/is (<= (:value (first evolved)) (:value (first pop))) (str parallel? ": the best point is not worse"))))))

(t/deftest evolve-complex
  (let [[complex ev] (points-of p/sphere [[1.0 1.0] [2.0 0.0] [-1.0 2.0] [0.0 -3.0] [3.0 3.0]])
        before (mapv (comp vec :x) complex)
        settings {:subcomplex-size 3 :evolution-steps 20 :evaluate ev
                  :lo (double-array [-4.0 -4.0]) :hi (double-array [4.0 4.0])}
        evolved (#'sut/evolve-complex complex settings (r/rng :jdk 1))]
    (t/is (= before (mapv (comp vec :x) complex)) "the coordinates of the input are not changed")
    (t/is (vector? evolved))
    (t/is (= 5 (count evolved)))
    (t/is (apply <= (values-of evolved)) "sorted by value")
    (t/is (<= (:value (first evolved)) (:value (first complex))) "the best point is not worse")
    (t/is (every? (fn [{[x y] :x}] (and (<= -4.0 x 4.0) (<= -4.0 y 4.0))) evolved) "in the bounds")
    (t/is (every? #(= (:value %) (p/sphere (:x %))) evolved) "values are values of the function")
    (t/testing "a sub-complex of the size of the complex and one step"
      (let [one-step (#'sut/evolve-complex complex (assoc settings :subcomplex-size 5 :evolution-steps 1) (r/rng :jdk 1))]
        (t/is (= 5 (count one-step)))
        (t/is (apply <= (values-of one-step)))))
    (t/testing "a sub-complex of two points, the smallest one"
      (let [two (#'sut/evolve-complex complex (assoc settings :subcomplex-size 2) (r/rng :jdk 1))]
        (t/is (= 5 (count two)))
        (t/is (apply <= (values-of two)))))
    (t/testing "steps use one to three evaluations"
      (let [counting (#'sut/evaluator p/sphere 1.0 1000000)
            steps 30
            _ (#'sut/evolve-complex complex (assoc settings :evaluate counting :evolution-steps steps) (r/rng :jdk 2))]
        (t/is (<= steps (counting) (* 3 steps)))))
    (t/testing "the reflection is clipped to the bounds"
      (let [[cx cev] (points-of p/sphere [[0.9 0.9] [1.0 1.0] [-1.0 -1.0] [0.95 -0.9] [0.0 0.0]])
            tight {:subcomplex-size 2 :evolution-steps 50 :evaluate cev :lo (double-array [-1.0 -1.0]) :hi (double-array [1.0 1.0])}
            out (#'sut/evolve-complex cx tight (r/rng :jdk 3))]
        (t/is (every? (fn [{[x y] :x}] (and (<= -1.0 x 1.0) (<= -1.0 y 1.0))) out))))
    (t/testing "a complex of equal points stays valid"
      (let [[same sev] (points-of p/sphere (repeat 5 [1.0 1.0]))
            out (#'sut/evolve-complex same (assoc settings :evaluate sev) (r/rng :jdk 4))]
        (t/is (= 5 (count out)))
        (t/is (every? #(= 2.0 (:value %)) (take 1 out)))))))

;; hot path: the same results as the previous implementations

(t/deftest hot-path-equivalence
  (t/testing "the triangular index equals the version with floor"
    (let [rng (r/rng :jdk 7)]
      (doseq [m [1 2 3 7 41 100 2001]
              :let [us (concat [0.0 0.5 0.9999999999999999] (repeatedly 150000 #(r/drandom rng)))]]
        (t/is (every? (fn [u] (= (triangular-index-reference m u) (#'sut/triangular-index m u))) us) (str m)))))
  (t/testing "the centroid equals v/average-vectors exactly, also with signed zeros"
    (let [rng (r/rng :jdk 3)
          ev (#'sut/evaluator (fn [_] 0.0) 1.0 1000000)
          random-x (fn [n] (double-array (repeatedly n #(* 10.0 (- (r/drandom rng) 0.5)))))]
      (doseq [size [2 3 5 11]
              n [1 2 5 20]
              q (range 2 (inc size))
              seed (range 1 4)
              :let [complex (mapv (fn [_] (ev (random-x n))) (range size))
                    complex (if (= seed 3) (assoc complex 0 (ev (double-array (take n (cycle [-0.0 0.0 1.0]))))) complex)
                    parents (#'sut/select-parents size q (r/rng :jdk seed))
                    expected (centroid-reference complex parents)]]
        (t/is (Arrays/equals ^doubles (#'sut/centroid complex parents q) ^doubles expected) (str size " " n " " q " " seed)))))
  (t/testing "the random point equals the previous one and draws the same numbers"
    (doseq [n [1 2 5 20]
            seed (range 1 101)
            :let [lo (double-array (take n (cycle [-1.0 0.0 2.5])))
                  hi (double-array (take n (cycle [1.0 3.0 2.6])))
                  rng-a (r/rng :jdk seed)
                  rng-b (r/rng :jdk seed)]]
      (t/is (Arrays/equals ^doubles (#'sut/random-point lo hi rng-a) ^doubles (random-point-reference lo hi rng-b)) (str n " " seed))
      (t/is (= (r/drandom rng-a) (r/drandom rng-b)) (str n " " seed ": the same number of draws")))))

;; evaluations skipped for a candidate equal to the worst point

(t/deftest skipped-evaluations
  (t/testing "in a complex collapsed in one point only the random point is evaluated, one evaluation per step"
    (doseq [subcomplex-size [2 3 5]
            :let [[same _] (points-of p/sphere (repeat 5 [1.0 1.0]))
                  counting (#'sut/evaluator p/sphere 1.0 1000000)
                  steps 30
                  out (#'sut/evolve-complex same {:subcomplex-size subcomplex-size :evolution-steps steps :evaluate counting
                                                  :lo (double-array [-4.0 -4.0]) :hi (double-array [4.0 4.0])}
                                            (r/rng :jdk 4))]]
      (t/is (= steps (counting)) (str subcomplex-size))
      (t/is (= 5 (count out)) (str subcomplex-size))))
  (t/testing "an optimum in a corner of the bounds needs fewer evaluations"
    ;; the previous implementation took 16451 evaluations for these seeds
    (let [runs (for [seed (range 1 11)]
                 (run p/sphere seed {:bounds (p/sphere-bounds 2) :goal :maximize :stats? true}))]
      (t/is (<= (reduce + (map :evaluations runs)) 12000))
      (t/is (every? #(m/delta-eq (* 2 5.12 5.12) (:value %) 1.0e-3) runs)))))

;; termination helpers

(t/deftest termination-helpers
  (t/testing "range of the population"
    (let [settings {:lo (double-array [0.0 0.0]) :hi (double-array [10.0 100.0])}
          pop (fn [& pts] (mapv (fn [pt] {:x (double-array pt) :value 0.0}) pts))
          range-of (fn [& pts] (#'sut/population-range (apply pop pts) settings))]
      (t/is (= 0.0 (range-of [1 1])))
      (t/is (= 0.0 (range-of [1 1] [1 1] [1 1])))
      (t/is (m/delta-eq 0.1 (range-of [1 1] [2 1])))
      (t/is (m/delta-eq 0.5 (range-of [1 1] [2 1] [1 51])) "the largest of the dimensions")
      (t/is (m/delta-eq 1.0 (range-of [0 0] [10 100])))
      (t/is (= 0.0 (#'sut/population-range (pop [0.0 Double/MIN_VALUE]) {:lo (double-array [0.0 0.0]) :hi (double-array [Double/MIN_VALUE 1.0])}))
            "a tiny range of the bounds does not overflow")))
  (t/testing "improvement over the last loops"
    (let [conv? (fn [bests loops tol] (#'sut/improvement-converged? bests loops tol))]
      (t/is (false? (conv? [1.0] 1 0.1)) "no loop yet")
      (t/is (false? (conv? [1.0 1.0] 2 0.1)) "fewer loops than needed")
      (t/is (true? (conv? [1.0 1.0] 1 0.1)))
      (t/is (true? (conv? [1.0 1.0 1.0] 2 0.0)))
      (t/is (false? (conv? [2.0 1.0] 1 0.1)))
      (t/is (true? (conv? [1.1 1.0] 1 0.11)) "a difference within the relative tolerance")
      (t/is (false? (conv? [1.1 1.0] 1 0.09)))
      (t/is (true? (conv? [1.0e-11 0.0] 1 1.0e-5)) "near zero the allowed difference is the square of the tolerance")
      (t/is (false? (conv? [1.0e-9 0.0] 1 1.0e-5)))
      (t/is (false? (conv? [1.0 0.0] 1 1.0e-5)))
      (t/is (true? (conv? [0.0 0.0] 1 0.0)))
      (t/is (false? (conv? [##Inf ##Inf] 1 0.1)) "non-finite values never converge")
      (t/is (false? (conv? [1.0 ##Inf] 1 0.1)))
      (t/is (false? (conv? [##-Inf ##-Inf] 1 0.1)))
      (t/is (false? (conv? [##NaN ##NaN] 1 0.1)))
      (t/is (true? (conv? [5.0 3.0 1.0 1.0] 1 0.0)) "only the last loops count")
      (t/is (false? (conv? [5.0 3.0 1.0 1.0] 2 0.0)))))
  (t/testing "stalled"
    (t/is (true? (#'sut/stalled? [2.0 2.0])))
    (t/is (true? (#'sut/stalled? [1.0 2.0])))
    (t/is (false? (#'sut/stalled? [2.0 1.0])))
    (t/is (true? (#'sut/stalled? [##Inf ##Inf])))))

;; recovery of a collapsed population

(def ^:private sphere3 (fn [x] (reduce + (map #(* % %) x))))
(def ^:private unit-cube-3 {:lo (double-array [0.0 0.0 0.0]) :hi (double-array [1.0 1.0 1.0])})

(t/deftest lost-directions
  (let [unit-norm (fn [d] (Math/sqrt (reduce + (map #(* % %) d))))
        lost (fn [pop settings] (#'sut/lost-directions pop settings))]
    (t/testing "points on a line: two directions are lost, both orthogonal to the line"
      (let [[pop _] (points-of sphere3 (for [t (range 0.0 1.0 0.05)] [t t t]))
            dirs (lost pop unit-cube-3)
            along (v/normalize (v/vec3 1 1 1))]
        (t/is (= 2 (count dirs)))
        (t/is (every? #(m/delta-eq 1.0 (unit-norm %) 1.0e-9) dirs))
        (t/is (every? #(m/delta-eq 0.0 (reduce + (map * % (vec along))) 1.0e-9) dirs))))
    (t/testing "points on a plane: one direction"
      (let [[pop _] (points-of sphere3 (for [s (range 0.0 1.0 0.2) t (range 0.0 1.0 0.2)] [s t (+ s t)]))]
        (t/is (= 1 (count (lost pop unit-cube-3))))))
    (t/testing "random points of the cube: none"
      (let [rng (r/rng :jdk 1)
            [pop _] (points-of sphere3 (repeatedly 60 #(vec (repeatedly 3 (fn [] (r/drandom rng))))))]
        (t/is (empty? (lost pop unit-cube-3)))))
    (t/testing "identical points have no variance: nothing is lost"
      (let [[pop _] (points-of sphere3 (repeat 10 [0.5 0.5 0.5]))]
        (t/is (empty? (lost pop unit-cube-3)))))
    (t/testing "five and more dimensions (the eigenvectors of the decomposition for 5 and more dimensions)"
      (let [n 6
            [pop _] (points-of (fn [x] (reduce + x)) (for [t (range 0.0 1.0 0.05)] (vec (repeat n t))))]
        (t/is (= (dec n) (count (lost pop {:lo (double-array (repeat n 0.0)) :hi (double-array (repeat n 1.0))}))))))
    (t/testing "the scale of the bounds does not matter"
      (let [[pop _] (points-of sphere3 (for [t (range 0.0 100.0 5.0)] [t (* 2.0 t) 3.0]))]
        (t/is (= 2 (count (lost pop {:lo (double-array [0.0 0.0 0.0]) :hi (double-array [100.0 200.0 3.0])}))))))
    (t/testing "directions are double arrays"
      (let [[pop _] (points-of sphere3 (for [t (range 0.0 1.0 0.05)] [t t t]))]
        (t/is (every? #(instance? (Class/forName "[D") %) (lost pop unit-cube-3)))))))

(defn- lost-directions-reference
  "The previous implementation of the lost directions: the covariance matrix from `stats/covariance-matrix` of boxed columns."
  [population {:keys [lo hi]}]
  (let [extent (v/sub hi lo)
        rows (map #(mapv m// (v/sub (:x %) lo) extent) population)
        eigen (mat/eigen-decomposition (mat/rows->mat (stats/covariance-matrix (apply mapv vector rows)))
                                       {:eigenvectors-scaling :raw})
        values (vec (mat/decomposition-component eigen :real-eigenvalues))
        largest (stats/maximum values)]
    (if (pos? largest)
      (into []
            (comp (keep-indexed (fn [i direction] (when (< (values i) (* 1.0e-3 largest)) direction)))
                  (map v/vec->array))
            (mat/decomposition-component eigen :eigenvectors))
      [])))

(defn- projector
  "The matrix of the orthogonal projection on the span of the directions; it does not depend on the basis chosen in the span."
  [directions n]
  (for [i (range n)]
    (for [j (range n)]
      (reduce + (map #(* (nth % i) (nth % j)) directions)))))

(defn- max-abs-difference [a b]
  (apply max 0.0 (map (fn [row-a row-b] (apply max 0.0 (map #(Math/abs (double (- %1 %2))) row-a row-b))) a b)))

(t/deftest lost-directions-equivalence
  ;; a multiple eigenvalue has no unique eigenvectors, so the lost subspaces are compared (projectors), and the directions
  ;; themselves when only one is lost
  (let [rng (r/rng :jdk 2)
        random-cube (fn [n size] (for [_ (range size)] (vec (repeatedly n #(r/drandom rng)))))
        populations {"random cube 3D" [3 (random-cube 3 60)]
                     "random cube 5D" [5 (random-cube 5 55)]
                     "random cube 20D" [20 (random-cube 20 205)]
                     "line 3D" [3 (for [t (range 0.0 1.0 0.05)] [t t t])]
                     "line 6D" [6 (for [t (range 0.0 1.0 0.05)] (vec (repeat 6 t)))]
                     "plane 3D" [3 (for [s (range 0.0 1.0 0.2) t (range 0.0 1.0 0.2)] [s t (+ s t)])]
                     "plane in 20D" [20 (for [_ (range 60)] (let [s (r/drandom rng) t (r/drandom rng)] (vec (take 20 (cycle [s t (* 0.5 (+ s t))])))))]
                     "identical points" [3 (repeat 10 [0.5 0.5 0.5])]
                     "two points" [3 [[0.1 0.1 0.1] [0.9 0.5 0.2]]]}]
    (doseq [[label [n coords]] populations
            :let [[pop _] (points-of (fn [_] 0.0) coords)
                  bounds {:lo (double-array (repeat n 0.0)) :hi (double-array (repeat n 1.0))}
                  old (lost-directions-reference pop bounds)
                  new (#'sut/lost-directions pop bounds)]]
      (t/is (= (count old) (count new)) label)
      (t/is (> 1.0e-9 (max-abs-difference (projector (map vec old) n) (projector (map vec new) n))) (str label ": the same lost subspace"))
      (when (= 1 (count old))
        (let [d-old (vec (first old)) d-new (vec (first new))]
          (t/is (or (v/delta-eq d-old d-new 1.0e-9) (v/delta-eq d-old (v/mult d-new -1.0) 1.0e-9)) (str label ": the same direction up to the sign")))))
    (t/testing "bounds of different scales"
      (let [[pop _] (points-of (fn [_] 0.0) (for [t (range 0.0 100.0 5.0)] [t (* 2.0 t) 3.0]))
            bounds {:lo (double-array [0.0 0.0 0.0]) :hi (double-array [100.0 200.0 3.0])}]
        (t/is (= (count (lost-directions-reference pop bounds)) (count (#'sut/lost-directions pop bounds))))))
    (t/testing "the covariance matrix equals stats/covariance-matrix"
      (let [coords (random-cube 4 40)
            [pop _] (points-of (fn [_] 0.0) coords)
            lo (double-array (repeat 4 0.0))
            hi (double-array (repeat 4 1.0))]
        (t/is (> 1.0e-12 (max-abs-difference (map vec (#'sut/unit-cube-covariance pop lo hi))
                                              (stats/covariance-matrix (apply mapv vector (map :x pop))))))))))

(t/deftest recover-dimensions
  (let [line (for [t (range 0.0 1.0 0.05)] [t t t])
        [pop _] (points-of sphere3 line)
        recover (fn [pop seed]
                  (let [ev (#'sut/evaluator sphere3 1.0 1000000)]
                    [(#'sut/recover-dimensions pop (assoc unit-cube-3 :evaluate ev) (r/rng :jdk seed)) ev]))
        [recovered ev] (recover pop 1)
        best-before (vec (:x (first pop)))]
    (t/is (vector? recovered))
    (t/is (= (count pop) (count recovered)))
    (t/is (<= (:value (first recovered)) (:value (first pop))) "the best value is not worse")
    (t/is (some #(= best-before (vec (:x %))) recovered) "the best point stays in the population")
    (t/is (apply <= (values-of recovered)) "sorted")
    (t/is (pos? (count (filter (fn [{[x y] :x}] (> (Math/abs (double (- x y))) 1.0e-9)) recovered))) "points left the line")
    (t/is (= 4 (ev)) "moved points are evaluated and counted: 2 lost directions times 2 points (10% of 20)")
    (t/is (every? (fn [{[x y z] :x}] (and (<= 0.0 x 1.0) (<= 0.0 y 1.0) (<= 0.0 z 1.0))) recovered) "in the bounds")
    (t/is (every? #(= (:value %) (sphere3 (:x %))) recovered) "values are values of the function")
    (t/testing "the moved points are different ones, the best is never moved"
      (doseq [seed (range 1 21)
              :let [[out _] (recover pop seed)
                    unchanged (count (filter (set pop) out))]]
        (t/is (<= 16 unchanged 18) (str "seed " seed))
        (t/is (= best-before (vec (:x (first out)))) (str "seed " seed))))
    (t/testing "a population of two and three points: only the second and third can move"
      (doseq [size [2 3]
              :let [[small _] (points-of sphere3 (for [t (take size [0.1 0.5 0.9])] [t t t]))
                    [out ev] (recover small 1)]]
        (t/is (= size (count out)) (str size))
        (t/is (some #(= (vec (:x (first small))) (vec (:x %))) out) (str size ": the best stays"))
        (t/is (= 2 (ev)) (str size ": two lost directions, one point each"))))
    (t/testing "nothing is lost: the population stays as it is"
      (let [rng (r/rng :jdk 1)
            [random-pop _] (points-of sphere3 (repeatedly 60 #(vec (repeatedly 3 (fn [] (r/drandom rng))))))
            [out ev] (recover random-pop 2)]
        (t/is (= random-pop out))
        (t/is (zero? (ev)))))))

(t/deftest pca-recovery-runs
  (let [f p/rosenbrock]
    (t/testing "the value of the result is the value of the function at its point"
      (doseq [seed (range 1 11)
              n [2 3 4 6]
              :let [[pt val] (run f seed {:bounds (p/rosenbrock-bounds n) :pca-recovery? true})]]
        (t/is (= val (f pt)) (str seed " " n))))
    (t/testing "stats of runs with recovery"
      (doseq [seed (range 1 4)
              :let [s (run f seed (merge {:bounds (p/rosenbrock-bounds 3) :stats? true :pca-recovery? true} (never-stop 6)))]]
        (t/is (= 6 (:iterations s)))
        (t/is (= :max-iterations (:status s)))
        (t/is (= (:value s) (f (:point s))))))
    (t/testing "one dimension and too small populations skip the recovery"
      (t/is (some? (run (fn [[x]] (* x x)) {:bounds [[-1 1]] :pca-recovery? true})))
      (t/is (some? (run p/rosenbrock {:bounds (p/rosenbrock-bounds 4) :pca-recovery? true :complexes 1 :complex-size 3 :subcomplex-size 2 :max-iters 5})))
      (t/is (some? (run p/rosenbrock {:bounds (p/rosenbrock-bounds 4) :pca-recovery? true :complexes 1 :complex-size 4 :max-iters 5}))))
    (t/testing "with parallel evolution and reduction"
      (let [[pt val] (run f {:bounds (p/rosenbrock-bounds 3) :pca-recovery? true :parallel? true :complexes 4 :min-complexes 2})]
        (t/is (= val (f pt)))))))

;; parallel mode

(t/deftest parallel-mode
  (doseq [complexes [1 2 5 9]
          :let [[pt val] (run p/himmelblau {:bounds hb :parallel? true :complexes complexes})]]
    (t/is (< val 1.0e-3) (str complexes))
    (t/is (= val (p/himmelblau pt)) (str complexes)))
  (t/testing "the objective is called from many threads"
    (let [threads (atom #{})
          f (fn [x] (swap! threads conj (Thread/currentThread)) (p/himmelblau x))]
      (run f {:bounds hb :parallel? true :complexes 8 :max-iters 3})
      (t/is (> (count @threads) 1))))
  (t/testing "sequential mode uses the calling thread only"
    (let [threads (atom #{})
          f (fn [x] (swap! threads conj (Thread/currentThread)) (p/himmelblau x))]
      (run f {:bounds hb :complexes 8 :max-iters 3})
      (t/is (= #{(Thread/currentThread)} @threads)))))
