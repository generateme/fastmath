(ns fastmath.optimization.sceua-test
  (:require [fastmath.optimization.sceua :as sut]
            [fastmath.optimization.problems :as p]
            [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.test :as t])
  (:import [java.util.concurrent.atomic AtomicLong]))

;; References: Himmelblau's function has four minima with value 0: (3, 2), (-2.805118, 3.131313),
;; (-3.779310, -3.283186), (3.584428, -1.848126) (Wikipedia). Sphere, Rosenbrock (minimum 0 at (1, ..., 1)) and Ackley
;; (minimum 0 at the origin) are described at https://www.sfu.ca/~ssurjano/optimization.html. problem02 (sin x + sin(10x/3)
;; on [2.7, 7.5]) has its global minimum -1.899599 at 5.145735 (https://infinity77.net/global_optimization/test_functions_1d.html).
;;
;; The optimizer is stochastic, so all runs use a seeded generator. Tolerances of the values below are the worst
;; absolute value over 30 seeds with :stop-loops 25 times at least 10, except Ackley (2D) where the worst of 30 seeds,
;; 1.5e-3, is limited by :stop-range (1e-3 of the range 65.5 is 0.065) and the tolerance is 1e-2. With the default
;; :stop-loops (7) about 10% of runs on Himmelblau stop before the minimum, so the convergence tests set :stop-loops 25.

(def ^:private himmelblau-minima [[3.0 2.0] [-2.805118 3.131313] [-3.779310 -3.283186] [3.584428 -1.848126]])
(def ^:private hb (p/himmelblau-bounds))
(def ^:private ackley (p/->ackley))

(defn- ex-data-of [thunk]
  (try (thunk) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))

(defn- near-minimum? [pt] (boolean (some #(v/delta-eq (vec pt) % 1.0e-3) himmelblau-minima)))

(defn- run
  "Runs the optimizer with a seeded generator and with 25 loops of patience."
  ([f opts] (run f 1 opts))
  ([f seed opts] (sut/sceua f (merge {:rng (r/rng :jdk seed) :stop-loops 25} opts))))

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
      (t/is (m/delta-eq 5.145735 (first pt) 1.0e-2) (str "seed " seed))
      (t/is (m/delta-eq -1.899599 val 1.0e-4) (str "seed " seed)))))

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
  (let [ev (#'sut/evaluator p/sphere 1.0 (AtomicLong. 0) 100000)
        lo (double-array [-1.0 -2.0])
        hi (double-array [1.0 2.0])
        build (fn [initial size jitter] (#'sut/initial-population ev lo hi initial size jitter (r/rng :jdk 1)))]
    (doseq [size [1 2 25 100]
            initial [nil (double-array [0.5 0.5])]
            jitter [0.0 0.25 1.0]
            :let [pop (build initial size jitter)
                  label (str size " " (some? initial) " " jitter)]]
      (t/is (= size (alength pop)) label)
      (t/is (every? #(= 3 (alength ^doubles %)) pop) label)
      (t/is (apply <= (map #(aget ^doubles % 2) pop)) (str label ": sorted by value"))
      (t/is (every? (fn [^doubles pt] (and (<= -1.0 (aget pt 0) 1.0) (<= -2.0 (aget pt 1) 2.0))) pop) (str label ": in the bounds"))
      (t/is (every? (fn [^doubles pt] (= (aget pt 2) (p/sphere [(aget pt 0) (aget pt 1)]))) pop) (str label ": values"))
      (when initial
        (t/is (some (fn [^doubles pt] (and (= 0.5 (aget pt 0)) (= 0.5 (aget pt 1)))) pop) (str label ": the initial point is there"))))
    (t/testing "the generator decides the points"
      (t/is (= (map vec (build nil 10 0.25)) (map vec (build nil 10 0.25))))
      (t/is (not= (map vec (build nil 10 0.25)) (map vec (#'sut/initial-population ev lo hi nil 10 0.25 (r/rng :jdk 2))))))
    (t/testing "more than 14 dimensions use a Sobol sequence"
      (doseq [n [14 15 16]
              :let [lo (double-array (repeat n 0.0))
                    hi (double-array (repeat n 1.0))
                    pop (#'sut/initial-population (#'sut/evaluator (fn [_] 0.0) 1.0 (AtomicLong. 0) 1000) lo hi nil 10 0.25 (r/rng :jdk 1))]]
        (t/is (= 10 (alength pop)) (str n))
        (t/is (every? #(= (inc n) (alength ^doubles %)) pop) (str n))))))

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
    (let [pop (object-array (range 15))
          reduced (fn [ngs minimum] (#'sut/reduce-complexes pop ngs minimum 5))]
      (t/is (= [10 2] [(alength ^objects (first (reduced 3 2))) (second (reduced 3 2))]))
      (t/is (= (seq (range 10)) (seq (first (reduced 3 2)))) "the worst points are the last ones")
      (t/is (= [15 3] [(alength ^objects (first (reduced 3 3))) (second (reduced 3 3))]))
      (t/is (= [10 1] [(alength ^objects (first (reduced 2 1))) (second (reduced 2 1))]) "one complex of the size is removed")
      (t/is (= [15 1] [(alength ^objects (first (#'sut/reduce-complexes pop 1 1 15))) (second (#'sut/reduce-complexes pop 1 1 15))]))))
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
        (let [res (fn [] (sut/sceua f (assoc opts :rng (make) :parallel? parallel? :stop-loops 25)))]
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
  (t/testing "the function gets a new vector for every call"
    (let [seen (atom [])]
      (run (fn [x] (swap! seen conj x) (p/himmelblau x)) {:bounds hb :max-iters 2})
      (t/is (every? vector? @seen))
      (t/is (every? #(= 2 (count %)) @seen))
      (t/is (= (count @seen) (count (set (map #(System/identityHashCode %) @seen)))))))
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

(defn- make-complex
  "A sorted complex of points of the 2D sphere function."
  [coords]
  (let [ev (#'sut/evaluator p/sphere 1.0 (AtomicLong. 0) 1000000)
        pop (object-array (map #(#'sut/make-point ev (double-array %)) coords))]
    (java.util.Arrays/sort pop @#'sut/by-value)
    [pop ev]))

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
          counts (frequencies (repeatedly n #(aget ^longs (#'sut/select-parents m 1 rng) 0)))]
      (doseq [i (range m)]
        (t/is (m/delta-eq (/ (* 2.0 (- m i)) (* m (inc m))) (/ (get counts i 0) (double n)) 0.02) (str i)))))
  (t/testing "different indices in the ascending order, for any size of the sub-complex"
    (let [rng (r/rng :jdk 5)]
      (doseq [m [2 3 7 20]
              q (range 1 (inc m))
              _ (range 20)
              :let [idx (vec (#'sut/select-parents m q rng))]]
        (t/is (= q (count idx)) (str m " " q))
        (t/is (apply < idx) (str m " " q " " idx))
        (t/is (every? #(<= 0 % (dec m)) idx) (str m " " q " " idx)))))
  (t/testing "the whole complex is selected when the sub-complex is the complex"
    (t/is (= [0 1 2 3 4] (vec (#'sut/select-parents 5 5 (r/rng :jdk 1)))))))

(t/deftest evolve-complex
  (let [[complex ev] (make-complex [[1.0 1.0] [2.0 0.0] [-1.0 2.0] [0.0 -3.0] [3.0 3.0]])
        before (mapv vec complex)
        settings {:subcomplex-size 3 :evolution-steps 20 :evaluate ev
                  :lo (double-array [-4.0 -4.0]) :hi (double-array [4.0 4.0])}
        evolved (#'sut/evolve-complex complex settings (r/rng :jdk 1))]
    (t/is (= before (mapv vec complex)) "the input is not changed")
    (t/is (not (identical? complex evolved)))
    (t/is (= 5 (alength evolved)))
    (t/is (apply <= (map #(aget ^doubles % 2) evolved)) "sorted by value")
    (t/is (<= (aget ^doubles (aget evolved 0) 2) (aget ^doubles (aget complex 0) 2)) "the best point is not worse")
    (t/is (every? (fn [^doubles pt] (<= -4.0 (aget pt 0) 4.0)) evolved) "in the bounds")
    (t/is (every? (fn [^doubles pt] (= (aget pt 2) (p/sphere [(aget pt 0) (aget pt 1)]))) evolved) "values are values of the function")
    (t/testing "a sub-complex of the size of the complex and one step"
      (let [one-step (#'sut/evolve-complex complex (assoc settings :subcomplex-size 5 :evolution-steps 1) (r/rng :jdk 1))]
        (t/is (= 5 (alength one-step)))
        (t/is (apply <= (map #(aget ^doubles % 2) one-step)))))
    (t/testing "steps use one to three evaluations"
      (let [counter (AtomicLong. 0)
            ev (#'sut/evaluator p/sphere 1.0 counter 1000000)
            steps 30
            _ (#'sut/evolve-complex complex (assoc settings :evaluate ev :evolution-steps steps) (r/rng :jdk 2))]
        (t/is (<= steps (.get counter) (* 3 steps)))))
    (t/testing "the reflection is clipped to the bounds"
      (let [[cx cev] (make-complex [[0.9 0.9] [1.0 1.0] [-1.0 -1.0] [0.95 -0.9] [0.0 0.0]])
            tight {:subcomplex-size 2 :evolution-steps 50 :evaluate cev :lo (double-array [-1.0 -1.0]) :hi (double-array [1.0 1.0])}
            out (#'sut/evolve-complex cx tight (r/rng :jdk 3))]
        (t/is (every? (fn [^doubles pt] (and (<= -1.0 (aget pt 0) 1.0) (<= -1.0 (aget pt 1) 1.0))) out))))))

;; termination helpers

(t/deftest termination-helpers
  (t/testing "range of the population"
    (let [lo (double-array [0.0 0.0])
          hi (double-array [10.0 100.0])
          pop (fn [& pts] (object-array (map #(double-array (conj (vec %) 0.0)) pts)))]
      (t/is (= 0.0 (#'sut/population-range (pop [1 1]) lo hi)))
      (t/is (= 0.0 (#'sut/population-range (pop [1 1] [1 1] [1 1]) lo hi)))
      (t/is (m/delta-eq 0.1 (#'sut/population-range (pop [1 1] [2 1]) lo hi)))
      (t/is (m/delta-eq 0.5 (#'sut/population-range (pop [1 1] [2 1] [1 51]) lo hi)) "the largest of the dimensions")
      (t/is (m/delta-eq 1.0 (#'sut/population-range (pop [0 0] [10 100]) lo hi)))))
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

(defn- population-of [coords f]
  (let [ev (#'sut/evaluator f 1.0 (AtomicLong. 0) 1000000)
        pop (object-array (map #(#'sut/make-point ev (double-array %)) coords))]
    (java.util.Arrays/sort pop @#'sut/by-value)
    [pop ev]))

(t/deftest lost-directions
  (let [lo (double-array [0.0 0.0 0.0])
        hi (double-array [1.0 1.0 1.0])
        sphere3 (fn [x] (reduce + (map #(* % %) x)))
        unit-norm (fn [d] (Math/sqrt (reduce + (map #(* % %) d))))]
    (t/testing "points on a line: two directions are lost, both orthogonal to the line"
      (let [[pop _] (population-of (for [t (range 0.0 1.0 0.05)] [t t t]) sphere3)
            dirs (#'sut/lost-directions pop lo hi)
            along (v/normalize (v/vec3 1 1 1))]
        (t/is (= 2 (count dirs)))
        (t/is (every? #(m/delta-eq 1.0 (unit-norm %) 1.0e-9) dirs))
        (t/is (every? #(m/delta-eq 0.0 (reduce + (map * % (vec along))) 1.0e-9) dirs))))
    (t/testing "points on a plane: one direction"
      (let [[pop _] (population-of (for [s (range 0.0 1.0 0.2) t (range 0.0 1.0 0.2)] [s t (+ s t)]) sphere3)]
        (t/is (= 1 (count (#'sut/lost-directions pop lo hi))))))
    (t/testing "random points of the cube: none"
      (let [rng (r/rng :jdk 1)
            [pop _] (population-of (repeatedly 60 #(vec (repeatedly 3 (fn [] (r/drandom rng))))) sphere3)]
        (t/is (empty? (#'sut/lost-directions pop lo hi)))))
    (t/testing "identical points have no variance: nothing is lost"
      (let [[pop _] (population-of (repeat 10 [0.5 0.5 0.5]) sphere3)]
        (t/is (empty? (#'sut/lost-directions pop lo hi)))))
    (t/testing "five and more dimensions (the eigenvectors of the decomposition for 5 and more dimensions)"
      (let [n 6
            [pop _] (population-of (for [t (range 0.0 1.0 0.05)] (vec (repeat n t))) (fn [x] (reduce + x)))]
        (t/is (= (dec n) (count (#'sut/lost-directions pop (double-array (repeat n 0.0)) (double-array (repeat n 1.0))))))))
    (t/testing "the scale of the bounds does not matter"
      (let [[pop _] (population-of (for [t (range 0.0 100.0 5.0)] [t (* 2.0 t) 3.0]) sphere3)]
        (t/is (= 2 (count (#'sut/lost-directions pop (double-array [0.0 0.0 0.0]) (double-array [100.0 200.0 3.0])))))))))

(t/deftest recover-dimensions
  (let [lo (double-array [0.0 0.0 0.0])
        hi (double-array [1.0 1.0 1.0])
        sphere3 (fn [x] (reduce + (map #(* % %) x)))
        [pop _] (population-of (for [t (range 0.0 1.0 0.05)] [t t t]) sphere3)
        before (mapv vec pop)
        best-before (vec (aget pop 0))
        counter (AtomicLong. 0)
        ev (#'sut/evaluator sphere3 1.0 counter 1000000)
        recovered (#'sut/recover-dimensions! pop (r/rng :jdk 1) ev lo hi)]
    (t/is (identical? pop recovered) "the population is changed in place")
    (t/is (= (count before) (alength recovered)))
    (t/is (<= (aget ^doubles (aget recovered 0) 3) (get best-before 3)) "the best value is not worse")
    (t/is (some #(= best-before (vec %)) recovered) "the best point stays in the population")
    (t/is (apply <= (map #(aget ^doubles % 3) recovered)) "sorted")
    (t/is (pos? (count (filter (fn [^doubles pt] (> (Math/abs (- (aget pt 0) (aget pt 1))) 1.0e-9)) recovered))) "points left the line")
    (t/is (= 4 (.get counter)) "moved points are evaluated and counted: 2 lost directions times 2 points (10% of 20)")
    (t/is (every? (fn [^doubles pt] (and (<= 0.0 (aget pt 0) 1.0) (<= 0.0 (aget pt 1) 1.0) (<= 0.0 (aget pt 2) 1.0))) recovered) "in the bounds")
    (t/is (every? (fn [^doubles pt] (= (aget pt 3) (sphere3 [(aget pt 0) (aget pt 1) (aget pt 2)]))) recovered) "values are values of the function"))
  (t/testing "nothing is lost: the population stays as it is"
    (let [rng (r/rng :jdk 1)
          sphere3 (fn [x] (reduce + (map #(* % %) x)))
          [pop ev] (population-of (repeatedly 60 #(vec (repeatedly 3 (fn [] (r/drandom rng))))) sphere3)
          before (mapv vec pop)]
      (#'sut/recover-dimensions! pop (r/rng :jdk 2) ev (double-array [0.0 0.0 0.0]) (double-array [1.0 1.0 1.0]))
      (t/is (= before (mapv vec pop)))))
  (t/testing "sampling of indices"
    (doseq [[from to k] [[1 2 1] [1 10 9] [1 10 1] [0 5 5] [3 4 1]]
            :let [idx (vec (#'sut/sample-indices (r/rng :jdk 1) from to k))]]
      (t/is (= k (count idx)) (str from to k))
      (t/is (= k (count (set idx))) (str from to k))
      (t/is (every? #(and (<= from %) (< % to)) idx) (str from to k)))))

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
