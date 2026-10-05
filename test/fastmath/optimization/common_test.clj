(ns fastmath.optimization.common-test
  (:require [fastmath.optimization.common :as sut]
            [clojure.test :as t]
            [clojure.string :as str])
  (:import [org.apache.commons.math3.optim.nonlinear.scalar GoalType]))

(defn- ex-data-of
  "Returns ex-data of the ExceptionInfo thrown by `thunk` or nil when nothing (or something else) is thrown."
  [thunk]
  (try (thunk) nil
       (catch clojure.lang.ExceptionInfo e (ex-data e))))

(t/deftest throw-unknown
  (let [d (ex-data-of #(sut/throw-unknown :goal :foo #{:a :b}))]
    (t/is (= {:goal :foo :allowed #{:a :b}} d)))
  (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #":goal.*:foo.*allowed"
                          (sut/throw-unknown :goal :foo [:a])))
  ;; nil value and nil allowed are reported, not swallowed
  (t/is (= {:x nil :allowed nil} (ex-data-of #(sut/throw-unknown :x nil nil)))))

(t/deftest throw-invalid-option
  (t/is (= {:option :jitter :value 2 :reason "too big"} (ex-data-of #(sut/throw-invalid-option :jitter 2 "too big"))))
  (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #":jitter.*too big" (sut/throw-invalid-option :jitter 2 "too big")))
  (t/is (= {:option :x :value nil :reason nil} (ex-data-of #(sut/throw-invalid-option :x nil nil)))))

(t/deftest positive-integer-option
  (t/testing "missing and nil give the default"
    (t/is (= 5 (sut/positive-integer-option {} :n 5)))
    (t/is (= 5 (sut/positive-integer-option {:n nil} :n 5)))
    (t/is (= 5 (sut/positive-integer-option nil :n 5))))
  (t/testing "integers, also the limit values"
    (t/are [v] (= (long v) (sut/positive-integer-option {:n v} :n 5))
      1 2 (int 3) (short 4) (byte 5) Long/MAX_VALUE))
  (t/testing "everything else throws"
    (doseq [bad [0 -1 Long/MIN_VALUE 1.5 2.0 1/2 3N "3" ##NaN true :a [1] {} 0.0]
            :let [d (ex-data-of #(sut/positive-integer-option {:n bad} :n 5))]]
      (t/is (= :n (:option d)) (pr-str bad))
      (t/is (= bad (:value d)) (pr-str bad))
      (t/is (string? (:reason d)) (pr-str bad)))))

(t/deftest nonnegative-number-option
  (t/testing "missing and nil give the default as a double"
    (t/is (= 0.5 (sut/nonnegative-number-option {} :x 0.5)))
    (t/is (= 1.0 (sut/nonnegative-number-option {:x nil} :x 1)))
    (t/is (= 0.0 (sut/nonnegative-number-option nil :x 0))))
  (t/testing "finite numbers not less than zero, also the limit values"
    (t/are [v] (= (double v) (sut/nonnegative-number-option {:x v} :x 5))
      0 0.0 1 1/2 1.0e300 Double/MAX_VALUE Double/MIN_VALUE 3N (float 0.5)))
  (t/testing "negative, non-finite and not numbers throw"
    (doseq [bad [-1 -1.0e-300 -0.5 ##NaN ##Inf ##-Inf "a" :a [1] true {}]
            :let [d (ex-data-of #(sut/nonnegative-number-option {:x bad} :x 5))]]
      (t/is (= :x (:option d)) (pr-str bad))
      (t/is (= (pr-str bad) (pr-str (:value d))) (pr-str bad))
      (t/is (string? (:reason d)) (pr-str bad)))))

(t/deftest parse-goal
  (t/are [in out] (= out (sut/parse-goal in))
    nil :minimize
    :minimize :minimize
    :maximize :maximize)
  (doseq [bad [:min :max "minimize" 0 false :MINIMIZE []]]
    (t/is (= {:goal bad :allowed #{:minimize :maximize}} (ex-data-of #(sut/parse-goal bad)))
          (str "should reject " (pr-str bad)))))

(t/deftest goal-type
  (t/is (= GoalType/MINIMIZE (sut/goal-type nil)))
  (t/is (= GoalType/MINIMIZE (sut/goal-type :minimize)))
  (t/is (= GoalType/MAXIMIZE (sut/goal-type :maximize)))
  (t/is (= {:goal :foo :allowed #{:minimize :maximize}} (ex-data-of #(sut/goal-type :foo)))))

(t/deftest resolve-vector-arg
  ;; nil -> default per method, only brent differs
  (doseq [m [:lbfgsb :bobyqa :cmaes :sceua :nelder-mead :multidirectional-simplex :powell :gradient :non-linear-gradient]]
    (t/is (true? (sut/resolve-vector-arg? m nil)) (str m)))
  (t/is (false? (sut/resolve-vector-arg? :brent nil)))
  ;; explicit values win, for every method including brent
  (doseq [m [:brent :lbfgsb :powell]]
    (t/is (true? (sut/resolve-vector-arg? m true)))
    (t/is (false? (sut/resolve-vector-arg? m false))))
  ;; always a real boolean
  (t/is (true? (sut/resolve-vector-arg? :brent :yes)))
  (t/is (false? (sut/resolve-vector-arg? :lbfgsb false))))

(t/deftest ->vector-fn
  (let [f (fn [x y] (+ x y))
        g (fn [v] (reduce + v))]
    (t/is (identical? g (sut/->vector-fn g true)))
    (t/is (= 3 ((sut/->vector-fn f false) [1 2])))
    ;; any sequential or array-like input, also lazy seq
    (t/is (= 3 ((sut/->vector-fn f false) (map inc [0 1]))))
    (t/is (== 3 ((sut/->vector-fn f false) (double-array [1 2]))))
    ;; zero and one argument boundaries
    (t/is (= :none ((sut/->vector-fn (fn [] :none) false) [])))
    (t/is (= 5 ((sut/->vector-fn (fn [x] x) false) [5])))
    ;; arity mismatch is a caller error and surfaces
    (t/is (thrown? clojure.lang.ArityException ((sut/->vector-fn f false) [1])))))

(t/deftest limits
  (t/is (= 10000 (.getMaxEval (sut/max-eval nil))))
  (t/is (= 10000 (.getMaxIter (sut/max-iter nil))))
  (t/is (= 5 (.getMaxEval (sut/max-eval 5))))
  (t/is (= 1 (.getMaxIter (sut/max-iter 1))))
  (t/is (= Integer/MAX_VALUE (.getMaxEval (sut/max-eval Integer/MAX_VALUE)))))

(defn- reason-of
  "Reason string of the bounds exception thrown for given arguments, nil when valid. Also checks the ex-data shape."
  [method bounds initial]
  (let [d (ex-data-of #(sut/normalize-bounds method bounds initial))]
    (when d
      (t/is (= #{:method :bounds :initial :reason} (set (keys d))) "ex-data keys")
      (t/is (= method (:method d)))
      (t/is (= bounds (:bounds d)))
      (t/is (= initial (:initial d)))
      (:reason d))))

(defn- rejected? [method bounds initial & [fragment]]
  (let [r (reason-of method bounds initial)]
    (and (string? r) (or (nil? fragment) (str/includes? r fragment)))))

(def ^:private all-methods [:brent :bobyqa :cmaes :sceua :nelder-mead :multidirectional-simplex :lbfgsb :powell :gradient :non-linear-gradient])

(t/deftest normalize-bounds-canonical-form
  ;; integers and ratios become doubles; lazy seqs and vectors are accepted
  (t/is (= [[0.0 1.0] [-2.0 0.5]] (sut/normalize-bounds :bobyqa [[0 1] [-2 1/2]] nil)))
  (t/is (= [[-5.0 10.0] [-5.0 10.0] [-5.0 10.0]] (sut/normalize-bounds :lbfgsb (repeat 3 [-5 10]) nil)))
  (t/is (= [[-5.0 10.0]] (sut/normalize-bounds :lbfgsb '((-5 10)) nil)))
  ;; result is a vector of vectors of doubles
  (let [r (sut/normalize-bounds :cmaes [[0 1] [2 3]] nil)]
    (t/is (vector? r))
    (t/is (every? vector? r))
    (t/is (every? double? (flatten r))))
  ;; flat 1d form
  (doseq [m [:brent :lbfgsb :bobyqa :cmaes :powell]]
    (t/is (= [[0.0 10.0]] (sut/normalize-bounds m [0 10] nil)) (str m)))
  ;; a flat pair is only 1d: a 2-element sequence of pairs is two dimensions
  (t/is (= [[0.0 1.0] [2.0 3.0]] (sut/normalize-bounds :lbfgsb [[0 1] [2 3]] nil))))

(t/deftest normalize-bounds-structure
  (doseq [m all-methods
          bad [[] '() [[]] [[1]] [[1 2 3]] [[nil 1]] [[1 nil]] [[nil nil]] [["a" "b"]] [[0 1] nil] [[0 1] [2]]
               5 "ab" :a {:lo 0 :hi 1} [1] [1 2 3]]]
    (t/is (rejected? m bad [0 0]) (str m " should reject " (pr-str bad) " but gave " (pr-str (try (sut/normalize-bounds m bad [0 0]) (catch Exception e (ex-message e))))))
    (t/is (rejected? m bad nil) (str m " (no initial) should reject " (pr-str bad)))))

(t/deftest normalize-bounds-nan
  (doseq [m all-methods
          bad [[[##NaN 1]] [[0 ##NaN]] [[##NaN ##NaN]] [[0 1] [##NaN 1]]]
          :when (not (and (= :brent m) (= 2 (count bad))))] ;; brent rejects two pairs earlier
    (t/is (rejected? m bad nil "NaN") (str m " " (pr-str bad)))))

(t/deftest normalize-bounds-order
  (doseq [m all-methods]
    (t/is (rejected? m [[5 -5]] [0] "greater") (str m " lo>hi"))
    (t/is (rejected? m [[##Inf ##Inf]] [0] "+Inf") (str m " lo=+Inf"))
    (t/is (rejected? m [[##-Inf ##-Inf]] [0] "-Inf") (str m " hi=-Inf")))
  ;; every dimension is checked, not only the first
  (doseq [m (remove #{:brent} all-methods)]
    (t/is (rejected? m [[0 1] [3 2]] [0 0] "greater") (str m " lo>hi in 2nd dimension")))
  ;; lo = hi is a degenerate but legal range where the method accepts it
  (doseq [m [:bobyqa :cmaes :lbfgsb :powell :gradient :non-linear-gradient]]
    (t/is (= [[1.0 1.0] [0.0 2.0]] (sut/normalize-bounds m [[1 1] [0 2]] [1 0])) (str m)))
  (doseq [m [:brent :sceua :nelder-mead :multidirectional-simplex]]
    (t/is (rejected? m [[1 1]] [1] "less than") (str m " lo=hi"))))

(t/deftest normalize-bounds-brent
  (t/is (= [[2.7 7.5]] (sut/normalize-bounds :brent [[2.7 7.5]] nil)))
  (t/is (= [[2.7 7.5]] (sut/normalize-bounds :brent [2.7 7.5] 5.0)))
  (t/is (= [[2.7 7.5]] (sut/normalize-bounds :brent [2.7 7.5] [5.0])))
  (t/is (rejected? :brent nil nil "required"))
  ;; extra pairs were silently ignored before
  (t/is (rejected? :brent [[2.7 7.5] [1 2]] nil "exactly one"))
  (t/is (rejected? :brent [[2.7 7.5] [1 2]] [5.0 1.5] "exactly one"))
  (t/is (rejected? :brent [[0 ##Inf]] 1 "finite"))
  (t/is (rejected? :brent [[##-Inf 0]] 1 "finite"))
  (t/is (rejected? :brent [[##-Inf ##Inf]] nil "finite"))
  (t/is (rejected? :brent [[0 1]] [0.1 0.2] "differs")))

(t/deftest normalize-bounds-required-and-finite
  (doseq [m [:bobyqa :cmaes :sceua :lbfgsb]]
    (t/is (rejected? m nil nil "required") (str m))
    (t/is (rejected? m nil [0 0] "required") (str m)))
  (doseq [m [:bobyqa :cmaes :sceua]]
    (t/is (rejected? m [[0 ##Inf]] [0] "finite") (str m))
    (t/is (rejected? m [[##-Inf 0]] [0] "finite") (str m))
    (t/is (rejected? m [[##-Inf ##Inf]] nil "finite") (str m))))

(t/deftest normalize-bounds-sceua
  (t/is (= [[0.0 1.0] [-2.0 0.5]] (sut/normalize-bounds :sceua [[0 1] [-2 1/2]] nil)))
  (t/is (= [[0.0 1.0] [-2.0 0.5]] (sut/normalize-bounds :sceua [[0 1] [-2 1/2]] [0.5 0])))
  ;; one dimension: a flat pair and a bare number as the initial point
  (t/is (= [[0.0 10.0]] (sut/normalize-bounds :sceua [0 10] nil)))
  (t/is (= [[0.0 10.0]] (sut/normalize-bounds :sceua [[0 10]] 5)))
  ;; the range of every dimension has to be positive, not only the first one
  (t/is (rejected? :sceua [[0 1] [2 2]] nil "less than"))
  (t/is (rejected? :sceua [[2 2]] nil "less than"))
  ;; a tiny range is a legal range
  (t/is (= [[0.0 4.9E-324]] (sut/normalize-bounds :sceua [[0 Double/MIN_VALUE]] nil)))
  (t/is (rejected? :sceua [[0 1] [0 1]] [0.5] "differs")))

(t/deftest normalize-bounds-simplex
  (doseq [m [:nelder-mead :multidirectional-simplex]]
    ;; bounds optional: they only size the simplex
    (t/is (nil? (sut/normalize-bounds m nil nil)) (str m))
    (t/is (nil? (sut/normalize-bounds m nil [1 1])) (str m))
    (t/is (= [[0.0 1.0]] (sut/normalize-bounds m [[0 1]] [0.5])) (str m))
    (t/is (rejected? m [[0 ##Inf]] [0] "finite") (str m))))

(t/deftest normalize-bounds-lbfgsb-infinite
  ;; infinite (one-sided and two-sided) bounds are fine when the initial point is given
  (t/is (= [[##-Inf ##Inf] [0.0 1.0]] (sut/normalize-bounds :lbfgsb [[##-Inf ##Inf] [0 1]] [0 0])))
  (t/is (= [[0.0 ##Inf]] (sut/normalize-bounds :lbfgsb [[0 ##Inf]] [1])))
  (t/is (= [[##-Inf 0.0]] (sut/normalize-bounds :lbfgsb [[##-Inf 0]] (double-array [-1]))))
  ;; ... and an error when it is not: no midpoint can be defined
  (t/is (rejected? :lbfgsb [[##-Inf ##Inf]] nil "initial"))
  (t/is (rejected? :lbfgsb [[0 ##Inf]] nil "initial"))
  (t/is (rejected? :lbfgsb [[0 1] [##-Inf 5]] nil "initial"))
  ;; finite bounds do not need the initial point
  (t/is (= [[0.0 1.0]] (sut/normalize-bounds :lbfgsb [[0 1]] nil))))

(t/deftest normalize-bounds-unconstrained-methods
  (doseq [m [:powell :gradient :non-linear-gradient]]
    (t/is (nil? (sut/normalize-bounds m nil nil)) (str m))
    (t/is (nil? (sut/normalize-bounds m nil [1 2 3])) (str m))
    ;; bounds are ignored by the algorithm, but when given they are checked
    (t/is (= [[##-Inf ##Inf]] (sut/normalize-bounds m [[##-Inf ##Inf]] [0])) (str m))
    (t/is (rejected? m [[##-Inf ##Inf]] nil "initial") (str m))))

(t/deftest normalize-bounds-dimensions
  ;; the number of bounds has to match the initial point (powell and gradient silently accepted any length before)
  (doseq [m (remove #{:brent} all-methods)]
    (t/is (rejected? m [[0 1] [0 1]] [0 0 0] "differs") (str m))
    (t/is (rejected? m [[0 1] [0 1]] [0] "differs") (str m))
    (t/is (rejected? m [[0 1] [0 1]] [] "differs") (str m))
    (t/is (= [[0.0 1.0] [0.0 1.0]] (sut/normalize-bounds m [[0 1] [0 1]] [0.5 0.5])) (str m))
    (t/is (= [[0.0 1.0] [0.0 1.0]] (sut/normalize-bounds m [[0 1] [0 1]] (double-array [0.5 0.5]))) (str m))
    (t/is (= [[0.0 1.0] [0.0 1.0]] (sut/normalize-bounds m [[0 1] [0 1]] (list 0.5 0.5))) (str m)))
  ;; a bare number is a 1d initial point, not a length-1 sequence of anything else
  (doseq [m [:lbfgsb :powell :bobyqa]]
    (t/is (= [[0.0 1.0]] (sut/normalize-bounds m [[0 1]] 0.5)) (str m))
    (t/is (rejected? m [[0 1] [0 1]] 0.5 "differs") (str m))))

(t/deftest normalize-bounds-unknown-method
  (doseq [m [:foo nil "brent" :bfgs :l-bfgs-b]]
    (let [d (ex-data-of #(sut/normalize-bounds m [[0 1]] nil))]
      (t/is (= m (:method d)) (str m))
      (t/is (= (set all-methods) (:allowed d)) (str m)))))

(t/deftest multivariate-optimize
  ;; thin glue: result shape for both stats? values, on a trivial one-evaluation optimizer
  (let [opt (org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer. 1.0e-10 1.0e-10)
        f (org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction.
           (reify org.apache.commons.math3.analysis.MultivariateFunction
             (value [_ x] (let [a (aget ^doubles x 0)] (* (- a 2.0) (- a 2.0))))))
        data [(sut/goal-type :minimize) (sut/max-eval nil) (sut/max-iter nil) f
              (org.apache.commons.math3.optim.InitialGuess. (double-array [0.0]))]
        [p v] (sut/multivariate-optimize opt data false)
        s (sut/multivariate-optimize (org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer. 1.0e-10 1.0e-10)
                                     data true)]
    (t/is (vector? p))
    (t/is (< (Math/abs (- 2.0 (double (first p)))) 1.0e-5))
    (t/is (< v 1.0e-9))
    (t/is (= #{:point :value :evaluations :iterations} (set (keys s))))
    (t/is (pos? (:evaluations s)))))
