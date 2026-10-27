(ns fastmath.solver.common-test
  (:require [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.solver.common :as c])
  (:import [clojure.lang ExceptionInfo]
           [fastmath.solver.common Problem]))

(t/deftest user-exceptions-are-not-failures
  (let [f (fn [_] (throw (ex-info "mine" {:status :nan})))]
    (t/is (= "mine" (try (c/solve-bracketed (fn [_] 0.0) false f 0.0 1.0 {:info? true}) (catch ExceptionInfo e (ex-message e))))
          "an ex-info of the user function is rethrown also with :info?"))
  (t/is (= :nan (:status (c/solve-bracketed (fn [_] 0.0) false (fn [_] ##NaN) 0.0 1.0 {:info? true}))) "a solver failure is returned"))

(defn- thrown-data
  "Data of the `ex-info` thrown by `thunk`, or `::no-throw`."
  [thunk]
  (try (thunk) ::no-throw (catch ExceptionInfo e (ex-data e))))

(defn- halving
  "A minimal core: bisection down to the tolerance, used to test the driver."
  [^Problem p]
  (let [ev (.ev p)]
    (loop [a (.a p) fa (.fa p) b (.b p)]
      (let [m (+ (* 0.5 a) (* 0.5 b))]
        (if (<= (- b a) (* 2.0 (c/tolerance (.atol p) (.rtol p) m)))
          m
          (let [fm (.eval ev m)]
            (cond (zero? fm) m
                  (= (m/signum fm) (m/signum fa)) (recur m fm b)
                  :else (recur a fa m))))))))

(defn- solve [f lo hi & [opts]] (c/solve-bracketed halving false f lo hi opts))

(t/deftest parse-options-test
  (let [o (c/parse-options nil)]
    (t/is (= 1200 (:max-evals o)) "default evaluation limit")
    (t/is (== c/default-relative-accuracy (:rtol o)))
    (t/is (== c/default-absolute-accuracy (:atol o)))
    (t/is (m/nan? (:start o)) "no initial value is NaN")
    (t/is (false? (:info? o))))
  (let [o (c/parse-options {:absolute-accuracy 0 :relative-accuracy 1e-3 :max-iters 7 :initial-value 2 :info? 1})]
    (t/is (== 0.0 (:atol o)) "zero accuracy is allowed")
    (t/is (== 1e-3 (:rtol o)))
    (t/is (= 7 (:max-evals o)))
    (t/is (== 2.0 (:start o)))
    (t/is (true? (:info? o)) "truthiness"))
  (doseq [[k v] [[:absolute-accuracy -1.0] [:absolute-accuracy ##NaN] [:absolute-accuracy ##Inf] [:absolute-accuracy "1e-3"]
                 [:relative-accuracy -1e-16] [:relative-accuracy ##NaN]
                 [:max-iters 0] [:max-iters -5] [:max-iters 50.0] [:max-iters 1/2] [:max-iters "50"] [:max-iters 10000000000000000000N]
                 [:initial-value ##NaN] [:initial-value ##Inf] [:initial-value "1"]]]
    (let [d (thrown-data #(c/parse-options {k v}))]
      (t/is (= k (:option d)) (str k " " (pr-str v) " is rejected"))
      (t/is (not (contains? d :status)))))
  (t/is (= :options (:option (thrown-data #(c/parse-options [:max-iters 5])))) "options have to be a map")
  (t/is (= 5 (:max-evals (c/parse-options {:max-iters (int 5)}))) "int accepted")
  (t/is (= 5 (:max-evals (c/parse-options {:max-iters (short 5)}))) "short accepted")
  (t/is (= 1200 (:max-evals (c/parse-options {:max-iters nil}))) "nil means default (no NPE)")
  (t/is (= 1 (:max-evals (c/parse-options {:max-iters 1}))) "the smallest limit"))

(t/deftest tolerance-test
  (t/is (== (+ 1e-6 (* 1e-3 10.0)) (c/tolerance 1e-6 1e-3 10.0)) "abs + rel |x|")
  (t/is (== (+ 1e-6 (* 1e-3 10.0)) (c/tolerance 1e-6 1e-3 -10.0)) "symmetric")
  (t/is (== (* 2.0 (m/ulp 1.0)) (c/tolerance 0.0 0.0 1.0)) "floor of two ulps")
  (t/is (pos? (c/tolerance 0.0 0.0 0.0)) "positive at zero")
  (t/is (== (* 2.0 (m/ulp 1e300)) (c/tolerance 0.0 0.0 1e300)) "huge x"))

(t/deftest evaluator-test
  (let [ev (c/evaluator (fn [x] (* 2.0 x)) 3)]
    (t/is (== 2.0 (.eval ev 1.0)))
    (t/is (== 0 (.eval ev 0.0)))
    (t/is (= 2 (.evaluations ev)))
    (t/is (== ##Inf (.eval (c/evaluator (fn [_] ##Inf) 3) 1.0)) "infinity is returned as it is")
    (t/is (== 3.0 (.eval (c/evaluator (fn [_] 3) 3) 1.0)) "long is coerced")
    (t/is (== 0.5 (.eval (c/evaluator (fn [_] 1/2) 3) 1.0)) "ratio is coerced"))
  (let [ev (c/evaluator (fn [x] x) 2)]
    (.eval ev 1.0)
    (.eval ev 1.0)
    (let [d (thrown-data #(.eval ev 1.0))]
      (t/is (= :max-evals (:status d)))
      (t/is (= 2 (:evaluations d)))
      (t/is (= 2 (:limit d)))
      (t/is (= 2 (.evaluations ev)) "no evaluation is made after the limit")))
  (let [ev (c/evaluator (fn [x] x) 1)]
    (t/is (== 1.0 (.eval ev 1.0)) "the limit 1 allows one evaluation")
    (t/is (= :max-evals (:status (thrown-data #(.eval ev 1.0))))))
  (let [d (thrown-data #(.eval (c/evaluator (fn [_] ##NaN) 5) 1.5))]
    (t/is (= :nan (:status d)))
    (t/is (= 1.5 (:x d))))
  (doseq [bad [nil "1" :a [1]]]
    (t/is (= :invalid-result (:status (thrown-data #(.eval (c/evaluator (fn [_] bad) 5) 0.0))))
          (str (pr-str bad) " is not a number"))))

(t/deftest solve-bracketed-test
  ;; the roots
  (t/is (m/delta-eq 1.0 (solve (fn [x] (- x 1.0)) 0.0 3.0) 1e-12))
  (t/is (m/delta-eq 1.0 (solve (fn [x] (- x 1.0)) 3.0 0.0) 1e-12) "reversed bounds are swapped")
  (t/is (m/delta-eq 1.0 (solve (fn [x] (- x 1)) 0 3) 1e-12) "long bounds")
  (t/is (m/delta-eq 1.0 (solve (fn [x] (- x 1)) 0 3.0M) 1e-12) "BigDecimal bound")
  ;; endpoint roots cost at most two evaluations
  (doseq [[lo hi r] [[1.0 3.0 1.0] [0.0 1.0 1.0] [-1.0 0.0 0.0] [0.0 2.0 0.0]]]
    (let [info (solve (fn [x] (- x r)) lo hi {:info? true})]
      (t/is (== r (:root info)) (str "root at an end of [" lo ", " hi "]"))
      (t/is (<= (:evaluations info) 2))
      (t/is (= :converged (:status info)))))
  ;; equal bounds
  (t/is (== 2.0 (solve (fn [x] (- x 2.0)) 2.0 2.0)) "equal bounds at the root")
  (t/is (= :no-bracket (:status (thrown-data #(solve (fn [x] (- x 2.0)) 1.0 1.0)))) "equal bounds, not a root")
  ;; no sign change
  (let [d (thrown-data #(solve (fn [x] (+ 1.0 (* x x))) -1.0 2.0))]
    (t/is (= :no-bracket (:status d)))
    (t/is (= 2 (:evaluations d)))
    (t/is (= [-1.0 2.0] [(:a d) (:b d)])))
  ;; invalid arguments
  (doseq [[lo hi which] [[##Inf 1.0 :lower-bound] [0.0 ##Inf :upper-bound] [##NaN 1.0 :lower-bound] [0.0 ##NaN :upper-bound]
                         [nil 1.0 :lower-bound] [0.0 "1" :upper-bound] [##-Inf 1.0 :lower-bound]]]
    (t/is (= which (:option (thrown-data #(solve (fn [x] x) lo hi)))) (str "bounds " lo " " hi)))
  (t/is (= :f (:option (thrown-data #(solve nil 0.0 1.0)))) "nil function")
  (t/is (= :f (:option (thrown-data #(solve 5 0.0 1.0)))) "a number is not a function")
  (t/is (= :max-iters (:option (thrown-data #(solve (fn [x] x) -1.0 1.0 {:max-iters 0})))))
  ;; the tolerance is used
  (t/is (> (m/abs (- 1.0 (solve (fn [x] (- x 1.0)) 0.0 3.7 {:absolute-accuracy 1e-2 :relative-accuracy 0.0}))) 1e-12)
        "a loose tolerance stops early")
  ;; budget and info
  (let [info (solve (fn [x] (- x 1.0)) 0.0 3.0 {:max-iters 5 :info? true})]
    (t/is (= :max-evals (:status info)))
    (t/is (m/nan? (:root info)))
    (t/is (= 5 (:evaluations info))))
  (t/is (= :max-evals (:status (thrown-data #(solve (fn [x] (- x 1.0)) 0.0 3.0 {:max-iters 5})))))
  (t/is (= :nan (:status (solve (fn [x] (if (> x 0.5) ##NaN (- x 1.0))) 0.0 3.0 {:info? true}))))
  (t/is (= :invalid-result (:status (solve (fn [x] (if (> x 0.5) nil (- x 1.0))) 0.0 3.0 {:info? true}))))
  (t/is (= :no-bracket (:status (solve (fn [_] 1.0) 0.0 3.0 {:info? true}))))
  (t/is (= :absolute-accuracy (:option (thrown-data #(solve (fn [x] x) -1.0 1.0 {:info? true :max-iters 0 :absolute-accuracy -1}))))
        "invalid options throw even with :info?")
  ;; limit 1 cannot even evaluate both ends
  (t/is (= :max-evals (:status (solve (fn [x] (- x 1.0)) 0.0 3.0 {:max-iters 1 :info? true}))))
  ;; the first evaluation is the lower end, so a root at the lower end needs one
  (t/is (== 0.0 (solve (fn [x] x) 0.0 3.0 {:max-iters 1})) "one evaluation is enough for a root at the lower end")
  ;; initial value
  (t/is (m/delta-eq 1.0 (solve (fn [x] (- x 1.0)) 0.0 3.0 {:initial-value 50.0}) 1e-12) "ignored when the solver does not use it")
  (t/is (= :initial-value (:option (thrown-data #(c/solve-bracketed halving true (fn [x] (- x 1.0)) 0.0 3.0 {:initial-value 3.0})))))
  (t/is (= :initial-value (:option (thrown-data #(c/solve-bracketed halving true (fn [x] (- x 1.0)) 0.0 3.0 {:initial-value 0.0})))))
  (t/is (= :initial-value (:option (thrown-data #(c/solve-bracketed halving true (fn [x] (- x 1.0)) 0.0 3.0 {:initial-value -1.0})))))
  (t/is (m/delta-eq 1.0 (c/solve-bracketed halving true (fn [x] (- x 1.0)) 0.0 3.0 {:initial-value 2.0}) 1e-12)))

(defn- problem [ev lo hi atol rtol] (Problem. ev lo 0.0 hi 0.0 atol rtol ##NaN))

(t/deftest certified-test
  (let [f (fn [x] (- x 1.0))
        ev (c/evaluator f 100)]
    (t/is (c/certified? ev 1.0 0.0 (problem ev 0.0 3.0 1e-9 0.0)) "zero value")
    (t/is (c/certified? ev 1.0000000001 1e-10 (problem ev 0.0 3.0 1e-9 0.0)) "sign change at 2 tol")
    (t/is (not (c/certified? ev 1.5 0.5 (problem ev 0.0 3.0 1e-9 0.0))) "far from the root"))
  (let [ev (c/evaluator (fn [x] (* (- x 1.0) (- x 1.0))) 100)]
    (t/is (not (c/certified? ev 1.0000001 1e-14 (problem ev 0.0 3.0 1e-9 0.0))) "even multiplicity has no sign change"))
  (let [ev (c/evaluator (fn [x] (- x 1.0)) 100)]
    (t/is (not (c/certified? ev 3.0 2.0 (problem ev 0.0 3.0 1e-9 0.0))) "at the end of the interval the outer probe is skipped, the inner one has the same sign")))
