(ns fastmath.solver-test
  (:require [fastmath.core :as m]
            [clojure.test :as t]
            [fastmath.solver :as solver])
  (:import [clojure.lang ExceptionInfo]))

(defn- thrown-data [thunk] (try (thunk) ::no-throw (catch ExceptionInfo e (ex-data e))))

(def f (fn [^double x] (- (m/cos x) x)))
(def root 0.7390851332151607)

(t/deftest find-root-dispatch
  (t/is (= (set (keys solver/solvers))
           #{:brent :bisection :illinois :pegasus :regula-falsi :ridders :muller :muller2 :secant :toms748 :newton :halley}))
  (doseq [s (keys solver/solvers)]
    (t/is (<= (m/abs (- (solver/find-root f 0.0 1.0 {:solver s}) root)) 1e-13) (str s " through find-root"))
    (t/is (== ((get solver/solvers s) f 0.0 1.0) (solver/find-root f 0.0 1.0 {:solver s})) (str s " same result as the direct call")))
  (t/is (== (solver/find-root f 0.0 1.0) (solver/find-root f 0.0 1.0 {:solver :brent})) "the default is brent")
  (t/is (== (solver/find-root f 0.0 1.0) (solver/find-root f 0.0 1.0 {:solver nil})) "nil solver is the default")
  (t/is (= :solver (:option (thrown-data #(solver/find-root f 0.0 1.0 {:solver :newtons})))) "unknown solver throws")
  (t/is (= :solver (:option (thrown-data #(solver/find-root f 0.0 1.0 {:solver "brent"})))))
  (t/is (= :solver (:option (thrown-data #(solver/find-root f 0.0 1.0 {:solver 5})))))
  (t/is (= :no-bracket (:status (thrown-data #(solver/find-root (fn [x] (+ 1.0 (* x x))) -1.0 2.0)))))
  (t/is (= {:root 0.5 :status :converged} (select-keys (solver/find-root (fn [x] (- x 0.5)) 0.0 1.0 {:info? true :solver :bisection}) [:root :status])))
  (t/is (= root (solver/find-root f 0.0 1.0 {:df (fn [x] (- (- (m/sin x)) 1.0)) :solver :newton :absolute-accuracy 1e-16})) "newton through find-root with a derivative"))

#_{:clj-kondo/ignore [:unresolved-var]}
(t/deftest exports
  (t/is (= solver/brent @(resolve 'fastmath.solver.brent/brent)))
  (doseq [s [#'solver/brent #'solver/bisection #'solver/toms748 #'solver/newton #'solver/halley #'solver/fixed-point]]
    (t/is (string? (:doc (meta s))) (str s " has the docstring"))
    (t/is (seq (:arglists (meta s))) (str s " has arglists"))))

(t/deftest find-bracket-finite
  (let [seen (atom ##Inf)
        d (thrown-data #(solver/find-bracket (fn [x] (swap! seen min (- (m/abs x))) (+ 1.0 (* x x))) 0.0 {:max-iters 5000}))]
    (t/is (= :no-bracket (:status d)) "the search ends at the finite limits")
    (t/is (> @seen ##-Inf) "f is never called at an infinite point"))
  (t/is (= :upper (:option (thrown-data #(solver/find-bracket identity 0.0 {:upper -1.0})))) "the violated limit is named")
  (t/is (= :lower (:option (thrown-data #(solver/find-bracket identity 0.0 {:lower 1.0}))))))

(t/deftest find-bracket
  (let [[a b] (solver/find-bracket (fn [x] (- x 100.0)) 0.0)]
    (t/is (< a 100.0 b) "grows until the sign changes")
    (t/is (< a b)))
  (let [[a b] (solver/find-bracket (fn [x] (- x 3.0)) 0.0 {:step 1.0 :factor 2.0})]
    (t/is (<= a 3.0 b)))
  (t/is (= [-1.0 1.0] (solver/find-bracket (fn [x] (- x 1.0)) 0.0 {:step 1.0})) "a zero value at an end counts")
  (t/is (= [-1.0 1.0] (solver/find-bracket (fn [x] x) 0.0 {:step 1.0})) "a sign change at the first interval")
  (t/is (every? #(<= (m/abs (- % 1e5)) 1e6) (solver/find-bracket (fn [x] (- x 1e5)) 0.0)))
  (t/is (= :max-evals (:status (thrown-data #(solver/find-bracket (fn [x] (+ 1.0 (* x x))) 0.0)))) "no sign change within the evaluation limit")
  (t/is (= :no-bracket (:status (thrown-data #(solver/find-bracket (fn [x] (- x 100.0)) 0.0 {:lower -10.0 :upper 10.0})))) "limits")
  (t/is (= :max-evals (:status (thrown-data #(solver/find-bracket (fn [x] (- x 100.0)) 0.0 {:upper 10.0})))) "one limit: the other end runs away")
  (t/is (= :nan (:status (thrown-data #(solver/find-bracket (fn [_] ##NaN) 0.0)))))
  (t/is (= :x0 (:option (thrown-data #(solver/find-bracket (fn [x] x) ##NaN)))))
  (t/is (= :step (:option (thrown-data #(solver/find-bracket (fn [x] x) 0.0 {:step 0.0})))))
  (t/is (= :factor (:option (thrown-data #(solver/find-bracket (fn [x] x) 0.0 {:factor 1.0})))))
  (t/is (= :max-iters (:option (thrown-data #(solver/find-bracket (fn [x] x) 0.0 {:max-iters 0})))))
  (t/is (= :f (:option (thrown-data #(solver/find-bracket nil 0.0)))))
  ;; the point is outside the limits
  (t/is (some? (:option (thrown-data #(solver/find-bracket (fn [x] x) 0.0 {:lower 1.0}))))))

(t/deftest find-roots
  (let [rs (solver/find-roots (fn [x] (m/sin x)) 0.5 10.0)]
    (t/is (= 3 (count rs)) "pi, 2pi, 3pi")
    (t/is (every? true? (map #(<= (m/abs (- %1 %2)) 1e-14) rs [m/PI (* 2 m/PI) (* 3 m/PI)])))
    (t/is (= rs (sort rs))))
  (t/is (= 10 (count (solver/find-roots (fn [x] (m/sin x)) 0.1 31.5))) "sin x: pi ... 10 pi")
  (t/is (= [] (solver/find-roots (fn [x] (+ 1.0 (* x x))) -5.0 5.0)) "no roots")
  (t/is (= [0.0] (solver/find-roots (fn [x] x) -1.0 1.0)) "a grid point is a root")
  (t/is (= [1.0] (solver/find-roots (fn [x] (- x 1.0)) 0.0 1.0)) "a root at the upper end")
  (t/is (= [0.0] (solver/find-roots (fn [x] x) 0.0 1.0)) "a root at the lower end")
  (let [rs (solver/find-roots (fn [x] (* (+ x 2.0) (- x 0.5) (- x 3.0))) -5.0 5.0)]
    (t/is (= 3 (count rs)))
    (t/is (every? true? (map #(<= (m/abs (- %1 %2)) 1e-14) rs [-2.0 0.5 3.0]))))
  ;; a pole is a sign change but not a root
  (t/is (= [] (solver/find-roots (fn [x] (/ 1.0 x)) -1.0 2.0)) "1/x")
  (t/is (= 3 (count (solver/find-roots (fn [x] (m/tan x)) -4.0 4.0))) "tan x: roots -pi, 0, pi, not the poles")
  ;; close roots need a finer grid
  (let [g (fn [x] (* (- x 1.00001) (- x 1.00011)))]
    (t/is (= [] (solver/find-roots g 0.0 2.0 {:samples 10})) "two roots in one part are missed")
    (t/is (= 2 (count (solver/find-roots g 0.0 2.0 {:samples 100000}))) "a finer grid resolves them"))
  ;; options are passed to the solver
  (t/is (= 1 (count (solver/find-roots (fn [x] (- x 0.3)) 0.0 1.0 {:solver :bisection :absolute-accuracy 1e-3}))))
  (t/is (= :solver (:option (thrown-data #(solver/find-roots (fn [x] x) -1.0 2.0 {:solver :nope})))))
  (t/is (= :samples (:option (thrown-data #(solver/find-roots (fn [x] x) -1.0 2.0 {:samples 0})))))
  (t/is (= :samples (:option (thrown-data #(solver/find-roots (fn [x] x) -1.0 2.0 {:samples 10.5})))))
  (t/is (= :upper-bound (:option (thrown-data #(solver/find-roots (fn [x] x) 1.0 1.0)))) "empty interval")
  (t/is (= :upper-bound (:option (thrown-data #(solver/find-roots (fn [x] x) 2.0 1.0)))) "reversed interval")
  (t/is (= :lower-bound (:option (thrown-data #(solver/find-roots (fn [x] x) ##-Inf 1.0)))))
  (t/is (= :f (:option (thrown-data #(solver/find-roots nil 0.0 1.0)))))
  ;; NaN regions are skipped
  (t/is (= 1 (count (solver/find-roots (fn [x] (if (< x -0.5) ##NaN (- x 0.3))) -1.0 1.0)))))
