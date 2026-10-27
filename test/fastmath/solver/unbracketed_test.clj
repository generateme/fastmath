(ns fastmath.solver.unbracketed-test
  "The methods without a bracket (`secant`, `muller2`) fail with `:diverged`, they never return a point which is not a root (FR-9)."
  (:require [fastmath.core :as m]
            [clojure.test :as t]
            [fastmath.solver.secant :as secant]
            [fastmath.solver.muller2 :as muller2])
  (:import [clojure.lang ExceptionInfo]))

(defn- outcome
  "The root, or the status of the failure."
  [solve f lo hi & [opts]]
  (try (solve f lo hi opts) (catch ExceptionInfo e (:status (ex-data e)))))

(defn- right?
  "True when the outcome is a failure status, or a point where the function has a sign change within 1e-6."
  [f r]
  (or (keyword? r)
      (not (== (m/signum (double (f (- r 1e-6)))) (m/signum (double (f (+ r 1e-6))))))))

(t/deftest never-a-non-root
  (doseq [[nm solve] {:secant secant/secant :muller2 muller2/muller2}
          [fname f lo hi] [["(x-1)^2 (x-3) on [2,5]" (fn [x] (* (- x 1.0) (- x 1.0) (- x 3.0))) 2.0 5.0]
                           ["exp x - 1e10 on [0,100]" (fn [x] (- (m/exp x) 1e10)) 0.0 100.0]
                           ["1/x on [-1,2]" (fn [x] (/ 1.0 x)) -1.0 2.0]
                           ["sin x on [0.5,10]" (fn [x] (m/sin x)) 0.5 10.0]
                           ["sign(x-0.3)" (fn [x] (m/signum (- x 0.3))) 0.0 1.0]
                           ["tanh steep" (fn [x] (m/tanh (* 100.0 (- x 0.3)))) 0.0 1.0]
                           ["sin x on [3,4]" (fn [x] (m/sin x)) 3.0 4.0]]]
    (let [r (outcome solve f lo hi)]
      (t/is (right? f r) (str nm " " fname " gave " r))
      (when-not (keyword? r)
        (t/is (<= lo r hi) (str nm " " fname " result inside the interval"))))))

(t/deftest runaway-is-diverged
  (t/is (= :diverged (outcome secant/secant (fn [x] (m/tanh (* 100.0 (- x 0.3)))) 0.0 1.0)) "equal consecutive values")
  (t/is (= :diverged (outcome secant/secant (fn [x] (- (m/exp x) 1e10)) 0.0 100.0)) "the iteration leaves the domain of the function or does not converge to a root")
  (t/is (= :diverged (outcome muller2/muller2 (fn [x] (- (m/exp x) 1e10)) 0.0 100.0)) "an iterate leaves the interval"))

(t/deftest results-agree-with-the-bracketing-solvers
  (doseq [[nm solve] {:secant secant/secant :muller2 muller2/muller2}
          [fname f lo hi root] [["cos x - x" (fn [x] (- (m/cos x) x)) 0.0 1.0 0.7390851332151607]
                                ["x^3 - 2 on [0,3]" (fn [x] (- (* x x x) 2.0)) 0.0 3.0 1.2599210498948732]
                                ["atan x - 1" (fn [x] (- (m/atan x) 1.0)) 0.0 10.0 1.5574077246549023]]]
    (t/is (<= (m/abs (- (solve f lo hi) root)) 1e-13) (str nm " " fname))))
