(ns fastmath.solver.iteration-test
  "Newton, Halley, fixed point and the evaluation counts of TOMS748."
  (:require [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.solver.newton :as newton]
            [fastmath.solver.halley :as halley]
            [fastmath.solver.fixed-point :as fp]
            [fastmath.solver.brent :as brent]
            [fastmath.solver.toms748 :as toms748])
  (:import [clojure.lang ExceptionInfo]))

(defn- thrown-data [thunk] (try (thunk) ::no-throw (catch ExceptionInfo e (ex-data e))))

(def f (fn [^double x] (- (m/cos x) x)))
(def df (fn [^double x] (- (- (m/sin x)) 1.0)))
(def d2f (fn [^double x] (- (m/cos x))))
(def root 0.7390851332151607)

(t/deftest newton-halley
  (doseq [[nm solve] {:newton newton/newton :halley halley/halley}]
    (let [i (solve f 0.0 1.0 {:df df :d2f d2f :info? true})]
      (t/is (<= (m/abs (- (:root i) root)) 1e-15) (str nm " analytic derivatives"))
      (t/is (<= (:evaluations i) 12) (str nm " fast convergence, evaluations " (:evaluations i))))
    (let [i (solve f 0.0 1.0 {:info? true})]
      (t/is (<= (m/abs (- (:root i) root)) 1e-13) (str nm " finite difference derivatives"))
      (t/is (= :converged (:status i))))
    (t/is (<= (m/abs (- (solve f 1.0 0.0 {:df df :d2f d2f}) root)) 1e-15) (str nm " reversed bounds"))
    (t/is (<= (m/abs (- (solve f 0.0 1.0 {:df df :d2f d2f :initial-value 0.9}) root)) 1e-15) (str nm " start"))
    (t/is (== 1.0 (solve (fn [x] (- x 1.0)) 1.0 3.0 {:df (fn [_] 1.0) :d2f (fn [_] 0.0)})) (str nm " a root at the end")))
  (t/is (< (:evaluations (newton/newton f 0.0 1.0 {:df df :info? true}))
           (:evaluations (brent/brent f 0.0 1.0 {:info? true}))) "Newton needs fewer evaluations than Brent on a smooth function"))

(t/deftest safeguard
  ;; a zero derivative at the start, a derivative that sends Newton out of the bracket, a derivative which is NaN
  (let [g (fn [^double x] (- (* x x x) 2.0 (* 0.0 x)))]
    (t/is (<= (m/abs (- (newton/newton g -3.0 3.0 {:df (fn [^double x] (* 3.0 x x)) :initial-value 0.0}) (m/cbrt 2.0))) 1e-14) "zero derivative at the start")
    (t/is (<= (m/abs (- (newton/newton g 0.0 3.0 {:df (fn [_] ##NaN)}) (m/cbrt 2.0))) 1e-14) "NaN derivative is a bisection")
    (t/is (<= (m/abs (- (newton/newton g 0.0 3.0 {:df (fn [_] ##Inf)}) (m/cbrt 2.0))) 1e-14) "infinite derivative")
    (t/is (<= (m/abs (- (halley/halley g 0.0 3.0 {:df (fn [_] 0.0) :d2f (fn [_] 0.0)}) (m/cbrt 2.0))) 1e-14) "zero derivatives"))
  ;; arctan x: Newton diverges from x0 = 2 without a bracket, here it is safeguarded
  (t/is (<= (m/abs (newton/newton (fn [^double x] (m/atan x)) -10.0 20.0 {:df (fn [^double x] (/ 1.0 (+ 1.0 (* x x)))) :initial-value 15.0})) 1e-12) "atan, far start")
  ;; a flat multiple root: x^3 converges linearly but converges
  (t/is (<= (m/abs (newton/newton (fn [^double x] (* x x x)) -1.0 2.0 {:df (fn [^double x] (* 3.0 x x)) :absolute-accuracy 1e-10})) 1e-9) "triple root"))

(t/deftest newton-options
  (t/is (= :df (:option (thrown-data #(newton/newton f 0.0 1.0 {:df 5})))))
  (t/is (= :d2f (:option (thrown-data #(halley/halley f 0.0 1.0 {:d2f "x"})))))
  (t/is (= :invalid-result (:status (thrown-data #(newton/newton f 0.0 1.0 {:df (fn [_] nil)})))))
  (t/is (= :no-bracket (:status (thrown-data #(newton/newton (fn [x] (+ 1.0 (* x x))) -1.0 1.0 {:df (fn [x] (* 2.0 x))})))))
  (t/is (= :initial-value (:option (thrown-data #(newton/newton f 0.0 1.0 {:initial-value 1.0})))))
  (t/is (= :max-evals (:status (thrown-data #(newton/newton f 0.0 1.0 {:max-iters 3 :absolute-accuracy 0.0 :relative-accuracy 0.0}))))
        "the finite differences count to the limit")
  (t/is (= :max-evals (:status (newton/newton f 0.0 1.0 {:max-iters 3 :absolute-accuracy 0.0 :relative-accuracy 0.0 :info? true})))))

(t/deftest fixed-point
  (let [g (fn [^double x] (m/cos x))]
    (t/is (<= (m/abs (- (fp/fixed-point g 1.0) root)) 1e-13) "cos x = x")
    (t/is (<= (m/abs (- (fp/fixed-point g 1.0 {:steffensen? true}) root)) 1e-15) "Steffensen")
    (let [plain (fp/fixed-point g 1.0 {:info? true})
          st (fp/fixed-point g 1.0 {:steffensen? true :info? true})]
      (t/is (< (:evaluations st) (:evaluations plain)) "Steffensen needs fewer evaluations")))
  (t/is (== 2.0 (fp/fixed-point (fn [_] 2.0) 5.0)) "a constant")
  (t/is (m/delta-eq 2.0 (fp/fixed-point (fn [^double x] (m/sqrt (+ 2.0 x))) 0.0) 1e-12) "sqrt(2 + x) = x")
  (t/is (= 1.0 (fp/fixed-point (fn [x] x) 1.0)) "every point is a fixed point of the identity")
  (t/is (= :max-evals (:status (thrown-data #(fp/fixed-point (fn [x] (- 1.0 x)) 0.0)))) "an oscillation never converges")
  (t/is (= :diverged (:status (thrown-data #(fp/fixed-point (fn [x] (* 10.0 x)) 1.0 {:max-iters 100000})))) "growth to infinity")
  (t/is (= :nan (:status (thrown-data #(fp/fixed-point (fn [_] ##NaN) 1.0)))))
  (t/is (= :diverged (:status (fp/fixed-point (fn [x] (* 10.0 x)) 1.0 {:info? true :max-iters 100000}))))
  (t/is (= :x0 (:option (thrown-data #(fp/fixed-point (fn [x] x) ##NaN)))))
  (t/is (= :g (:option (thrown-data #(fp/fixed-point nil 1.0)))))
  (t/is (= :max-iters (:option (thrown-data #(fp/fixed-point (fn [x] x) 1.0 {:max-iters 0}))))))

;; evaluation counts of scipy.optimize.toms748 (utils/fastmath/dev/generate_solver_toms748_counts.py), xtol 2e-12, rtol 8.9e-16
(t/deftest toms748-counts
  (let [fs {"cos x - x" [(fn [x] (- (m/cos x) x)) 0.0 1.0 8]
            "x^3 - 2" [(fn [x] (- (* x x x) 2.0)) 0.0 3.0 12]
            "tanh steep" [(fn [x] (m/tanh (* 100.0 (- x 0.3)))) 0.0 1.0 14]
            "(x-1)^3" [(fn [x] (let [d (- x 1.0)] (* d d d))) 0.0 3.0 89]
            "exp x - 1e10" [(fn [x] (- (m/exp x) 1e10)) 0.0 100.0 25]
            "atan x - 1" [(fn [x] (- (m/atan x) 1.0)) 0.0 10.0 13]
            "x^3 + x" [(fn [x] (+ (* x x x) x)) -1.0 2.0 12]
            "1/x - 2" [(fn [x] (- (/ 1.0 x) 2.0)) 0.1 5.0 14]
            "(x-1)^2 (x-3)" [(fn [x] (* (- x 1.0) (- x 1.0) (- x 3.0))) 2.0 5.0 11]
            "x^5 - x - 1" [(fn [x] (- (m/pow x 5) x 1.0)) 1.0 2.0 13]}]
    (doseq [[nm [g lo hi scipy-count]] fs]
      (let [i (toms748/toms748 g lo hi {:absolute-accuracy 2e-12 :relative-accuracy 8.9e-16 :info? true})]
        (t/is (<= (:evaluations i) (* 1.25 scipy-count)) (str nm ": " (:evaluations i) " evaluations, scipy " scipy-count))))))
