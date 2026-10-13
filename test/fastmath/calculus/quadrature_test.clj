(ns fastmath.calculus.quadrature-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.calculus.quadrature :as sut]))

;; Reference values: QuadGK.jl, `utils/fastmath/dev/generate_gk_reference.jl`.

(def ^:private integrands
  {"poly" (fn ^double [^double x] (m/sq x))
   "const" (fn ^double [^double _] 1.0)
   "sqrtabs" (fn ^double [^double x] (m/sqrt (m/abs x)))
   "abssin20" (fn ^double [^double x] (m/abs (m/sin (m/* 20.0 x))))
   "invsqrt" (fn ^double [^double x] (m// (m/sqrt x)))
   "logx" (fn ^double [^double x] (m/log x))
   "expinf" (fn ^double [^double x] (m/exp (m/- x)))
   "halfexp" (fn ^double [^double x] (m// (m/exp (m/- x)) (m/sqrt x)))
   "cauchy" (fn ^double [^double x] (m// (m/inc (m/sq x))))
   "osc" (fn ^double [^double x] (m/sin (m/* 50.0 x)))
   "runge" (fn ^double [^double x] (m// (m/inc (m/* 25.0 x x))))
   "gauss" (fn ^double [^double x] (m/exp (m/- (m/sq x))))
   "peak" (fn ^double [^double x] (m/exp (m/* -1.0e6 x x)))})

(def ^:private gk-reference (delay (edn/read-string (slurp (io/resource "calculus/gk_reference.edn")))))
(def ^:private kronrod-reference (delay (edn/read-string (slurp (io/resource "calculus/kronrod_reference.edn")))))

(defn- counted
  "Returns `[f' calls]` where `calls` is an atom counting the calls of `f'`."
  [f]
  (let [calls (atom 0)]
    [(fn ^double [^double x] (swap! calls inc) (f x)) calls]))

(defn- ex-data-of [thunk]
  (try (thunk) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))

(t/deftest kronrod-rules
  ;; measured max difference to QuadGK.kronrod: 5.6e-16 (nodes and weights, orders 1..50)
  (doseq [{:keys [order x w gw]} @kronrod-reference]
    (let [rule (sut/kronrod order)]
      (t/testing (str "order " order)
        (t/is (= (count x) (count (:x rule))))
        (t/is (= (count w) (count (:w rule))))
        (t/is (= (count gw) (count (:gw rule))))
        (t/is (every? true? (map #(m/delta-eq %1 %2 1.0e-15) x (:x rule))))
        (t/is (every? true? (map #(m/delta-eq %1 %2 1.0e-15) w (:w rule))))
        (t/is (every? true? (map #(m/delta-eq %1 %2 1.0e-15) gw (:gw rule))))
        (t/is (= (m/inc (m/* 4 order)) (m/dec (:gk-evals rule))) "gk-evals counts the evaluations of two segments")))))

(t/deftest kronrod-order-one
  (let [{:keys [x w gw]} (sut/kronrod 1)]
    (t/is (every? true? (map #(m/delta-eq %1 %2 1.0e-15) [(m/- (m/sqrt 0.6)) 0.0] x)))
    (t/is (every? true? (map #(m/delta-eq %1 %2 1.0e-15) [(m// 5.0 9.0) (m// 8.0 9.0)] w)))
    (t/is (every? true? (map #(m/delta-eq %1 %2 1.0e-15) [2.0] gw))))
  (t/testing "the 3-point rule integrates polynomials of degree 5 exactly"
    (t/is (m/delta-eq (m// 1.0 3.0) (sut/gk-quadrature (integrands "poly") 0.0 1.0 {:integration-points 1}) 1.0e-15))
    (t/is (m/delta-eq (m// 1.0 6.0) (sut/gk-quadrature (fn [^double x] (m/pow x 5.0)) 0.0 1.0 {:integration-points 1}) 1.0e-15))))

(t/deftest reference-cases
  ;; measured: :result within 1e-15 relative; :error differs by cancellation noise of K-G (at most ~4 eps * |I|)
  (doseq [{:keys [id lower upper rel abs order points initdiv result error evaluations]} @gk-reference]
    (let [[f calls] (counted (integrands id))
          res (sut/gk-quadrature f lower upper (cond-> {:rel rel :abs abs :integration-points order :initdiv initdiv
                                                        :max-iters 100000 :info? true}
                                                 (seq points) (assoc :points points)))
          label (str id " [" lower "," upper "] order " order " points " points " initdiv " initdiv)]
      (t/testing label
        (t/is (m/delta-eq result (:result res) 1.0e-12 1.0e-12) "result")
        (t/is (m/delta-eq error (:error res) (m/* 1.0e-13 (m/max 1.0 (m/abs result))) 1.0e-5) "error")
        (t/is (= evaluations (:evaluations res)) "evaluations of the reference")
        (t/is (= @calls (:evaluations res)) "evaluations are the calls of the integrand")
        (t/is (false? (:fail? res)))))))

(t/deftest limits-and-failure
  (let [hard (integrands "abssin20")
        easy (integrands "poly")]
    (t/testing "converged calls report :fail? false"
      (t/is (= false (:fail? (sut/gk-quadrature easy 0.0 1.0 {:info? true}))))
      (t/is (contains? (sut/gk-quadrature easy 0.0 1.0 {:info? true}) :fail?)))
    (t/testing ":max-iters bounds the rounds including the initial one"
      (doseq [mi [1 2 3 10]]
        (let [res (sut/gk-quadrature hard -1.0 1.0 {:rel 1.0e-12 :max-iters mi :info? true})]
          (t/is (= mi (:iterations res)))
          (t/is (= :max-iters (:fail? res)))
          (t/is (= (m/* 15 (m/inc (m/* 2 (m/dec mi)))) (:evaluations res))))))
    (t/testing "convergence in the last allowed round wins over the limit"
      (t/is (= false (:fail? (sut/gk-quadrature easy 0.0 1.0 {:max-iters 1 :info? true})))))
    (t/testing ":max-evals is never exceeded, unless the initial round alone is above it"
      (doseq [me [100 200 1000]]
        (let [res (sut/gk-quadrature hard -1.0 1.0 {:rel 1.0e-12 :max-evals me :info? true})]
          (t/is (<= (:evaluations res) me))
          (t/is (> (+ (:evaluations res) 30) me) "no more room for another split")
          (t/is (= :max-evals (:fail? res)))))
      (let [res (sut/gk-quadrature hard -1.0 1.0 {:rel 1.0e-12 :max-evals 10 :initdiv 4 :info? true})]
        (t/is (= 60 (:evaluations res)))
        (t/is (= :max-evals (:fail? res)))))
    (t/testing "a failed call still returns the best estimate"
      (let [res (sut/gk-quadrature hard -1.0 1.0 {:rel 1.0e-14 :max-iters 40 :info? true})]
        (t/is (m/delta-eq 1.2591918 (:result res) 1.0e-3))
        (t/is (= :max-iters (:fail? res)))))
    (t/testing "without :info? a number is returned"
      (t/is (number? (sut/gk-quadrature hard -1.0 1.0 {:max-iters 2}))))
    (t/testing "a NaN integrand never converges and is refined until the limit"
      (let [res (sut/gk-quadrature (fn ^double [^double x] (if (m/< x 0.3) ##NaN 1.0)) 0.0 1.0 {:max-iters 6 :info? true})]
        (t/is (m/nan? (:result res)))
        (t/is (= :max-iters (:fail? res)))))))

(t/deftest tied-errors
  (t/testing "segments with equal (zero) errors are all kept"
    (doseq [k [2 3 8 17]]
      (let [res (sut/gk-quadrature (integrands "poly") 0.0 1.0 {:initdiv k :info? true})]
        (t/is (= k (:subdivisions res)))
        (t/is (m/delta-eq (m// 1.0 3.0) (:result res) 1.0e-15)))))
  (t/testing "a symmetric integrand has equal errors in both halves"
    (t/is (m/delta-eq m/PI (sut/gk-quadrature (integrands "cauchy") ##-Inf ##Inf {:rel 1.0e-12}) 1.0e-12))
    (t/is (m/delta-eq 1.7724538509055159 (sut/gk-quadrature (integrands "gauss") -8.0 8.0 {:rel 1.0e-12}) 1.0e-12))))

(t/deftest bounds
  (let [f (integrands "expinf")]
    (t/testing "reversed bounds negate the result"
      (t/is (m/delta-eq -0.5 (sut/gk-quadrature (fn ^double [^double x] x) 1.0 0.0) 1.0e-15))
      (let [a (sut/gk-quadrature f 0.0 2.0 {:info? true})
            b (sut/gk-quadrature f 2.0 0.0 {:info? true})]
        (t/is (= (m/- (:result a)) (:result b)))
        (t/is (= (:error a) (:error b)) "the error estimate is not negated")
        (t/is (= (:evaluations a) (:evaluations b)))))
    (t/testing "infinite bounds in reverse order"
      (t/is (m/delta-eq -1.0 (sut/gk-quadrature f ##Inf 0.0 {:rel 1.0e-12}) 1.0e-11))
      (t/is (m/delta-eq -1.0 (sut/gk-quadrature (fn ^double [^double x] (m/exp x)) 0.0 ##-Inf {:rel 1.0e-12}) 1.0e-11))
      (t/is (m/delta-eq (m/- m/PI) (sut/gk-quadrature (integrands "cauchy") ##Inf ##-Inf {:rel 1.0e-12}) 1.0e-11)))
    (t/testing "empty interval"
      (t/is (= 0.0 (sut/gk-quadrature f 1.5 1.5)))
      (t/is (= 0.0 (:result (sut/gk-quadrature f 0.0 0.0 {:info? true})))))
    (t/testing "NaN bounds and two equal infinite bounds throw ex-info"
      (doseq [[lower upper] [[##NaN 1.0] [0.0 ##NaN] [##Inf ##Inf] [##-Inf ##-Inf]]]
        (t/is (some? (ex-data-of #(sut/gk-quadrature f lower upper))))))))

(t/deftest break-points
  (let [f (integrands "poly")]
    (t/testing "points split the interval and are sorted, deduplicated and clipped"
      (doseq [[points n] [[[0.5 1.0 2.0] 4]
                          [[2.0 0.5 1.0] 4]
                          [[0.5 0.5 1.0] 3]
                          [[-1.0 0.0 3.0 5.0] 1]
                          [[1.5] 2]
                          [[] 1]
                          [nil 1]]]
        (let [res (sut/gk-quadrature f 0.0 3.0 {:points points :info? true})]
          (t/is (= n (:subdivisions res)) (pr-str points))
          (t/is (m/delta-eq 9.0 (:result res) 1.0e-14)))))
    (t/testing "points and :initdiv combine"
      (t/is (= 6 (:subdivisions (sut/gk-quadrature f 0.0 3.0 {:points [1.5] :initdiv 3 :info? true})))))
    (t/testing "a point at the singularity gives the exact answer in one round"
      (let [res (sut/gk-quadrature (integrands "sqrtabs") -1.0 1.0 {:points [0.0] :rel 1.0e-12 :info? true})]
        (t/is (m/delta-eq (m// 4.0 3.0) (:result res) 1.0e-8))
        (t/is (< (:evaluations res) 945))))
    (t/testing "points of improper integrals are given in the original variable"
      (doseq [points [[1.0] [0.5 2.0 10.0] [30.0 1.0e-3]]]
        (t/is (m/delta-eq 1.0 (sut/gk-quadrature (integrands "expinf") 0.0 ##Inf {:points points :rel 1.0e-12}) 1.0e-11)))
      (t/is (m/delta-eq m/PI (sut/gk-quadrature (integrands "cauchy") ##-Inf ##Inf {:points [-5.0 0.0 3.0] :rel 1.0e-12}) 1.0e-11))
      (t/is (m/delta-eq 1.0 (sut/gk-quadrature (fn ^double [^double x] (m/exp x)) ##-Inf 0.0 {:points [-1.0 -10.0] :rel 1.0e-12}) 1.0e-11)))
    (t/testing "non numeric points throw ex-info"
      (t/is (= {:option :points :value [:a] :expected "sequence of numbers"}
               (ex-data-of #(sut/gk-quadrature f 0.0 1.0 {:points [:a]})))))))

(t/deftest option-validation
  (doseq [bad [0 -1 2.5 7.0 "7" :a nil ##NaN]]
    (t/testing (str ":integration-points " (pr-str bad))
      (let [data (ex-data-of #(sut/gk-quadrature (integrands "poly") 0.0 1.0 {:integration-points bad}))]
        (t/is (= :integration-points (:option data)))
        (t/is (= bad (:value data) )))))
  (doseq [bad [0 -2 1.5 "a"]]
    (t/testing (str ":initdiv " (pr-str bad))
      (t/is (= :initdiv (:option (ex-data-of #(sut/gk-quadrature (integrands "poly") 0.0 1.0 {:initdiv bad})))))))
  (t/testing "every positive integer order works"
    (doseq [n [1 2 3 4 5 9 12 30]]
      (t/is (m/delta-eq (m/- 1.0 (m// m/E)) (sut/gk-quadrature (integrands "expinf") 0.0 1.0 {:integration-points n :rel 1.0e-12}) 1.0e-11)))))
