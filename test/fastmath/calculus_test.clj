(ns fastmath.calculus-test
  "Public API of `fastmath.calculus`. The engines are tested in `fastmath.calculus.{quadrature,cubature,vegas,finite,divided,samples,common}-test`."
  (:require [fastmath.calculus :as sut]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.special :as special])
  (:import [org.apache.commons.math3.analysis.integration SimpsonIntegrator]))

(defn- ex-data-of [thunk]
  (try (thunk) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))

(t/deftest public-api
  (t/is (= '#{vegas cubature integrate integrate-samples extrapolate fx->gx+h fx->gx-h derivative f' f'' f'''
              gradient jacobian hessian hessian-diagonal laplacian divided divided-all}
           (set (keys (ns-publics 'fastmath.calculus)))))
  (t/testing "every public var is documented"
    (doseq [[sym var] (ns-publics 'fastmath.calculus)]
      (t/is (string? (:doc (meta var))) (str sym))
      (t/is (seq (:arglists (meta var))) (str sym)))))

;; univariate integration

(defn- exp-minus ^double [^double x] (m/exp (m/- x)))
(defn- cauchy ^double [^double x] (m// (m/inc (m/sq x))))

(t/deftest univariate-gk
  (t/are [f a b res] (m/delta-eq res (sut/integrate f a b {:integrator :gauss-kronrod}))
    (fn [^double x] (m/exp x)) -1 0 (m/- 1.0 (m// m/E))
    (fn [^double x] (m/fpow (m/cos x) 200)) 0 1 0.088512
    (fn [^double x] (m// 1.0 (m/inc (m/sq x)))) 0 1 m/QUARTER_PI
    (fn [^double x] (m/- (m// 2.0 (m/sqrt x))
                        (m/* 5.0 (m/sqrt x))
                        (m// (m/* x (m/sqrt x))))) 1 4 (m// -61.0 3.0)
    (fn [^double x] (->> (range 1 100)
                        (map (fn [^long v]
                               (m/cos (m// x v))))
                        (reduce m/* (m// (m/sin (m/* 4.0 x)) x)))) 0 ##Inf m/HALF_PI
    (fn [^double x] (m// (m/ln x) (m/inc (m/sq x)))) 0 ##Inf 0.0
    (fn [^double x] (m// (m/cosh x))) ##-Inf ##Inf m/PI
    (partial r/pdf (r/distribution :laplace)) ##-Inf ##Inf 1.0))

(t/deftest integrators
  (t/testing "every integrator on proper, improper and reversed intervals"
    (doseq [integrator [:gauss-kronrod :gauss-legendre :romberg :simpson :midpoint :trapezoid]
            :let [tol (case integrator (:midpoint :trapezoid) 1.0e-5 :simpson 1.0e-7 1.0e-8)
                  opts {:integrator integrator :rel (case integrator :midpoint 1.0e-5 :trapezoid 1.0e-6 1.0e-9) :abs 1.0e-12 :max-evals 5000000}
                  improper-tol (case integrator (:midpoint :trapezoid) 1.0e-3 1.0e-4)]]
      (t/testing (name integrator)
        (t/is (m/delta-eq (m/- 1.0 (m// m/E)) (sut/integrate exp-minus 0.0 1.0 opts) tol))
        (t/is (m/delta-eq (m/- (m/- 1.0 (m// m/E))) (sut/integrate exp-minus 1.0 0.0 opts) tol) "reversed")
        (t/is (m/delta-eq 1.0 (sut/integrate exp-minus 0.0 ##Inf opts) improper-tol) "half line")
        (t/is (m/delta-eq -1.0 (sut/integrate exp-minus ##Inf 0.0 opts) improper-tol) "reversed half line")
        (t/is (m/delta-eq m/PI (sut/integrate cauchy ##-Inf ##Inf opts) improper-tol) "whole line")
        (t/is (m/delta-eq (m/- m/PI) (sut/integrate cauchy ##Inf ##-Inf opts) improper-tol) "reversed whole line")
        (t/is (= 0.0 (sut/integrate exp-minus 0.5 0.5 opts)) "empty interval"))))
  (t/testing "default bounds and arities"
    (t/is (m/delta-eq (m/- 1.0 (m// m/E)) (sut/integrate exp-minus)))
    (t/is (m/delta-eq (m/- 1.0 (m// m/E)) (sut/integrate exp-minus 0 1))))
  (t/testing "an integrator instance"
    (let [res (sut/integrate exp-minus 0.0 1.0 {:integrator (SimpsonIntegrator.) :info? true})]
      (t/is (m/delta-eq (m/- 1.0 (m// m/E)) (:result res) 1.0e-6))
      (t/is (false? (:fail? res))))))

(t/deftest result-maps
  (t/testing "Gauss-Kronrod"
    (let [res (sut/integrate exp-minus 0.0 1.0 {:info? true})]
      (t/is (= #{:result :error :evaluations :iterations :subdivisions :fail?} (set (keys res))))
      (t/is (false? (:fail? res)))))
  (t/testing "the other integrators"
    (doseq [integrator [:gauss-legendre :romberg :simpson :midpoint :trapezoid]]
      (let [res (sut/integrate exp-minus 0.0 1.0 {:integrator integrator :info? true})]
        (t/is (= #{:result :iterations :evaluations :integrator :fail?} (set (keys res))) (name integrator))
        (t/is (= integrator (:integrator res)))
        (t/is (false? (:fail? res)))))))

(t/deftest failure-contract
  (let [slow (fn ^double [^double x] (m/abs (m/sin (m/* 50.3 x))))]
    (t/testing "Gauss-Kronrod returns the best estimate and the reason"
      (let [res (sut/integrate slow 0.0 1.0 {:rel 1.0e-14 :abs 0.0 :max-iters 3 :info? true})]
        (t/is (= :max-iters (:fail? res)))
        (t/is (m/delta-eq (m// 2.0 m/PI) (:result res) 0.2)))
      (let [res (sut/integrate slow 0.0 1.0 {:rel 1.0e-14 :abs 0.0 :max-evals 100 :info? true})]
        (t/is (= :max-evals (:fail? res)))
        (t/is (<= (:evaluations res) 100))))
    (t/testing "the other integrators return NaN and the reason"
      (doseq [integrator [:romberg :simpson :midpoint :trapezoid :gauss-legendre]]
        (let [res (sut/integrate slow 0.0 1.0 {:integrator integrator :rel 1.0e-14 :abs 0.0 :max-evals 50 :info? true})]
          (t/is (= :max-evals (:fail? res)) (name integrator))
          (t/is (m/nan? (:result res)) (name integrator)))
        (let [res (sut/integrate slow 0.0 1.0 {:integrator integrator :rel 1.0e-14 :abs 0.0 :max-iters 5 :info? true})]
          (t/is (= :max-iters (:fail? res)) (name integrator))
          (t/is (m/nan? (:result res)) (name integrator)))))
    (t/testing "a slowly decaying improper integral is stopped by the default evaluation limit"
      (let [res (sut/integrate (fn ^double [^double x] (m/pow x -1.5)) 1.0 ##Inf {:integrator :romberg :info? true})]
        (t/is (= :max-evals (:fail? res)))
        (t/is (<= (:evaluations res) 10000000))))))

(t/deftest peak-guard
  (let [peak (fn ^double [^double x] (m/exp (m/* -1.0e6 x x)))
        exact (m// m/SQRTPI 1000.0)]
    (t/testing "the first quadrature points miss a narrow peak, which the default absolute error accepts"
      (t/is (m/< (sut/integrate peak -1.0 1.0) (m/* 1.0e-6 exact))))
    (t/testing ":peak-guard sets the default absolute error to 0"
      (t/is (m/delta-eq exact (sut/integrate peak -1.0 1.0 {:peak-guard true}) (m/* 1.0e-6 exact)))
      (t/is (m/delta-eq exact (sut/integrate peak -1.0 1.0 {:peak-guard true :rel 1.0e-10}) (m/* 1.0e-8 exact))))
    (t/testing "an explicit :abs wins"
      (t/is (m/< (sut/integrate peak -1.0 1.0 {:peak-guard true :abs 1.0e-3}) (m/* 1.0e-6 exact))))
    (t/testing "ignored by the other integrators"
      (t/is (m/delta-eq (m/- 1.0 (m// m/E)) (sut/integrate exp-minus 0.0 1.0 {:integrator :romberg :peak-guard true}) 1.0e-8)))))

(t/deftest break-points-and-options
  (t/testing ":points and :initdiv through integrate"
    (let [f (fn ^double [^double x] (m/sqrt (m/abs x)))]
      (t/is (m/delta-eq (m// 4.0 3.0) (sut/integrate f -1.0 1.0 {:points [0.0] :rel 1.0e-12}) 1.0e-8))
      (t/is (= 3 (:subdivisions (sut/integrate f -1.0 1.0 {:points [-0.5 0.5] :max-iters 1 :info? true}))))
      (t/is (= 8 (:subdivisions (sut/integrate f -1.0 1.0 {:initdiv 8 :max-iters 1 :info? true}))))
      (t/is (m/delta-eq 1.0 (sut/integrate exp-minus 0.0 ##Inf {:points [1.0 5.0] :rel 1.0e-12}) 1.0e-10))))
  (t/testing ":integration-points is the order of the Gauss-Kronrod rule"
    (t/is (= 15 (:evaluations (sut/integrate exp-minus 0.0 1.0 {:info? true}))))
    (let [res (sut/integrate exp-minus 0.0 1.0 {:info? true :integration-points 1})]
      (t/is (zero? (mod (:evaluations res) 3)) "the 3 point rule")
      (t/is (m/delta-eq (m/- 1.0 (m// m/E)) (:result res) 1.0e-8)))
    (t/is (= 43 (:evaluations (sut/integrate exp-minus 0.0 1.0 {:info? true :integration-points 21}))))))

(t/deftest integrate-validation
  (t/testing "ex-info with the option and the value"
    (doseq [[options option value] [[{:integration-points 0} :integration-points 0]
                                    [{:integration-points 2.5} :integration-points 2.5]
                                    [{:integrator :gauss-legendre :integration-points 2.5} :integration-points 2.5]
                                    [{:initdiv 0} :initdiv 0]
                                    [{:points [:a]} :points [:a]]
                                    [{:integrator :foo} :integrator :foo]]]
      (let [data (ex-data-of #(sut/integrate exp-minus 0.0 1.0 options))]
        (t/is (= option (:option data)) (pr-str options))
        (t/is (= value (:value data)) (pr-str options)))))
  (t/testing "NaN bounds and two equal infinite bounds"
    (doseq [[lower upper] [[##NaN 1.0] [0.0 ##NaN] [##Inf ##Inf] [##-Inf ##-Inf]]
            integrator [:gauss-kronrod :romberg]]
      (t/is (some? (ex-data-of #(sut/integrate exp-minus lower upper {:integrator integrator})))))))

;; sampled data (details in samples-test)

(t/deftest integrate-samples
  (t/is (m/delta-eq 8.0 (sut/integrate-samples [0 1 2 3 4] 1.0 {:method :simpson})) "linear data")
  (t/is (m/delta-eq 16.0 (sut/integrate-samples [0 2 4 6 8] 1.0)) "arity with a spacing")
  (t/is (m/delta-eq 4.0 (sut/integrate-samples [1 2 3] [0 1 2] {:method :trapezoid})) "abscissae")
  (t/is (= [0.0 1.5 4.0] (sut/integrate-samples [1 2 3] 1.0 {:cumulative? true})) "cumulative")
  (t/is (= 0.0 (sut/integrate-samples [])))
  (t/is (m/delta-eq 4.0 (sut/integrate-samples [1.0 1.0 1.0 1.0 1.0]))))

;; multivariate integration (details in cubature-test and vegas-test)

(t/deftest multivariate
  (let [f (fn [[^double x ^double y]] (m/* (m/sin x) (m/cos y)))
        exact (m/* (m/- 1.0 (m/cos 1.0)) (m/sin 2.0))]
    (t/is (m/delta-eq exact (sut/cubature f [0 0] [1 2]) 1.0e-6))
    (t/is (m/delta-eq (m/- exact) (sut/cubature f [1 0] [0 2]) 1.0e-6) "one reversed dimension")
    (t/is (map? (sut/cubature f [0 0] [1 2] {:info? true})))
    (t/is (m/delta-eq exact (sut/vegas f [0 0] [1 2] {:rng (r/rng :jdk 5) :nevals 4000 :max-iters 5 :rel 1.0e-4}) 2.0e-3))
    (t/is (m/delta-eq (m/- exact) (sut/vegas f [1 0] [0 2] {:rng (r/rng :jdk 5) :nevals 4000 :max-iters 5 :rel 1.0e-4}) 2.0e-3))
    (t/is (m/delta-eq 1.0 (sut/cubature (fn [[^double x]] (m/exp (m/- x))) [0] [##Inf]) 1.0e-6) "one dimension")))

;; https://docs.scipy.org/doc/scipy/tutorial/integrate.html
(t/deftest scipy-example
  (t/is (m/delta-eq (sut/integrate (fn [^double x] (special/bessel-J 2.5 x)) 0.0 4.5)
                    1.117817938088701 1.0e-11))
  (t/is (v/delta-eq (map (fn [^double x]
                           (sut/integrate (fn [^double t] (m// (m/exp (m/* -1.0 x t))
                                                              (m/pow t 3.0))) 1.0 ##Inf))
                         (range 1.0 4.0 0.5))
                    [0.1097,  0.0567,  0.0301,  0.0163,  0.0089,  0.0049] 1.0e-4)))

;; derivatives (details in finite-test)

(defn xsinx ^double [^double x] (m/* x (m/sin x)))
(defn sinx+xcosx ^double [^double x] (m/+ (m/sin x) (m/* x (m/cos x))))
(defn twocosx-xsinx ^double [^double x] (m/- (m/* 2.0 (m/cos x))
                                          (xsinx x)))
(defn -threesinx-xcosx ^double [^double x] (m/- (m/* -3.0 (m/sin x))
                                             (m/* x (m/cos x))))
(defn xsinx-4cosx ^double [^double x] (m/- (xsinx x)
                                        (m/* 4.0 (m/cos x))))

(t/deftest derivative-tests
  (t/are [f f'] (->> (repeatedly 1000 r/grand)
                     (map (juxt (sut/derivative f 1 {:acc 4}) f'))
                     (every? (fn [[^double a ^double b]]
                               (m/delta-eq a b))))
    m/sq (fn [^double x] (m/* 2.0 x))
    (fn [^double x] (m/atan x)) (fn [^double x] (m// (m/inc (m/sq x))))

    (fn [^double x] (m/exp (m/sin (m/sq x))))
    (fn [^double x] (m/* 2.0 x (m/exp (m/sin (m/sq x))) (m/cos (m/sq x))))

    (fn [^double x] (m// (m/cos x) (m/inc (m/* x x))))
    (fn [^double x] (let [d (m/inc (m/* x x))]
                     (m/- (m// (m/+ (m/* d (m/sin x)) (m/* 2.0 x (m/cos x))) (m/sq d))))))
  (t/are [d n] (->> (repeatedly 1000 r/grand)
                    (map (juxt (sut/derivative xsinx n {:acc 6}) d))
                    (every? (fn [[^double a ^double b]]
                              (m/delta-eq a b 1.0e-3))))
    xsinx 0
    sinx+xcosx 1
    twocosx-xsinx 2
    -threesinx-xcosx 3
    xsinx-4cosx 4)
  (t/testing "shortcuts"
    (t/is (m/delta-eq (m/cos 0.7) ((sut/f' m/sin) 0.7) 1.0e-8))
    (t/is (m/delta-eq (m/- (m/sin 0.7)) ((sut/f'' m/sin) 0.7) 1.0e-6))
    (t/is (m/delta-eq (m/- (m/cos 0.7)) ((sut/f''' m/sin) 0.7) 1.0e-4))))

(t/deftest gradient-hessian-tests
  (t/are [f point res] (every? identity (map (fn [h v] (v/delta-eq h v 1.0e-3))
                                             ((sut/hessian f {:h 1.0e-3}) point)
                                             res))
    (fn [[^double x ^double y]] (m/* (m/sin x) (m/cos y)))
    [m/HALF_PI m/PI] [[1 0] [0 1]]

    (fn [[^double x ^double y ^double z]] (m/+ (m/sq x) (m/sq y) (m/* x y) z (m/* x y z)))
    [1 1 1] [[2 2 1] [2 2 1] [1 1 0]]

    (fn [[^double x ^double y]] (m/+ (m/cb x) (m/* -2.0 x y) (m/* -1.0 (m/sq (m/cb y)))))
    [1 2] [[6 -2] [-2 -480]])
  (t/testing "the other multivariate derivatives are exported with vectors as results"
    (let [f (fn [[^double x ^double y]] (m/+ (m/sq x) (m/* 3.0 y y)))]
      (t/is (v/delta-eq [2.0 12.0] ((sut/gradient f) [1.0 2.0]) 1.0e-8))
      (t/is (vector? ((sut/gradient f) [1.0 2.0])))
      (t/is (v/delta-eq [2.0 6.0] ((sut/hessian-diagonal f) [1.0 2.0]) 1.0e-6))
      (t/is (vector? ((sut/hessian-diagonal f) [1.0 2.0])))
      (t/is (m/delta-eq 8.0 ((sut/laplacian f) [1.0 2.0]) 1.0e-6))
      (t/is (every? true? (map #(v/delta-eq %1 %2 1.0e-6)
                               [[3.0 2.0] [1.0 1.0]]
                               ((sut/jacobian (fn [[^double x ^double y]] [(m/* x y) (m/+ x y)])) [2.0 3.0])))))))

;; extrapolation

(t/deftest richardson-extrapolation
  (t/is (m/delta-eq ((sut/extrapolate (sut/fx->gx+h (fn [^double x] (m// (m/sin x) x)))
                                      {:abs 0 :rel 1.0e-10
                                       :power 2 :init-h 1.0}) 0.0) 1.0 1.0e-15))
  (t/is (m/delta-eq ((sut/extrapolate (sut/fx->gx-h (fn [^double x] (m// (m/sin x) x)))
                                      {:abs 0 :rel 1.0e-10
                                       :power 2 :init-h 1.0}) 0.0) 1.0 1.0e-15))
  (t/are [d n] (->> (repeatedly 1000 r/grand)
                    (map (juxt (sut/derivative xsinx n {:acc 4 :extrapolate? true}) d))
                    (every? (fn [[^double a ^double b]]
                              (m/delta-eq a b 1.0e-5))))
    xsinx 0
    sinx+xcosx 1
    twocosx-xsinx 2
    -threesinx-xcosx 3
    xsinx-4cosx 4))

;; divided differences

(t/deftest divided-differences
  (t/testing "Simple case, logmean as 1/f[dt1,dt2]"
    (let [dt1 30 dt2 20]
      (t/is (m/delta-eq 24.663 (m// (sut/divided m/log [dt1 dt2])) 1.0e-4))))
  (t/testing "wikipedia example (Newton polynomial article)"
    (let [dd (vec (sut/divided-all m/tan [-3/2 -3/4 0 3/4 3/2]))]
      (t/are [id res] (v/delta-eq res (dd id) 1.0e-4)
        0 [0.0]
        1 [4.83484 4.83484]
        2 [-10.8784 0.0 10.8784]
        3 [17.5597 1.24213 1.24213 17.5597]
        4 [-14.1014 -0.931596 0.0 0.931596 14.1014]))
    (t/is (m/zero? (sut/divided m/tan [-3/2 -3/4 0 3/4 3/2])))
    (let [dd (vec (sut/divided-all [6 9 2 5] [0 1 2 3]))]
      (t/are [id res] (v/delta-eq res (dd id) 1.0e-10)
        0 [10/3]
        1 [-5 5]
        2 [3 -7 3]
        3 [6 9 2 5])))
  (t/testing "repeated nodes"
    (t/is (m/delta-eq (m/- (m/exp 1.0) 2.0) (sut/divided m/exp [0 0 1] {:derivatives [m/exp]}) 1.0e-12))
    (t/is (m/delta-eq 0.5 (sut/divided m/exp [0 0 0] {:derivatives [m/exp m/exp]}) 1.0e-12))))
