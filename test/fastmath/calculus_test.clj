(ns fastmath.calculus-test
  (:require [fastmath.calculus :as sut]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.random :as r]))

;; https://github.com/JohannesBuchner/cuba/blob/master/demo/demo-c.c
;; + wolfram alpha results

(def fun1-res 0.664669679781377268445214)
(defn fun1 [[^double x ^double y ^double z]] (m/* (m/sin x) (m/cos y) (m/exp z)))

(def fun1+-res 3.019450507398802024611853)
(defn fun1+ [[^double x ^double y ^double z]] (m/+ (m/sin x) (m/cos y) (m/exp z)))

(def fun2-res 5.26852)
(defn fun2 [[^double x ^double y ^double z]] (m/* (m// (m/+ (m/sq (m/+ x y)) 0.003))
                                               (m/cos y) (m/exp z)))

(def fun3-res 0.307807)
(defn fun3 [[^double x ^double y ^double z]] (m// (m/- 3.75 (m/cos (m/* m/PI x))
                                                    (m/cos (m/* m/PI y))
                                                    (m/cos (m/* m/PI z)))))

(def fun4-res 0.877314)
(defn fun4 [v] (m/abs (m/- (v/magsq v) 0.125)))

(def fun5-res 0.416538385886638169609660243601001780010941391101)
(defn fun5 [v] (m/exp (m/- (v/magsq v))))

(def fun6-res 1.202056902995100878845769442835673986997509412384)
(defn fun6 [[^double x ^double y ^double z]] (m// 1.0 (m/- 1.0 (m/* x y z) -1.0e-10)))

(def fun7-res 0.709615885636934)
(defn fun7 [[^double x ^double y ^double z]] (m/sqrt (m/abs (m/- x y z))))

(def fun8-res 0.891213)
(defn fun8 [[^double x ^double y ^double z]] (m/exp (m/- (m/* x y z))))

(def fun9-res 0.0801856)
(defn fun9 [[^double x ^double y ^double z]] (m// (m/sq x)
                                               (m/+ 5.0 (m/cos (m/+ 1.0 x y z)))))

#_(def fun10-res (m/+ 2.29149 0.104757))
#_(defn fun10 [[^double x ^double y ^double z]] (if (m/> x 0.5)
                                               (m// (m/+ (m/sqrt (m/* x y z)) 1.0e-5))
                                               (m/sqrt (m/* x y z))))

(def fun11-res 0.5235987755982988)
(defn fun11 [v] (if (m/< (v/magsq v) 1.0) 1.0 0.0))

(def normal-pdf (partial r/pdf (r/distribution :multi-normal {:means [0 0 0]
                                                            :covariances [[1 0 0]
                                                                          [0 0.5 0]
                                                                          [0 0 2]]})))

(defn gaussian [[^double _phi ^double r]] (m/* r (m/exp (m/* -0.5 r r))))

;; https://tutorial.math.lamar.edu/classes/calciii/TISphericalCoords.aspx
(defn polar-3d
  [[^double ro ^double _theta ^double phi]]
  (m/* 16.0 ro ro ro (m/sin phi) (m/cos phi)))

(t/deftest multivariate-integrals-vegas+
  (t/are [f res] (m/delta-eq res (:result (sut/vegas f [0 0 0] [1 1 1]
                                                     {:niters 20 :warmup 5 :res 1.0e-5 :abs 1.0e-5
                                                      :nstrats 20 :nevals 15000})) 1.0e-3)
    fun1+ fun1+-res
    fun1 fun1-res
    fun2 fun2-res
    fun3 fun3-res
    fun4 fun4-res
    fun5 fun5-res
    fun6 fun6-res
    fun7 fun7-res
    fun8 fun8-res
    fun9 fun9-res
    fun11 fun11-res)
  (t/are [lower upper res] (m/delta-eq res (:result (sut/vegas normal-pdf lower upper
                                                               {:niters 20 :warmup 5 :res 1.0e-4 :abs 1.0e-4
                                                                :nstrats 20 :nevals 10000})) 1.0e-3)
    [##-Inf ##-Inf ##-Inf] [##Inf ##Inf ##Inf] 1.0
    [##-Inf ##-Inf ##-Inf] [0 0 0] 0.125
    [0 0 0] [##Inf ##Inf ##Inf] 0.125
    [##-Inf 0 ##-Inf] [0 ##Inf ##Inf] 0.25)
  (t/is (m/delta-eq m/TWO_PI (:result (sut/vegas gaussian [0.0 0.0] [m/TWO_PI ##Inf]
                                                 {:max-iters 20 :warmup 5
                                                  :random-sequence :sobol})) 1.0e-3))
  (t/is (m/delta-eq (m/* 4.0 m/PI)
                    (:result (sut/vegas polar-3d [0 0 0] [1.0 m/TWO_PI m/HALF_PI]
                                        {:max-iters 20 :warmup 5
                                         :nstrats 20 :nevals 20000})) 1.0e-3)))

(t/deftest multivariate-integrals-cubature
  (t/are [f res] (m/delta-eq res (:result (sut/cubature f [0 0 0] [1 1 1]
                                                        {:max-iters 1000
                                                         :initdiv 18})) 1.0e-3)
    fun1+ fun1+-res
    fun1 fun1-res
    fun2 fun2-res
    fun3 fun3-res
    fun4 fun4-res
    fun5 fun5-res
    fun6 fun6-res
    fun7 fun7-res
    fun8 fun8-res
    fun9 fun9-res
    fun11 fun11-res)
  (t/are [lower upper res] (m/delta-eq res (:result (sut/cubature normal-pdf lower upper
                                                                  {:initdiv 15})) 1.0e-3)
    [##-Inf ##-Inf ##-Inf] [##Inf ##Inf ##Inf] 1.0
    [##-Inf ##-Inf ##-Inf] [0 0 0] 0.125
    [0 0 0] [##Inf ##Inf ##Inf] 0.125
    [##-Inf 0 ##-Inf] [0 ##Inf ##Inf] 0.25)
  (t/is (m/delta-eq m/TWO_PI (:result (sut/cubature gaussian [0.0 0.0] [m/TWO_PI ##Inf]))))
  (t/is (m/delta-eq (m/* 4.0 m/PI)
                    (:result (sut/cubature polar-3d [0 0 0] [1.0 m/TWO_PI m/HALF_PI] {:initdiv 15})))))


;;

(t/deftest univariate-gk
  (t/are [f a b res] (m/delta-eq res (:result (sut/integrate f a b {:integrator :gauss-kronrod})))
    (fn [^double x] (m/exp x)) -1 0 (m/- 1.0 (m// m/E))
    (fn [^double x] (m/fpow (m/cos x) 200)) 0 1 0.088512
    (fn [^double x] (m// 1.0 (m/inc (m/sq x)))) 0 1 m/QUARTER_PI
    (fn [^double x] (m/- (m// 2.0 (m/sqrt x))
                        (m/* 5.0 (m/sqrt x))
                        (m// (m/* x (m/sqrt x))))) 1 4 -61/3
    (fn [^double x] (->> (range 1 100)
                        (map (fn [^long v]
                               (m/cos (m// x v))))
                        (reduce m/fast* (m// (m/sin (m/* 4.0 x)) x)))) 0 ##Inf m/HALF_PI
    (fn [^double x] (m// (m/ln x) (m/inc (m/sq x)))) 0 ##Inf 0.0
    (fn [^double x] (m// (m/cosh x))) ##-Inf ##Inf m/PI
    (partial r/pdf (r/distribution :laplace)) ##-Inf ##Inf 1.0))

;;

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
                     (m/- (m// (m/+ (m/* d (m/sin x)) (m/* 2.0 x (m/cos x))) (m/sq d)))))))

(t/deftest hessian-tests
  (t/are [f point res] (every? identity (map (fn [h v] (v/delta-eq h v 1.0e-3))
                                             ((sut/hessian f {:h 1.0e-3}) point)
                                             res))
    (fn [[^double x ^double y]] (m/* (m/sin x) (m/cos y)))
    [m/HALF_PI m/PI] [[1 0] [0 1]]

    (fn [[^double x ^double y ^double z]] (m/+ (m/sq x) (m/sq y) (m/* x y) z (m/* x y z)))
    [1 1 1] [[2 2 1] [2 2 1] [1 1 0]]

    (fn [[^double x ^double y]] (m/+ (m/cb x) (m/* -2.0 x y) (m/* -1.0 (m/sq (m/cb y)))))
    [1 2] [[6 -2] [-2 -480]]))
