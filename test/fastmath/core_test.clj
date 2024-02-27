(ns fastmath.core-test
  (:require [fastmath.core :as m]
            [clojure.test :as t]
            [fastmath.vector :as v]))

(m/use-primitive-operators)

(t/deftest angles
  (t/is (= 180.0 (m/degrees m/PI)))
  (t/is (= m/PI (m/radians 180.0))))

(t/deftest frac
  (t/is (= 0.342225 (m/frac 3.342225)))
  (t/is (= 0.342225 (m/frac -3.342225)))
  (t/is (= 0.342225 (m/sfrac 3.342225)))
  (t/is (= -0.342225 (m/sfrac -3.342225))))

(t/deftest round-up-down
  (t/is (> (m/next-double 4.44) 4.44))
  (t/is (< (m/prev-double 4.44) 4.44))
  (t/is (= (m/round-up-pow2 1023) 1024))
  (t/is (= (m/round-up-pow2 1024) 1024))
  (t/is (= (m/round-up-pow2 1025) 2048)))

(t/deftest sgn
  (t/is (= -1.0 (m/signum -2)))
  (t/is (= 0.0 (m/signum 0)))
  (t/is (= 1.0 (m/signum 2)))
  (t/is (= -1.0 (m/sgn -2)))
  (t/is (= 1.0 (m/sgn 0)))
  (t/is (= 1.0 (m/sgn 2))))

(t/deftest norm
  (t/is (= (m/constrain -2 -1 1) -1))
  (t/is (= (m/constrain 0 -1 1) 0))
  (t/is (= (m/constrain 2 -1 1) 1))
  (t/is (= (m/norm 2 0 10) 0.2))
  (t/is (= (m/norm 2 0 10 0 100) 20.0)))

(t/deftest floating-points
  (t/is (= (m/double-exponent 2.0) 1))
  (t/is (= (m/double-exponent 0.5) -1))
  (t/is (= (m/double-significand 3.0) (Long/parseLong "1000000000000000000000000000000000000000000000000000" 2)))
  (t/is (= (m/double-significand 7.0) (Long/parseLong "1100000000000000000000000000000000000000000000000000" 2))))

;;


;; Abramowitz and Stegun p.511
(t/deftest kummers-m
  (t/are [a b x res] (m/delta-eq res (m/kummers-M a b x))
    0.3 0.2 -0.1 0.8578490
    -0.1 0.2 0.1 (* 0.8578490 (m/pow m/E 0.1))
    17 16 1 2.8881744
    -1 16 -1 1.0625
    -1.3 0.2 0.1 0.3582123
    -1.3 1.2 0.1 0.8924108
    -0.3 1.2 0.1 0.9745952
    1 1 1 m/E
    2 2 2 (m/exp 2)
    0.3 0.4 0.5 (/ 1.724128 (/ 0.7 0.6))
    ;; p.533
    -1 1 9 -8
    -1 0.6 9 -14
    0 1 9 1
    1 1 9 (m/exp 9))
  (t/are [a b x res acc] (m/delta-eq res (m/kummers-M a b x) acc)
    0.9 0.1 10 1227235 1
    -52.5 0.1 1 -16.34 0.2)
  (t/is (m/pos-inf? (m/kummers-M 1 0 1)))
  (t/is (m/neg-inf? (m/kummers-M -1 0 1)))
  (t/is (m/pos-inf? (m/kummers-M -1 0 -1)))
  (t/is (m/neg-inf? (m/kummers-M 1 0 -1)))
  (t/is (m/nan? (m/kummers-M -1 -1 2)))) ;; example 4

(t/deftest whittaker-m
  (t/is (m/delta-eq 1.10622 (m/whittaker-M 0 -0.4 1) 1.0e-5)))

;;

(defn lp5 [^double x] (/ (m/mevalpoly x 120.0 -600.0 600.0 -200.0 25.0 -1.0) 120.0))

(t/deftest laguerre
  (t/is (== 1.0 (m/laguerre-polynomials 0 1)))
  (t/is (== 1.0 (m/laguerre-polynomials 0 2 1)))
  (t/is (== -1.0 (m/laguerre-polynomials 1 2)))
  (t/is (== -0.5 (m/laguerre-polynomials 1 0.5 2)))
  (t/is (== (lp5 1) (m/laguerre-polynomials 5 1)))
  (t/is (== (lp5 5) (m/laguerre-polynomials 5 5)))
  (t/is (m/delta-eq (lp5 -1.234) (m/laguerre-polynomials 5 -1.234))))

(m/unuse-primitive-operators)

;;

(t/deftest besselk
  (t/are [order ress] (v/delta-eq ress (mapv (partial m/bessel-k-half order) [0.5 1 1.33 2.5 5]))
    1 [1.075047603 0.461068504 0.287423621 0.065065943 0.003776613]
    3 [3.225142810 0.922137009 0.503531608 0.091092320 0.004531936]
    5 [20.425904466  3.227479531  1.423209202  0.174376728  0.006495775]
    7 [207.48418748  17.05953466   5.85394214   0.43984578   0.01102771]
    9 [2.925204529e+03 1.226442222e+02 3.223343101e+01 1.405944900e+00 2.193457048e-02]))
