(ns fastmath.core-test
  (:require [fastmath.core :as m]
            [clojure.test :refer :all]))

(m/use-primitive-operators)

(deftest angles
  (is (= 180.0 (m/degrees m/PI)))
  (is (= m/PI (m/radians 180.0))))

(deftest frac
  (is (= 0.342225 (m/frac 3.342225)))
  (is (= 0.342225 (m/frac -3.342225)))
  (is (= 0.342225 (m/sfrac 3.342225)))
  (is (= -0.342225 (m/sfrac -3.342225))))

(deftest round-up-down
  (is (> (m/next-double 4.44) 4.44))
  (is (< (m/prev-double 4.44) 4.44))
  (is (= (m/round-up-pow2 1023) 1024))
  (is (= (m/round-up-pow2 1024) 1024))
  (is (= (m/round-up-pow2 1025) 2048)))

(deftest sgn
  (is (= -1.0 (m/signum -2)))
  (is (= 0.0 (m/signum 0)))
  (is (= 1.0 (m/signum 2)))
  (is (= -1.0 (m/sgn -2)))
  (is (= 1.0 (m/sgn 0)))
  (is (= 1.0 (m/sgn 2))))

(deftest norm
  (is (= (m/constrain -2 -1 1) -1))
  (is (= (m/constrain 0 -1 1) 0))
  (is (= (m/constrain 2 -1 1) 1))
  (is (= (m/norm 2 0 10) 0.2))
  (is (= (m/norm 2 0 10 0 100) 20.0)))

(deftest floating-points
  (is (= (m/double-exponent 2.0) 1))
  (is (= (m/double-exponent 0.5) -1))
  (is (= (m/double-significand 3.0) (Long/parseLong "1000000000000000000000000000000000000000000000000000" 2)))
  (is (= (m/double-significand 7.0) (Long/parseLong "1100000000000000000000000000000000000000000000000000" 2))))
