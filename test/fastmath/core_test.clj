(ns fastmath.core-test
  (:require [fastmath.core :as m]
            [clojure.test :as t]
            [fastmath.vector :as v]))

(m/use-primitive-operators)

;; primitive ops

(let [numbers (repeatedly 10000 #(* 100.0 (- (rand) 0.5)))]
  (t/deftest add
    (t/is (int? (+ 1 2 3)))
    (t/is (double? (+ 1 2 3.0)))
    (t/is (== 6 (+ 1 2 3)))
    (t/is (== 6.0 (+ 1 2 3.0)))
    (t/is (== 6 (m/long-add 1 2 3)))
    (t/is (== 6 (m/long-add 1 2 3.0)))
    (t/is (== (reduce clojure.core/+ numbers)
              (reduce + numbers)))
    (t/is (== (apply clojure.core/+ numbers)
              (apply + numbers)))
    (t/is (== (reduce clojure.core/+ (map long numbers))
              (reduce m/long-add numbers)))))

(let [numbers (repeatedly 10000 #(* 100.0 (- (rand) 0.5)))]
  (t/deftest sub
    (t/is (int? (- 1)))
    (t/is (double? (- 1.0)))
    (t/is (== -1.0 (- 1.0)))
    (t/is (== -1 (- 1)))
    (t/is (int? (- 1 2 3)))
    (t/is (double? (- 1 2 3.0)))
    (t/is (== -4 (- 1 2 3)))
    (t/is (== -4.0 (- 1 2 3.0)))
    (t/is (== -4 (m/long-sub 1 2 3)))
    (t/is (== -4 (m/long-sub 1 2 3.0)))
    (t/is (== (reduce clojure.core/- numbers)
              (reduce - numbers)))
    (t/is (== (apply clojure.core/- numbers)
              (apply - numbers)))
    (t/is (== (reduce clojure.core/- (map long numbers))
              (reduce m/long-sub numbers)))))



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


(m/unuse-primitive-operators)

;;

