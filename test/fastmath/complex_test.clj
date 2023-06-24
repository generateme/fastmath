(ns fastmath.complex-test
  (:refer-clojure :exclude [abs])
  (:require [fastmath.complex :refer :all]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.test :refer :all]))

;; expected reults verified with Wolfram Alpha

(deftest complex-operations
  (is (= (v/vec2 1.0 1.0) (complex 1.0 1.0)))
  (is (= (complex 1.0 1.0) (add I ONE)))
  (is (= (m/sqrt 2.0) (abs (complex 1.0 1.0))))
  (is (= (complex 1.0 1.0) (sub ONE I-)))
  (is (= -90.0 (m/degrees (arg I-))))
  (is (= (complex 0.0 -1.0) (conjugate I)))
  (is (= (complex 0.44 0.08) (div (complex 1 2) (complex 3 4))))
  (is (= (complex 0.12 -0.16) (reciprocal (complex 3 4))))
  (is (= (complex -5.0 10.0) (mult (complex 1 2) (complex 3 4))))
  (is (= (complex -1.0 2.0) (neg (complex 1.0 -2.0))))
  (is (= (complex -3.0 4.0) (sq (complex 1.0 2.0))))
  (is (= (complex 25.0 0.0) (sq (complex 5.0 0.0))))
  (is (= (complex (m/sqrt 2.0) 0.0) (sqrt (complex 2.0 0.0))))
  (is (= (complex 1.0 2.0) (v/approx (sqrt (complex -3.0 4.0)))))
  (is (= (complex 4.08 -2.94) (v/approx (sqrt1z (complex 3 4))))))

(deftest complex-trig
  (is (= (complex 2.03 -3.05) (v/approx (cos (complex 1 2)))))
  (is (= (complex 3.17 1.96) (v/approx (sin (complex 1 2)))))
  (is (= (complex -0.64 1.07) (v/approx (cosh (complex 1 2)))))
  (is (= (complex -0.49 1.4) (v/approx (sinh (complex 1 2)))))
  (is (= (complex 0.03 1.01) (v/approx (tan (complex 1 2)))))
  (is (= (complex 1.17 -0.24) (v/approx (tanh (complex 1 2)))))
  (is (= (complex 0.15 0.23) (v/approx (sec (complex 1 2)))))
  (is (= (complex 0.23 -0.14) (v/approx (csc (complex 1 2)))))
  (is (= (complex 1.14 -1.53) (v/approx (acos (complex 1 2)))))
  (is (= (complex 0.43 1.53) (v/approx (asin (complex 1 2)))))
  (is (= (complex 1.34 0.4) (v/approx (atan (complex 1 2))))))

(deftest log-exp
  (is (= ZERO (v/approx (add (exp (complex 0 m/PI)) ONE))))
  (is (= (complex -1.13 2.47) (v/approx (exp (complex 1 2)))))
  (is (= ONE (log (complex m/E 0.0))))
  (is (= (complex 0.8 1.11) (v/approx (log (complex 1 2)))))
  (is (= (complex -0.22 0.1) (v/approx (pow (complex 1 2) (complex 1 2))))))
