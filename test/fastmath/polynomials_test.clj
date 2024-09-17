(ns fastmath.polynomials-test
  (:require [fastmath.polynomials :as sut]
            [clojure.test :as t]
            [fastmath.core :as m]))

(defn lp5 [^double x] (/ (sut/mevalpoly x 120.0 -600.0 600.0 -200.0 25.0 -1.0) 120.0))

(t/deftest laguerre
  (t/is (== 1.0 (sut/eval-laguerre-L 0 1)))
  (t/is (== 1.0 (sut/eval-laguerre-L 0 2 1)))
  (t/is (== -1.0 (sut/eval-laguerre-L 1 2)))
  (t/is (== -0.5 (sut/eval-laguerre-L 1 0.5 2)))
  (t/is (== (lp5 1) (sut/eval-laguerre-L 5 1)))
  (t/is (== (lp5 5) (sut/eval-laguerre-L 5 5)))
  (t/is (m/delta-eq (lp5 -1.234) (sut/eval-laguerre-L 5 -1.234))))
