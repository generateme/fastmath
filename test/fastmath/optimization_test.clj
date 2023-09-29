(ns fastmath.optimization-test
  (:require [fastmath.optimization :as sut]
            [clojure.test :as t]
            [fastmath.vector :as v]
            [fastmath.core :as m]))

(t/deftest linear-programming
  (t/are [res target constrains opts] (let [[p v] res
                                            [rp rv] (sut/linear-optimization target constrains opts)]
                                        (and (v/delta-eq (vec rp) p)
                                             (m/delta-eq v rv)))
    ;; scipy
    [[10 -3] -22] [-1 4 0] [[-3 1] :<= 6
                            [-1 -2] :>= -4
                            [0 1] :>= -3] nil
    
    ;; http://people.brunel.ac.uk/~mastjjb/jeb/or/morelp.html
    [[45.0 6.25] 1.25] [1 1 -50] [[50 24] :leq 2400
                                  [30 33] :leq 2100
                                  [1 0] :geq 45
                                  [0 1] :geq 5] {:goal :maximize}
    
    ;; transportation problem https://dewwool.com/linear-programming-examples/
    [[0.0 0.0 10.0 30.0 20.0 20.0] 290.0] [5 4 6 3 2 5 0] '[[1 1 1 0 0 0] <= 50
                                                            [0 0 0 1 1 1] <= 70
                                                            [1 0 0 1 0 0] >= 30
                                                            [0 1 0 0 1 0] >= 20
                                                            [0 0 1 0 0 1] >= 30] {:non-negative? true}))
