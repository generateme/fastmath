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

;; sudoku solver

(defn sudoku-idx ^long [^long v ^long r ^long c] (+ (* 81 v) (* 9 r) c))

(defn sudoku-solver
  [input]
  (let [zeros (vec (repeat (* 9 9 9) 0))
        c1 (for [r (range 9)
                 c (range 9)]
             (reduce (fn [z v]
                       (assoc z (sudoku-idx v r c) 1)) zeros (range 9)))
        c2 (for [v (range 9)
                 c (range 9)]
             (reduce (fn [z r]
                       (assoc z (sudoku-idx v r c) 1)) zeros (range 9)))
        c3 (for [v (range 9)
                 r (range 9)]
             (reduce (fn [z c]
                       (assoc z (sudoku-idx v r c) 1)) zeros (range 9)))
        c4 (for [v (range 9)
                 p (range 3)
                 q (range 3)
                 :let [rc (for [r (range (* 3 p) (* 3 (inc p)))
                                c (range (* 3 q) (* 3 (inc q)))]
                            [r c])]]
             (reduce (fn [z [r c]]
                       (assoc z (sudoku-idx v r c) 1)) zeros rc))
        c5 (for [[v r c] input]
             (assoc zeros (sudoku-idx (dec v) r c) 1))
        all (reduce (fn [buff c]
                      (conj buff c :eq 1)) [] (mapcat identity [c1 c2 c3 c4 c5]))
        res (vec (first (sut/linear-optimization (conj zeros 0) all {:non-negative? true})))
        solution (vec (repeat 81 0))]
    (->> (for [r (range 9)
               c (range 9)
               v (range 9)
               :let [id (sudoku-idx v r c)]
               :when (m/one? (m/round (res id)))]
           [v r c])
         (reduce (fn [s [v r c]]
                   (assoc s (+ c (* r 9)) (inc v))) solution)
         (partition 9))))

(t/deftest sudoku
  ;; v - value 1-9
  ;; r - row 0-8
  ;; c - column 0-8
  ;;
  ;;                        v r c
  (t/is (= (sudoku-solver [[8 0 2]
                           [4 0 4]
                           [2 0 5]
                           [3 0 8]

                           [6 1 1]
                           [1 1 7]
                           
                           [3 2 1]
                           [7 2 5]

                           [5 3 5]
                           [3 3 6]

                           [4 4 1]
                           [9 4 3]
                           [1 4 5]
                           [7 4 7]

                           [5 5 2]
                           [2 5 3]

                           [3 6 3]
                           [4 6 7]

                           [5 7 1]
                           [6 7 7]

                           [7 8 0]
                           [6 8 3]
                           [1 8 4]
                           [9 8 6]])
           [[9 7 8 1 4 2 6 5 3]
            [2 6 4 5 9 3 7 1 8]
            [5 3 1 8 6 7 4 2 9]
            [1 2 7 4 8 5 3 9 6]
            [8 4 6 9 3 1 5 7 2]
            [3 9 5 2 7 6 1 8 4]
            [6 1 9 3 5 8 2 4 7]
            [4 5 3 7 2 9 8 6 1]
            [7 8 2 6 1 4 9 3 5]])))


