(ns fastmath.interpolation.barycentric-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.interpolation.barycentric :as sut]))

(set! *warn-on-reflection* true)

;; Oracles (utils/fastmath/dev/generate_interp1d_reference.py):
;; - :polynomial   mpmath Lagrange polynomial (60 digits) and scipy BarycentricInterpolator
;; - :barycentric  mpmath Floater-Hormann rational interpolant (50 digits), weights from the paper formula

(def ^:private ref-data
  (delay (edn/read-string (slurp (io/resource "interpolation/interp1d_reference.edn")))))

(defn- property
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (= :invalid-input (:type d)) (:property d))))))

(defn- ymax ^double [ys] (double (reduce max 0.0 (map #(Math/abs (double %)) ys))))

(defn- permute
  "Deterministic permutation of the indices: reversed and rotated by two (not the identity for n > 2)."
  [^long n]
  (map (fn [^long i] (mod (- (dec n) (+ i 2)) n)) (range n)))

(t/deftest polynomial-reference
  ;; measured vs mpmath, inside the domain: <= 1.3e-14 * max(ys) (equispaced 15 nodes, Neville 8.9e-16; scipy
  ;; BarycentricInterpolator differs from mpmath by 2.2e-14);
  ;; outside the domain: <= 6.0e-9 relative (the second barycentric form cancels when it extrapolates)
  (doseq [d (:polynomial @ref-data)]
    (let [f (sut/polynomial (:xs d) (:ys d))
          lo (apply min (:xs d))
          hi (apply max (:xs d))
          sc (ymax (:ys d))]
      (doseq [[q e s] (map vector (:qs d) (:expected d) (:scipy d))]
        (if (<= lo q hi)
          (do (t/is (m/delta-eq e (f q) (* 1e-13 sc) 1e-13) (str (:name d) " at " q " inside"))
              (t/is (m/delta-eq s (f q) (* 1e-12 sc) 1e-12) (str (:name d) " vs scipy at " q)))
          (t/is (m/delta-eq e (f q) 0.0 5e-8) (str (:name d) " at " q " outside")))))))

(t/deftest polynomial-edge-cases
  (t/testing "two points: the straight line"
    (let [f (sut/polynomial [0 2] [1 5])]
      (t/is (m/delta-eq 3.0 (f 1) 1e-15))
      (t/is (m/delta-eq 9.0 (f 4) 1e-14))))
  (t/testing "reproduces a polynomial of degree n-1 exactly (x^3 on 4 points), unsorted nodes"
    (let [f (sut/polynomial [2.0 -1.0 0.5 3.0] [8.0 -1.0 0.125 27.0])]
      (doseq [x [-2.0 0.0 0.3 1.7 2.9 5.0]
              :let [x (double x)]]
        (t/is (m/delta-eq (* x x x) (f x) 1e-12)))))
  (t/testing "equals barycentric with order n-1"
    (let [xs [0.0 0.3 1.1 1.5 2.7]
          ys [1.0 -2.0 0.5 3.0 -1.0]
          f (sut/polynomial xs ys)
          g (sut/barycentric xs ys {:order 4})]
      (doseq [x [-1.0 0.1 0.9 2.0 3.5]]
        (t/is (== (f x) (g x))))))
  (t/testing "exact node hit even with unsorted nodes returns the value of the node"
    (let [xs [3.0 1.0 2.0 0.0]
          ys [30.0 10.0 20.0 0.0]
          f (sut/polynomial xs ys)]
      (t/is (= ys (mapv f xs)))))
  (t/testing "many nodes on a tiny interval: the product of node distances underflows without scaling"
    ;; cheb100-tiny: 100 Chebyshev nodes on [0, 1e-4], the product of the distances is about 1e-450
    (let [d (last (:polynomial @ref-data))
          f (sut/polynomial (:xs d) (:ys d))]
      (t/is (= "cheb100-tiny" (:name d)))
      (doseq [[q e] (map vector (:qs d) (:expected d))]
        (t/is (m/delta-eq e (f q) 1e-13) (str "at " q)))))
  (t/testing "NaN query gives NaN"
    (t/is (m/nan? ((sut/polynomial [0 1 2] [0 1 4]) ##NaN)))))

(t/deftest polynomial-validation
  (t/is (= :too-few-points (property #(sut/polynomial [] []))))
  (t/is (= :too-few-points (property #(sut/polynomial [1] [1]))))
  (t/is (= :length-mismatch (property #(sut/polynomial [0 1 2] [0 1]))))
  (t/is (= :not-finite (property #(sut/polynomial [0 ##NaN 2] [0 1 2]))))
  (t/is (= :not-finite (property #(sut/polynomial [0 ##Inf 2] [0 1 2]))))
  (t/is (= :duplicate-nodes (property #(sut/polynomial [0 1 1 2] [0 1 2 3]))))
  (t/is (= :duplicate-nodes (property #(sut/polynomial [2 0 1 0] [0 1 2 3]))))
  (t/is (nil? (property #(sut/polynomial [2 0 1] [0 1 2])))))

(t/deftest barycentric-reference
  ;; measured vs mpmath Floater-Hormann: inside <= 4.6e-15 * max(ys), outside <= 2.9e-12 relative
  (doseq [d (:barycentric @ref-data)]
    (let [f (sut/barycentric (:xs d) (:ys d) {:order (:order d)})
          lo (apply min (:xs d))
          hi (apply max (:xs d))
          sc (ymax (:ys d))]
      (doseq [[q e] (map vector (:qs d) (:expected d))]
        (if (<= lo q hi)
          (t/is (m/delta-eq e (f q) (* 1e-12 sc) 1e-12) (str (:name d) " order " (:order d) " at " q))
          (t/is (m/delta-eq e (f q) 0.0 1e-10) (str (:name d) " order " (:order d) " at " q " outside")))))))

(t/deftest barycentric-unsorted
  ;; the exact node hit was a binary search: unsorted xs returned wrong values (0.07 and 0.005 in the audit)
  (doseq [d (:barycentric @ref-data)
          :when (> (count (:xs d)) 2)]
    (let [xs (:xs d)
          ys (:ys d)
          n (count xs)
          idx (permute n)
          pxs (mapv xs idx)
          pys (mapv ys idx)
          sorted (sut/barycentric xs ys {:order (:order d)})
          unsorted (sut/barycentric pxs pys {:order (:order d)})]
      (doseq [q (concat (:qs d) xs)]
        (t/is (m/delta-eq (sorted q) (unsorted q) 1e-12) (str (:name d) " order " (:order d) " at " q)))
      (t/testing "exact hit of every node"
        (t/is (= (vec pys) (mapv unsorted pxs)))))))

(t/deftest barycentric-options
  (let [xs [0.0 0.5 1.0 2.0 3.5]
        ys [1.0 3.0 2.0 5.0 4.0]]
    (t/testing "default order is 1"
      (let [f (sut/barycentric xs ys)
            g (sut/barycentric xs ys {:order 1})
            h (sut/barycentric xs ys {:order nil})]
        (doseq [x [0.2 1.7 2.9 -1.0 4.0]]
          (t/is (== (f x) (g x) (h x))))))
    (t/testing "a double with an integer value is accepted"
      (t/is (== ((sut/barycentric xs ys {:order 2}) 1.7) ((sut/barycentric xs ys {:order 2.0}) 1.7))))
    (t/testing "order 0 is Berrut: a rational function with alternating weights, through the data"
      (let [f (sut/barycentric xs ys {:order 0})]
        (t/is (= ys (mapv f xs)))
        ;; xs 0 1: weights +1 -1 -> r(0.5) = (y0/(0.5) - y1/(-0.5)) / (1/0.5 - 1/(-0.5)) = (y0 + y1) / 2
        (t/is (m/delta-eq 2.0 ((sut/barycentric [0.0 1.0] [1.0 3.0] {:order 0}) 0.5) 1e-15))))
    (t/testing "order n-1 equals the polynomial"
      (t/is (m/delta-eq ((sut/polynomial xs ys) 1.3) ((sut/barycentric xs ys {:order 4}) 1.3) 1e-15)))))

(t/deftest barycentric-validation
  (let [xs [0.0 1.0 2.0 3.0]
        ys [1.0 2.0 0.0 1.0]]
    (t/testing "order range (the assert was dropped with *assert* false and a float order was truncated)"
      (t/is (nil? (property #(sut/barycentric xs ys {:order 0}))))
      (t/is (nil? (property #(sut/barycentric xs ys {:order 3}))))
      (t/is (= :invalid-order (property #(sut/barycentric xs ys {:order 4}))))
      (t/is (= :invalid-order (property #(sut/barycentric xs ys {:order -1}))))
      (t/is (= :invalid-order (property #(sut/barycentric xs ys {:order 1.5}))))
      (t/is (= :invalid-order (property #(sut/barycentric xs ys {:order 2.0000001}))))
      (t/is (= :invalid-order (property #(sut/barycentric xs ys {:order ##NaN}))))
      (t/is (= :invalid-order (property #(sut/barycentric xs ys {:order ##Inf}))))
      (t/is (= :invalid-order (property #(sut/barycentric xs ys {:order "2"}))))
      (t/is (= :invalid-order (property #(sut/barycentric xs ys {:order :a}))))
      (t/is (= :invalid-order (property #(sut/barycentric xs ys {:order 100000000000})))))
    (t/testing "the order is checked even when asserts are disabled"
      (binding [*assert* false]
        (t/is (= :invalid-order (property #(sut/barycentric xs ys {:order 9}))))))
    (t/testing "minimal number of points is order + 1"
      (t/is (nil? (property #(sut/barycentric [1] [1] {:order 0}))))
      (t/is (= :invalid-order (property #(sut/barycentric [1] [1]))) "default order 1 needs two points")
      (t/is (= :invalid-order (property #(sut/barycentric [1 2] [1 2] {:order 2}))))
      (t/is (nil? (property #(sut/barycentric [1 2] [1 2] {:order 1}))))
      (t/is (= :too-few-points (property #(sut/barycentric [] [] {:order 0})))))
    (t/is (= :length-mismatch (property #(sut/barycentric [0 1 2] [0 1]))))
    (t/is (= :not-finite (property #(sut/barycentric [0 ##NaN 2] [0 1 2]))))
    (t/is (= :not-finite (property #(sut/barycentric [0 1 ##-Inf] [0 1 2]))))
    (t/testing "duplicated nodes gave partial NaN"
      (t/is (= :duplicate-nodes (property #(sut/barycentric [0 1 1 2] [0 1 2 3]))))
      (t/is (= :duplicate-nodes (property #(sut/barycentric [2 1 0 1] [0 1 2 3])))))
    (t/testing "unsorted xs are accepted"
      (t/is (nil? (property #(sut/barycentric [2 0 1 3] [0 1 2 3])))))
    (t/testing "ex-data names the method and the order"
      (let [d (try (sut/barycentric xs ys {:order 7}) (catch clojure.lang.ExceptionInfo e (ex-data e)))]
        (t/is (= :barycentric (:method d)))
        (t/is (= 7 (:order d)))))))

(t/deftest barycentric-docstring
  (let [doc (:doc (meta #'sut/barycentric))]
    (t/is (re-find #"Floater-Hormann" doc))
    (t/is (not (re-find #"Berrut and Trefethen" doc)))))
