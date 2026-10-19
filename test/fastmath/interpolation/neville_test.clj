(ns fastmath.interpolation.neville-test
  (:require [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.interpolation.acm :as sut]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- property
  "Returns [method property] of the invalid-input ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (= :invalid-input (:type d)) [(:method d) (:property d)])))))

(defn- fv
  ^double [f x]
  (double (f x)))

(defn- cube
  ^double [^double x]
  (m/- (m/* x x x) (m/* 2.0 x)))

(t/deftest neville-edge-cases
  (t/testing "points: empty, one, exactly two"
    (t/is (= [:neville :too-few-points] (property #(sut/neville [] []))))
    (t/is (= [:neville :too-few-points] (property #(sut/neville [1.0] [2.0]))))
    (let [f (sut/neville [0.0 2.0] [1.0 5.0])]
      (t/is (m/delta-eq 3.0 (fv f 1.0) 1.0e-13))
      (t/is (m/delta-eq 9.0 (fv f 4.0) 1.0e-13) "outside the domain the line continues")))
  (t/testing "lengths"
    (t/is (= [:neville :length-mismatch] (property #(sut/neville [0.0 1.0 2.0] [1.0 2.0]))))
    (t/is (= [:neville :length-mismatch] (property #(sut/neville [0.0 1.0] [1.0 2.0 3.0])))))
  (t/testing "NaN and infinite xs"
    (t/is (= [:neville :not-finite] (property #(sut/neville [0.0 ##NaN 2.0] [1.0 2.0 3.0]))))
    (t/is (= [:neville :not-finite] (property #(sut/neville [0.0 1.0 ##Inf] [1.0 2.0 3.0]))))
    (t/is (= [:neville :not-finite] (property #(sut/neville [##-Inf 1.0 2.0] [1.0 2.0 3.0])))))
  (t/testing "duplicate xs, also not adjacent (the order is not required)"
    (t/is (= [:neville :duplicate-points] (property #(sut/neville [0.0 1.0 1.0] [1.0 2.0 3.0]))))
    (t/is (= [:neville :duplicate-points] (property #(sut/neville [2.0 0.0 1.0 2.0] [1.0 2.0 3.0 4.0]))))
    (t/is (= [:neville :duplicate-points] (property #(sut/neville [1.0 1.0] [1.0 2.0])))))
  (t/testing "NaN in ys is allowed"
    (t/is (fn? (sut/neville [0.0 1.0 2.0] [1.0 ##NaN 3.0])))))

(t/deftest neville-unsorted-allowed
  (let [xs [3.0 -1.0 2.0 0.0 1.5]
        ys (map cube xs)
        sorted (sort-by first (map vector xs ys))
        f (sut/neville xs ys)
        fs (sut/neville (map first sorted) (map second sorted))]
    (doseq [x [-2.0 -1.0 0.3 1.5 2.9 3.0 5.0]]
      ;; degree 3 polynomial is reproduced by 5 points; observed error below 1e-12
      (t/is (m/delta-eq (cube x) (fv f x) 1.0e-11))
      (t/is (m/delta-eq (fv fs x) (fv f x) 1.0e-11)))))

(t/deftest divided-difference-edge-cases
  (t/testing "points: empty, one, exactly two"
    (t/is (= [:divided-difference :too-few-points] (property #(sut/divided-difference [] []))))
    (t/is (= [:divided-difference :too-few-points] (property #(sut/divided-difference [1.0] [2.0]))))
    (let [f (sut/divided-difference [0.0 2.0] [1.0 5.0])]
      (t/is (m/delta-eq 3.0 (fv f 1.0) 1.0e-13))
      (t/is (m/delta-eq 9.0 (fv f 4.0) 1.0e-13))))
  (t/testing "lengths"
    (t/is (= [:divided-difference :length-mismatch] (property #(sut/divided-difference [0.0 1.0 2.0] [1.0 2.0])))))
  (t/testing "NaN and infinite xs"
    (t/is (= [:divided-difference :not-finite] (property #(sut/divided-difference [0.0 ##NaN 2.0] [1.0 2.0 3.0]))))
    (t/is (= [:divided-difference :not-finite] (property #(sut/divided-difference [0.0 1.0 ##Inf] [1.0 2.0 3.0])))))
  (t/testing "sorted xs are required"
    (t/is (= [:divided-difference :not-increasing] (property #(sut/divided-difference [0.0 2.0 1.0] [1.0 2.0 3.0]))))
    (t/is (= [:divided-difference :not-increasing] (property #(sut/divided-difference [2.0 1.0 0.0] [1.0 2.0 3.0]))))
    (t/is (= [:divided-difference :not-increasing] (property #(sut/divided-difference [0.0 1.0 1.0] [1.0 2.0 3.0]))))))

(t/deftest polynomial-values
  ;; both are exact for a polynomial of degree n-1; observed error below 1e-12
  (let [xs [-1.0 0.0 0.5 2.0 3.0]
        ys (map cube xs)
        fn1 (sut/neville xs ys)
        fd (sut/divided-difference xs ys)]
    (doseq [x [-3.0 -1.0 -0.2 0.5 1.0 2.5 3.0 6.0]]
      (t/is (m/delta-eq (cube x) (fv fn1 x) 1.0e-11))
      (t/is (m/delta-eq (cube x) (fv fd x) 1.0e-11))))
  (t/testing "through the nodes, many points (the values keep the previous results)"
    (let [xs (mapv (fn [^double i] (m/* 0.5 i)) (range 12.0))
          ys (mapv (fn [^double x] (m/sin x)) xs)
          fn1 (sut/neville xs ys)
          fd (sut/divided-difference xs ys)]
      (doseq [[x y] (map vector xs ys)]
        (t/is (m/delta-eq y (fv fn1 x) 1.0e-9))
        (t/is (m/delta-eq y (fv fd x) 1.0e-9)))
      (t/is (m/delta-eq (fv fn1 2.75) (fv fd 2.75) 1.0e-9)))))
