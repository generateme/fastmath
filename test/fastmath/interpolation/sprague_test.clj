(ns fastmath.interpolation.sprague-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.interpolation.sprague :as sut]))

(set! *warn-on-reflection* true)

;; Oracle: colour-science SpragueInterpolator (CIE 167:2005 boundary rule); outside the domain the
;; same formula with the index of the first / last interval, see
;; utils/fastmath/dev/generate_interp1d_reference.py.

(def ^:private ref-data
  (delay (edn/read-string (slurp (io/resource "interpolation/interp1d_reference.edn")))))

(defn- property
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (= :invalid-input (:type d)) (:property d))))))

(defn- ymax ^double [ys] (double (reduce max 0.0 (map #(Math/abs (double %)) ys))))

(t/deftest sprague-reference
  ;; measured: 1.9e-15 * max(ys) inside, 4.0e-14 relative outside (5 steps away)
  (doseq [d (:sprague @ref-data)]
    (let [f (sut/sprague (:xs d) (:ys d))
          sc (ymax (:ys d))]
      (doseq [[q e] (map vector (:qs d) (:expected d))]
        (t/is (m/delta-eq e (f q) (* 1e-12 sc) 1e-12) (str (:name d) " at " q " inside")))
      (doseq [[q e] (map vector (:qs-outside d) (:expected-outside d))]
        (t/is (m/delta-eq e (f q) (* 1e-12 sc) 1e-11) (str (:name d) " at " q " outside"))))))

(t/deftest sprague-polynomial-exactness
  (let [xs (range 0 14)
        poly (fn [k] (fn [^double x] (Math/pow x (double k))))]
    (t/testing "degree 0 and 1 are reproduced on all intervals including the boundary ones"
      (doseq [k [0 1]
              :let [g (poly k)
                    f (sut/sprague xs (map g xs))]
              x [0.0 0.25 0.9 1.5 5.5 11.5 12.2 12.9 13.0]]
        (t/is (m/delta-eq (g x) (f x) 1e-12) (str "degree " k " at " x))))
    (t/testing "degree 4 is reproduced on the interior intervals (all six points are data)"
      (let [g (poly 4)
            f (sut/sprague xs (map g xs))]
        (doseq [x [2.0 2.3 4.5 7.77 10.0 10.99]]
          ;; measured relative error 2e-16
          (t/is (m/delta-eq (g x) (f x) 0.0 1e-12) (str "at " x)))))
    (t/testing "the boundary rule is second order only: x^2 is not reproduced on the boundary intervals"
      (let [g (poly 2)
            f (sut/sprague xs (map g xs))]
        ;; measured errors: 0.097 (first interval, max), 0.097 (last interval, max)
        (t/is (> (m/abs (- (double (g 0.5)) (double (f 0.5)))) 1e-2))
        (t/is (> (m/abs (- (double (g 12.5)) (double (f 12.5)))) 1e-2))
        (t/is (m/delta-eq (g 6.5) (f 6.5) 1e-11))))))

(t/deftest sprague-edge-cases
  (t/testing "minimal six points: passes through the data"
    (let [xs [0 1 2 3 4 5]
          ys [3.0 -1.0 2.0 2.5 0.0 7.0]
          f (sut/sprague xs ys)]
      (doseq [[x y] (map vector xs ys)]
        (t/is (m/delta-eq y (f x) 1e-12)))))
  (t/testing "knots of a longer series"
    (let [xs (map #(* 0.25 (double %)) (range 12))
          ys (map #(Math/sin %) xs)
          f (sut/sprague xs ys)]
      (doseq [[x y] (map vector xs ys)]
        (t/is (m/delta-eq y (f x) 1e-12)))))
  (t/testing "NaN query gives NaN"
    (t/is (m/nan? ((sut/sprague [0 1 2 3 4 5] [0 1 2 3 4 5]) ##NaN))))
  (t/testing "negative and shifted nodes give the same polynomial in the shifted variable"
    (let [ys [0.0 1.0 0.5 2.0 -1.0 3.0 1.5]
          f (sut/sprague (range 7) ys)
          g (sut/sprague (map #(- (double %) 16.0) (range 7)) ys)]
      (doseq [x [0.0 0.3 2.5 4.4 5.9 6.0]]
        (t/is (m/delta-eq (f x) (g (- (double x) 16.0)) 1e-12))))))

(t/deftest sprague-validation
  (t/testing "six points are required (one to five gave an ArrayIndexOutOfBoundsException)"
    (doseq [n [0 1 2 3 5]]
      (t/is (= :too-few-points (property #(sut/sprague (range n) (range n)))) (str n " points")))
    (t/is (nil? (property #(sut/sprague (range 6) (range 6))))))
  (t/is (= :length-mismatch (property #(sut/sprague (range 7) (range 6)))))
  (t/is (= :not-finite (property #(sut/sprague [0 1 ##NaN 3 4 5 6] (range 7)))))
  (t/is (= :not-finite (property #(sut/sprague [0 1 2 3 4 5 ##Inf] (range 7)))))
  (t/testing "unsorted and duplicated xs"
    (t/is (= :not-increasing (property #(sut/sprague [0 1 3 2 4 5 6] (range 7)))))
    (t/is (= :not-increasing (property #(sut/sprague [0 1 2 2 4 5 6] (range 7)))))
    (t/is (= :not-increasing (property #(sut/sprague (reverse (range 7)) (range 7)))))))

(t/deftest sprague-docstring
  (let [doc (:doc (meta #'sut/sprague))]
    (t/is (re-find #"CIE" doc))
    (t/is (re-find #"colour" doc))
    (t/is (re-find #"uniform" doc))
    (t/is (re-find #"second order" doc))
    (t/is (re-find #"at least 6" doc))
    (t/is (re-find #"extrapolation" doc))))
