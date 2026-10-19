(ns fastmath.interpolation.common-test
  (:require [clojure.test :as t]
            [fastmath.interpolation.common :as sut]
            [fastmath.matrix :as mat]
            [fastmath.vector :as v]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- property
  "Returns the :property of the ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (:method d)) (:property d))))))

(defn- da [& xs] (double-array xs))

(t/deftest validate-xy-valid
  (t/is (nil? (sut/validate-xy! :linear (da 0 1 2) (da 5 6 7) {:min-points 2})))
  (t/testing "boundary: exactly the minimal number of points"
    (t/is (nil? (sut/validate-xy! :linear (da 0 1) (da 5 6) {:min-points 2})))
    (t/is (nil? (sut/validate-xy! :step (da 3) (da 5) {}))))
  (t/testing "NaN and infinite ys are allowed"
    (t/is (nil? (sut/validate-xy! :linear (da 0 1 2) (da ##NaN ##Inf 7) {:min-points 2}))))
  (t/testing "unsorted xs are accepted when :sorted? is false"
    (t/is (nil? (sut/validate-xy! :neville (da 2 0 1) (da 5 6 7) {:min-points 2 :sorted? false})))))

(t/deftest validate-xy-length
  (t/is (= :length-mismatch (property #(sut/validate-xy! :linear (da 0 1 2) (da 5 6) {}))))
  (t/is (= :length-mismatch (property #(sut/validate-xy! :linear (da 0 1) (da 5 6 7) {}))))
  (t/is (= :length-mismatch (property #(sut/validate-xy! :linear (da) (da 1) {})))))

(t/deftest validate-xy-count
  (t/is (= :too-few-points (property #(sut/validate-xy! :linear (da 0) (da 5) {:min-points 2}))))
  (t/is (= :too-few-points (property #(sut/validate-xy! :linear (da) (da) {:min-points 2}))))
  (t/is (= :too-few-points (property #(sut/validate-xy! :akima (da 0 1 2 3) (da 0 1 2 3) {:min-points 5}))))
  (t/testing "default minimal count is 1"
    (t/is (= :too-few-points (property #(sut/validate-xy! :step (da) (da) {}))))
    (t/is (nil? (sut/validate-xy! :step (da 1) (da 1) {})))))

(t/deftest validate-xy-finite-and-sorted
  (t/is (= :not-finite (property #(sut/validate-xy! :linear (da 0 ##NaN 2) (da 1 2 3) {}))))
  (t/is (= :not-finite (property #(sut/validate-xy! :linear (da 0 1 ##Inf) (da 1 2 3) {}))))
  (t/is (= :not-finite (property #(sut/validate-xy! :neville (da 0 ##NaN 2) (da 1 2 3) {:sorted? false}))))
  (t/is (= :not-increasing (property #(sut/validate-xy! :linear (da 0 2 1) (da 1 2 3) {}))))
  (t/testing "duplicates are not strictly increasing"
    (t/is (= :not-increasing (property #(sut/validate-xy! :linear (da 0 1 1) (da 1 2 3) {}))))
    (t/is (= :not-increasing (property #(sut/validate-xy! :linear (da 1 1) (da 1 2) {})))))
  (t/testing "descending"
    (t/is (= :not-increasing (property #(sut/validate-xy! :linear (da 2 1 0) (da 1 2 3) {}))))))

(t/deftest invalid-input-data
  (let [e (try (sut/validate-xy! :cubic (da 0 1 1) (da 1 2 3) {})
               (catch clojure.lang.ExceptionInfo e e))
        d (ex-data e)]
    (t/is (= :invalid-input (:type d)))
    (t/is (= :cubic (:method d)))
    (t/is (= :not-increasing (:property d)))
    (t/is (= 2 (:index d)))
    (t/is (.startsWith ^String (ex-message e) "cubic: "))))

(t/deftest check-axis-grid
  (t/is (nil? (sut/check-axis! :bilinear :ys (da 0 1 2) {:min-points 2})))
  (t/is (= :not-increasing (property #(sut/check-axis! :bilinear :ys (da 0 0) {}))))
  (t/is (= :too-few-points (property #(sut/check-axis! :bicubic :xs (da 0) {:min-points 2})))))

(t/deftest solve-lu-or-svd
  (let [A (mat/rows->RealMatrix [[2.0 0.0] [0.0 4.0]])
        b (v/vec->RealVector [2.0 4.0])]
    (t/is (v/delta-eq [1.0 1.0] (sut/solve-LU-or-SVD A b))))
  (t/testing "singular matrix: minimum norm solution, nothing printed to stdout"
    (let [A (mat/rows->RealMatrix [[1.0 1.0] [1.0 1.0]])
          b (v/vec->RealVector [2.0 2.0])
          out (java.io.StringWriter.)
          res (binding [*out* out] (sut/solve-LU-or-SVD A b))]
      (t/is (v/delta-eq [1.0 1.0] res 1.0e-12))
      (t/is (= "" (str out))))))

(t/deftest binary-search-boundary
  (let [arr (double-array [0.0 1.0 2.0 3.0])]
    (t/is (= 0 (sut/binary-search arr -5.0)))
    (t/is (= 0 (sut/binary-search arr 0.0)))
    (t/is (= 1 (sut/binary-search arr 1.5)))
    (t/is (= 2 (sut/binary-search arr 3.0)))
    (t/is (= 2 (sut/binary-search arr 10.0)))))
