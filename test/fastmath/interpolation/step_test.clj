(ns fastmath.interpolation.step-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.interpolation.step :as sut]))

(set! *warn-on-reflection* true)

;; Oracles: scipy interp1d kinds 'next' (step-before), 'previous' (step-after) and 'nearest-up'
;; (step, ties go to the right knot), constant fill outside the domain; values are compared exactly,
;; see utils/fastmath/dev/generate_interp1d_reference.py.

(def ^:private ref-data
  (delay (edn/read-string (slurp (io/resource "interpolation/interp1d_reference.edn")))))

(defn- property
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (= :invalid-input (:type d)) (:property d))))))

(t/deftest step-reference
  (doseq [d (:step @ref-data)]
    (let [xs (:xs d)
          ys (:ys d)
          before (sut/step-before xs ys)
          after (sut/step-after xs ys)
          nearest (sut/step xs ys)]
      (doseq [[i q] (map-indexed vector (:qs d))]
        (t/is (== (double (nth (:step-before d) i)) (double (before q))) (str (:name d) " step-before at " q))
        (t/is (== (double (nth (:step-after d) i)) (double (after q))) (str (:name d) " step-after at " q))
        (t/is (== (double (nth (:step d) i)) (double (nearest q))) (str (:name d) " step at " q))))))

(t/deftest step-continuity
  (let [xs [0 1 3 4 6]
        ys [10 20 30 40 50]
        before (sut/step-before xs ys)
        after (sut/step-after xs ys)]
    (t/testing "step-before is left-continuous: the node takes its own value and the step is just after the node"
      (t/is (== 10.0 (before 0)))
      (t/is (== 20.0 (before 1)))
      (t/is (== 30.0 (before 1.0000001)))
      (t/is (== 20.0 (before 0.9999999)))
      (t/is (== 20.0 (before 0.5))
            "the value between nodes is the y of the nearest node on the RIGHT (smallest xs[i] >= x)")
      (t/is (== 50.0 (before 5.5))))
    (t/testing "step-after is right-continuous: the node takes its own value and the step is at the node"
      (t/is (== 20.0 (after 1)))
      (t/is (== 10.0 (after 0.9999999)))
      (t/is (== 20.0 (after 2.5))
            "the value between nodes is the y of the nearest node on the LEFT (largest xs[i] <= x)")
      (t/is (== 50.0 (after 6)))
      (t/is (== 40.0 (after 5.9999999))))
    (t/testing "outside the domain: first / last value"
      (t/is (== 10.0 (before -100) (after -100)))
      (t/is (== 50.0 (before 100) (after 100))))))

(t/deftest step-ties-and-point
  (let [xs [0 1 3 4 6]
        ys [10 20 30 40 50]]
    (t/testing "default :point 0.5: the nearest node, a tie goes to the right node"
      (let [f (sut/step xs ys)]
        (t/is (== 20.0 (f 1.9999)))
        (t/is (== 30.0 (f 2)))
        (t/is (== 30.0 (f 2.0001)))
        (t/is (== 10.0 (f 0.4999)))
        (t/is (== 20.0 (f 0.5)))
        (t/is (== 10.0 (f -5)))
        (t/is (== 50.0 (f 60)))))
    (t/testing ":point 0.0 equals step-before and :point 1.0 equals step-after"
      (let [f0 (sut/step xs ys {:point 0.0})
            f1 (sut/step xs ys {:point 1.0})
            before (sut/step-before xs ys)
            after (sut/step-after xs ys)]
        (doseq [q (concat xs [-1 0.5 1.5 2 2.5 3.5 4.5 5 5.9 7] (range -1.0 7.0 0.173))]
          (t/is (== (double (before q)) (double (f0 q))) (str "point 0 at " q))
          (t/is (== (double (after q)) (double (f1 q))) (str "point 1 at " q)))))
    (t/testing ":point 0.25 switches a quarter of the way"
      (let [f (sut/step [0 4] [1 2] {:point 0.25})]
        (t/is (== 1.0 (f 0.99)))
        (t/is (== 2.0 (f 1.0)))
        (t/is (== 2.0 (f 3.0)))))))

(t/deftest step-edge-cases
  (t/testing "one point (step-before and step-after): a constant"
    (doseq [ctor [sut/step-before sut/step-after]]
      (let [f (ctor [3] [7])]
        (t/is (every? #(== 7.0 (f %)) [-1e10 0 3 3.5 1e10])))))
  (t/testing "two points"
    (t/is (== 1.0 ((sut/step-before [0 1] [1 2]) 0.0)))
    (t/is (== 2.0 ((sut/step-before [0 1] [1 2]) 0.5)))
    (t/is (== 1.0 ((sut/step-after [0 1] [1 2]) 0.5)))
    (t/is (== 2.0 ((sut/step [0 1] [1 2]) 0.5)))
    (t/is (== 1.0 ((sut/step [0 1] [1 2]) 0.4999999))))
  (t/testing "negative and huge nodes"
    (let [f (sut/step-after [-1e300 0 1e300] [1 2 3])]
      (t/is (== 1.0 (f -1e299)))
      (t/is (== 2.0 (f 1e299)))
      (t/is (== 3.0 (f 1e300)))))
  (t/testing "NaN and infinity in ys are returned as data"
    (t/is (m/nan? ((sut/step-after [0 1 2] [0 ##NaN 2]) 1.5)))
    (t/is (== ##Inf ((sut/step-before [0 1 2] [0 ##Inf 2]) 0.5))))
  (t/testing "infinite query"
    (t/is (== 5.0 ((sut/step-before [0 1] [4 5]) ##Inf)))
    (t/is (== 4.0 ((sut/step-before [0 1] [4 5]) ##-Inf)))
    (t/is (== 5.0 ((sut/step-after [0 1] [4 5]) ##Inf)))
    (t/is (== 4.0 ((sut/step-after [0 1] [4 5]) ##-Inf)))))

(t/deftest step-validation
  (t/is (= :length-mismatch (property #(sut/step-before [0 1 2] [0 1]))))
  (t/is (= :length-mismatch (property #(sut/step-after [0 1 2] [0 1]))))
  (t/is (= :length-mismatch (property #(sut/step [0 1 2] [0 1]))))
  (t/testing "minimal number of points: 1 for step-before and step-after, 2 for step"
    (t/is (= :too-few-points (property #(sut/step-before [] []))))
    (t/is (= :too-few-points (property #(sut/step-after [] []))))
    (t/is (nil? (property #(sut/step-before [1] [1]))))
    (t/is (nil? (property #(sut/step-after [1] [1]))))
    (t/is (= :too-few-points (property #(sut/step [1] [1]))))
    (t/is (nil? (property #(sut/step [1 2] [1 2])))))
  (doseq [ctor [sut/step-before sut/step-after sut/step]]
    (t/is (= :not-finite (property #(ctor [0 ##NaN 2] [0 1 2]))))
    (t/is (= :not-finite (property #(ctor [0 1 ##Inf] [0 1 2]))))
    (t/testing "unsorted xs gave silently wrong values"
      (t/is (= :not-increasing (property #(ctor [0 2 1 3] [0 1 2 3]))))
      (t/is (= :not-increasing (property #(ctor [0 1 1 3] [0 1 2 3]))))
      (t/is (= :not-increasing (property #(ctor [3 2 1 0] [0 1 2 3])))))))

(t/deftest step-docstrings
  ;; the explanatory sentences were swapped: step-before returns the node on the RIGHT, step-after on the LEFT
  (t/is (re-find #"nearest node on the right" (:doc (meta #'sut/step-before))))
  (t/is (re-find #"nearest node on the left" (:doc (meta #'sut/step-after))))
  (t/is (re-find #"'next'" (:doc (meta #'sut/step-before))))
  (t/is (re-find #"'previous'" (:doc (meta #'sut/step-after))))
  (t/is (re-find #"nearest-up" (:doc (meta #'sut/step)))))
