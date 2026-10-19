(ns fastmath.interpolation.shepard-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.distance :as dist]
            [fastmath.vector :as v]
            [fastmath.interpolation.shepard :as sut]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; Reference values: utils/fastmath/dev/generate_scattered_reference.py (numpy inverse distance weighting
;; with the weights 1 / d^p, float64).
(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "interpolation/scattered_reference.edn")))))

(defn- property
  "Returns the :property of the invalid-input ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (= :shepard (:method d))) (:property d))))))

(defn- finite?
  [x]
  (m/valid-double? (double x)))

(defn- max-err
  ^double [f queries expected]
  (reduce (fn [^double acc [q e]]
            (m/max acc (m/abs (m/- (double (f q)) (double e)))))
          0.0 (map vector queries expected)))

(t/deftest numpy-reference
  (doseq [{:keys [name p xs ys queries values]} (:shepard @reference)
          :let [d1? (= name :d1)
                f (sut/shepard (if d1? (mapv first xs) xs) ys {:p p})
                qs (if d1? (mapv first queries) queries)]]
    ;; measured maximum 3.1e-15 over p in 0, 1, 2, 3.5 and 1-3 dimensions
    (t/is (< (max-err f qs values) 1e-14) (str name " p " p))))

(t/deftest default-power
  (let [{:keys [xs ys queries values]} (first (filter #(and (= :d2 (:name %)) (= 2.0 (:p %))) (:shepard @reference)))
        f (sut/shepard xs ys)]
    (t/is (< (max-err f queries values) 1e-14))
    (t/is (< (max-err (sut/shepard xs ys nil) queries values) 1e-14))
    (t/is (< (max-err (sut/shepard xs ys {:p nil}) queries values) 1e-14) "a nil option is a missing option")))

(t/deftest nodes-are-exact
  (doseq [{:keys [xs ys p name]} (:shepard @reference)
          :let [d1? (= name :d1)
                f (sut/shepard (if d1? (mapv first xs) xs) ys {:p p})]]
    (t/is (every? true? (map (fn [x y] (= (double y) (f (if d1? (first x) x)))) xs ys)) (str name " p " p))))

(t/deftest exact-hit-by-distance
  (t/testing "a long query on a double node"
    (let [f (sut/shepard [1.0 3.0] [2.0 4.0])]
      (t/is (= 2.0 (f 1)))
      (t/is (= 4.0 (f 3)))
      (t/is (= 2.0 (f 1.0)))))
  (t/testing "long vs double vectors, both ways"
    (let [fd (sut/shepard [[1.0 0.0] [3.0 3.0]] [2.0 4.0])
          fl (sut/shepard [[1 0] [3 3]] [2.0 4.0])]
      (t/is (= 2.0 (fd [1 0])))
      (t/is (= 2.0 (fl [1.0 0.0])))
      (t/is (= 4.0 (fd [3 3])))
      (t/is (= 4.0 (fl [3.0 3.0])))
      (t/is (= 2.0 (fd '(1 0))) "a list")))
  (t/testing "Vec2 vs vector, both ways"
    (let [fv (sut/shepard [[1.0 0.0] [3.0 3.0]] [2.0 4.0])
          f2 (sut/shepard [(v/vec2 1.0 0.0) (v/vec2 3.0 3.0)] [2.0 4.0])]
      (t/is (= 2.0 (fv (v/vec2 1.0 0.0))))
      (t/is (= 2.0 (f2 [1 0])))
      (t/is (= 4.0 (f2 (v/vec2 3 3))))
      (t/is (= (fv [1.5 0.5]) (f2 (v/vec2 1.5 0.5))))))
  (t/testing "an interior point is a weighted mean"
    (t/is (= 3.0 ((sut/shepard [0.0 2.0] [2.0 4.0]) 1.0)))))

(t/deftest tiny-distances
  (let [xs [[0.0 0.0] [1.0 0.0]]
        ys [10.0 20.0]]
    (doseq [p [0.0 1.0 2.0 3.0 3.5 10.0]
            e [1e-100 1e-150 1e-154 1e-155 1e-160 1e-200 1e-300 4.9e-324]
            :let [f (sut/shepard xs ys {:p p})
                  r (f [e 0.0])]]
      (t/is (finite? r) (str "p " p " distance " e))
      ;; the nearest node gets all the weight, but p = 0 gives the mean
      (t/is (m/delta-eq (if (m/zero? (double p)) 15.0 10.0) r 1e-12) (str "p " p " distance " e)))
    (t/testing "on the other axis and with the other node nearest"
      (let [f (sut/shepard xs ys {:p 3.0})]
        (t/is (m/delta-eq 10.0 (f [0.0 1e-150]) 1e-12))
        (t/is (m/delta-eq 20.0 (f [1.0 1e-150]) 1e-12))
        (t/is (m/delta-eq 20.0 (f [1.0000000000000002 1e-300]) 1e-12))))
    (t/testing "a point 1e-100 away"
      (t/is (m/delta-eq 10.0 ((sut/shepard xs ys) [1e-100 0.0]) 1e-12)))
    (t/testing "1-d"
      (let [f (sut/shepard [0.0 1.0] ys {:p 3.0})]
        (t/is (m/delta-eq 10.0 (f 1e-150) 1e-12))
        (t/is (m/delta-eq 10.0 (f 1e-300) 1e-12))
        (t/is (m/delta-eq 10.0 (f 5e-324) 1e-12))))
    (t/testing "a custom distance (the plain Euclidean one, squares underflow below 1e-162)"
      (let [f (sut/shepard xs ys {:distance dist/euclidean})]
        (doseq [e [1e-100 1e-154 1e-155 1e-160 1e-300]]
          (t/is (m/delta-eq 10.0 (f [e 0.0]) 1e-12) (str "distance " e)))
        (t/is (m/delta-eq 10.0 (f [0.0 0.0]) 1e-12))))))

(t/deftest far-away
  (let [xs [[0.0 0.0] [1.0 0.0]]
        ys [10.0 20.0]
        f (sut/shepard xs ys)]
    (t/is (m/delta-eq 15.0 (f [1e100 0.0]) 1e-12) "tends to the mean")
    (t/is (m/delta-eq 15.0 (f [1e200 1e200]) 1e-12))
    (t/is (m/delta-eq 15.0 (f [1e308 1e308]) 1e-12) "the squares of the differences would overflow")
    (t/is (m/delta-eq 15.0 (f [-1e300 0.0]) 1e-12))
    (t/testing "the distances which overflow to infinity"
      (t/is (m/nan? (f [1.5e308 1.5e308]))))
    (t/testing "far away, the huge scale of the nodes"
      (let [g (sut/shepard [[1e200 0.0] [3e200 0.0]] ys)]
        (t/is (= 10.0 (g [1e200 0.0])))
        (t/is (m/delta-eq 10.0 (g [1.0000000001e200 0.0]) 1e-9))
        (t/is (m/delta-eq 15.0 (g [2e200 5e199]) 1e-12))))))

(t/deftest duplicates-and-minimal-sizes
  (t/testing "the last of the repeated nodes wins"
    (let [f (sut/shepard [0 1 1] [0 1 5])]
      (t/is (= 5.0 (f 1.0)))
      (t/is (= 5.0 (f 1)))
      (t/is (= 0.0 (f 0)))
      (t/is (finite? (f 0.5))))
    (t/is (= 7.0 ((sut/shepard [[0 0] [1 1] [1.0 1.0] [1 1]] [0 3 4 7]) [1 1])))
    (t/is (= 0.0 ((sut/shepard [[0 0] [1 1] [1.0 1.0] [1 1]] [0 3 4 7]) [0 0.0]))))
  (t/testing "one point is a constant function"
    (let [f (sut/shepard [5] [7.0])]
      (t/is (= 7.0 (f 5)))
      (t/is (= 7.0 (f 100)))
      (t/is (= 7.0 (f -1e300)))
      (t/is (= 7.0 (f 5.0000001))))
    (t/is (= 7.0 ((sut/shepard [[1 2 3]] [7.0]) [4 5 6]))))
  (t/testing "two points, the midpoint is the mean for any p"
    (doseq [p [0.0 0.5 1.0 2.0 7.0]]
      (t/is (m/delta-eq 3.0 ((sut/shepard [0.0 2.0] [2.0 4.0] {:p p}) 1.0) 1e-14))))
  (t/testing "constant values stay constant"
    (let [f (sut/shepard [0.0 1.0 2.0 3.0] [4.0 4.0 4.0 4.0] {:p 3})]
      (t/is (every? #(m/delta-eq 4.0 (f %) 1e-14) [-5.0 0.3 1.5 2.9 1e5])))))

(t/deftest not-a-number
  (let [f (sut/shepard [[0.0 0.0] [1.0 0.0]] [10.0 20.0])]
    (t/is (m/nan? (f [##NaN 0.0])))
    (t/is (m/nan? (f [0.0 ##NaN]))))
  (t/is (m/nan? ((sut/shepard [0.0 1.0] [10.0 20.0]) ##NaN)))
  (t/testing "a NaN value makes a NaN result, except at the other nodes"
    (let [f (sut/shepard [0.0 1.0 2.0] [10.0 ##NaN 30.0])]
      (t/is (m/nan? (f 0.5)))
      (t/is (= 10.0 (f 0)))
      (t/is (= 30.0 (f 2)))
      (t/is (m/nan? (f 1))))))

(t/deftest query-forms
  (let [f1 (sut/shepard [0.0 1.0 3.0] [1.0 2.0 5.0])]
    (t/testing "1-d: a number or a one element sequence"
      (t/is (= (f1 0.7) (f1 [0.7])))
      (t/is (= (f1 0.7) (f1 (list 0.7))))
      (t/is (= :dimension-mismatch (property #(f1 (v/vec2 0.7 9.0)))))
      (t/is (= :dimension-mismatch (property #(f1 []))))))
  (t/testing "wrong dimension with the default distance"
    (let [f (sut/shepard [[0.0 0.0] [1.0 1.0]] [1.0 2.0])]
      (t/is (= :dimension-mismatch (property #(f [1.0 2.0 3.0]))))
      (t/is (= :dimension-mismatch (property #(f [1.0]))))
      (t/is (= :dimension-mismatch (property #(f 1.0)))))
    (t/is (= :dimension-mismatch (property #(sut/shepard [[0.0 0.0] [1.0]] [1.0 2.0]))))))

(t/deftest custom-distance
  (let [xs [[0.0 0.0] [1.0 2.0] [3.0 1.0]]
        ys [1.0 4.0 2.0]
        q [0.5 0.5]
        manhattan (fn [a b] (reduce + (map (fn [x y] (m/abs (- (double x) (double y)))) a b)))
        w (map (fn [x] (m// 1.0 (m/sq (double (manhattan q x))))) xs)
        expected (m// (double (reduce + (map * w ys))) (double (reduce + w)))]
    (t/is (m/delta-eq expected ((sut/shepard xs ys {:distance manhattan}) q) 1e-14))
    (t/is (= 4.0 ((sut/shepard xs ys {:distance manhattan}) [1 2])))
    (t/testing "the default distance is the Euclidean one"
      (t/is (m/delta-eq ((sut/shepard xs ys {:distance dist/euclidean}) q) ((sut/shepard xs ys) q) 1e-14)))
    (t/testing "1-d nodes with a custom distance"
      (t/is (m/delta-eq ((sut/shepard [0.0 1.0 3.0] [1.0 2.0 5.0] {:distance dist/euclidean-1d}) 0.7)
                        ((sut/shepard [0.0 1.0 3.0] [1.0 2.0 5.0]) 0.7) 1e-14)))))

(t/deftest many-points
  (let [n 3000
        rnd (java.util.Random. 12345)
        xs (vec (repeatedly n #(vector (.nextDouble rnd) (.nextDouble rnd) (.nextDouble rnd))))
        ys (vec (repeatedly n #(.nextGaussian rnd)))
        f (sut/shepard xs ys {:p 2.5})
        q [0.3 0.6 0.2]
        w (mapv (fn [x] (m// 1.0 (Math/pow (dist/euclidean x q) 2.5))) xs)
        expected (m// (double (reduce + (map * w ys))) (double (reduce + w)))]
    (t/is (m/delta-eq expected (f q) 1e-13))
    (t/is (= (nth ys 17) (f (nth xs 17))))))

(t/deftest invalid-input
  (t/testing "lengths and minimal size"
    (t/is (= :length-mismatch (property #(sut/shepard [0 1] [1.0]))))
    (t/is (= :length-mismatch (property #(sut/shepard [0 1] [1.0 2.0 3.0]))))
    (t/is (= :length-mismatch (property #(sut/shepard [[0 0] [1 1] [2 2] [3 3] [4 4]] [1 2 3 4])))
          "a node without a value is not silently dropped")
    (t/is (= :too-few-points (property #(sut/shepard [] []))))
    (t/is (= :too-few-points (property #(sut/shepard nil nil))))
    (t/is (nil? (property #(sut/shepard [0] [1.0])))))
  (t/testing "the power: the first invalid value and the boundary"
    (doseq [p [-1.0 -1e-300 -1 ##NaN ##Inf ##-Inf]]
      (t/is (= :p (property #(sut/shepard [0 1] [1.0 2.0] {:p p}))) (str "p " p)))
    (t/is (nil? (property #(sut/shepard [0 1] [1.0 2.0] {:p 0}))))
    (t/is (nil? (property #(sut/shepard [0 1] [1.0 2.0] {:p 0.0}))))
    (t/is (nil? (property #(sut/shepard [0 1] [1.0 2.0] {:p 1e-300}))))
    (t/is (nil? (property #(sut/shepard [0 1] [1.0 2.0] {:p 50})))))
  (t/testing "the exception data"
    (let [d (try (sut/shepard [0 1] [1.0 2.0] {:p -2.0}) (catch clojure.lang.ExceptionInfo e (ex-data e)))]
      (t/is (= {:type :invalid-input :method :shepard :property :p :option :p :value -2.0}
               (select-keys d [:type :method :property :option :value]))))))

(t/deftest large-power-is-stable
  ;; the weights are scaled by the nearest distance: 1/d^p would overflow for p 50 and d 1e-10
  (let [f (sut/shepard [0.0 1.0 2.0] [1.0 2.0 3.0] {:p 50.0})]
    (t/is (m/delta-eq 1.0 (f 1e-10) 1e-12))
    (t/is (m/delta-eq 2.0 (f 1.0000001) 1e-6))
    (t/is (finite? (f 0.5)))))
