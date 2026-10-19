(ns fastmath.interpolation.microsphere-test
  (:require [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.interpolation.microsphere :as sut]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- property
  "Returns the :property of the invalid-input ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (= :microsphere-projection (:method d))) (:property d))))))

(def ^:private grid-2d
  (vec (for [i (range 6) j (range 6)] [(m/* 0.2 (double i)) (m/* 0.2 (double j))])))

(def ^:private values-2d
  (mapv (fn [[^double x ^double y]] (m/+ (m/* x x) y)) grid-2d))

(def ^:private grid-3d
  (vec (for [i (range 4) j (range 4) k (range 4)] [(m/* 0.3 (double i)) (m/* 0.3 (double j)) (m/* 0.3 (double k))])))

(def ^:private values-3d
  (mapv (fn [[^double x ^double y ^double z]] (m/+ x (m/+ y z))) grid-3d))

(def ^:private queries-2d
  [[0.5 0.5] [0.33 0.61] [0.1 0.9] [0.77 0.23] [0.45 0.05]])

(def ^:private queries-3d
  [[0.4 0.4 0.4] [0.13 0.61 0.5] [0.7 0.2 0.35] [0.45 0.05 0.8]])

(t/deftest two-dimensions-are-deterministic
  (let [f1 (sut/microsphere-projection grid-2d values-2d)
        f2 (sut/microsphere-projection grid-2d values-2d)
        f3 (sut/microsphere-projection grid-2d values-2d {:rng (r/rng :isaac 1)})
        f4 (sut/microsphere-projection grid-2d values-2d {:rng (r/rng :isaac 2)})]
    ;; the old function drew a random sphere for every construction: two 2-d constructions differed by 0.26
    (t/is (= (mapv f1 queries-2d) (mapv f2 queries-2d)) "bitwise identical")
    (t/is (= (mapv f1 queries-2d) (mapv f3 queries-2d) (mapv f4 queries-2d)) ":rng is not used in two dimensions")
    (t/is (every? #(m/valid-double? (double %)) (map f1 queries-2d)))))

(t/deftest random-sphere-with-seeded-rng
  (let [mk (fn [seed] (sut/microsphere-projection grid-3d values-3d {:rng (r/rng :isaac seed)}))
        a (mapv (mk 1) queries-3d)
        b (mapv (mk 1) queries-3d)
        c (mapv (mk 2) queries-3d)]
    (t/is (= a b) "equal seeds give equal values")
    (t/is (not= a c) "different seeds give different values")
    (t/is (every? true? (map #(> (m/abs (m/- (double %1) (double %2))) 1e-6) a c)) "all the queries differ")
    (t/testing "other generators"
      (t/is (= (mapv (sut/microsphere-projection grid-3d values-3d {:rng (r/rng :mersenne 3)}) queries-3d)
               (mapv (sut/microsphere-projection grid-3d values-3d {:rng (r/rng :mersenne 3)}) queries-3d))))
    (t/testing "the generator is consumed"
      (let [rng (r/rng :isaac 1)
            f1 (sut/microsphere-projection grid-3d values-3d {:rng rng})
            f2 (sut/microsphere-projection grid-3d values-3d {:rng rng})]
        (t/is (not= (mapv f1 queries-3d) (mapv f2 queries-3d)))))
    (t/testing "without :rng two constructions differ"
      (let [f1 (sut/microsphere-projection grid-3d values-3d)
            f2 (sut/microsphere-projection grid-3d values-3d)]
        (t/is (not= (mapv f1 queries-3d) (mapv f2 queries-3d)))))
    (t/testing "an invalid generator"
      (t/is (thrown? clojure.lang.ExceptionInfo (sut/microsphere-projection grid-3d values-3d {:rng :isaac}))))))

(t/deftest one-dimension
  (let [xs [0.0 1.0 2.0 3.0 4.0 5.0]
        ys [0.0 1.0 4.0 9.0 16.0 25.0]
        mk (fn [pts seed] (sut/microsphere-projection pts ys {:rng (r/rng :isaac seed)}))]
    (t/testing "vectors of one element as training points and as queries"
      (let [fv (mk (mapv vector xs) 4)
            fnum (mk xs 4)]
        (t/is (= (fv [2.5]) (fv 2.5) (fnum [2.5]) (fnum 2.5)) "the same sphere gives the same value")
        (t/is (= (fv [2.5]) (fv (list 2.5))))))
    (t/testing "node values"
      (doseq [x xs
              :let [x (double x)]]
        (t/is (= (m/* x x) ((mk (mapv vector xs) 5) [x])))))
    (t/testing "wrong dimension of the query"
      (t/is (= :dimension-mismatch (property #((mk (mapv vector xs) 5) [1.0 2.0]))))
      (t/is (= :dimension-mismatch (property #((mk xs 5) [1.0 2.0]))))
      (t/is (= :dimension-mismatch (property #((mk xs 5) [])))))))

(t/deftest node-values-are-exact
  (t/testing "two dimensions"
    (let [f (sut/microsphere-projection grid-2d values-2d)]
      (t/is (every? true? (map (fn [x y] (= y (f x))) grid-2d values-2d)))
      (t/is (every? true? (map (fn [x y] (= y (f (first x) (second x)))) grid-2d values-2d)))))
  (t/testing "three dimensions"
    (let [f (sut/microsphere-projection grid-3d values-3d {:rng (r/rng :isaac 9)})]
      (t/is (every? true? (map (fn [x y] (= y (f x))) grid-3d values-3d)))))
  (t/testing "one dimension, numbers"
    (let [f (sut/microsphere-projection [0.0 1.0 2.0 3.0 4.0] [5.0 3.0 1.0 -1.0 2.0] {:rng (r/rng :isaac 9)})]
      (t/is (= [5.0 3.0 1.0 -1.0 2.0] (mapv f [0.0 1.0 2.0 3.0 4.0])))))
  (t/testing "very near a node (tolerance 1e-6) and the tolerance option"
    (let [f (sut/microsphere-projection grid-2d values-2d)]
      (doseq [e [1e-100 1e-154 1e-155 1e-300]
              :let [e (double e)]]
        (t/is (= 0.44000000000000006 (f [(m/+ 0.2 e) 0.4])) "below the resolution of a double")
        (t/is (m/delta-eq 0.0 (f [e e]) 1e-12))
        (t/is (= (nth values-2d 8) (f [(m/+ 0.2 1e-9) 0.4])) "inside the default tolerance"))
      (t/is (not= (nth values-2d 8) (f [0.2001 0.4])) "outside the tolerance the value is interpolated")
      (t/is (= (nth values-2d 8) ((sut/microsphere-projection grid-2d values-2d {:no-interpolation-tolerance 1e-2}) [0.2001 0.4]))))))

(t/deftest options
  (let [mk (fn [opts] (sut/microsphere-projection grid-2d values-2d opts))
        base (mk nil)
        q [0.33 0.61]]
    (t/testing "defaults: exponent 1, elements max(5, 4 + d, n / 5), max-dark-friction 0.9, dark-threshold 0.01, background 0, tolerance 1e-6"
      (t/is (= (base q) ((mk {:exponent 1.0}) q) ((mk {:elements 7}) q) ((mk {:max-dark-friction 0.9 :dark-threshold 0.01 :background 0.0 :no-interpolation-tolerance 1e-6}) q))
            "36 nodes: max(5, 6, 7) = 7 elements")
      (t/is (not= (base q) ((mk {:elements 8}) q)))
      (t/is (not= (base q) ((mk {:elements 6}) q))))
    (t/testing "the elements of other sizes"
      (let [xs (vec (for [i (range 12) j (range 12)] [(m/* 0.1 (double i)) (m/* 0.1 (double j))]))
            ys (mapv (fn [[^double x ^double y]] (m/+ x y)) xs)
            dflt (sut/microsphere-projection xs ys)]
        ;; 144 nodes: max(5, 6, 28) = 28
        (t/is (= (dflt q) ((sut/microsphere-projection xs ys {:elements 28}) q)))
        (t/is (not= (dflt q) ((sut/microsphere-projection xs ys {:elements 27}) q)))))
    (t/testing "the exponent"
      (t/is (not= (base q) ((mk {:exponent 2.0}) q)))
      (t/is (not= ((mk {:exponent 2.0}) q) ((mk {:exponent 0.5}) q)))
      (t/is (every? #(m/valid-double? (double %)) (map (mk {:exponent 0.0}) queries-2d)) "exponent 0 is valid"))
    (t/testing "the background where the sphere is dark"
      (let [far [100.0 100.0]]
        (t/is (= 0.0 (base far)))
        (t/is (= -99.0 ((mk {:max-dark-friction 0.1 :background -99.0}) far)))
        (t/is (not= -99.0 ((mk {:max-dark-friction 0.1 :background -99.0}) [0.5 0.5])) "not dark inside the data")
        (t/is (= -7.0 ((mk {:dark-threshold 1e9 :background -7.0}) [0.5 0.5])) "every facet is dark")
        (t/is (= 0.0 ((mk {:dark-threshold 1e9}) [0.5 0.5])))
        (t/is (= 0.0 (base [##NaN 0.5])) "a NaN query is dark: the background")))
    (t/testing "nil options are missing options"
      (t/is (= (mapv base queries-2d)
               (mapv (mk {:elements nil :exponent nil :max-dark-friction nil :dark-threshold nil :background nil
                          :no-interpolation-tolerance nil :rng nil}) queries-2d))))
    (t/testing "the defaults in three dimensions: 64 nodes, max(5, 7, 12) = 12 elements"
      (let [mk3 (fn [opts] (sut/microsphere-projection grid-3d values-3d (assoc opts :rng (r/rng :isaac 1))))]
        (t/is (= (mapv (mk3 {}) queries-3d) (mapv (mk3 {:elements 12}) queries-3d)))
        (t/is (not= (mapv (mk3 {}) queries-3d) (mapv (mk3 {:elements 13}) queries-3d)))))))

(t/deftest query-forms
  (let [f (sut/microsphere-projection grid-2d values-2d)
        fl (sut/microsphere-projection (mapv (fn [[x y]] [(Math/round (* 5.0 (double x))) (Math/round (* 5.0 (double y)))]) grid-2d) values-2d)
        fv (sut/microsphere-projection (mapv (fn [[x y]] (v/vec2 x y)) grid-2d) values-2d)]
    (t/is (= (f [0.5 0.5]) (f 0.5 0.5) (f (v/vec2 0.5 0.5)) (f '(0.5 0.5)) (fv [0.5 0.5]) (fv (v/vec2 0.5 0.5))))
    (t/testing "integer nodes and queries"
      (t/is (= (fl [2 3]) (fl [2.0 3.0]) (fl 2 3)))
      (t/is (= (nth values-2d 0) (fl [0 0]))))
    (t/testing "wrong dimension"
      (t/is (= :dimension-mismatch (property #(f [0.5]))))
      (t/is (= :dimension-mismatch (property #(f [0.5 0.5 0.5]))))
      (t/is (= :dimension-mismatch (property #(f 0.5)))))))

(t/deftest sizes
  (t/testing "one node: its value where the sphere is lit, the background where it is dark"
    (let [f (sut/microsphere-projection [[1.0 2.0]] [7.0])]
      (t/is (= 7.0 (f [1 2])))
      (t/is (m/delta-eq 7.0 (f [5 5]) 1e-12))
      (t/is (= 0.0 (f [-1e5 3.0]))))
    ;; a random one dimensional sphere has 5 facets of the direction +1 or -1, it can be dark on one side: seeded
    (t/is (= 7.0 ((sut/microsphere-projection [3.0] [7.0] {:rng (r/rng :isaac 1)}) 5)))
    (t/is (= 7.0 ((sut/microsphere-projection [[3.0]] [7.0] {:rng (r/rng :isaac 1)}) [3]))))
  (t/testing "two nodes, every dimension"
    (doseq [pts [[[0.0] [1.0]] [[0.0 0.0] [1.0 1.0]] [[0.0 0.0 0.0] [1.0 1.0 1.0]]]
            :let [f (sut/microsphere-projection pts [2.0 4.0] {:rng (r/rng :isaac 6)})
                  mid (vec (repeat (count (first pts)) 0.5))]]
      (t/is (= 2.0 (f (first pts))))
      (t/is (= 4.0 (f (second pts))))
      ;; the midpoint is exactly the mean only for the symmetric deterministic sphere of two dimensions
      (t/is (<= 2.0 (double (f mid)) 4.0) (str "between the values, dimension " (count mid)))
      (when (= 2 (count mid))
        (t/is (m/delta-eq 3.0 (f mid) 1e-9)))))
  (t/testing "many nodes"
    (let [rnd (java.util.Random. 3)
          xs (vec (repeatedly 600 #(vector (.nextDouble rnd) (.nextDouble rnd))))
          ys (mapv (fn [[^double x ^double y]] (m/+ (m/sin (m/* 3.0 x)) y)) xs)
          f (sut/microsphere-projection xs ys)]
      ;; 120 elements, a smooth function: the error inside the data is small
      (t/is (< (m/abs (m/- (double (f [0.5 0.5])) (m/+ (m/sin 1.5) 0.5))) 0.1))
      (t/is (= (nth ys 10) (f (nth xs 10)))))))

(t/deftest threads
  ;; the evaluation changes the state of the sphere: the calls are serialised, the result does not depend on the order
  (let [f (sut/microsphere-projection grid-3d values-3d {:rng (r/rng :isaac 11)})
        qs (vec (take 200 (cycle queries-3d)))
        expected (mapv f qs)]
    (t/is (= expected (vec (pmap f qs))))
    (t/is (= expected (mapv f (reverse (reverse qs)))))))

(t/deftest invalid-input
  (t/testing "lengths and minimal size"
    (t/is (= :length-mismatch (property #(sut/microsphere-projection [0.0 1.0] [1.0]))))
    (t/is (= :length-mismatch (property #(sut/microsphere-projection [[0.0 0.0] [1.0 1.0]] [1.0 2.0 3.0]))))
    (t/is (= :too-few-points (property #(sut/microsphere-projection [] []))))
    (t/is (= :too-few-points (property #(sut/microsphere-projection nil nil))))
    (t/is (nil? (property #(sut/microsphere-projection [0.0] [1.0])))))
  (t/testing "dimensions of the nodes"
    (t/is (= :dimension-mismatch (property #(sut/microsphere-projection [[0.0 0.0] [1.0]] [1.0 2.0])))))
  (t/testing "elements: the first invalid value and the boundary"
    (doseq [e [0 -1 ##NaN]]
      (t/is (= :elements (property #(sut/microsphere-projection grid-2d values-2d {:elements e}))) (str "elements " e)))
    (t/is (nil? (property #(sut/microsphere-projection grid-2d values-2d {:elements 1})))))
  (t/testing "exponent"
    (doseq [e [-1.0 -1e-300 ##NaN ##Inf]]
      (t/is (= :exponent (property #(sut/microsphere-projection grid-2d values-2d {:exponent e}))) (str "exponent " e)))
    (t/is (nil? (property #(sut/microsphere-projection grid-2d values-2d {:exponent 0.0})))))
  (t/testing "the exception data"
    (let [d (try (sut/microsphere-projection grid-2d values-2d {:exponent -2.0}) (catch clojure.lang.ExceptionInfo e (ex-data e)))]
      (t/is (= {:type :invalid-input :method :microsphere-projection :property :exponent :option :exponent :value -2.0}
               (select-keys d [:type :method :property :option :value]))))))
