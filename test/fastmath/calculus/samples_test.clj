(ns fastmath.calculus.samples-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.calculus.samples :as sut]))

;; Reference: scipy 1.18.1 (utils/fastmath/dev/generate_samples_reference.py).
;; Measured max absolute error against the reference over all cases:
;; trapezoid total 3.6e-15, simpson total 1.8e-15 (the sums are accumulated in a different order than
;; numpy's pairwise sum), cumulative variants 0.0 (bitwise identical).
;; Polynomial cases where Simpson is exact: max error 8.9e-16 against the exact rational value.

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "calculus/samples_reference.edn")))))

(def ^:private tol 1.0e-13) ; ~30x the measured max error 3.6e-15 (summation order)

(defn- close?
  ([a b] (close? a b tol))
  ([a b ^double eps]
   (if (number? a)
     (m/delta-eq a b eps)
     (and (= (count a) (count b))
          (every? true? (map #(m/delta-eq %1 %2 eps) a b))))))

(defn- ex-data-of
  [f]
  (try (f) nil (catch clojure.lang.ExceptionInfo e (ex-data e))))

(defn- run-case
  [c method cumulative?]
  (sut/integrate-samples (:ys c) (or (:x c) (:dx c)) {:method method :cumulative? cumulative?}))

(defn- check-group
  [group]
  (doseq [c (get @reference group)
          [method cumulative? k] [[:trapezoid false :trapezoid]
                                  [:simpson false :simpson]
                                  [:trapezoid true :trapezoid-cumulative]
                                  [:simpson true :simpson-cumulative]]]
    (let [r (run-case c method cumulative?)]
      (t/is (close? r (k c)) (str (:id c) " " method " cumulative=" cumulative?))
      (when cumulative?
        (t/is (vector? r))
        (t/is (== 0.0 (first r)))
        (t/is (= (count (:ys c)) (count r)))))))

(t/deftest scipy-reference
  (t/testing "uniform grid given as a spacing"
    (check-group :uniform-dx))
  (t/testing "uniform grid given as abscissae"
    (check-group :uniform-x))
  (t/testing "non-uniform abscissae"
    (check-group :nonuniform-x))
  (t/testing "negative spacing and decreasing abscissae"
    (check-group :negative))
  (t/testing "polynomials and hand cases"
    (check-group :known))
  (t/testing "reference covers all sample counts"
    (t/is (= #{2 3 4 5 6 7 10 11 50 101}
             (set (map (comp count :ys) (:uniform-dx @reference)))))))

(t/deftest exact-polynomials
  (doseq [c (:known @reference)
          :when (:exact c)]
    (t/is (close? (run-case c :simpson false) (:exact c))
          (str (:id c) " simpson"))
    (when (re-find #"^linear" (:id c))
      (t/is (close? (run-case c :trapezoid false) (:exact c))
            (str (:id c) " trapezoid is exact for linear data")))))

(t/deftest hand-values
  (t/testing "trapezoid"
    (t/is (== 4.0 (sut/integrate-samples [1 2 3])))
    (t/is (== 2.0 (sut/integrate-samples [1 2 3] 0.5)))
    (t/is (== 4.0 (sut/integrate-samples [1 2 3] [0 1 2])))
    (t/is (== 4.0 (sut/integrate-samples [1 2 3] nil)))
    (t/is (== 4.0 (sut/integrate-samples [1 2 3] 1 {:method :trapezoid})))
    (t/is (== 5.0 (sut/integrate-samples [1 2 3 4] [0 1 1 2]))))
  (t/testing "cumulative trapezoid, x with a repeated abscissa"
    (t/is (close? [0.0 3.0 3.0 7.0]
                  (sut/integrate-samples [1.0 2.0 5.0 3.0] [0 2 2 3] {:cumulative? true}))))
  (t/testing "simpson on three points is the parabola rule"
    (t/is (close? (/ 26.0 3.0) (sut/integrate-samples [1 4 9] 1 {:method :simpson})))
    ;; (x+1)^2 on [0, 2] sampled at 0, 1, 2 is integrated exactly
    (t/is (close? (/ 26.0 3.0) (sut/integrate-samples [1 4 9] [0 1 2] {:method :simpson}))))
  (t/testing "simpson, scipy docs"
    (t/is (close? 40.5 (sut/integrate-samples (range 10) (range 10) {:method :simpson})))
    (t/is (close? 1640.5 (sut/integrate-samples (map #(* % % %) (range 10)) (range 10)
                                                {:method :simpson}))))
  (t/testing "cumulative simpson with a negative spacing (scipy value)"
    (t/is (close? [0.0 (/ -2.0 3.0) (/ -7.0 3.0) (/ -109.0 24.0)]
                  (sut/integrate-samples [1 2 5 3] -0.5 {:method :simpson :cumulative? true}))))
  (t/testing "simpson on decreasing abscissae, even number of intervals (scipy value)"
    (t/is (close? -8.625 (sut/integrate-samples [1 2 5 3 7] [3 2 1.5 1 0] {:method :simpson})))))

(t/deftest edge-cases
  (doseq [method [:trapezoid :simpson]]
    (t/testing (str "no samples, " method)
      (t/is (== 0.0 (sut/integrate-samples [] 1.0 {:method method})))
      (t/is (== 0.0 (sut/integrate-samples nil nil {:method method})))
      (t/is (== 0.0 (sut/integrate-samples [] [] {:method method})))
      (t/is (instance? Double (sut/integrate-samples [] 1.0 {:method method})))
      (t/is (= [] (sut/integrate-samples [] 1.0 {:method method :cumulative? true})))
      (t/is (= [] (sut/integrate-samples (double-array 0) [] {:method method :cumulative? true}))))
    (t/testing (str "one sample, " method)
      (t/is (== 0.0 (sut/integrate-samples [5.0] 1.0 {:method method})))
      (t/is (== 0.0 (sut/integrate-samples [5.0] [3.0] {:method method})))
      (t/is (instance? Double (sut/integrate-samples [5.0] 1.0 {:method method})))
      (t/is (= [0.0] (sut/integrate-samples [5.0] 2.0 {:method method :cumulative? true}))))
    (t/testing (str "two samples are always the trapezoid rule, " method)
      (t/is (== 6.0 (sut/integrate-samples [1.0 5.0] 2.0 {:method method})))
      (t/is (== 12.0 (sut/integrate-samples [1.0 5.0] [0.0 4.0] {:method method})))
      (t/is (= [0.0 6.0] (sut/integrate-samples [1.0 5.0] 2.0 {:method method :cumulative? true})))
      (t/is (= [0.0 12.0] (sut/integrate-samples [1.0 5.0] [0.0 4.0]
                                                 {:method method :cumulative? true}))))
    (t/testing (str "zero and negative spacing, " method)
      (t/is (== 0.0 (sut/integrate-samples [1.0 2.0 5.0] 0.0 {:method method})))
      (t/is (== -4.0 (sut/integrate-samples [1 2 3] -1 {:method method})))
      (t/is (== -4.0 (sut/integrate-samples [1 2 3] [2 1 0] {:method method}))))
    (t/testing (str "constant data, " method)
      (t/is (close? 10.0 (sut/integrate-samples (repeat 11 2.0) 0.5 {:method method})))
      (t/is (close? 11.0 (sut/integrate-samples (repeat 12 2.0) 0.5 {:method method})))))
  (t/testing "simpson, even and odd number of samples agree on exact data"
    ;; x^2 on [0, 1]: Simpson is exact for odd N and, thanks to the correction, for even N too
    (doseq [n [3 4 5 6 7 8 9 10]]
      (let [xs (map #(/ (double %) (dec n)) (range n))
            ys (map #(* % %) xs)]
        (t/is (close? (/ 1.0 3.0) (sut/integrate-samples ys xs {:method :simpson})) (str "x " n))
        (t/is (close? (/ 1.0 3.0) (sut/integrate-samples ys (/ 1.0 (dec n)) {:method :simpson}))
              (str "dx " n)))))
  (t/testing "simpson with zero spacings is finite (scipy values)"
    (t/is (close? 9.333333333333334 (sut/integrate-samples [1 2 3 4 5] [0 1 1 2 3] {:method :simpson})))
    (t/is (close? 2.6666666666666665 (sut/integrate-samples [1 2 3 4] [0 1 1 2] {:method :simpson}))))
  (t/testing "reversing the abscissae flips the sign"
    (let [ys [0.3 -1.2 2.2 0.7 1.9 -0.4 0.8]
          xs [0.0 0.4 1.1 1.3 2.0 2.9 3.1]]
      (doseq [method [:trapezoid :simpson]]
        (t/is (close? (sut/integrate-samples ys xs {:method method})
                      (- (sut/integrate-samples (reverse ys) (reverse xs) {:method method}))))))))

(t/deftest cumulative-consistency
  (t/testing "the last element of the trapezoid cumulative integral is the total"
    (doseq [c (mapcat #(get @reference %) [:uniform-dx :uniform-x :nonuniform-x :negative])]
      (t/is (== (run-case c :trapezoid false) (last (run-case c :trapezoid true)))
            (:id c))))
  (t/testing "the last element of the Simpson cumulative integral is the total for odd number of samples"
    (doseq [c (mapcat #(get @reference %) [:uniform-dx :uniform-x :nonuniform-x :negative])
            :when (odd? (count (:ys c)))]
      (t/is (close? (run-case c :simpson false) (last (run-case c :simpson true)))
            (:id c))))
  (t/testing "a cumulative integral of a polynomial is exact where Simpson is"
    ;; x^2 on [0, 4] with a spacing of 1: the parabola through each triple is x^2
    (let [r (sut/integrate-samples [0 1 4 9 16] 1 {:method :simpson :cumulative? true})]
      (t/is (close? [0.0 (/ 1.0 3.0) (/ 8.0 3.0) 9.0 (/ 64.0 3.0)] r)))))

(t/deftest input-types
  (let [ys [0.5 1.5 -0.5 2.0 1.0]
        xs [0.0 0.5 1.5 1.75 3.0]
        expected-uniform (sut/integrate-samples ys 0.5 {:method :simpson})
        expected-x (sut/integrate-samples ys xs {:method :simpson})]
    (t/testing "types of samples"
      (t/is (== expected-uniform (sut/integrate-samples (double-array ys) 0.5 {:method :simpson})))
      (t/is (== expected-uniform (sut/integrate-samples (map identity ys) 0.5 {:method :simpson})))
      (t/is (== expected-uniform (sut/integrate-samples (list* ys) 0.5 {:method :simpson})))
      (t/is (== (sut/integrate-samples (map double (range 5)) 1)
                (sut/integrate-samples (range 5) 1))))
    (t/testing "types of abscissae"
      (t/is (== expected-x (sut/integrate-samples ys (double-array xs) {:method :simpson})))
      (t/is (== expected-x (sut/integrate-samples ys (map identity xs) {:method :simpson})))
      (t/is (== expected-x (sut/integrate-samples ys (list* xs) {:method :simpson})))))
  (t/testing "ratios and integers"
    (t/is (== 2.0 (sut/integrate-samples [1/2 1 3/2] 1)))
    (t/is (== 1.0 (sut/integrate-samples [1/2 1 3/2] 1/2)))
    (t/is (== 2.0 (sut/integrate-samples [1/2 1 3/2] [0 1/1 2])))
    (t/is (== 1.0 (sut/integrate-samples [1/2 1 3/2] [0 1/2 1])))
    (t/is (= [0.0 0.75 2.0] (sut/integrate-samples [1/2 1 3/2] 1 {:cumulative? true})))
    (t/is (close? (/ 26.0 3.0) (sut/integrate-samples [1 4 9] 1/1 {:method :simpson})))
    (t/is (== 4.0 (sut/integrate-samples (long-array [1 2 3]) 1))))
  (t/testing "the samples array is not modified"
    (let [ys (double-array [1.0 2.0 4.0 3.0])
          xs (double-array [0.0 1.0 3.0 4.0])]
      (sut/integrate-samples ys xs {:method :simpson :cumulative? true})
      (sut/integrate-samples ys 0.5 {:method :simpson})
      (t/is (= [1.0 2.0 4.0 3.0] (vec ys)))
      (t/is (= [0.0 1.0 3.0 4.0] (vec xs)))))
  (t/testing "results are doubles"
    (t/is (instance? Double (sut/integrate-samples [1 2 3] 1)))
    (t/is (every? #(instance? Double %)
                  (sut/integrate-samples [1 2 3] 1 {:cumulative? true :method :simpson})))))

(t/deftest validation
  (t/testing "mismatched lengths"
    (t/is (= {:ys-count 3 :x-count 2} (ex-data-of #(sut/integrate-samples [1 2 3] [0 1]))))
    (t/is (= {:ys-count 3 :x-count 4} (ex-data-of #(sut/integrate-samples [1 2 3] [0 1 2 3]))))
    (t/is (= {:ys-count 0 :x-count 1} (ex-data-of #(sut/integrate-samples [] [0]))))
    (t/is (= {:ys-count 1 :x-count 0} (ex-data-of #(sut/integrate-samples [1] []))))
    (t/is (= {:ys-count 3 :x-count 2}
             (ex-data-of #(sut/integrate-samples [1 2 3] (double-array [0 1])
                                                 {:method :simpson :cumulative? true}))))
    (t/is (thrown? clojure.lang.ExceptionInfo (sut/integrate-samples [1 2 3] [0 1]))))
  (t/testing "unknown method"
    (t/is (= {:option :method :value :romberg :expected #{:trapezoid :simpson}}
             (ex-data-of #(sut/integrate-samples [1 2 3] 1 {:method :romberg}))))
    (t/is (= {:option :method :value nil :expected #{:trapezoid :simpson}}
             (ex-data-of #(sut/integrate-samples [1 2 3] 1 {:method nil}))))
    (t/is (= {:option :method :value "simpson" :expected #{:trapezoid :simpson}}
             (ex-data-of #(sut/integrate-samples [] 1 {:method "simpson" :cumulative? true}))))
    (t/is (thrown? clojure.lang.ExceptionInfo (sut/integrate-samples [1 2 3] 1 {:method :boole}))))
  (t/testing "samples and spacing of a wrong type"
    (t/is (= :ys (:option (ex-data-of #(sut/integrate-samples 5 1)))))
    (t/is (= :ys (:option (ex-data-of #(sut/integrate-samples "abc" 1)))))
    (t/is (= :x-or-dx (:option (ex-data-of #(sut/integrate-samples [1 2 3] "abc")))))
    (t/is (= :x-or-dx (:option (ex-data-of #(sut/integrate-samples [1 2 3] {:a 1}))))))
  (t/testing "valid input does not throw"
    (t/is (nil? (ex-data-of #(sut/integrate-samples [1 2 3] 1 {:method :simpson}))))
    (t/is (nil? (ex-data-of #(sut/integrate-samples [1 2 3] 1 nil))))
    (t/is (nil? (ex-data-of #(sut/integrate-samples [1 2 3] 1 {}))))))
