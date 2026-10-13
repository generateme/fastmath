(ns fastmath.calculus.divided-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.calculus.divided :as sut]))

(set! *warn-on-reflection* true)

;; helpers

(defn- thrown-data
  [thunk]
  (try (thunk) nil
       (catch clojure.lang.ExceptionInfo e (ex-data e))))

(defn- invalid?
  [thunk option]
  (let [d (thrown-data thunk)]
    (and (= option (:option d))
         (contains? d :value)
         (string? (:expected d)))))

(defn- rel-err
  ^double [^double a ^double b]
  (if (zero? b) (Math/abs a) (/ (Math/abs (- a b)) (Math/abs b))))

(defn- mixed-err
  "Absolute error when the expected value is below 1, relative otherwise."
  ^double [^double a ^double b]
  (/ (Math/abs (- a b)) (Math/max 1.0 (Math/abs b))))

;; functions of the reference by id (defined identically in generate_divided_reference.py)

(def ^:private functions
  {:exp (fn ^double [^double x] (Math/exp x))
   :sin (fn ^double [^double x] (Math/sin x))
   :cos (fn ^double [^double x] (Math/cos x))
   :log (fn ^double [^double x] (Math/log x))
   :poly (fn ^double [^double x] (+ (Math/pow x 5.0) (* -2.0 (Math/pow x 3.0)) x))})

(defn- derivative-fn
  "The k-th derivative of the function `id`."
  [id ^long k]
  (case id
    :exp (fn ^double [^double x] (Math/exp x))
    :sin (fn ^double [^double x] (Math/sin (+ x (* k m/HALF_PI))))
    :cos (fn ^double [^double x] (Math/cos (+ x (* k m/HALF_PI))))
    :log (fn ^double [^double x] (/ (* (if (odd? k) 1.0 -1.0) (m/factorial (dec k))) (Math/pow x k)))
    :poly (fn ^double [^double x]
            (case k
              1 (+ (* 5.0 (Math/pow x 4.0)) (* -6.0 x x) 1.0)
              2 (+ (* 20.0 (Math/pow x 3.0)) (* -12.0 x))
              3 (+ (* 60.0 x x) -12.0)
              4 (* 120.0 x)
              5 120.0
              0.0))))

(defn- derivatives
  [id]
  (mapv #(derivative-fn id %) (range 1 8)))

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "calculus/divided_reference.edn")))))

(defn- max-run
  "Length of the longest run of equal nodes."
  [xs]
  (apply max (map count (partition-by identity xs))))

;; ------------------------------------------------------------------ Hermite divided differences

(t/deftest hermite-exact-derivatives-vs-sympy
  (t/is (= 15 (count @reference)))
  (t/testing "the whole table of every case, exact derivatives: measured max error (relative, floor 1) 1.5e-15"
    (doseq [{:keys [f xs levels]} @reference
            :let [res (sut/divided-all (functions f) xs {:derivatives (derivatives f)})
                  expected (reverse levels)]]
      (t/is (= (count expected) (count res)) (str f xs))
      (doseq [[level-res level-ref] (map vector res expected)]
        (t/is (= (count level-ref) (count level-res)) (str f xs))
        (t/is (every? #(<= ^double % 1.0e-14) (map mixed-err level-res level-ref)) (str f xs)))))
  (t/testing "the highest level equals divided"
    (doseq [{:keys [f xs levels]} @reference
            :let [top (ffirst (reverse levels))]]
      (t/is (<= (mixed-err (sut/divided (functions f) xs {:derivatives (derivatives f)}) top) 1.0e-14) (str f xs)))))

(t/deftest hermite-spec-cases
  (let [exp-ref (first (filter #(and (= :exp (:f %)) (= [0.0 0.0 1.0 1.0] (:xs %))) @reference))
        sin-ref (first (filter #(and (= :sin (:f %)) (= [0.0 0.0 0.0 1.0 1.0] (:xs %))) @reference))]
    (t/testing "exp at [0 0 1 1]: relative error to 1e-14 (measured 5.9e-16)"
      (t/is (< (rel-err (sut/divided (functions :exp) [0 0 1 1] {:derivatives (derivatives :exp)})
                        (first (last (:levels exp-ref))))
               1.0e-14))
      (t/is (every? true? (map (fn [a b] (< (rel-err a b) 1.0e-14))
                               (flatten (sut/divided-all (functions :exp) [0 0 1 1] {:derivatives (derivatives :exp)}))
                               (flatten (reverse (:levels exp-ref)))))))
    (t/testing "sin at [0 0 0 1 1]: relative error to 1e-14 (measured 6.8e-15)"
      (t/is (every? true? (map (fn [a b] (< (rel-err a b) 1.0e-14))
                               (flatten (sut/divided-all (functions :sin) [0 0 0 1 1] {:derivatives (derivatives :sin)}))
                               (flatten (reverse (:levels sin-ref)))))))))

(t/deftest hermite-closed-forms
  (t/testing "repeated node: f^(k)(x) / k!"
    (let [d (derivatives :exp)]
      (t/is (= (Math/exp 0.5) (sut/divided (functions :exp) [0.5] {:derivatives d})))
      (t/is (= (Math/exp 0.5) (sut/divided (functions :exp) [0.5 0.5] {:derivatives d})))
      (t/is (= (/ (Math/exp 0.5) 2.0) (sut/divided (functions :exp) [0.5 0.5 0.5] {:derivatives d})))
      (t/is (< (rel-err (sut/divided (functions :exp) [0.5 0.5 0.5 0.5 0.5] {:derivatives d})
                        (/ (Math/exp 0.5) 24.0))
               1.0e-15))))
  (t/testing "the leading coefficient of a polynomial of degree d is the divided difference of d + 1 nodes, repeated or not"
    (let [cube (fn ^double [^double x] (* x x x))
          d [(fn ^double [^double x] (* 3.0 x x)) (fn ^double [^double x] (* 6.0 x)) (constantly 6.0)]]
      (doseq [xs [[0 1 2 3] [0 0 1 1] [2 2 2 2] [-1 -1 -1 4] [1 1 5 5] [0 0 3 7]]]
        (t/is (< (Math/abs (- 1.0 ^double (sut/divided cube xs {:derivatives d}))) 1.0e-14) (str xs)))
      (t/testing "and a lower degree gives zero"
        (t/is (< (Math/abs ^double (sut/divided cube [0 0 1 1 2] {:derivatives (conj d (constantly 0.0))})) 1.0e-14)))))
  (t/testing "order of the runs and of the nodes of a run does not matter"
    (let [d (derivatives :sin)
          f (functions :sin)
          a (sut/divided f [0 0 1 1 1 2] {:derivatives d})]
      (t/is (< (rel-err (sut/divided f [2 1 1 1 0 0] {:derivatives d}) a) 1.0e-14))
      (t/is (< (rel-err (sut/divided f [1 1 1 0 0 2] {:derivatives d}) a) 1.0e-14))))
  (t/testing "-0.0 and 0.0 are the same node"
    (t/is (= 1.0 (sut/divided (functions :exp) [0.0 -0.0] {:derivatives (derivatives :exp)})))
    (t/is (invalid? #(sut/divided (functions :exp) [0.0 1.0 -0.0] {:derivatives (derivatives :exp)}) :xs)))
  (t/testing "nodes as ratios and longs"
    (t/is (< (rel-err (sut/divided (functions :exp) [1/2 1/2 1] {:derivatives (derivatives :exp)})
                      (sut/divided (functions :exp) [0.5 0.5 1.0] {:derivatives (derivatives :exp)}))
             1.0e-15))))

(t/deftest hermite-extra-derivatives-are-ignored
  (t/is (= (sut/divided (functions :exp) [0 0 1] {:derivatives (derivatives :exp)})
           (sut/divided (functions :exp) [0 0 1] {:derivatives [(derivative-fn :exp 1)]}))))

(t/deftest hermite-fallback-finite-differences
  (t/testing "without :derivatives: relative error below 1e-6 for runs of up to 3 nodes (k <= 2), measured max 3e-11"
    (doseq [{:keys [f xs levels]} @reference
            :when (<= (max-run xs) 3)
            :let [top (ffirst (reverse levels))
                  res (sut/divided (functions f) xs)]]
      (t/is (<= (mixed-err res top) 1.0e-6) (str f xs))
      (t/is (<= (rel-err res top) 1.0e-6) (str f xs))))
  (t/testing "the whole table of such a case"
    (let [{:keys [f xs levels]} (first @reference)]
      (t/is (every? true? (map (fn [a b] (<= (mixed-err a b) 1.0e-8))
                               (flatten (sut/divided-all (functions f) xs))
                               (flatten (reverse levels)))))))
  (t/testing "longer runs lose digits (k = 3, 5): relative error below 1e-5, measured 3.5e-8 for k = 5"
    (doseq [{:keys [f xs levels]} @reference
            :when (> (max-run xs) 3)
            :let [top (ffirst (reverse levels))]]
      (t/is (<= (mixed-err (sut/divided (functions f) xs) top) 1.0e-5) (str f xs))))
  (t/testing "nodes close to the border of the domain (log at 1e-3: the starting step 0.5 would leave the domain)"
    (t/is (< (rel-err (sut/divided (functions :log) [1.0e-3 1.0e-3]) 1000.0) 1.0e-6)))
  (t/testing "huge nodes"
    (t/is (< (rel-err (sut/divided (functions :log) [1.0e8 1.0e8]) 1.0e-8) 1.0e-6)))
  (t/testing "var and partial functions are functions"
    (t/is (< (rel-err (sut/divided #'identity [3.0 3.0]) 1.0) 1.0e-8))
    (t/is (< (rel-err (sut/divided (partial * 2.0) [3.0 3.0]) 2.0) 1.0e-8))))

;; ------------------------------------------------------------------ boundary cases

(t/deftest empty-and-single-node
  (t/testing "empty nodes: ex-info"
    (t/is (invalid? #(sut/divided m/log []) :xs))
    (t/is (invalid? #(sut/divided m/log nil) :xs))
    (t/is (invalid? #(sut/divided m/log '()) :xs))
    (t/is (invalid? #(sut/divided-all m/log []) :xs))
    (t/is (invalid? #(sut/divided-all m/log nil {:derivatives []}) :xs)))
  (t/testing "single node gives f(x0)"
    (t/is (= (Math/log 2.5) (sut/divided m/log [2.5])))
    (t/is (= (Math/log 2.5) (sut/divided m/log '(2.5))))
    (t/is (= [[(Math/log 2.5)]] (map vec (sut/divided-all m/log [2.5]))))
    (t/is (= 20.0 (sut/divided [10 20 30] [1])))
    (t/is (= -3.0 (sut/divided (fn [x] (- x)) [3]))))
  (t/testing "single node with :derivatives"
    (t/is (= (Math/exp 1.0) (sut/divided (functions :exp) [1.0] {:derivatives (derivatives :exp)})))))

(t/deftest non-adjacent-repeated-nodes
  (t/is (invalid? #(sut/divided m/log [1 2 1]) :xs))
  (t/is (invalid? #(sut/divided m/log [1 1 2 1] {:derivatives [m/log]}) :xs))
  (t/is (invalid? #(sut/divided-all m/log [3 4 5 4 6]) :xs))
  (t/is (invalid? #(sut/divided-all m/log [1.0 2.0 2.0 3.0 1.0] {:derivatives [identity]}) :xs))
  (t/testing "the value in ex-data is the nodes"
    (t/is (= [1 2 1] (:value (thrown-data #(sut/divided m/log [1 2 1])))))))

(t/deftest non-function-f
  (t/testing "vector and map lookups work with distinct nodes (as before)"
    (t/is (= 1.0 (sut/divided {0 1 1 3 2 7} [0 1 2])))
    (t/is (= 2.0 (sut/divided {0 1 1 3} [0 1])))
    (t/is (m/delta-eq (/ 10.0 3.0) (sut/divided [6 9 2 5] [0 1 2 3])))
    (t/is (= [[1.0]] (map vec (sut/divided-all {5 1.0} [5])))))
  (t/testing "repeated nodes without :derivatives: ex-info"
    (t/is (invalid? #(sut/divided [1 2 3] [0 0 1]) :f))
    (t/is (invalid? #(sut/divided {0 1 1 3} [0 0 1]) :f))
    (t/is (invalid? #(sut/divided-all [1 2 3] [1 1]) :f))
    (t/is (invalid? #(sut/divided-all :kw [1 1]) :f)))
  (t/testing "repeated nodes with :derivatives work for a vector or a map"
    (t/is (= 5.0 (sut/divided [3 7] [0 0] {:derivatives [(constantly 5.0)]})))
    (t/is (= 2.0 (sut/divided {0 1.0 1 3.0} [0 1] {:derivatives []}))))
  (t/testing "not callable f"
    (t/is (invalid? #(sut/divided 42 [1 2]) :f))
    (t/is (invalid? #(sut/divided nil [1 2]) :f))
    (t/is (invalid? #(sut/divided-all "f" [1 2]) :f))))

(t/deftest invalid-derivatives
  (t/testing "too few derivatives for the longest run"
    (t/is (invalid? #(sut/divided m/log [1 1] {:derivatives []}) :derivatives))
    (t/is (invalid? #(sut/divided m/log [1 1 1] {:derivatives [identity]}) :derivatives))
    (t/is (invalid? #(sut/divided-all m/log [1 1 1 1] {:derivatives [identity identity]}) :derivatives)))
  (t/testing "not a sequence"
    (t/is (invalid? #(sut/divided m/log [1 1] {:derivatives identity}) :derivatives))
    (t/is (invalid? #(sut/divided m/log [1 1] {:derivatives 3}) :derivatives)))
  (t/testing ":derivatives are not needed (and not checked) for distinct nodes"
    (t/is (= (sut/divided m/log [1 2 4]) (sut/divided m/log [1 2 4] {:derivatives []})))))

;; ------------------------------------------------------------------ shape, order, regressions

(t/deftest divided-all-shape
  (let [res (sut/divided-all (functions :exp) [0 1 2 3 4])]
    (t/testing "a list of lists, the highest level first"
      (t/is (list? res))
      (t/is (every? seq? res))
      (t/is (= [1 2 3 4 5] (map count res)))
      (t/is (every? #(every? double? %) res)))
    (t/testing "the last element is the function values, the first one is the value of divided"
      (t/is (= (map #(Math/exp %) [0 1 2 3 4]) (last res)))
      (t/is (= (sut/divided (functions :exp) [0 1 2 3 4]) (ffirst res))))
    (t/testing "levels are the divided differences of the neighbours"
      (t/is (= (map #(/ (- ^double %2 ^double %1) 1.0) (last res) (rest (last res))) (nth res 3))))))

(t/deftest divided-regressions
  (t/testing "distinct nodes, as before"
    (t/is (m/delta-eq 24.663 (m// (sut/divided m/log [30 20])) 1.0e-4))
    (t/is (m/delta-eq (double 10/3) (sut/divided [6 9 2 5] [0 1 2 3])))
    (t/is (m/zero? (sut/divided m/tan [-3/2 -3/4 0 3/4 3/2]))))
  (t/testing "the nodes need not be sorted and the divided difference is symmetric"
    (t/is (< (rel-err (sut/divided m/log [3 1 2 5]) (sut/divided m/log [1 2 3 5])) 1.0e-14))
    (t/is (< (rel-err (sut/divided m/log [5 3 2 1]) (sut/divided m/log [1 2 3 5])) 1.0e-14)))
  (t/testing "doubles in a lazy seq or an array of nodes"
    (t/is (= (sut/divided m/log [1.0 2.0 4.0]) (sut/divided m/log (map double [1 2 4]))))
    (t/is (= (sut/divided m/log [1.0 2.0 4.0]) (sut/divided m/log (double-array [1.0 2.0 4.0])))))
  (t/testing "close nodes are still a smooth function of the nodes"
    (t/is (< (rel-err (sut/divided (functions :exp) [0.0 1.0e-8]) 1.0) 1.0e-7)))
  (t/testing "negative and large nodes"
    (t/is (< (rel-err (sut/divided (fn ^double [^double x] (* x x)) [-1.0e4 3.0e4]) 2.0e4) 1.0e-12))
    (t/is (< (rel-err (sut/divided (fn ^double [^double x] (* x x)) [-1.0e8 1.0e8 3.0e8]) 1.0) 1.0e-6))))
