(ns fastmath.interpolation.kriging-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.distance :as dist]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.kernel.variogram :as vgram]
            [fastmath.interpolation.kriging :as sut]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; Reference values: utils/fastmath/dev/generate_kriging_reference.R (gstat::krige with an explicit vgm,
;; ordinary kriging z ~ 1 and universal kriging z ~ x + y, the measurement error `Err`).
(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "interpolation/kriging_reference.edn")))))

(def ^:private models
  {:spherical vgram/spherical
   :exponential vgram/exponential
   :gaussian vgram/gaussian})

(defn- property
  "Returns the :property of the invalid-input ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (= :kriging (:method d))) (:property d))))))

(defn- max-err
  ^double [f queries expected]
  (reduce (fn [^double acc [q e]]
            (m/max acc (m/abs (m/- (double (f q)) (double e)))))
          0.0 (map vector queries expected)))

(defn- max-diff
  ^double [f g queries]
  (max-err f queries (map g queries)))

(defn- linear-terms
  [x]
  (if (number? x) [1.0 x] (into [1.0] x)))

(defn- model-of
  [{:keys [model nugget psill range]}]
  ((models model) {:nugget nugget :psill psill :range range}))

(t/deftest gstat-ordinary-and-universal
  (doseq [{:keys [name drift err xs ys queries values] :as c} (:cases @reference)
          :let [f (sut/kriging xs ys (model-of c)
                               (cond-> {}
                                 (= drift :linear) (assoc :polynomial-terms linear-terms)
                                 (pos? (double err)) (assoc :error (m/sqrt (double err)))))]]
    ;; measured maximum 2.9e-15 over the 14 cases
    (t/is (< (max-err f queries values) 1e-11) (str name))))

(t/deftest gstat-error-per-node
  (doseq [{:keys [name errors xs ys queries values] :as c} (:error-sequence @reference)
          :let [f (sut/kriging xs ys (model-of c) {:error errors})]]
    ;; measured maximum 1.8e-15, the system is solved in R with the gstat semivariogram
    (t/is (< (max-err f queries values) 1e-11) (str name))
    (t/testing "a sequence of the length of the nodes, as a list, a vector and a double array"
      (t/is (= (mapv f queries) (mapv (sut/kriging xs ys (model-of c) {:error (map identity errors)}) queries)))
      (t/is (= (mapv f queries) (mapv (sut/kriging xs ys (model-of c) {:error (double-array errors)}) queries))))))

(t/deftest error-is-a-standard-deviation
  ;; the variance error^2 is subtracted from the diagonal: a scalar error 0.3 equals the sequence of 0.3
  ;; and gstat Err = 0.09 (the case d2-gaussian-ordinary-err uses Err 0.04 for the error 0.2)
  (let [{:keys [xs ys queries] :as c} (first (:cases @reference))
        mk (fn [e] (sut/kriging xs ys (model-of c) {:error e}))]
    (t/is (< (max-diff (mk 0.3) (mk (repeat (count xs) 0.3)) queries) 1e-12))
    (t/is (< (max-diff (mk 0.3) (mk -0.3) queries) 1e-12) "only the square matters")
    (t/is (= (mapv (mk 0.0) queries) (mapv (mk nil) queries) (mapv (mk []) queries) (mapv (mk '()) queries)) "no error")
    (t/is (> (max-diff (mk 0.3) (mk 0.0) queries) 1e-3) "the error changes the result")
    (t/testing "the interpolation of the nodes: exact without error, smoothing with the error"
      (t/is (< (max-err (mk 0.0) xs ys) 1e-9))
      (t/is (> (max-err (mk 0.3) xs ys) 1e-3)))))

(t/deftest seeded-default-variogram
  (let [{:keys [xs ys queries]} (first (filter #(= :d2-gaussian-ordinary (:name %)) (:cases @reference)))
        mk (fn [seed opts] (sut/kriging xs ys nil (assoc opts :rng (r/rng :isaac seed))))]
    (t/testing "the same seed gives bitwise identical interpolators"
      (t/is (= (mapv (mk 42 {}) queries) (mapv (mk 42 {}) queries)))
      (t/is (= (mapv (mk 1 {:N 20}) queries) (mapv (mk 1 {:N 20}) queries)))
      (t/is (= (mapv (mk 7 {:polynomial-terms linear-terms}) queries)
               (mapv (mk 7 {:polynomial-terms linear-terms}) queries))))
    (t/testing "the result is reasonable: a smooth fit of the data"
      (t/is (< (max-err (mk 42 {}) xs ys) 1e-9) "passes through the nodes")
      (t/is (every? #(m/valid-double? (double %)) (map (mk 42 {}) queries))))
    (t/testing "1-d: scalars and one element vectors give the same function"
      (let [{:keys [xs ys queries]} (first (filter #(= 1 (:dim %)) (:cases @reference)))
            f1 (sut/kriging xs ys nil {:rng (r/rng :isaac 7)})
            f2 (sut/kriging (mapv vector xs) ys nil {:rng (r/rng :isaac 7)})]
        (t/is (= (mapv f1 queries) (mapv f2 queries) (mapv f1 (map vector queries))))))))

(t/deftest default-variogram-forwards-options
  (let [{:keys [xs ys]} (first (filter #(= :d2-gaussian-ordinary (:name %)) (:cases @reference)))
        seen (atom {})
        rng (r/rng :isaac 5)
        dd (fn [a b] (m/* 100.0 (dist/euclidean a b)))]
    (with-redefs [vgram/empirical (fn [& args] (swap! seen assoc :empirical args) [{:h 1.0 :gamma 1.0 :n 1}])
                  vgram/fit (fn [& args] (swap! seen assoc :fit args) (vgram/gaussian {:nugget 0.0 :psill 1.0 :range 3.0}))]
      (sut/kriging xs ys nil {:rng rng :N 11 :distance dd})
      (let [[emp-xs emp-ys emp-opts] (:empirical @seen)
            [_ model fit-opts] (:fit @seen)]
        (t/is (= (count xs) (count emp-xs) (count emp-ys)))
        (t/is (identical? dd (:distance emp-opts)) ":distance goes to empirical")
        (t/is (= :gaussian model))
        (t/is (identical? rng (:rng fit-opts)) ":rng goes to fit")
        (t/is (= 11 (:N fit-opts)) ":N goes to fit"))
      (reset! seen {})
      (sut/kriging xs ys nil)
      (t/testing "no rng and no N when not given"
        (t/is (= {} (nth (:fit @seen) 2)))))))

(t/deftest scaled-distance-default-path
  (let [{:keys [xs ys queries]} (first (filter #(= :d2-gaussian-ordinary (:name %)) (:cases @reference)))
        dd (fn [a b] (m/* 100.0 (dist/euclidean a b)))
        via-default (sut/kriging xs ys nil {:distance dd :rng (r/rng :isaac 42)})
        explicit-vg (vgram/fit (vgram/empirical xs ys {:distance dd}) :gaussian {:rng (r/rng :isaac 42)})
        explicit (sut/kriging xs ys explicit-vg {:distance dd})
        unscaled (sut/kriging xs ys nil {:rng (r/rng :isaac 42)})]
    ;; the old implementation ignored :distance in the estimation of the variogram: 1.33 apart
    (t/is (< (max-diff via-default explicit queries) 1e-11))
    ;; scaling the distance does not change the function: the fitted range scales too (measured 1.7e-8, the fit is random)
    (t/is (< (max-diff via-default unscaled queries) 1e-6))))

(t/deftest query-forms
  (let [c (first (filter #(= :d2-gaussian-ordinary (:name %)) (:cases @reference)))
        {:keys [xs ys]} c
        vg (model-of c)
        f (sut/kriging xs ys vg)
        q [3.3 4.4]]
    (t/testing "long nodes, Vec2 nodes, Vec2 and list queries"
      (let [fl (sut/kriging (mapv #(mapv (comp double long) %) xs) ys vg)
            f2 (sut/kriging (mapv (fn [[x y]] (v/vec2 x y)) xs) ys vg)]
        (t/is (< (m/abs (m/- (double (f q)) (double (fl q)))) 1.0) "just a different data set (rounded nodes)")
        (t/is (= (f q) (f2 q) (f (v/vec2 3.3 4.4)) (f (list 3.3 4.4)) (f2 (v/vec2 3.3 4.4))))))
    (t/testing "integer nodes and queries equal the double ones"
      (let [fi (sut/kriging [[0 0] [3 0] [0 4] [5 5]] [1.0 2.0 3.0 4.0] vg)
            fd (sut/kriging [[0.0 0.0] [3.0 0.0] [0.0 4.0] [5.0 5.0]] [1.0 2.0 3.0 4.0] vg)]
        (t/is (= (fi [1 2]) (fd [1.0 2.0]) (fi [1.0 2.0]) (fd [1 2])))
        (t/is (m/delta-eq 3.0 (fi [0 4]) 1e-12))))
    (t/testing "a wrong dimension of the query"
      (t/is (= :dimension-mismatch (property #(f [1.0 2.0 3.0]))))
      (t/is (= :dimension-mismatch (property #(f [1.0]))))
      (t/is (= :dimension-mismatch (property #(f 1.0)))))
    (t/is (= :dimension-mismatch (property #(sut/kriging [[0.0 0.0] [1.0]] [1.0 2.0] vg))))))

(t/deftest nodes-and-edge-positions
  (let [vg (vgram/spherical {:nugget 0.0 :psill 1.0 :range 3.0})]
    (t/testing "one point is a constant function"
      (let [f (sut/kriging [2.0] [3.0] vg)]
        (t/is (= 3.0 (f 2.0)))
        (t/is (m/delta-eq 3.0 (f 100.0) 1e-12))
        (t/is (m/delta-eq 3.0 (f -1e300) 1e-12))))
    (t/testing "two points, the midpoint of a symmetric configuration is the mean"
      (let [f (sut/kriging [0.0 2.0] [1.0 5.0] vg)]
        (t/is (m/delta-eq 3.0 (f 1.0) 1e-12))
        (t/is (m/delta-eq 1.0 (f 0.0) 1e-12))
        (t/is (m/delta-eq 5.0 (f 2.0) 1e-12))))
    (t/testing "very near a node (no nugget) and far away"
      (let [f (sut/kriging [0.0 1.0 2.0 4.0] [1.0 5.0 2.0 3.0] vg)]
        (doseq [e [1e-100 1e-154 1e-155 1e-300]
                :let [e (double e)]]
          (t/is (m/delta-eq 5.0 (f (m/+ 1.0 e)) 1e-12) "the offset is below the resolution of the double")
          (t/is (m/delta-eq 1.0 (f e) 1e-12))
          (t/is (m/delta-eq 1.0 (f (m/- e)) 1e-12)))
        (t/is (m/valid-double? (f 1e100)))
        (t/is (m/valid-double? (f -1e300)))))
    (t/testing "a nugget makes the interpolant discontinuous at the nodes, the nodes are still exact"
      (let [f (sut/kriging [0.0 1.0 2.0] [1.0 5.0 2.0] (vgram/spherical {:nugget 0.5 :psill 1.0 :range 3.0}))]
        (t/is (m/delta-eq 5.0 (f 1.0) 1e-12))
        (t/is (> (m/abs (m/- (double (f 1.000001)) 5.0)) 1e-3))))
    (t/testing "NaN query"
      (t/is (m/nan? ((sut/kriging [0.0 1.0 2.0] [1.0 5.0 2.0] vg) ##NaN))))
    (t/testing "many points"
      (let [rnd (java.util.Random. 5)
            xs (vec (repeatedly 150 #(vector (m/* 10.0 (.nextDouble rnd)) (m/* 10.0 (.nextDouble rnd)))))
            ys (mapv (fn [[^double x ^double y]] (m/+ (m/sin x) (m/cos y))) xs)
            f (sut/kriging xs ys (vgram/gaussian {:nugget 0.0 :psill 1.0 :range 2.0}) {:error 1e-3})]
        (t/is (< (max-err f xs ys) 0.05) "smoothed by the tiny error")
        (t/is (m/delta-eq (m/+ (m/sin 5.0) (m/cos 5.0)) (f [5.0 5.0]) 0.1))))))

(t/deftest repeated-nodes
  (let [vg (vgram/spherical {:nugget 0.0 :psill 1.0 :range 3.0})
        out (java.io.StringWriter.)
        f (binding [*out* out]
            (sut/kriging [0.0 1.0 1.0 2.0] [0.0 1.0 3.0 4.0] vg))]
    (t/is (= "" (str out)) "no output on stdout")
    (t/is (m/delta-eq 2.0 (f 1.0) 1e-9) "the mean of the repeated values")
    (t/is (m/valid-double? (f 0.5)))))

(t/deftest scales
  (let [xs [[0.0 0.0] [1.0 0.0] [0.0 1.0] [1.0 1.0] [0.5 0.5] [0.2 0.8]]
        ys [1.0 2.0 3.0 5.0 4.0 2.5]
        qs [[0.3 0.3] [0.7 0.4] [0.1 0.9]]
        base (sut/kriging xs ys (vgram/gaussian {:nugget 0.0 :psill 1.0 :range 1.5}))]
    (doseq [c [1e-6 1e6 1e100]
            :let [c (double c)
                  sc (fn [p] (mapv (fn [^double x] (m/* c x)) p))
                  f (sut/kriging (map sc xs) ys (vgram/gaussian {:nugget 0.0 :psill 1.0 :range (m/* c 1.5)}))]]
      (t/is (< (max-err f (map sc qs) (map base qs)) 1e-7) (str "scale " c)))
    (t/testing "the values scale with the sill and shift with the mean"
      (let [f (sut/kriging xs (map #(m/+ 1000.0 (m/* 10.0 (double %))) ys) (vgram/gaussian {:nugget 0.0 :psill 100.0 :range 1.5}))]
        (t/is (< (max-err f qs (map #(m/+ 1000.0 (m/* 10.0 (double (base %)))) qs)) 1e-6))))))

(t/deftest drift-options
  (let [{:keys [xs ys queries] :as c} (first (filter #(= :d2-gaussian-linear (:name %)) (:cases @reference)))
        vg (model-of c)
        ordinary (sut/kriging xs ys vg)]
    (t/is (= (mapv ordinary queries) (mapv (sut/kriging xs ys vg {:polynomial-terms (constantly [1.0])}) queries))
          "the default drift is a constant")
    (t/is (= (mapv ordinary queries) (mapv (sut/kriging xs ys vg nil) queries)))
    (t/is (= (mapv ordinary queries) (mapv (sut/kriging xs ys vg {:error nil :distance nil :rng nil :N nil}) queries))
          "nil options are missing options")
    (t/testing "a linear drift reproduces a plane"
      (let [plane (fn [[x y]] (m/+ 1.0 (m/+ (m/* 2.0 (double x)) (m/* -3.0 (double y)))))
            f (sut/kriging xs (mapv plane xs) vg {:polynomial-terms linear-terms})]
        (t/is (< (max-err f queries (map plane queries)) 1e-9))))
    (t/testing "no drift: the plain variogram system"
      (t/is (every? #(m/valid-double? (double %)) (map (sut/kriging xs ys vg {:polynomial-terms nil}) queries))))))

(t/deftest invalid-input
  (let [vg (vgram/spherical {:nugget 0.0 :psill 1.0 :range 3.0})]
    (t/testing "lengths and minimal sizes"
      (t/is (= :length-mismatch (property #(sut/kriging [0.0 1.0] [1.0] vg))))
      (t/is (= :length-mismatch (property #(sut/kriging [0.0 1.0] [1.0 2.0 3.0] vg))))
      (t/is (= :length-mismatch (property #(sut/kriging [0.0 1.0] [1.0]))))
      (t/is (= :too-few-points (property #(sut/kriging [] [] vg))))
      (t/is (= :too-few-points (property #(sut/kriging nil nil vg))))
      (t/is (nil? (property #(sut/kriging [0.0] [1.0] vg))))
      (t/testing "the default variogram needs two points"
        (t/is (= :too-few-points (property #(sut/kriging [0.0] [1.0]))))
        (t/is (= :too-few-points (property #(sut/kriging [] []))))))
    (t/testing "the error: a wrong length gives a clear message"
      (let [e (try (sut/kriging [0.0 1.0 2.0] [1.0 2.0 3.0] vg {:error [0.1 0.1]}) (catch clojure.lang.ExceptionInfo e e))]
        (t/is (= :error (:property (ex-data e))))
        (t/is (= :kriging (:method (ex-data e))))
        (t/is (re-find #"expected 3, got 2" (ex-message e))))
      (t/is (= :error (property #(sut/kriging [0.0 1.0 2.0] [1.0 2.0 3.0] vg {:error [0.1 0.1 0.1 0.1]}))))
      (t/is (= :error (property #(sut/kriging [0.0 1.0 2.0] [1.0 2.0 3.0] vg {:error [0.1 ##NaN 0.1]}))))
      (t/is (= :error (property #(sut/kriging [0.0 1.0 2.0] [1.0 2.0 3.0] vg {:error ##NaN}))))
      (t/is (= :error (property #(sut/kriging [0.0 1.0 2.0] [1.0 2.0 3.0] vg {:error ##Inf}))))
      (t/is (= :error (property #(sut/kriging [0.0 1.0 2.0] [1.0 2.0 3.0] vg {:error "0.1"}))))
      (t/is (nil? (property #(sut/kriging [0.0 1.0 2.0] [1.0 2.0 3.0] vg {:error [0.1 0.1 0.1]}))))
      (t/is (nil? (property #(sut/kriging [0.0 1.0 2.0] [1.0 2.0 3.0] vg {:error []})))))))
