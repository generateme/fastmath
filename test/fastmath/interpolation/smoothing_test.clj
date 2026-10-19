(ns fastmath.interpolation.smoothing-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.interpolation.cubic :as cubic]
            [fastmath.interpolation.ssj :as sut])
  (:import [umontreal.ssj.functionfit SmoothingCubicSpline]))

(set! *warn-on-reflection* true)

;; scipy 1.18.1 make_smoothing_spline (lam = (1 - rho) / rho) and the numpy least squares line for
;; rho = 0, see utils/fastmath/dev/generate_bspline_bicubic_reference.py
(def ^:private reference
  (delay (:smoothing (edn/read-string (slurp (io/resource "interpolation/bspline_bicubic_reference.edn"))))))

(defn- property
  "Returns the :property of the ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (= :cubic-smoothing (:method d))) (:property d))))))

(defn- opts
  [{:keys [weights rho]}]
  (cond-> {:rho rho} weights (assoc :weights weights)))

(t/deftest scipy-make-smoothing-spline
  ;; inside the data range; measured max abs error 4.8e-12 (nonuni17, weights, rho 0.001, the worst
  ;; conditioned case), 6e-13 unweighted, 1.9e-15 for rho 1 and the rho 0 line
  (let [{:keys [datasets cases]} @reference]
    (t/is (= 33 (count cases)))
    (doseq [{:keys [dataset rho queries expected] :as c} cases
            :let [{:keys [xs ys]} (get datasets dataset)
                  f (sut/cubic-smoothing xs ys (opts c))
                  lo (double (first xs))
                  hi (double (peek xs))]]
      (t/testing (str dataset " rho " rho (when (:weights c) " weighted"))
        (t/is (every? true? (map (fn [q e] (or (< q lo) (> q hi) (m/delta-eq (f q) e 1.0e-11)))
                                 queries expected)))))))

(t/deftest extrapolates-linearly-with-the-end-tangent
  ;; scipy gives the end value and the end derivative; outside the spline is the tangent line
  (let [{:keys [datasets cases]} @reference]
    (doseq [{:keys [dataset rho d0 d1 v0 v1] :as c} cases
            :let [{:keys [xs ys]} (get datasets dataset)
                  f (sut/cubic-smoothing xs ys (opts c))
                  lo (double (first xs))
                  hi (double (peek xs))]
            d [1.0e-3 0.5 1.0 7.0 100.0]]
      (t/testing (str dataset " rho " rho " d " d)
        (t/is (m/delta-eq (- v0 (* d d0)) (f (- lo d)) 1.0e-9 1.0e-10))
        (t/is (m/delta-eq (+ v1 (* d d1)) (f (+ hi d)) 1.0e-9 1.0e-10))))))

(t/deftest equals-ssj-directly
  ;; valid input keeps the previous values: the validation does not change the spline
  (let [{:keys [datasets cases]} @reference]
    (doseq [{:keys [dataset rho queries weights] :as c} cases
            :let [{:keys [xs ys]} (get datasets dataset)
                  f (sut/cubic-smoothing xs ys (opts c))
                  ^SmoothingCubicSpline s (if weights
                                            (SmoothingCubicSpline. (double-array xs) (double-array ys) (double-array weights) (double rho))
                                            (SmoothingCubicSpline. (double-array xs) (double-array ys) (double rho)))]]
      (t/is (every? true? (map (fn [q] (m/delta-eq (f q) (.evaluate s (double q)) 1.0e-15)) queries))
            (str dataset " rho " rho)))))

(t/deftest rho-limits
  (let [{:keys [xs ys]} (get (:datasets @reference) "nonuni17")
        interp (cubic/cubic xs ys)]
    (t/testing "rho 1 is the natural interpolating cubic spline"
      (let [f (sut/cubic-smoothing xs ys {:rho 1.0})
            g (sut/cubic-smoothing xs ys)]
        (doseq [q (concat xs (map #(+ 0.3 (* 0.25 (double %))) (range 18)))]
          (t/is (m/delta-eq (interp q) (f q) 1.0e-12))
          (t/is (= (f q) (g q)) "rho 1 is the default"))))
    (t/testing "rho 0 is accepted and gives the linear regression line"
      (let [f (sut/cubic-smoothing xs ys {:rho 0.0})
            n (double (count xs))
            mx (/ (reduce + xs) n)
            my (/ (reduce + ys) n)
            slope (/ (reduce + (map #(* (- %1 mx) (- %2 my)) xs ys))
                     (reduce + (map #(m/sq (- % mx)) xs)))]
        (doseq [q [-5.0 0.0 1.5 3.3 9.0]]
          ;; measured error 2.4e-12 at -5, 1.3e-12 at 9 (the line is evaluated far from the data)
          (t/is (m/delta-eq (+ my (* slope (- q mx))) (f q) 1.0e-10) (str q)))))
    (t/testing "tiny positive rho approaches the line"
      (let [f (sut/cubic-smoothing xs ys {:rho 0.0})
            g (sut/cubic-smoothing xs ys {:rho 1.0e-12})]
        (t/is (m/delta-eq (f 2.2) (g 2.2) 1.0e-6))))
    (t/testing "nil rho and nil weights are missing options"
      (t/is (= ((sut/cubic-smoothing xs ys) 1.7) ((sut/cubic-smoothing xs ys {:rho nil :weights nil}) 1.7))))))

(t/deftest minimal-number-of-points
  (t/testing "two points: the line through both, for every rho; linear outside"
    (doseq [rho [0.0 0.3 1.0]
            :let [f (sut/cubic-smoothing [0.0 1.0] [1.0 3.0] {:rho rho})]]
      (doseq [[x e] [[0.0 1.0] [0.5 2.0] [1.0 3.0] [3.0 7.0] [-2.0 -3.0]]]
        (t/is (m/delta-eq e (f x) 1.0e-12) (str "rho " rho " x " x)))))
  (t/testing "three points rho 0: least squares line; rho 1: natural spline"
    (t/is (m/delta-eq 1.75 ((sut/cubic-smoothing [0.0 1.0 2.0] [1.0 3.0 2.0] {:rho 0.0}) 0.5) 1.0e-12))
    (t/is (m/delta-eq 2.28125 ((sut/cubic-smoothing [0.0 1.0 2.0] [1.0 3.0 2.0]) 0.5) 1.0e-12))))

(t/deftest weights
  (let [{:keys [xs ys]} (get (:datasets @reference) "nonuni17")
        f (sut/cubic-smoothing xs ys {:rho 0.5})
        g (sut/cubic-smoothing xs ys {:rho 0.5 :weights (repeat (count xs) 1.0)})]
    (t/is (m/delta-eq (f 2.1) (g 2.1) 1.0e-12) "unit weights are the default")
    ;; a large weight pulls the spline to that data point
    (let [heavy (sut/cubic-smoothing xs ys {:rho 0.5 :weights (assoc (vec (repeat (count xs) 1.0)) 5 1.0e6)})
          res (fn [g] (m/abs (- (double (g (xs 5))) (double (ys 5)))))]
      (t/is (< (res heavy) (* 0.1 (res f))))
      (t/is (> (res f) 1.0e-2)))))

(t/deftest nan-and-scales
  (let [xs [0.0 1.0 2.0 3.0 4.0]
        ys [1.0 3.0 2.0 4.0 1.0]]
    (t/is (m/nan? ((sut/cubic-smoothing xs [1.0 ##NaN 2.0 4.0 1.0] {:rho 0.5}) 0.5)) "NaN in ys")
    (t/is (m/nan? ((sut/cubic-smoothing xs ys {:rho 0.5}) ##NaN)) "NaN query")
    (doseq [c [1.0e-9 1.0e9]
            q [-1.0 0.5 2.5 5.0]]
      (t/is (m/delta-eq (* c (double ((sut/cubic-smoothing xs ys {:rho 0.5}) q)))
                        ((sut/cubic-smoothing xs (map #(* c (double %)) ys) {:rho 0.5}) q)
                        1.0e-12 1.0e-12)
            (str "ys scale " c))
      (t/is (m/delta-eq ((sut/cubic-smoothing xs ys {:rho 1.0}) q)
                        ((sut/cubic-smoothing (map #(* c (double %)) xs) ys {:rho 1.0}) (* c (double q)))
                        1.0e-9 1.0e-9)
            (str "xs scale " c)))))

(t/deftest invalid-input
  (let [xs [0.0 1.0 2.0 3.0]
        ys [1.0 3.0 2.0 4.0]]
    (t/is (nil? (property #(sut/cubic-smoothing xs ys))))
    (t/testing "points"
      (t/is (= :too-few-points (property #(sut/cubic-smoothing [] []))))
      (t/is (= :too-few-points (property #(sut/cubic-smoothing [1.0] [2.0]))))
      (t/is (nil? (property #(sut/cubic-smoothing [1.0 2.0] [2.0 3.0])))))
    (t/testing "sorted, finite xs"
      (t/is (= :not-increasing (property #(sut/cubic-smoothing [0.0 2.0 1.0 3.0] ys))))
      (t/is (= :not-increasing (property #(sut/cubic-smoothing [0.0 1.0 1.0 3.0] ys))))
      (t/is (= :not-increasing (property #(sut/cubic-smoothing [3.0 2.0 1.0 0.0] ys))))
      (t/is (= :not-finite (property #(sut/cubic-smoothing [0.0 ##NaN 2.0 3.0] ys))))
      (t/is (= :not-finite (property #(sut/cubic-smoothing [0.0 1.0 2.0 ##Inf] ys)))))
    (t/testing "lengths"
      (t/is (= :length-mismatch (property #(sut/cubic-smoothing xs [1.0 2.0 3.0]))))
      (t/is (= :length-mismatch (property #(sut/cubic-smoothing xs (conj ys 1.0)))))
      (t/is (= :length-mismatch (property #(sut/cubic-smoothing xs ys {:weights [1.0 1.0 1.0]}))))
      (t/is (= :length-mismatch (property #(sut/cubic-smoothing xs ys {:weights [1.0 1.0 1.0 1.0 1.0]})))))
    (t/testing "rho in [0, 1]"
      (doseq [rho [-0.1 -1.0e-12 1.0000001 2.0 ##NaN ##Inf ##-Inf]]
        (t/is (= :rho (property #(sut/cubic-smoothing xs ys {:rho rho}))) (str rho)))
      (doseq [rho [0.0 1.0e-300 0.5 0.9999999999 1.0 -0.0]]
        (t/is (nil? (property #(sut/cubic-smoothing xs ys {:rho rho}))) (str rho))))
    (t/testing "weights positive and finite"
      (doseq [w [0.0 -1.0 ##NaN ##Inf]]
        (t/is (= :weights (property #(sut/cubic-smoothing xs ys {:weights [1.0 w 1.0 1.0]}))) (str w))))
    (t/testing "ex-data"
      (let [d (try (sut/cubic-smoothing xs ys {:rho 1.5}) nil
                   (catch clojure.lang.ExceptionInfo e (ex-data e)))]
        (t/is (= {:type :invalid-input :method :cubic-smoothing :property :rho :option :rho :value 1.5 :expected "[0, 1]"} d))))))
