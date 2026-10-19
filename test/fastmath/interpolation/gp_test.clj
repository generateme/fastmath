(ns fastmath.interpolation.gp-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.kernel :as k]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.interpolation.gp :as sut]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; Reference values: utils/fastmath/dev/generate_scattered_reference.py (scikit-learn GaussianProcessRegressor,
;; ConstantKernel(kscale, fixed) * RBF(sigma, fixed), alpha = noise, optimizer = None; with :normalize? the
;; values are normalised with the sample standard deviation first).
(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "interpolation/scattered_reference.edn")))))

(defn- property
  "Returns the :property of the invalid-input ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (= :gp (:method d))) (:property d))))))

(defn- gauss
  "The gaussian covariance kernel with the given sigma."
  [^double sigma]
  (k/kernel :gaussian {:sigma sigma}))

(defn- points
  "The reference points as vectors."
  [xs]
  (mapv #(if (number? %) [%] (vec %)) xs))

(defn- model
  [{:keys [sigma kscale noise normalize? xs ys]}]
  (sut/gaussian-process xs ys {:kernel (gauss sigma) :kscale kscale :noise noise :normalize? normalize?}))

(defn- max-rel-err
  "Maximum of the errors relative to max(1, abs(expected))."
  ^double [actual expected]
  (reduce (fn [^double acc [a e]]
            (let [e (double e)]
              (m/max acc (m// (m/abs (m/- (double a) e)) (m/max 1.0 (m/abs e))))))
          0.0 (map vector actual expected)))

(t/deftest sklearn-mean-std
  (doseq [{:keys [name queries mean std] :as c} (:gp @reference)
          :let [g (model c)
                all (sut/predict-all g queries true)
                single (map #(sut/predict g % true) queries)]]
    ;; measured maximum relative error: mean 5.8e-13, std 7.7e-11 (uniform-400, conditioning 1e8)
    (t/is (< (max-rel-err (map first all) mean) 1e-9) (str name " predict-all mean"))
    (t/is (< (max-rel-err (map second all) std) 1e-9) (str name " predict-all std"))
    (t/is (< (max-rel-err (map first single) mean) 1e-9) (str name " predict mean"))
    (t/is (< (max-rel-err (map second single) std) 1e-9) (str name " predict std"))
    (t/is (< (max-rel-err (sut/predict-all g queries) mean) 1e-9) (str name " predict-all without std"))
    (t/is (= (map first all) (sut/predict-all g queries)) (str name " the same means with and without std"))))

(t/deftest sklearn-covariance
  (doseq [{:keys [name queries cov] :as c} (:gp @reference)
          :let [g (model c)
                kv (#'sut/cross-cov g (points queries))
                ours (#'sut/posterior-cov g (points queries) kv)]]
    ;; measured maximum relative error 6.4e-12 (normalized-2d)
    (t/is (< (max-rel-err (mapcat (fn [i] (.getRow ^org.apache.commons.math3.linear.RealMatrix ours (int i))) (range (count queries))) (apply concat cov)) 1e-9)
          (str name " covariance"))))

(t/deftest sklearn-log-marginal-likelihood
  (doseq [{:keys [name L] :as c} (:gp @reference)
          :let [g (model c)
                lazy (sut/gaussian-process (:xs c) (:ys c) {:kernel (gauss (:sigma c)) :kscale (:kscale c) :noise (:noise c)
                                                           :normalize? (:normalize? c) :L? false})]]
    ;; measured 7e-8 for the 400 points (L 3224.54, the condition number of the covariance is 1e8), 2e-12 and less for the others
    (t/is (< (m/abs (m/- (sut/L g) (double L))) 1e-6) (str name))
    (t/is (m/nan? (:L lazy)) "L? false stores NaN")
    (t/is (= (sut/L g) (sut/L lazy)) "and computes it on demand"))
  (t/testing "the determinant of the 400-point covariance underflows, the stable form is finite"
    (let [c (first (filter #(= :uniform-400 (:name %)) (:gp @reference)))]
      (t/is (m/valid-double? (sut/L (model c))))
      (t/is (< 3000.0 (sut/L (model c)) 3500.0))))
  (t/testing "40 points, sigma 0.01, kscale 1e20: the determinant overflows to infinity"
    (let [c (first (filter #(= :narrow-40 (:name %)) (:gp @reference)))
          l (sut/L (model c))]
      (t/is (m/valid-double? l))
      (t/is (< -1000.0 l -900.0))))
  (t/testing "an exactly solvable case: one point"
    ;; y = 2, K = 1 + noise: L = -1/2 (y^2 / K + log K + log 2 pi)
    (let [g (sut/gaussian-process [0.0] [2.0] {:noise 0.5})]
      (t/is (m/delta-eq (* -0.5 (+ (/ 4.0 1.5) (m/log 1.5) (m/log m/TWO_PI))) (sut/L g) 1e-12)))))

(t/deftest predict-all-has-no-nan
  (t/testing "queries at the training points, kscale 1e6 and noise 1e-10: the variance is a tiny negative number"
    (let [rnd (java.util.Random. 11)]
      (dotimes [_ 12]
        (let [xs (vec (repeatedly 10 #(.nextDouble rnd)))
              g (sut/gaussian-process xs (map #(m/sin (m/* 5.0 (double %))) xs)
                                      {:kernel (gauss 0.3) :kscale 1e6 :noise 1e-10})
              res (sut/predict-all g xs true)]
          (t/is (every? (fn [[mu sd]] (and (m/valid-double? mu) (m/valid-double? sd) (m/>= (double sd) 0.0))) res))
          (t/is (every? #(m/valid-double? (double (second (sut/predict g % true)))) xs))))))
  (t/testing "random configurations: sigma, kscale, noise, normalisation, queries at the nodes and on a grid"
    (let [rnd (java.util.Random. 7)
          pick (fn [coll] (nth coll (.nextInt rnd (count coll))))
          results (doall
                   (for [_ (range 200)
                         :let [n (m/+ 5 (.nextInt rnd 40))
                               xs (vec (repeatedly n #(.nextDouble rnd)))
                               sigma (pick [0.01 0.1 0.3 1.0])
                               kscale (pick [1.0 1e6 1e-3])
                               noise (pick [1e-8 1e-10 1e-12])
                               normalize? (.nextBoolean rnd)
                               g (try (sut/gaussian-process xs (map #(m/sin (m/* 5.0 (double %))) xs)
                                                            {:kernel (gauss sigma) :kscale kscale :noise noise
                                                             :normalize? normalize?})
                                      (catch clojure.lang.ExceptionInfo _ nil))]
                         :when g]
                     (sut/predict-all g (concat xs (range 0.0 1.0 0.05)) true)))]
      ;; measured: 183 of the 200 constructed (the others are numerically singular: kscale 1e6 and noise 1e-10 or less), no NaN
      (t/is (> (count results) 150))
      (t/is (every? (fn [res] (every? (fn [[mu sd]] (and (m/valid-double? mu) (m/valid-double? sd))) res)) results)))))

(t/deftest normalize-options
  (let [xs [0.0 0.5 1.0 1.5 2.0]]
    (t/testing "constant values with :normalize? true: the standard deviation is 1, predictions equal the constant"
      (doseq [c [3.0 0.0 -7.5 1e10]
              :let [c (double c)
                    g (sut/gaussian-process xs (repeat 5 c) {:normalize? true :kernel (gauss 0.7)})
                    plain (sut/gaussian-process xs (repeat 5 c) {:kernel (gauss 0.7)})]]
        (t/is (m/delta-eq c (sut/predict g 0.7) 1e-9 ) (str "constant " c))
        (t/is (m/delta-eq c (sut/predict g 100.0) 1e-9))
        (t/is (every? #(m/delta-eq c % 1e-9) (sut/predict-all g [0.1 0.9 3.0])))
        (t/is (every? #(and (m/valid-double? (first %)) (m/valid-double? (second %))) (sut/predict-all g [0.1 0.9 3.0] true)))
        (t/is (m/valid-double? (sut/L g)))
        (t/is (= 1.0 (:ystddev g)))
        (t/is (m/delta-eq c (:ymean g) (m/max 1e-9 (m/* 1e-15 (m/abs (double c))))))
        (t/is (m/valid-double? (second (sut/predict g 0.7 true))))
        (t/testing "the standard deviation comes from the kernel only: the same as for the not normalised constant"
          (t/is (m/delta-eq (second (sut/predict plain 0.7 true)) (second (sut/predict g 0.7 true)) 1e-9)))))
    (t/testing "a single point is never normalised"
      (let [g (sut/gaussian-process [1.0] [5.0] {:normalize? true})]
        (t/is (= 0.0 (:ymean g)))
        (t/is (= 1.0 (:ystddev g)))
        (t/is (m/delta-eq 5.0 (sut/predict g 1.0) 1e-6))))
    (t/testing "the sample standard deviation (n - 1) is used"
      (let [ys [1.0 2.0 4.0 7.0 11.0]
            g (sut/gaussian-process xs ys {:normalize? true})]
        (t/is (m/delta-eq 5.0 (:ymean g) 1e-12))
        (t/is (m/delta-eq (m/sqrt 16.5) (:ystddev g) 1e-12) "variance (16 + 9 + 1 + 4 + 36) / 4 = 16.5")))
    (t/testing "two points with different values, and nearly constant values"
      (let [g (sut/gaussian-process [0.0 1.0] [1.0 3.0] {:normalize? true :kernel (gauss 0.5)})]
        (t/is (m/delta-eq 1.0 (sut/predict g 0.0) 1e-6))
        (t/is (m/delta-eq 3.0 (sut/predict g 1.0) 1e-6)))
      (let [g (sut/gaussian-process xs [1.0 1.0 1.0 1.0 1.0000001] {:normalize? true :kernel (gauss 0.7)})]
        (t/is (every? m/valid-double? (sut/predict-all g [0.2 0.8 5.0])))))))

(t/deftest default-kernel-and-options
  (t/testing "the default kernel is the gaussian with sigma 1"
    ;; one point at 0 with y = 1: the mean at 1 is k(1, 0) / (1 + noise) = exp(-1/2) / (1 + 1e-8)
    (let [g (sut/gaussian-process [0.0] [1.0])]
      (t/is (m/delta-eq (m// (m/exp -0.5) (m/+ 1.0 1e-8)) (sut/predict g 1.0) 1e-14))
      (t/is (= (sut/predict g 0.3) (sut/predict (sut/gaussian-process [0.0] [1.0] {:kernel (k/kernel :gaussian)}) 0.3)))
      (t/is (= (sut/predict g 0.3) (sut/predict (sut/gaussian-process [0.0] [1.0] {:kernel (k/kernel :gaussian {:sigma 1.0})}) 0.3)))
      (t/is (not= (sut/predict g 0.3) (sut/predict (sut/gaussian-process [0.0] [1.0] {:kernel (gauss 0.1)}) 0.3)))))
  (t/testing "kscale multiplies the kernel: the mean does not change without noise, the standard deviation scales"
    (let [mk (fn [^double ks] (sut/gaussian-process [0.0 1.0 2.5] [1.0 0.0 2.0] {:kernel (gauss 0.8) :kscale ks :noise 1e-12}))
          [mu1 sd1] (sut/predict (mk 1.0) 1.7 true)
          [mu9 sd9] (sut/predict (mk 9.0) 1.7 true)]
      (t/is (m/delta-eq (double mu1) (double mu9) 1e-7))
      (t/is (m/delta-eq (m/* 3.0 (double sd1)) (double sd9) 1e-7))))
  (t/testing "noise is a variance added to the diagonal: a larger noise smooths the data"
    (let [mk (fn [nz] (sut/gaussian-process [0.0 1.0 2.0] [1.0 0.0 1.0] {:noise nz :kernel (gauss 0.8)}))]
      (t/is (m/delta-eq 0.0 (sut/predict (mk 1e-8) 1.0) 1e-6))
      (t/is (> (sut/predict (mk 1.0) 1.0) 0.1))
      (t/is (> (sut/predict (mk 10.0) 1.0) 0.0))
      (t/is (< (sut/predict (mk 10.0) 1.0) (sut/predict (mk 1.0) 1.0)) "the shrinkage towards the zero mean")))
  (t/testing "nil options are missing options"
    (let [g1 (sut/gaussian-process [0.0 1.0] [1.0 2.0])
          g2 (sut/gaussian-process [0.0 1.0] [1.0 2.0] {:kscale nil :kernel nil :noise nil :normalize? nil :L? nil})]
      (t/is (= (sut/predict g1 0.4) (sut/predict g2 0.4)))
      (t/is (= (sut/L g1) (sut/L g2))))))

(t/deftest gp-function
  (let [xs [0.0 1.0 2.0 3.0]
        ys [1.0 3.0 2.0 5.0]
        f (sut/gp xs ys (gauss 0.8))
        g (sut/gaussian-process xs ys {:kernel (gauss 0.8)})]
    (t/is (= (sut/predict g 1.5) (f 1.5)))
    (t/is (m/delta-eq 3.0 (f 1.0) 1e-6))
    (t/is (= (f 1.5) ((sut/gp xs ys (gauss 0.8) nil) 1.5)))
    (t/is (= (f 1.5) ((sut/gp xs ys (gauss 0.8) {:noise 1e-8}) 1.5)))
    (t/is (= (f 1.5) ((sut/gp xs ys (gauss 0.8) {:L? true}) 1.5)) "gp does not compute L")
    (t/is (= (sut/predict (sut/gaussian-process xs ys) 1.5) ((sut/gp xs ys) 1.5)) "the default kernel")
    (t/is (= (sut/predict (sut/gaussian-process xs ys) 1.5) ((sut/gp xs ys nil) 1.5)))
    (t/is (= (sut/predict (sut/gaussian-process xs ys) 1.5) ((sut/gp xs ys nil nil) 1.5)))))

(t/deftest record-as-function
  (let [xs [0.0 1.0 2.0 3.0]
        g (sut/gaussian-process xs [1.0 3.0 2.0 5.0] {:kernel (gauss 0.8)})]
    (t/is (= (sut/predict g 1.5) (g 1.5)))
    (t/is (= (sut/predict g 1.5) (g [1.5])))
    (t/is (= (sut/predict g 1.5 true) (g 1.5 true)))
    (t/is (= (sut/predict-all g [0.5 1.5]) (g [[0.5] [1.5]])))
    (t/is (= (sut/predict-all g [0.5 1.5] true) (g [[0.5] [1.5]] true)))))

(t/deftest query-forms-and-dimensions
  (let [xs [[0.0 0.0] [1.0 0.0] [0.0 1.0] [1.0 1.0] [0.5 0.5]]
        ys [1.0 2.0 3.0 5.0 4.0]
        g (sut/gaussian-process xs ys {:kernel (gauss 0.8)})]
    (t/testing "the query as a vector, a list, a Vec2, longs, a sequence of them"
      (let [q (sut/predict g [0.3 0.7])]
        (t/is (= q (sut/predict g '(0.3 0.7)) (sut/predict g (v/vec2 0.3 0.7)) (g (v/vec2 0.3 0.7))))
        (t/is (= (sut/predict g [1 0]) (sut/predict g [1.0 0.0])))
        (t/is (= [q] (sut/predict-all g [[0.3 0.7]]) (sut/predict-all g [(v/vec2 0.3 0.7)]) (g [(v/vec2 0.3 0.7)])))))
    (t/testing "the training points as Vec2, lists and longs"
      (let [q [0.3 0.7]
            base (sut/predict g q)]
        (t/is (= base (sut/predict (sut/gaussian-process (map (fn [[x y]] (v/vec2 x y)) xs) ys {:kernel (gauss 0.8)}) q)))
        (t/is (= base (sut/predict (sut/gaussian-process (map seq xs) ys {:kernel (gauss 0.8)}) q)))
        (t/is (= (sut/predict (sut/gaussian-process [0 1 2] [1.0 2.0 0.5]) 0.5)
                 (sut/predict (sut/gaussian-process [0.0 1.0 2.0] [1.0 2.0 0.5]) 0.5)
                 (sut/predict (sut/gaussian-process [[0] [1] [2]] [1 2 0.5]) [0.5])))))
    (t/testing "the query dimension has to equal the training dimension"
      (t/is (= :dimension-mismatch (property #(sut/predict g [0.5]))))
      (t/is (= :dimension-mismatch (property #(sut/predict g 0.5))))
      (t/is (= :dimension-mismatch (property #(sut/predict g [0.5 0.5 0.5]))))
      (t/is (= :dimension-mismatch (property #(sut/predict g [0.5] true))))
      (t/is (= :dimension-mismatch (property #(sut/predict-all g [[0.5 0.5] [0.5]]))))
      (t/is (= :dimension-mismatch (property #(sut/predict-all g [[0.5 0.5 0.5]]))))
      (t/is (= :dimension-mismatch (property #(sut/predict-all g [0.5 0.6]))))
      (t/is (= :dimension-mismatch (property #(g [0.5]))))
      (t/is (= :dimension-mismatch (property #(sut/posterior-samples g [[0.5]]))))
      (t/is (= :dimension-mismatch (property #(sut/prior-samples g [[0.5 0.5 0.5]]))))
      (t/is (nil? (property #(sut/predict g [0.5 0.5])))))
    (t/testing "one dimension"
      (let [g1 (sut/gaussian-process [0.0 1.0 2.0] [1.0 2.0 0.5])]
        (t/is (= :dimension-mismatch (property #(sut/predict g1 [0.5 0.5]))))
        (t/is (= :dimension-mismatch (property #(sut/predict-all g1 [[0.5 0.5]]))))
        (t/is (= (sut/predict g1 0.5) (sut/predict g1 [0.5])))
        (t/is (= (sut/predict-all g1 [0.5 0.7]) (sut/predict-all g1 [[0.5] [0.7]])))))))

(t/deftest invalid-input
  (t/testing "ragged training points"
    (t/is (= :dimension-mismatch (property #(sut/gaussian-process [[0.0 0.0] [1.0]] [1.0 2.0]))))
    (t/is (= :dimension-mismatch (property #(sut/gaussian-process [[0.0] [1.0 1.0] [2.0]] [1.0 2.0 3.0]))))
    (t/is (= :dimension-mismatch (property #(sut/gaussian-process [0.0 [1.0 1.0]] [1.0 2.0])))))
  (t/testing "lengths and minimal size"
    (t/is (= :length-mismatch (property #(sut/gaussian-process [0.0 1.0] [1.0]))))
    (t/is (= :length-mismatch (property #(sut/gaussian-process [0.0 1.0] [1.0 2.0 3.0]))))
    (t/is (= :too-few-points (property #(sut/gaussian-process [] []))))
    (t/is (= :too-few-points (property #(sut/gaussian-process nil nil))))
    (t/is (nil? (property #(sut/gaussian-process [0.0] [1.0])))))
  (t/testing "kscale and noise: the first invalid values and the boundary"
    (doseq [ks [0.0 -1.0 ##NaN ##Inf]]
      (t/is (= :kscale (property #(sut/gaussian-process [0.0 1.0] [1.0 2.0] {:kscale ks}))) (str "kscale " ks)))
    (doseq [nz [-1e-300 -1.0 ##NaN ##Inf]]
      (t/is (= :noise (property #(sut/gaussian-process [0.0 1.0] [1.0 2.0] {:noise nz}))) (str "noise " nz)))
    (t/is (nil? (property #(sut/gaussian-process [0.0 1.0] [1.0 2.0] {:kscale 1e-300}))))
    (t/is (nil? (property #(sut/gaussian-process [0.0 1.0] [1.0 2.0] {:noise 0.0})))))
  (t/testing "a covariance matrix which is not positive definite gives a readable exception"
    (t/is (= :not-positive-definite (property #(sut/gaussian-process [0.0 0.0] [1.0 2.0] {:noise 0.0})))))
  (t/testing "repeated points with a tiny noise (the former absolute threshold 1e-10 refused 1e-11)"
    (let [g (sut/gaussian-process [0.0 0.0 1.0] [1.0 1.0 2.0] {:noise 1e-12 :kernel (gauss 0.5)})]
      (t/is (m/delta-eq 1.0 (sut/predict g 0.0) 1e-4))
      (t/is (m/valid-double? (sut/L g))))))

;; the samples

(defn- mean
  ^double [xs]
  (m// (double (reduce + xs)) (double (count xs))))

(defn- sd
  ^double [xs]
  (let [mu (mean xs)]
    (m/sqrt (m// (double (reduce (fn [^double acc ^double x] (m/+ acc (m/sq (m/- x mu)))) 0.0 xs))
                 (m/dec (count xs))))))

(t/deftest posterior-samples-statistics
  ;; 5000 seeded samples at a point: the mean within 4 standard errors of the mean of predict and the standard
  ;; deviation within 4 standard errors of the one of predict (the standard error of a sample sd is sd / sqrt(2 (N - 1)))
  (doseq [[nm ys normalize?] [["normalised, ystddev about 20" [190.0 215.0 240.0 205.0 180.0 230.0 250.0 195.0] true]
                              ["normalised, ystddev about 0.3" [0.1 0.5 0.9 0.4 0.2 0.8 1.0 0.35] true]
                              ["not normalised" [1.0 3.0 2.0 5.0 4.0 3.5 2.5 0.5] false]]
          :let [xs [0.0 0.4 0.9 1.3 1.8 2.2 2.7 3.1]
                g (sut/gaussian-process xs ys {:kernel (gauss 0.9) :kscale 2.0 :noise 1e-6 :normalize? normalize?})
                n 5000]
          x [0.65 2.45 5.0]
          :let [[mu s] (sut/predict g x true)
                mu (double mu)
                s (double s)
                rng (r/rng :isaac 1234)
                samples (vec (repeatedly n #(first (sut/posterior-samples g [x] false rng))))
                se-mean (m// s (m/sqrt (double n)))
                se-sd (m// s (m/sqrt (m/* 2.0 (m/dec n))))]]
    (t/is (< (m/abs (m/- (mean samples) mu)) (m/* 4.0 se-mean)) (str nm " mean at " x))
    (t/is (< (m/abs (m/- (sd samples) s)) (m/* 4.0 se-sd)) (str nm " sd at " x))))

(t/deftest posterior-samples-covariance
  (let [xs [0.0 0.4 0.9 1.3 1.8 2.2]
        ys [190.0 215.0 240.0 205.0 180.0 230.0]
        g (sut/gaussian-process xs ys {:kernel (gauss 0.9) :noise 1e-6 :normalize? true})
        qs [0.65 0.85 2.0]
        rng (r/rng :isaac 99)
        n 6000
        samples (vec (repeatedly n #(vec (sut/posterior-samples g qs false rng))))
        cov (#'sut/posterior-cov g (points qs) (#'sut/cross-cov g (points qs)))]
    (t/testing "the sample covariance equals the posterior covariance"
      (doseq [i (range 3)
              j (range i 3)
              :let [ci (map #(nth % i) samples)
                    cj (map #(nth % j) samples)
                    mi (mean ci)
                    mj (mean cj)
                    sample-cov (m// (double (reduce + (map (fn [^double a ^double b] (m/* (m/- a mi) (m/- b mj))) ci cj)))
                                    (double (dec n)))
                    expected (.getEntry ^org.apache.commons.math3.linear.RealMatrix cov (int i) (int j))
                    ;; the standard error of a sample covariance: sqrt((s_ii s_jj + s_ij^2) / N)
                    se (m/sqrt (m// (m/+ (m/* (.getEntry ^org.apache.commons.math3.linear.RealMatrix cov (int i) (int i))
                                              (.getEntry ^org.apache.commons.math3.linear.RealMatrix cov (int j) (int j)))
                                         (m/sq expected))
                                    (double n)))]]
        (t/is (< (m/abs (m/- sample-cov expected)) (m/* 4.0 se)) (str "entry " i j))))))

(t/deftest posterior-samples-options
  (let [xs [0.0 0.4 0.9 1.3 1.8]
        g (sut/gaussian-process xs [1.0 3.0 2.0 5.0 4.0] {:kernel (gauss 0.9)})
        qs [0.2 0.7 1.1]]
    (t/testing "a seeded rng gives the same samples, another seed other samples"
      (t/is (= (sut/posterior-samples g qs false (r/rng :isaac 5)) (sut/posterior-samples g qs false (r/rng :isaac 5))))
      (t/is (not= (sut/posterior-samples g qs false (r/rng :isaac 5)) (sut/posterior-samples g qs false (r/rng :isaac 6))))
      (t/is (= (sut/posterior-samples g qs true (r/rng :isaac 5)) (sut/posterior-samples g qs true (r/rng :isaac 5)))))
    (t/testing "the default arities use the default generator"
      (let [draw (fn [f] (with-redefs [r/default-rng (r/rng :isaac 77)] (f)))
            a (draw #(sut/posterior-samples g qs))]
        (t/is (= a (draw #(sut/posterior-samples g qs false)) (draw #(sut/posterior-samples g qs false nil))))
        (t/is (not= a (draw #(sut/posterior-samples g qs false (r/rng :isaac 78))))))
      (t/is (= 3 (count (sut/posterior-samples g qs)))))
    (t/testing "with the standard deviation: pairs of the sample and the posterior sd of predict"
      (let [res (sut/posterior-samples g qs true (r/rng :isaac 5))
            plain (sut/posterior-samples g qs false (r/rng :isaac 5))]
        (t/is (= plain (map first res)))
        (t/is (< (max-rel-err (map second res) (map #(second (sut/predict g % true)) qs)) 1e-9))))
    (t/testing "an invalid generator"
      (t/is (thrown? clojure.lang.ExceptionInfo (sut/posterior-samples g qs false :isaac)))
      (t/is (thrown? clojure.lang.ExceptionInfo (sut/prior-samples g qs 42))))))

(t/deftest posterior-samples-dense-grids
  (let [xs (vec (range 0.0 5.0 0.5))
        ys (mapv #(m/sin (double %)) xs)]
    (doseq [[nm kernel] [["gaussian 0.3" (gauss 0.3)] ["gaussian 0.1" (gauss 0.1)] ["matern 5/2" (k/kernel :matern-52)]]
            normalize? [false true]
            :let [g (sut/gaussian-process xs ys {:kernel kernel :normalize? normalize? :noise 1e-8})]]
      ;; the former implementation threw NonSymmetricMatrixException from 10 to 30 grid points
      (doseq [n [10 30 100 200 400]
              :let [grid (vec (m/slice-range 0.0 5.0 (long n)))
                    s (sut/posterior-samples g grid false (r/rng :isaac n))]]
        (t/is (= n (count s)) (str nm " grid " n))
        (t/is (every? #(m/valid-double? (double %)) s) (str nm " grid " n))))
    (t/testing "at the training points the samples equal the data (noise 1e-8) also for the noise 1e-12"
      (doseq [nz [1e-8 1e-12]
              :let [g (sut/gaussian-process xs ys {:kernel (gauss 0.5) :noise nz :normalize? true})
                    s (sut/posterior-samples g xs false (r/rng :isaac 3))]]
        (t/is (< (max-rel-err s ys) 1e-3) (str "noise " nz))))
    (t/testing "a single point, a repeated point and a far away point"
      (let [g (sut/gaussian-process xs ys {:kernel (gauss 0.5)})]
        (t/is (= 1 (count (sut/posterior-samples g [2.2]))))
        (t/is (= 3 (count (sut/posterior-samples g [2.2 2.2 2.2]))))
        (t/is (every? #(m/valid-double? (double %)) (sut/posterior-samples g [1e6 -1e6 2.2])))
        (t/is (empty? (sut/posterior-samples g [])))))))

(t/deftest prior-samples-statistics
  (doseq [[nm ys normalize?] [["normalised" [190.0 215.0 240.0 205.0 180.0 230.0] true]
                              ["not normalised" [1.0 3.0 2.0 5.0 4.0 3.5] false]]
          :let [xs [0.0 0.4 0.9 1.3 1.8 2.2]
                g (sut/gaussian-process xs ys {:kernel (gauss 0.9) :kscale 2.0 :noise 0.25 :normalize? normalize?})
                n 5000
                rng (r/rng :isaac 4321)
                samples (vec (repeatedly n #(first (sut/prior-samples g [0.7] rng))))
                ;; the prior at a point: the mean ymean, the variance ystddev^2 (kscale + noise)
                expected-sd (m/* (double (:ystddev g)) (m/sqrt 2.25))
                se-mean (m// expected-sd (m/sqrt (double n)))
                se-sd (m// expected-sd (m/sqrt (m/* 2.0 (m/dec n))))]]
    (t/is (< (m/abs (m/- (mean samples) (double (:ymean g)))) (m/* 4.0 se-mean)) (str nm " mean"))
    (t/is (< (m/abs (m/- (sd samples) expected-sd)) (m/* 4.0 se-sd)) (str nm " sd")))
  (t/testing "options: seeded rng, default arity, dense grids, empty input"
    (let [g (sut/gaussian-process [0.0 1.0 2.0] [1.0 2.0 0.5] {:kernel (gauss 0.5)})
          grid (vec (range 0.0 3.0 0.01))]
      (t/is (= (sut/prior-samples g [0.1 0.2] (r/rng :isaac 8)) (sut/prior-samples g [0.1 0.2] (r/rng :isaac 8))))
      (t/is (= 2 (count (sut/prior-samples g [0.1 0.2]))))
      (t/is (= (count grid) (count (sut/prior-samples g grid (r/rng :isaac 8)))))
      (t/is (empty? (sut/prior-samples g []))))))

(t/deftest predict-contract
  (let [xs [0.0 1.0 2.0 3.0]
        ys [1.0 3.0 2.0 5.0]
        g (sut/gaussian-process xs ys {:kernel (gauss 0.8)})]
    (t/testing "returns a double, and a Vec2 with the standard deviation"
      (t/is (instance? Double (sut/predict g 1.5)))
      (t/is (instance? fastmath.vector.Vec2 (sut/predict g 1.5 true)))
      (t/is (seq? (sut/predict-all g [1.5 2.5])))
      (t/is (= 2 (count (sut/predict-all g [1.5 2.5] true)))))
    (t/testing "at a training point the standard deviation is about zero, far away it is the prior one"
      (t/is (< (second (sut/predict g 1.0 true)) 1e-3))
      (t/is (m/delta-eq 1.0 (second (sut/predict g 1e3 true)) 1e-9))
      (t/is (m/delta-eq 0.0 (sut/predict g 1e3) 1e-12) "the zero mean far away")
      (let [gn (sut/gaussian-process xs ys {:kernel (gauss 0.8) :normalize? true})]
        (t/is (m/delta-eq (:ystddev gn) (second (sut/predict gn 1e3 true)) 1e-9) "the prior sd in the original units")
        (t/is (m/delta-eq (:ymean gn) (sut/predict gn 1e3) 1e-9))))
    (t/testing "very near a training point, and the huge scale"
      (doseq [e [1e-100 1e-154 1e-155 1e-300]]
        (t/is (m/delta-eq 3.0 (sut/predict g (m/+ 1.0 (double e))) 1e-6)))
      (t/is (m/delta-eq 1.0 (sut/predict g 1e-300) 1e-6))
      (t/is (m/valid-double? (sut/predict g 1e300)))
      (t/is (m/nan? (sut/predict g ##NaN))))
    (t/testing "empty query"
      (t/is (empty? (sut/predict-all g []))))))
