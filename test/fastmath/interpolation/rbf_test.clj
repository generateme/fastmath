(ns fastmath.interpolation.rbf-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.kernel :as k]
            [fastmath.vector :as v]
            [fastmath.interpolation.rbf :as sut]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; Reference values: utils/fastmath/dev/generate_scattered_reference.py (scipy RBFInterpolator for the
;; interpolation and the ridge `smoothing`, numpy for the normal equations).
(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "interpolation/scattered_reference.edn")))))

(defn- property
  "Returns the :property of the invalid-input ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (= :rbf (:method d))) (:property d))))))

(defn- max-err
  ^double [f queries expected]
  (reduce (fn [^double acc [q e]]
            (m/max acc (m/abs (m/- (double (f q)) (double e)))))
          0.0 (map vector queries expected)))

(defn- max-diff
  "Maximum absolute difference of two functions on the queries."
  ^double [f g queries]
  (max-err f queries (map g queries)))

;; arithmetic helpers with primitive hints
(defn- plane
  ^double [[^double x ^double y]]
  (m/+ 1.0 (m/+ (m/* 2.0 x) (m/* -3.0 y))))

(defn- scale-point
  [^double c p]
  (mapv (fn [^double x] (m/* c x)) p))

(defn- linear-terms
  [x]
  (into [1.0] x))

(defn- poly-terms
  [degree]
  (case (long degree)
    -1 nil
    0 (constantly [1.0])
    1 linear-terms))

(t/deftest scipy-interpolation
  (doseq [{:keys [name kernel params degree xs ys queries cases]} (:rbf @reference)
          :let [f (sut/rbf xs ys (k/rbf kernel params) (cond-> {} (poly-terms degree) (assoc :polynomial-terms (poly-terms degree))))
                {:keys [values node-values]} (first (filter #(m/zero? (double (:smoothing %))) cases))]]
    ;; measured maximum 1.1e-12 (gaussian and multiquadric 2-d, conditioning limited), others 6e-14 and below
    (t/is (< (max-err f queries values) 1e-10) (str name " 500 random points"))
    (t/is (< (max-err f xs node-values) 1e-10) (str name " at the nodes"))
    (t/is (< (max-err f xs ys) 1e-9) (str name " passes through the data"))))

(t/deftest scipy-ridge-smoothing
  (doseq [{:keys [name kernel params degree xs ys queries cases]} (:rbf @reference)
          {:keys [smoothing values]} cases
          :let [f (sut/rbf xs ys (k/rbf kernel params)
                           (cond-> {:lambda smoothing :regularization :ridge}
                             (poly-terms degree) (assoc :polynomial-terms (poly-terms degree))))]]
    ;; measured maximum 4.7e-15 for smoothing 1e-3 and 1e-1, 1.1e-12 for 0
    (t/is (< (max-err f queries values) 1e-10) (str name " smoothing " smoothing))))

(t/deftest normal-equations-regularization
  (doseq [{:keys [scale kscale lambda degree xs ys queries values]} (:rbf-normal @reference)
          :let [f (sut/rbf xs ys (k/rbf :gaussian {:scale scale})
                           (cond-> {:lambda lambda :kscale kscale}
                             (poly-terms degree) (assoc :polynomial-terms (poly-terms degree))))
                g (sut/rbf xs ys (k/rbf :gaussian {:scale scale})
                           (cond-> {:lambda lambda :kscale kscale :regularization :normal}
                             (poly-terms degree) (assoc :polynomial-terms (poly-terms degree))))]]
    ;; the normal equations square the condition number: measured maximum 1.3e-11 (kscale 10, lambda 1e-4), others 6e-13
    (t/is (< (max-err f queries values) 1e-9) (str "degree " degree " kscale " kscale " lambda " lambda))
    (t/is (= (mapv f queries) (mapv g queries)) ":normal is the default")))

(t/deftest default-unchanged
  ;; golden values of the implementation before the regularization option (HEAD 64adc9f), the default is :normal
  (let [xs [[0.0 0.0] [1.0 0.0] [0.0 1.0] [1.0 1.0] [0.5 0.5] [0.2 0.8]]
        ys [1.0 2.0 3.0 5.0 4.0 2.5]
        qs [[0.3 0.3] [0.7 0.4] [2.0 2.0]]
        check (fn [expected f] (t/is (< (max-err f qs expected) 1e-12)))]
    (check [2.685306923846209 4.580624468715261 0.06957712386321736]
           (sut/rbf xs ys (k/rbf :gaussian {:scale 0.7}) {:lambda 1e-3 :kscale 2.0}))
    (check [2.9385247542664343 4.4128295264286805 8.760425122219436]
           (sut/rbf xs ys (k/rbf :gaussian {:scale 0.7}) {:lambda 1e-3 :kscale 2.0 :polynomial-terms linear-terms}))
    (check [2.901372762219837 4.239307942087227 7.070280356553653]
           (sut/rbf xs ys (k/rbf :radial-powers-3) {:polynomial-terms linear-terms}))
    (check [2.6789699028837406 3.68733422236705 8.307499601353037]
           (sut/rbf xs ys (k/rbf :thin-plate) {:polynomial-terms linear-terms :lambda 0.01}))
    (let [f (sut/rbf [0.0 1.0 2.0 3.0 4.0] [0.0 1.0 0.0 1.0 0.0])]
      (t/is (< (max-err f [0.5 1.5 2.5 10.0] [0.6886627520835753 0.4554551555412274 0.4554551555412273 -1.3232321200512878E-16])
               1e-12)))))

(t/deftest kscale-and-lambda
  (let [xs [[0.0 0.0] [1.0 0.0] [0.0 1.0] [1.0 1.0] [0.5 0.5] [0.2 0.8] [0.9 0.3]]
        ys [1.0 2.0 3.0 5.0 4.0 2.5 1.5]
        qs [[0.3 0.3] [0.7 0.4] [0.1 0.9] [0.5 0.2]]
        kern (k/rbf :gaussian {:scale 0.8})
        mk (fn [opts] (sut/rbf xs ys kern opts))]
    (t/testing "kscale cancels without regularization"
      ;; measured 2.1e-14 for the factors 0.1, 10 and 1000
      (doseq [ks [0.1 10.0 1000.0]]
        (t/is (< (max-diff (mk {:kscale ks}) (mk {}) qs) 1e-10) (str "kscale " ks))))
    (t/testing ":normal without polynomial: lambda acts as lambda / kscale^2"
      (doseq [ks [0.5 3.0 20.0]
              lam [1e-4 1e-2]
              :let [ks (double ks) lam (double lam)]]
        (t/is (< (max-diff (mk {:kscale ks :lambda lam})
                           (mk {:lambda (m// lam (m/* ks ks))}) qs)
                 1e-10) ; measured 5.3e-13
              (str "kscale " ks " lambda " lam))))
    (t/testing ":normal with polynomial terms: the polynomial block is not scaled, the relation does not hold"
      (let [pt {:polynomial-terms linear-terms}]
        (t/is (> (max-diff (mk (assoc pt :kscale 3.0 :lambda 1e-2))
                           (mk (assoc pt :lambda (/ 1e-2 9.0))) qs)
                 1e-4))))
    (t/testing ":ridge: lambda acts as lambda / kscale, also with polynomial terms"
      (doseq [ks [0.5 3.0 20.0]
              lam [1e-4 1e-2]
              pt [{} {:polynomial-terms linear-terms}]
              :let [ks (double ks) lam (double lam)]]
        (t/is (< (max-diff (mk (assoc pt :kscale ks :lambda lam :regularization :ridge))
                           (mk (assoc pt :lambda (m// lam ks) :regularization :ridge)) qs)
                 1e-10) ; measured 2.8e-14
              (str "kscale " ks " lambda " lam " " pt))))))

(t/deftest regularization-options
  (let [xs [0.0 1.0 2.0 3.0 4.0 5.0]
        ys [0.0 1.0 0.5 2.0 1.5 3.0]
        kern (k/rbf :gaussian {:scale 1.5})
        interp (sut/rbf xs ys kern)]
    (t/testing "lambda 0 is the interpolation for both"
      (t/is (= (mapv interp [0.5 2.5 7.0])
               (mapv (sut/rbf xs ys kern {:lambda 0.0 :regularization :ridge}) [0.5 2.5 7.0])))
      (t/is (= (mapv interp [0.5 2.5 7.0])
               (mapv (sut/rbf xs ys kern {:lambda 0 :regularization :normal}) [0.5 2.5 7.0]))))
    (t/testing "the ridge smooths: the residual at the nodes grows with lambda, the huge lambda gives the zero function"
      (let [res (fn [lam] (max-err (sut/rbf xs ys kern {:lambda lam :regularization :ridge}) xs ys))]
        (t/is (< (res 0.0) 1e-8))
        (t/is (< (res 0.0) (res 1e-4) (res 1e-2) (res 1.0)))
        (t/is (< (max-err (sut/rbf xs ys kern {:lambda 1e12 :regularization :ridge}) xs (repeat 0.0)) 1e-9))))
    (t/testing "nil options are missing options"
      (t/is (= (mapv interp [0.5 2.5])
               (mapv (sut/rbf xs ys kern {:kscale nil :lambda nil :regularization nil :distance nil :polynomial-terms nil}) [0.5 2.5])))
      (t/is (= (mapv interp [0.5 2.5]) (mapv (sut/rbf xs ys kern nil) [0.5 2.5]))))
    (t/testing "the ridge keeps the polynomial block unregularized: a plane stays a plane for a huge lambda"
      (let [pxs [[0.0 0.0] [1.0 0.0] [0.0 1.0] [1.0 1.0] [0.5 0.3] [0.2 0.9]]
            pys (mapv plane pxs)
            f (sut/rbf pxs pys (k/rbf :thin-plate) {:lambda 1e6 :regularization :ridge :polynomial-terms linear-terms})
            qs [[0.3 0.3] [0.9 0.1] [2.0 -1.0]]]
        (t/is (< (max-err f qs (map plane qs)) 1e-6))))))

(t/deftest polynomial-reproduction
  (let [rnd (java.util.Random. 7)
        xs (vec (repeatedly 25 #(vector (.nextDouble rnd) (.nextDouble rnd))))
        qs (vec (repeatedly 30 #(vector (m/- (m/* 3.0 (.nextDouble rnd)) 1.0) (m/- (m/* 3.0 (.nextDouble rnd)) 1.0))))]
    ;; conditionally positive definite kernels with a polynomial of the matching degree reproduce it (also far outside the data):
    ;; a plane for the thin plate and the cubic kernel (degree 1), a constant for the linear and the multiquadric kernel (degree 0)
    (doseq [[nm kern deg] [["thin-plate" (k/rbf :thin-plate) 1]
                           ["cubic" (k/rbf :radial-powers-3) 1]
                           ["linear" (k/rbf :radial-powers {:beta 1.0}) 0]
                           ["multiquadric" (k/rbf :generalized-multiquadratic {:beta 0.5}) 0]]
            :let [target (if (m/== 1 (long deg)) plane (constantly 3.7))
                  f (sut/rbf xs (mapv target xs) kern {:polynomial-terms (poly-terms deg)})]]
      (t/is (< (max-err f qs (map target qs)) 1e-8) nm))))

(t/deftest cpd-kernel-without-polynomial-differs
  ;; documented: a conditionally positive definite kernel without the polynomial is not the usual interpolant
  (let [xs [0.0 1.0 2.0 3.0 4.0]
        ys [0.0 1.0 4.0 9.0 16.0]
        with (sut/rbf xs ys (k/rbf :radial-powers-3) {:polynomial-terms (fn [x] [1.0 x])})
        without (sut/rbf xs ys (k/rbf :radial-powers-3))]
    (t/is (< (max-err with xs ys) 1e-9))
    (t/is (< (max-err without xs ys) 1e-9) "both pass through the data")
    (t/is (> (m/abs (- (double (with 6.0)) (double (without 6.0)))) 1e-3) "but differ elsewhere")))

(t/deftest minimal-sizes
  (t/testing "one point"
    (let [f (sut/rbf [2.0] [3.0])]
      (t/is (= 3.0 (f 2.0)))
      (t/is (m/delta-eq 3.0 (f 2.0000001) 1e-6))
      (t/is (= 0.0 (f 100.0)) "the gaussian vanishes far away"))
    (let [f (sut/rbf [[1.0 2.0]] [3.0] (k/rbf :radial-powers-3) {:polynomial-terms (constantly [1.0])})]
      (t/is (= 3.0 (f [1.0 2.0])))))
  (t/testing "two points, a constant polynomial"
    (let [f (sut/rbf [0.0 1.0] [2.0 4.0] (k/rbf :radial-powers {:beta 1.0}) {:polynomial-terms (constantly [1.0])})]
      (t/is (m/delta-eq 2.0 (f 0.0) 1e-12))
      (t/is (m/delta-eq 4.0 (f 1.0) 1e-12))
      (t/is (m/delta-eq 3.0 (f 0.5) 1e-12) "the linear radial power with a constant is the linear interpolation between two points"))))

(t/deftest query-forms
  (let [xs [[1.0 0.0] [3.0 3.0] [0.0 2.0] [2.0 1.0]]
        ys [2.0 4.0 1.0 3.0]
        kern (k/rbf :gaussian {:scale 2.0})
        f (sut/rbf xs ys kern)]
    (t/testing "long vs double, Vec2 and lists, in the nodes and in the queries"
      (let [fl (sut/rbf [[1 0] [3 3] [0 2] [2 1]] ys kern)
            fv (sut/rbf (mapv (fn [[x y]] (v/vec2 x y)) xs) ys kern)
            q [1.3 0.7]]
        (t/is (m/delta-eq 2.0 (f [1 0]) 1e-9))
        (t/is (m/delta-eq 2.0 (fl [1.0 0.0]) 1e-9))
        (t/is (m/delta-eq 2.0 (fv [1 0]) 1e-9))
        (t/is (= (f q) (fl q) (fv q) (f (v/vec2 1.3 0.7)) (f (list 1.3 0.7)) (fv (v/vec2 1.3 0.7))))))
    (t/testing "one dimension: numbers, long nodes, one element query"
      (let [f1 (sut/rbf [0 1 2 3] [0.0 1.0 0.0 1.0])
            f2 (sut/rbf [0.0 1.0 2.0 3.0] [0 1 0 1])]
        (t/is (= (f1 1.5) (f2 1.5) (f1 [1.5])))))
    (t/testing "a wrong dimension of the query"
      (t/is (= :dimension-mismatch (property #(f [1.0 2.0 3.0]))))
      (t/is (= :dimension-mismatch (property #(f [1.0]))))
      (t/is (= :dimension-mismatch (property #(f 1.0)))))
    (t/is (= :dimension-mismatch (property #(sut/rbf [[0.0 0.0] [1.0]] [1.0 2.0]))))))

(t/deftest custom-distance-and-scales
  (let [xs [0.0 0.7 1.5 2.0 3.1]
        ys [1.0 0.0 2.0 1.0 3.0]
        base (sut/rbf xs ys (k/rbf :gaussian {:scale 1.0}))
        qs [0.3 1.0 2.5 4.0]]
    (t/testing "a distance doubled with a scale doubled is the same function"
      (t/is (< (max-diff (sut/rbf xs ys (k/rbf :gaussian {:scale 2.0}) {:distance (fn [a b] (m/* 2.0 (m/abs (m/- (double a) (double b)))))})
                         base qs)
               1e-9)))
    (t/testing "the same function with the default and an explicit distance"
      (t/is (< (max-diff (sut/rbf xs ys (k/rbf :gaussian) {:distance (fn [a b] (m/abs (- (double a) (double b))))}) base qs) 1e-12))))
  (t/testing "huge and tiny scales of the coordinates with the matching kernel scale give the same function"
    (let [xs [[0.0 0.0] [1.0 0.0] [0.0 1.0] [1.0 1.0] [0.5 0.5] [0.2 0.8]]
          ys [1.0 2.0 3.0 5.0 4.0 2.5]
          qs [[0.3 0.3] [0.7 0.4] [0.1 0.9]]
          base (sut/rbf xs ys (k/rbf :gaussian {:scale 0.7}))]
      (doseq [c [1e-8 1e8 1e150]
              :let [c (double c)
                    scaled (sut/rbf (map #(scale-point c %) xs) ys (k/rbf :gaussian {:scale (m/* c 0.7)}))]]
        (t/is (< (max-err scaled (map #(scale-point c %) qs) (map base qs)) 1e-8) (str "scale " c)))))
  (t/testing "very near a node and far away"
    (let [f (sut/rbf [0.0 1.0 2.0] [1.0 5.0 2.0] (k/rbf :gaussian {:scale 0.5}))]
      (doseq [e [1e-100 1e-154 1e-155 1e-300]]
        (t/is (m/delta-eq 1.0 (f e) 1e-12) "a node at zero")
        (t/is (m/delta-eq 1.0 (f (m/- (double e))) 1e-12)))
      (t/is (= 0.0 (f 1e150)))
      (t/is (= 0.0 (f -1e300))))
    (let [f (sut/rbf [[0.0 0.0] [1.0 0.0] [0.0 1.0] [1.0 1.0]] [1.0 5.0 2.0 3.0] (k/rbf :thin-plate) {:polynomial-terms linear-terms})]
      (t/is (m/valid-double? (f [1e100 1e100])))
      (t/is (m/delta-eq 5.0 (f [1.0 1e-300]) 1e-12))
      (t/is (m/delta-eq 5.0 (f [1.0 1e-155]) 1e-12)))))

(t/deftest duplicates-and-nan
  (t/testing "repeated nodes: the SVD fallback gives the mean of the values, silently"
    (let [out (java.io.StringWriter.)
          err (java.io.StringWriter.)
          f (binding [*out* out *err* err]
              (sut/rbf [0.0 1.0 1.0 2.0] [0.0 1.0 3.0 4.0]))]
      (t/is (= "" (str out)) "no output on stdout")
      (t/is (= "" (str err)))
      (t/is (m/delta-eq 2.0 (f 1.0) 1e-6))
      (t/is (m/delta-eq 0.0 (f 0.0) 1e-6))
      (t/is (m/valid-double? (f 0.5)))))
  (t/testing "a repeated node with the same value is harmless"
    (let [f (sut/rbf [0.0 1.0 1.0 2.0] [0.0 1.0 1.0 4.0])]
      (t/is (m/delta-eq 1.0 (f 1.0) 1e-6))))
  (t/testing "NaN"
    (t/is (m/nan? ((sut/rbf [0.0 1.0 2.0] [0.0 1.0 0.5]) ##NaN)))
    (t/is (m/nan? ((sut/rbf [0.0 1.0 2.0] [0.0 ##NaN 0.5]) 0.5)))))

(t/deftest many-points
  (let [n 200
        rnd (java.util.Random. 99)
        xs (vec (repeatedly n #(vector (.nextDouble rnd) (.nextDouble rnd))))
        ys (mapv (fn [[^double x ^double y]] (m/+ (m/sin (m/* 3.0 x)) (m/* x y))) xs)
        f (sut/rbf xs ys (k/rbf :thin-plate) {:polynomial-terms linear-terms})]
    (t/is (< (max-err f xs ys) 1e-9))
    (t/is (m/delta-eq (m/+ (m/sin 1.2) 0.24) (f [0.4 0.6]) 0.02) "an interpolation of a smooth function")))

(t/deftest invalid-input
  (t/testing "lengths and minimal size"
    (t/is (= :length-mismatch (property #(sut/rbf [0 1] [1.0]))))
    (t/is (= :length-mismatch (property #(sut/rbf [0 1] [1.0 2.0 3.0]))))
    (t/is (= :too-few-points (property #(sut/rbf [] []))))
    (t/is (= :too-few-points (property #(sut/rbf nil nil))))
    (t/is (nil? (property #(sut/rbf [0] [1.0])))))
  (t/testing "lambda: the first invalid values and the boundary"
    (doseq [lam [-1.0 -1e-300 -1 ##NaN ##Inf ##-Inf]]
      (t/is (= :lambda (property #(sut/rbf [0 1] [1.0 2.0] nil {:lambda lam}))) (str "lambda " lam))
      (t/is (= :lambda (property #(sut/rbf [0 1] [1.0 2.0] nil {:lambda lam :regularization :ridge}))) (str "ridge lambda " lam)))
    (t/is (nil? (property #(sut/rbf [0 1] [1.0 2.0] nil {:lambda 0}))))
    (t/is (nil? (property #(sut/rbf [0 1] [1.0 2.0] nil {:lambda 0.0}))))
    (t/is (nil? (property #(sut/rbf [0 1] [1.0 2.0] nil {:lambda 1e-300 :regularization :ridge})))))
  (t/testing "regularization"
    (t/is (= :regularization (property #(sut/rbf [0 1] [1.0 2.0] nil {:regularization :tikhonov}))))
    (t/is (= :regularization (property #(sut/rbf [0 1] [1.0 2.0] nil {:regularization "ridge"}))))
    (t/is (nil? (property #(sut/rbf [0 1] [1.0 2.0] nil {:regularization :ridge}))))
    (t/is (nil? (property #(sut/rbf [0 1] [1.0 2.0] nil {:regularization :normal})))))
  (t/testing "the exception data"
    (let [d (try (sut/rbf [0 1] [1.0 2.0] nil {:lambda -2.0}) (catch clojure.lang.ExceptionInfo e (ex-data e)))]
      (t/is (= {:type :invalid-input :method :rbf :property :lambda :option :lambda :value -2.0}
               (select-keys d [:type :method :property :option :value]))))))
