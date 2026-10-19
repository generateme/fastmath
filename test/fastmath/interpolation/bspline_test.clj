(ns fastmath.interpolation.bspline-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.interpolation.bspline :as sut]))

(set! *warn-on-reflection* true)

;; scipy 1.18.1 make_interp_spline reference, see
;; utils/fastmath/dev/generate_bspline_bicubic_reference.py
(def ^:private reference
  (delay (:bspline (edn/read-string (slurp (io/resource "interpolation/bspline_bicubic_reference.edn"))))))

(defn- property
  "Returns the :property of the ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (= :b-spline (:method d))) (:property d))))))

(defn- agree?
  "Absolute or relative closeness (the larger of both tolerances)."
  [a b tol]
  (m/delta-eq (double a) (double b) tol tol))

(defn- max-err
  "Maximal combined error between two sequences."
  ^double [as bs]
  (double (reduce max 0.0 (map (fn [a b]
                                 (let [a (double a)
                                       b (double b)]
                                   (/ (m/abs (- a b)) (max 1.0 (m/abs b)))))
                               as bs))))

(t/deftest scipy-make-interp-spline
  ;; 9 datasets x degrees 1..5; queries: on every knot, mid points, random interior,
  ;; 1 and 10 percent of the span outside on both sides (end polynomial continued).
  ;; Measured: max combined (absolute or relative) error 7.1e-13 (jump9, degree 5, the worst
  ;; conditioned collocation matrix, at 10 percent outside; value 233, relative 3e-14),
  ;; all other cases below 1e-13.
  (let [{:keys [datasets cases]} @reference]
    (t/is (= 45 (count cases)))
    (t/is (= 9 (count datasets)))
    (doseq [{:keys [dataset degree expected]} cases
            :let [{:keys [xs ys queries]} (get datasets dataset)
                  f (sut/b-spline xs ys {:degree degree})
                  err (max-err (map f queries) expected)]]
      (t/testing (str dataset " degree " degree)
        (t/is (< err 1.0e-12) (str "error " err))))))

(t/deftest scipy-default-degree-is-3
  (let [{:keys [datasets cases]} @reference]
    (doseq [{:keys [dataset degree expected]} cases
            :when (= 3 degree)
            :let [{:keys [xs ys queries]} (get datasets dataset)]]
      (t/is (< (max-err (map (sut/b-spline xs ys) queries) expected) 1.0e-12))
      (t/is (< (max-err (map (sut/b-spline xs ys nil) queries) expected) 1.0e-12))
      (t/is (< (max-err (map (sut/b-spline xs ys {}) queries) expected) 1.0e-12))
      (t/is (< (max-err (map (sut/b-spline xs ys {:degree nil :knots nil}) queries) expected) 1.0e-12)))))

(t/deftest passes-through-all-nodes
  ;; measured: at most 6.7e-16 relative to the data scale
  (let [{:keys [datasets cases]} @reference]
    (doseq [{:keys [dataset degree]} cases
            :let [{:keys [xs ys]} (get datasets dataset)
                  scale (double (reduce max 1.0 (map m/abs ys)))
                  f (sut/b-spline xs ys {:degree degree})]]
      (t/testing (str dataset " degree " degree)
        (t/is (every? true? (map (fn [x y] (m/delta-eq (f x) y (* 1.0e-13 scale))) xs ys)))))))

(t/deftest scipy-knots-passed-explicitly
  ;; the default knots of scipy given as :knots reproduce the default result exactly
  (let [{:keys [datasets cases]} @reference]
    (doseq [{:keys [dataset degree knots]} cases
            :let [{:keys [xs ys queries]} (get datasets dataset)
                  f (sut/b-spline xs ys {:degree degree})
                  g (sut/b-spline xs ys {:degree degree :knots knots})]]
      (t/is (< (max-err (map f queries) (map g queries)) 1.0e-13) (str dataset " " degree)))))

(t/deftest custom-knots-scipy
  ;; clamped knots with chosen interior knots, and uniform knots spanning the data
  ;; measured: 1.1e-14 and 1.4e-14
  (doseq [{:keys [name degree xs ys knots queries expected]} (:custom-knots @reference)
          :let [f (sut/b-spline xs ys {:degree degree :knots knots})]]
    (t/is (< (max-err (map f queries) expected) 1.0e-12) name)
    (t/is (every? true? (map (fn [x y] (m/delta-eq (f x) y 1.0e-13)) xs ys)) name)))

(t/deftest custom-knots-validation
  (let [xs [0.0 1.0 2.0 3.0 4.0]
        ys [1.0 2.0 0.0 2.0 1.0]
        ok [0.0 0.0 0.0 0.0 2.0 4.0 4.0 4.0 4.0]]
    (t/is (nil? (property #(sut/b-spline xs ys {:knots ok}))))
    (t/testing "wrong number of knots"
      (t/is (= :knots-count (property #(sut/b-spline xs ys {:knots (butlast ok)}))))
      (t/is (= :knots-count (property #(sut/b-spline xs ys {:knots (conj ok 4.0)}))))
      (t/is (= :knots-count (property #(sut/b-spline xs ys {:knots []}))))
      (t/is (= :knots-count (property #(sut/b-spline xs ys {:degree 2 :knots ok})))))
    (t/testing "knots not sorted, not finite"
      (t/is (= :knots-not-sorted (property #(sut/b-spline xs ys {:knots [0.0 0.0 0.0 0.0 4.0 2.0 4.0 4.0 4.0]}))))
      (t/is (= :not-finite (property #(sut/b-spline xs ys {:knots [0.0 0.0 0.0 0.0 ##NaN 4.0 4.0 4.0 4.0]}))))
      (t/is (= :not-finite (property #(sut/b-spline xs ys {:knots [0.0 0.0 0.0 0.0 2.0 4.0 4.0 4.0 ##Inf]})))))
    (t/testing "nodes outside of the knot interval"
      (t/is (= :knots-domain (property #(sut/b-spline xs ys {:knots [0.5 0.5 0.5 0.5 2.0 4.0 4.0 4.0 4.0]}))))
      (t/is (= :knots-domain (property #(sut/b-spline xs ys {:knots [0.0 0.0 0.0 0.0 2.0 3.5 3.5 3.5 3.5]}))))
      (t/is (= :knots-domain (property #(sut/b-spline xs ys {:knots [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0]})))))
    (t/testing "Schoenberg-Whitney violated: the collocation matrix is singular"
      (let [xs6 [0.0 1.0 2.0 3.0 4.0 5.0]
            ys6 [1.0 2.0 0.0 2.0 1.0 0.5]]
        (t/is (nil? (property #(sut/b-spline xs6 ys6 {:knots [0.0 0.0 0.0 0.0 2.0 3.0 5.0 5.0 5.0 5.0]}))))
        ;; x1 = 1 is not below the knot t5 = 0.5
        (t/is (= :schoenberg-whitney (property #(sut/b-spline xs6 ys6 {:knots [0.0 0.0 0.0 0.0 0.2 0.5 5.0 5.0 5.0 5.0]}))))
        ;; x4 = 4 is not above the knot t4 = 4.5
        (t/is (= :schoenberg-whitney (property #(sut/b-spline xs6 ys6 {:knots [0.0 0.0 0.0 0.0 4.5 4.8 5.0 5.0 5.0 5.0]}))))))))

(defn- poly
  "Evaluates polynomial given by ascending coefficients."
  ^double [cs ^double x]
  (double (reduce (fn [^double acc c] (+ (* acc x) (double c))) 0.0 (reverse cs))))

(t/deftest polynomial-reproduction
  ;; a spline of degree k contains all polynomials of degree k: any number of points >= k+1,
  ;; also outside of the data range (the end polynomial is continued)
  (doseq [k (range 1 6)
          n [(inc k) (+ k 2) (+ k 5) 20]
          :let [cs (map #(/ (double (inc %)) (double (+ 2 (* 2 %)))) (range (inc k)))
                xs (mapv #(+ -1.5 (* 0.37 (double %)) (* 0.05 (m/sin (double (* 3 %))))) (range n))
                ys (mapv #(poly cs %) xs)
                f (sut/b-spline xs ys {:degree k})
                span (- (double (peek xs)) (double (first xs)))
                qs (concat [(- (first xs) 0.7)] (map #(+ (first xs) (* span (/ % 40.0))) (range 41))
                           [(+ (peek xs) 0.4)])]]
    (t/is (every? true? (map (fn [q] (agree? (f q) (poly cs q) 1.0e-9)) qs))
          (str "degree " k " points " n))))

(t/deftest minimal-count-is-the-interpolating-polynomial
  ;; n = degree + 1: single polynomial through all points (Lagrange form as oracle)
  (doseq [k (range 1 6)
          :let [xs (mapv #(+ (* 0.9 (double %)) (* 0.1 (m/sq (double %)))) (range (inc k)))
                ys (mapv #(m/sin (* 1.7 (double %))) xs)
                lagrange (fn [^double x]
                           (reduce + (map-indexed
                                      (fn [i yi]
                                        (* (double yi)
                                           (reduce * (map-indexed (fn [j xj] (if (= i j) 1.0 (/ (- x (double xj)) (- (double (xs i)) (double xj))))) xs))))
                                      ys)))
                f (sut/b-spline xs ys {:degree k})]]
    (doseq [q [-1.0 0.0 0.33 1.1 2.0 3.7 (double (peek xs)) 9.0]]
      (t/is (m/delta-eq (f q) (lagrange q) 1.0e-9 1.0e-11) (str "degree " k " at " q)))))

(t/deftest degree-one-is-piecewise-linear
  (let [xs [0.0 1.0 3.0 4.0]
        ys [2.0 -1.0 5.0 4.0]
        f (sut/b-spline xs ys {:degree 1})]
    (doseq [[x e] [[0.0 2.0] [0.5 0.5] [1.0 -1.0] [2.0 2.0] [3.0 5.0] [3.5 4.5] [4.0 4.0]
                   ;; end segments are continued linearly
                   [-1.0 5.0] [5.0 3.0]]]
      (t/is (m/delta-eq e (f x) 1.0e-14) (str "x " x)))))

(t/deftest left-of-the-domain-follows-the-spline
  ;; I-03: the former implementation returned the last y left of the domain
  (let [xs [0.0 1.0 2.0 3.0 4.0 5.0]
        ys (mapv #(- (* % % %) (* 2.0 %)) xs)
        f (sut/b-spline xs ys)]
    (t/is (m/delta-eq (- (* -1.65 -1.65 -1.65) (* 2.0 -1.65)) (f -1.65) 1.0e-12))
    (t/is (m/delta-eq (- (* 7.5 7.5 7.5) (* 2.0 7.5)) (f 7.5) 1.0e-10))
    (t/is (not (m/delta-eq (last ys) (f -1.65) 1.0)))
    (t/is (m/delta-eq (f -1.0e-9) (f 0.0) 1.0e-8) "continuous at the left end")
    (t/is (m/delta-eq (f (+ 5.0 1.0e-9)) (f 5.0) 1.0e-6) "continuous at the right end (slope 73)")))

(t/deftest smoothness-at-nodes
  ;; degree k: continuous derivatives up to k-1 at the simple knots (finite difference check)
  (let [xs [0.0 0.8 1.7 2.1 3.5 4.0 5.2 6.0 7.5]
        ys [0.3 -1.0 0.5 2.0 1.1 -0.4 0.9 0.2 -1.5]]
    (doseq [k [2 3 4 5]
            :let [f (sut/b-spline xs ys {:degree k})
                  h 1.0e-5]
            x (subvec xs (inc (quot k 2)) (- (count xs) (inc (quot k 2))))]
      (t/is (m/delta-eq (f (- x 1.0e-9)) (f (+ x 1.0e-9)) 1.0e-7) (str "value k " k " x " x))
      (t/is (m/delta-eq (/ (- (f x) (f (- x h))) h) (/ (- (f (+ x h)) (f x)) h) 1.0e-3)
            (str "slope k " k " x " x)))))

(t/deftest not-a-knot-cubic
  ;; odd degree 3 default knots: the first two and the last two intervals are one polynomial
  (let [xs [0.0 1.0 2.5 3.0 4.2 5.0]
        ys [1.0 3.0 -1.0 2.0 0.5 1.5]
        f (sut/b-spline xs ys)
        ;; Lagrange cubic through four samples of the first piece [x0, x2]
        qs [0.2 0.9 1.7 2.3]
        piece-values (mapv f qs)
        lag (fn [^double x] (reduce + (map-indexed
                                       (fn [i yi]
                                         (* (double yi)
                                            (reduce * (map-indexed (fn [j qj] (if (= i j) 1.0 (/ (- x (double qj)) (- (double (qs i)) (double qj))))) qs))))
                                       piece-values)))]
    (doseq [q [0.0 0.5 1.0 1.4 2.0 2.4]]
      (t/is (m/delta-eq (f q) (lag q) 1.0e-11) (str "x " q)))
    ;; x_2 = 2.5 is a knot: the same cubic does not describe the next interval
    (t/is (not (m/delta-eq (f 2.8) (lag 2.8) 1.0e-3)))))

(t/deftest invalid-input
  (t/testing "too few points: degree + 1 are needed"
    (t/is (= :too-few-points (property #(sut/b-spline [] []))))
    (t/is (= :too-few-points (property #(sut/b-spline [1.0] [2.0]))))
    (t/is (= :too-few-points (property #(sut/b-spline [0.0 1.0 2.0] [1.0 2.0 3.0]))))
    (t/is (= :too-few-points (property #(sut/b-spline [0.0 1.0] [1.0 2.0] {:degree 2}))))
    (t/is (= :too-few-points (property #(sut/b-spline (range 5) (range 5) {:degree 5}))))
    (t/is (= :too-few-points (property #(sut/b-spline (range 5) (range 5) {:degree 1.0e30}))))
    (t/is (nil? (property #(sut/b-spline [0.0 1.0 2.0 3.0] [1.0 2.0 3.0 5.0]))))
    (t/is (nil? (property #(sut/b-spline [0.0 1.0] [1.0 2.0] {:degree 1})))))
  (t/testing "xs not strictly increasing or not finite"
    (t/is (= :not-increasing (property #(sut/b-spline [0.0 2.0 1.0 3.0] [1.0 2.0 3.0 4.0]))))
    (t/is (= :not-increasing (property #(sut/b-spline [0.0 1.0 1.0 3.0] [1.0 2.0 3.0 4.0]))))
    (t/is (= :not-increasing (property #(sut/b-spline [3.0 2.0 1.0 0.0] [1.0 2.0 3.0 4.0]))))
    (t/is (= :not-finite (property #(sut/b-spline [0.0 ##NaN 2.0 3.0] [1.0 2.0 3.0 4.0]))))
    (t/is (= :not-finite (property #(sut/b-spline [0.0 1.0 2.0 ##Inf] [1.0 2.0 3.0 4.0])))))
  (t/testing "different lengths"
    (t/is (= :length-mismatch (property #(sut/b-spline [0.0 1.0 2.0 3.0] [1.0 2.0 3.0]))))
    (t/is (= :length-mismatch (property #(sut/b-spline [0.0 1.0 2.0 3.0] [1.0 2.0 3.0 4.0 5.0])))))
  (t/testing "degree is an integer >= 1"
    (doseq [d [0 -1 -3 0.5 2.5 ##NaN ##Inf ##-Inf]]
      (t/is (= :degree (property #(sut/b-spline [0.0 1.0 2.0 3.0] [1.0 2.0 3.0 4.0] {:degree d}))) (str d)))
    (t/is (nil? (property #(sut/b-spline [0.0 1.0 2.0 3.0] [1.0 2.0 3.0 4.0] {:degree 3.0}))))
    (t/is (nil? (property #(sut/b-spline [0.0 1.0 2.0 3.0] [1.0 2.0 3.0 4.0] {:degree 1})))))
  (t/testing "ex-data of the exception"
    (let [d (try (sut/b-spline [0.0 1.0 1.0] [1.0 2.0 3.0] {:degree 1}) nil
                 (catch clojure.lang.ExceptionInfo e (ex-data e)))]
      (t/is (= {:type :invalid-input :method :b-spline :property :not-increasing :axis :xs :index 2 :value 1.0} d)))))

(t/deftest nan-and-infinite-values
  (let [xs [0.0 1.0 2.0 3.0 4.0]
        f (sut/b-spline xs [1.0 ##NaN 2.0 0.0 1.0])]
    (t/is (m/nan? (f 2.5)) "NaN in ys")
    (t/is (m/nan? (f 0.0)))
    (t/is (m/nan? ((sut/b-spline xs [1.0 ##Inf 2.0 0.0 1.0]) 2.5))))
  (let [f (sut/b-spline [0.0 1.0 2.0 3.0 4.0] [1.0 0.0 2.0 0.0 1.0])]
    (t/is (m/nan? (f ##NaN)) "NaN query")
    (t/is (number? (f 1.0e300)) "huge argument: some number, no exception")
    (t/is (not (m/nan? (f 1.0e10))))))

(t/deftest huge-and-tiny-scales
  ;; scaling of values and of nodes is a symmetry of the interpolation
  (let [xs [0.0 0.7 1.1 2.0 3.2 3.9 5.0]
        ys [0.3 -1.0 0.5 2.0 1.1 -0.4 0.9]
        base (sut/b-spline xs ys)]
    (doseq [c [1.0e-9 1.0e9]
            :let [fy (sut/b-spline xs (map #(* c (double %)) ys))
                  fx (sut/b-spline (map #(* c (double %)) xs) ys)]
            q [-0.5 0.35 1.1 2.7 4.4 5.0 6.0]]
      (t/is (m/delta-eq (* c (double (base q))) (fy q) 1.0e-12 1.0e-12) (str "ys scale " c " q " q))
      (t/is (m/delta-eq (base q) (fx (* c (double q))) 1.0e-9 1.0e-9) (str "xs scale " c " q " q)))))

(t/deftest many-points
  ;; banded solve: 20000 points, nodes reproduced, no dense matrix
  ;; measured max errors at the nodes below 1e-12, at the mid points 1.25e-7 (degree 1, h^2/8 f''),
  ;; 2.6e-9 (degree 3) and 1.0e-9 (degree 5, rounding amplification of the conditioning)
  (let [n 20000
        xs (mapv #(* 0.001 (double %)) (range n))
        ys (mapv #(m/sin %) xs)]
    (doseq [k [1 3 5]
            :let [f (sut/b-spline xs ys {:degree k})
                  node-err (double (reduce max 0.0 (map (fn [i] (m/abs (- (double (f (xs i))) (double (ys i)))))
                                                        (range 0 n 97))))
                  mid-err (double (reduce max 0.0 (map (fn [i] (let [x (* 0.001 (+ 0.5 (double i)))]
                                                                 (m/abs (- (double (f x)) (m/sin x)))))
                                                       (range 0 (dec n) 101))))]]
      (t/is (< node-err 1.0e-11) (str "degree " k))
      (t/is (< mid-err 2.0e-7) (str "degree " k)))))

(t/deftest evaluation-is-independent-of-order
  ;; the same function can be called from many threads
  (let [xs (mapv double (range 12))
        f (sut/b-spline xs (mapv #(m/sin %) xs) {:degree 4})
        qs (map #(+ 0.5 (* 0.01 (double %))) (range 1000))
        expected (mapv f qs)]
    (t/is (= expected (vec (pmap f qs))))))
