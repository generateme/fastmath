(ns fastmath.interpolation.grid-test
  (:require [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [fastmath.core :as m]
            [fastmath.interpolation.linear :as il]
            [fastmath.interpolation.cubic :as ic]))

(set! *warn-on-reflection* true)

;; Oracles (utils/fastmath/dev/generate_interp1d_reference.py):
;; - :bilinear  scipy RegularGridInterpolator(method='linear', bounds_error=False, fill_value=None)
;; - :cubic-2d  two passes of scipy CubicSpline(bc_type='natural'): along ys through every row of vss,
;;              then along xs

(def ^:private ref-data
  (delay (edn/read-string (slurp (io/resource "interpolation/interp1d_reference.edn")))))

(defn- property
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (= :invalid-input (:type d)) (:property d))))))

(defn- vmax ^double [vss] (double (reduce max 0.0 (map #(Math/abs (double %)) (flatten vss)))))

(t/deftest bilinear-reference
  ;; measured: 1.0e-15 relative
  (doseq [d (:bilinear @ref-data)]
    (let [f (il/bilinear (:xs d) (:ys d) (:vss d))
          sc (vmax (:vss d))]
      (doseq [[[x y] e] (map vector (:pts d) (:expected d))]
        (t/is (m/delta-eq e (f x y) (* 1e-13 sc) 1e-13) (str (:name d) " at " [x y]))
        (t/is (== (f x y) (f [x y])) "both arities give the same value")))))

(t/deftest cubic-2d-reference
  ;; measured: 1.4e-14 relative
  (doseq [d (:cubic-2d @ref-data)]
    (let [f (ic/cubic-2d (:xs d) (:ys d) (:vss d))
          sc (vmax (:vss d))]
      (doseq [[[x y] e] (map vector (:pts d) (:expected d))]
        (t/is (m/delta-eq e (f x y) (* 1e-12 sc) 1e-12) (str (:name d) " at " [x y]))
        (t/is (== (f x y) (f [x y])) "both arities give the same value")))))

(t/deftest grid-edge-cases
  (let [xs [0.0 1.0 3.0]
        ys [-1.0 2.0]
        vss [[1.0 2.0] [4.0 -3.0] [0.5 0.25]]]
    (doseq [ctor [il/bilinear ic/cubic-2d]]
      (let [f (ctor xs ys vss)]
        (t/testing "nodes return the data"
          (doseq [[i x] (map-indexed vector xs)
                  [j y] (map-indexed vector ys)]
            (t/is (m/delta-eq (get-in vss [i j]) (f x y) 1e-13))))))
    (t/testing "bilinear is exact for a bilinear function everywhere (also outside)"
      (let [g (fn ^double [^double x ^double y] (+ 1.0 (* 2.0 x) (* -3.0 y) (* 0.5 x y)))
            f (il/bilinear xs ys (map (fn [x] (map #(g x %) ys)) xs))]
        (doseq [x [-5.0 -0.1 0.0 0.7 2.0 3.0 9.0]
                y [-7.0 -1.0 0.3 2.0 5.0]]
          (t/is (m/delta-eq (g x y) (f x y) 1e-11)))))
    (t/testing "cubic-2d reproduces a bilinear function (natural splines are exact for linear data)"
      (let [g (fn ^double [^double x ^double y] (+ 1.0 (* 2.0 x) (* -3.0 y) (* 0.5 x y)))
            f (ic/cubic-2d xs ys (map (fn [x] (map #(g x %) ys)) xs))]
        (doseq [x [-5.0 -0.1 0.0 0.7 2.0 3.0 9.0]
                y [-7.0 -1.0 0.3 2.0 5.0]]
          (t/is (m/delta-eq (g x y) (f x y) 1e-11)))))
    (t/testing "NaN query gives NaN"
      (t/is (m/nan? ((il/bilinear xs ys vss) ##NaN 0.5)))
      (t/is (m/nan? ((ic/cubic-2d xs ys vss) 0.5 ##NaN))))
    (t/testing "the minimal 2x2 grid"
      (let [f (il/bilinear [0 1] [0 1] [[0 1] [2 3]])
            g (ic/cubic-2d [0 1] [0 1] [[0 1] [2 3]])]
        (t/is (m/delta-eq 1.5 (f 0.5 0.5) 1e-15))
        (t/is (m/delta-eq 1.5 (g 0.5 0.5) 1e-15))
        (t/is (m/delta-eq 5.5 (f 2 1.5) 1e-14))
        (t/is (m/delta-eq 5.5 (g 2 1.5) 1e-14))))
    (t/testing "sequence types: arrays, lists, double-array of double-arrays"
      (let [f (il/bilinear (double-array xs) (list -1.0 2.0) (into-array (map double-array vss)))]
        (t/is (m/delta-eq (f 0.5 0.5) ((il/bilinear xs ys vss) 0.5 0.5) 1e-15))))))

(t/deftest grid-validation
  (let [xs [0.0 1.0 3.0]
        ys [-1.0 2.0]
        vss [[1.0 2.0] [4.0 -3.0] [0.5 0.25]]]
    (doseq [[label ctor] [[:bilinear il/bilinear] [:cubic-2d ic/cubic-2d]]]
      (t/testing (name label)
        (t/is (nil? (property #(ctor xs ys vss))))
        (t/testing "axes: too few points"
          (t/is (= :too-few-points (property #(ctor [0.0] ys [[1.0 2.0]]))))
          (t/is (= :too-few-points (property #(ctor xs [0.0] [[1.0] [2.0] [3.0]]))))
          (t/is (= :too-few-points (property #(ctor [] ys []))))
          (t/is (= :too-few-points (property #(ctor xs [] [[] [] []])))))
        (t/testing "axes: unsorted xs or ys gave silently wrong values (1.80 vs 1.38, 0.658 vs 1.361)"
          (t/is (= :not-increasing (property #(ctor [0.0 3.0 1.0] ys vss))))
          (t/is (= :not-increasing (property #(ctor xs [2.0 -1.0] vss))))
          (t/is (= :not-increasing (property #(ctor [3.0 1.0 0.0] ys vss)))))
        (t/testing "axes: duplicates gave NaN"
          (t/is (= :not-increasing (property #(ctor [0.0 1.0 1.0] ys vss))))
          (t/is (= :not-increasing (property #(ctor xs [2.0 2.0] vss)))))
        (t/testing "axes: not finite"
          (t/is (= :not-finite (property #(ctor [0.0 ##NaN 3.0] ys vss))))
          (t/is (= :not-finite (property #(ctor xs [-1.0 ##Inf] vss)))))
        (t/testing "shape of vss: rows = count xs, every row has count ys values"
          (t/is (= :shape-mismatch (property #(ctor xs ys [[1.0 2.0] [4.0 -3.0]]))) "missing row")
          (t/is (= :shape-mismatch (property #(ctor xs ys (conj vss [1.0 2.0])))) "extra row")
          (t/is (= :shape-mismatch (property #(ctor xs ys [[1.0 2.0] [4.0] [0.5 0.25]]))) "short row")
          (t/is (= :shape-mismatch (property #(ctor xs ys [[1.0 2.0] [4.0 -3.0 7.0] [0.5 0.25]]))) "long row")
          (t/is (= :shape-mismatch (property #(ctor xs ys [[1.0 2.0 3.0] [4.0 -3.0 7.0] [0.5 0.25 1.0]]))) "transposed")
          (t/is (= :shape-mismatch (property #(ctor xs ys []))) "empty")
          (t/is (= :shape-mismatch (property #(ctor xs ys nil))) "nil"))
        (t/testing "ex-data names the method"
          (t/is (= label (:method (try (ctor xs ys [[1.0 2.0]]) (catch clojure.lang.ExceptionInfo e (ex-data e)))))))))))

(t/deftest grid-docstrings
  (t/is (re-find #"extrapolates linearly" (:doc (meta #'il/bilinear))))
  (t/is (re-find #"RegularGridInterpolator" (:doc (meta #'il/bilinear))))
  (t/is (re-find #"natural" (:doc (meta #'ic/cubic-2d))))
  (t/is (re-find #"not-a-knot" (:doc (meta #'ic/cubic-2d))))
  (t/is (re-find #"RectBivariateSpline" (:doc (meta #'ic/cubic-2d))))
  (t/is (re-find #"(?i)extrapolat" (:doc (meta #'ic/cubic-2d)))))
