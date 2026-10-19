(ns fastmath.interpolation.bicubic-test
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :as t]
            [fastmath.core :as m]
            [fastmath.interpolation.bicubic :as sut]
            [fastmath.interpolation.cubic :as cubic]))

(set! *warn-on-reflection* true)

;; Hermite patch clone (pure Python) with the same derivative rule, see
;; utils/fastmath/dev/generate_bspline_bicubic_reference.py
(def ^:private reference
  (delay (:bicubic (edn/read-string (slurp (io/resource "interpolation/bspline_bicubic_reference.edn"))))))

(defn- property
  "Returns the :property of the ex-info thrown by `f`, nil when it does not throw."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e
         (let [d (ex-data e)]
           (when (and (= :invalid-input (:type d)) (= :bicubic (:method d))) (:property d))))))

(defn- target
  "Test function sin(1.3 x) cos(0.9 y) + 0.2 x y"
  ^double [^double x ^double y]
  (+ (* (m/sin (* 1.3 x)) (m/cos (* 0.9 y))) (* 0.2 x y)))

(defn- values
  [f xs ys]
  (mapv (fn [x] (mapv (fn [y] (f x y)) ys)) xs))

(defn- linspace
  [^double a ^double b ^long n]
  (mapv (fn [i] (+ a (* (- b a) (/ (double i) (double (dec n)))))) (range n)))

(def ^:private grids
  {:uni9x9 [(linspace 0.0 4.0 9) (linspace 0.0 4.0 9)]
   :uni5x5 [(linspace 0.0 4.0 5) (linspace 0.0 4.0 5)]
   :nonuni7x6 [[0.0 0.35 0.9 1.7 2.4 3.3 4.0] [0.0 0.6 1.3 2.5 3.1 4.0]]})

(defn- cell-errors
  "Maximal errors of `f` against the target over 12 x 12 sample points of every cell: `[ring inner]`,
  ring = the cells touching the border of the grid."
  [f xs ys]
  (let [nx (count xs)
        ny (count ys)
        errs (for [i (range (dec nx))
                   j (range (dec ny))
                   a (range 12)
                   b (range 12)
                   :let [x (+ (double (xs i)) (* (/ a 11.0) (- (double (xs (inc i))) (double (xs i)))))
                         y (+ (double (ys j)) (* (/ b 11.0) (- (double (ys (inc j))) (double (ys j)))))]]
               [(or (zero? i) (zero? j) (= i (- nx 2)) (= j (- ny 2)))
                (m/abs (- (double (f x y)) (target x y)))])
        mx (fn [ring?] (reduce max 0.0 (map second (filter #(= ring? (first %)) errs))))]
    [(mx true) (mx false)]))

(t/deftest hermite-patch-clone
  ;; independent pure python Hermite patch with the same derivative rule; measured max difference
  ;; 6e-14 (uni9x9), 2e-15, 1.5e-14, 2e-16, 7e-16; points inside, outside (extrapolation) and on nodes
  (doseq [{:keys [name xs ys vss points expected]} @reference
          :let [f (sut/bicubic xs ys vss)]]
    (t/testing name
      (t/is (every? true? (map (fn [[x y] e] (m/delta-eq (f x y) e 1.0e-12 1.0e-12)) points expected)))
      (t/is (every? true? (map (fn [[x y] e] (m/delta-eq (f [x y]) e 1.0e-12 1.0e-12)) points expected))
            "vector argument"))))

(t/deftest outer-ring-accuracy
  ;; I-02: Apache bicubic had zero boundary derivatives. Measured max error over 12x12 samples per cell
  ;; (bicubic ring / cubic-2d ring / Apache ring):
  ;;   uni9x9   0.0145 / 0.0288 / 0.127
  ;;   uni5x5   0.0989 / 0.130  / 0.302
  ;;   nonuni7x6 0.0629 / 0.0692 / 0.230
  ;; inner cells (the same as Apache): 5.37e-3, 5.04e-2, 4.18e-2 (cubic-2d 6.8e-3, 2.1e-2, 1.3e-2)
  (doseq [[k [xs ys]] grids
          :let [vss (values target xs ys)
                [ring inner] (cell-errors (sut/bicubic xs ys vss) xs ys)
                [cring _] (cell-errors (cubic/cubic-2d xs ys vss) xs ys)
                inner-limit ({:uni9x9 5.4e-3 :uni5x5 5.1e-2 :nonuni7x6 4.2e-2} k)]]
    (t/testing (name k)
      (t/is (<= ring cring) (str "ring " ring " cubic-2d " cring))
      (t/is (<= inner inner-limit) (str "inner " inner)))))

(t/deftest passes-through-the-nodes
  (doseq [[k [xs ys]] grids
          :let [f (sut/bicubic xs ys (values target xs ys))]
          x xs
          y ys]
    (t/is (m/delta-eq (target x y) (f x y) 1.0e-13) (str (name k) " " x " " y))))

(defn- bilinear-fn
  ^double [^double x ^double y]
  (+ 0.7 (* -1.3 x) (* 2.1 y) (* 0.45 x y)))

(t/deftest exact-for-bilinear-data
  ;; a + b x + c y + d x y is reproduced on every cell: edge cells, outside of the grid, every shape
  (doseq [[xs ys] [[[0.0 1.5] [0.2 2.0]]
                   [[0.0 1.5] [0.2 0.9 2.0 2.1 5.0]]
                   [[0.0 0.7 2.0] [0.0 0.5 1.5 2.2]]
                   [[-3.0 -2.9 0.1 0.15 4.0 10.0] [1.0 1.001 7.0]]
                   (grids :nonuni7x6)
                   (grids :uni9x9)]
          :let [f (sut/bicubic xs ys (values bilinear-fn xs ys))
                ;; query range: half of the first/last cell outside of the grid on both sides
                lo-x (- (double (first xs)) (* 0.5 (- (double (xs 1)) (double (first xs)))))
                hi-x (+ (double (peek xs)) (* 0.5 (- (double (peek xs)) (double (xs (- (count xs) 2))))))
                lo-y (- (double (first ys)) (* 0.5 (- (double (ys 1)) (double (first ys)))))
                hi-y (+ (double (peek ys)) (* 0.5 (- (double (peek ys)) (double (ys (- (count ys) 2))))))]
          px (linspace lo-x hi-x 12)
          py (linspace lo-y hi-y 12)]
    (t/is (m/delta-eq (bilinear-fn px py) (f px py) 1.0e-11 1.0e-11)
          (str (count xs) "x" (count ys) " at " px " " py))))

(t/deftest two-by-two-grid-is-bilinear
  ;; documented: with two nodes per axis the derivatives are the secants and the patch is the bilinear interpolation
  (let [xs [1.0 3.5] ys [-2.0 0.5]
        v00 0.3 v01 -1.2 v10 2.2 v11 0.9
        f (sut/bicubic xs ys [[v00 v01] [v10 v11]])
        bil (fn [^double x ^double y]
              (let [u (/ (- x 1.0) 2.5) v (/ (+ y 2.0) 2.5)]
                (+ (* (- 1.0 u) (- 1.0 v) v00) (* (- 1.0 u) v v01) (* u (- 1.0 v) v10) (* u v v11))))]
    (doseq [x [1.0 1.4 2.2 3.5 0.0 5.0] y [-2.0 -1.1 0.5 -3.0 1.7]]
      (t/is (m/delta-eq (bil x y) (f x y) 1.0e-13) (str x " " y)))))

(t/deftest c1-continuity
  ;; value and first derivatives continuous across the cell borders (one-sided differences, h = 1e-6;
  ;; the jump is O(h f'') = 1e-6 at most, tolerance 1e-4)
  (doseq [[k [xs ys]] grids
          :let [f (sut/bicubic xs ys (values target xs ys))
                h 1.0e-6
                slope-x (fn [x y] [(/ (- (f x y) (f (- x h) y)) h) (/ (- (f (+ x h) y) (f x y)) h)])
                slope-y (fn [x y] [(/ (- (f x y) (f x (- y h))) h) (/ (- (f x (+ y h)) (f x y)) h)])]]
    (doseq [x (rest (butlast xs))
            y [0.123 0.97 2.2 3.71 (second ys)]]
      (t/is (m/delta-eq (f (- x 1.0e-9) y) (f (+ x 1.0e-9) y) 1.0e-8) (str (name k) " value across x " x))
      (let [[l r] (slope-x x y)]
        (t/is (m/delta-eq l r 1.0e-4) (str (name k) " fx across x " x " at y " y))))
    (doseq [y (rest (butlast ys))
            x [0.123 0.97 2.2 3.71 (second xs)]]
      (t/is (m/delta-eq (f x (- y 1.0e-9)) (f x (+ y 1.0e-9)) 1.0e-8) (str (name k) " value across y " y))
      (let [[l r] (slope-y x y)]
        (t/is (m/delta-eq l r 1.0e-4) (str (name k) " fy across y " y " at x " x))))))

(t/deftest second-derivative-jumps
  ;; the surface is C1 and not C2: the second derivative has a jump at interior nodes (documented)
  (let [[xs ys] (grids :uni9x9)
        f (sut/bicubic xs ys (values target xs ys))
        h 1.0e-3
        y 1.3
        x (double (xs 3))
        d2 (fn [x0] (/ (+ (f (- x0 h) y) (f (+ x0 h) y) (* -2.0 (f x0 y))) (* h h)))]
    (t/is (> (m/abs (- (double (d2 (- x h))) (double (d2 (+ x h))))) 1.0e-3))))

(t/deftest symmetry-of-axes
  (let [[xs ys] (grids :nonuni7x6)
        vss (values target xs ys)
        vst (apply mapv vector vss)
        f (sut/bicubic xs ys vss)
        g (sut/bicubic ys xs vst)]
    (doseq [x [-0.5 0.0 0.31 1.7 2.9 4.0 4.5] y [-0.2 0.0 0.8 2.5 3.3 4.0 4.9]]
      (t/is (m/delta-eq (f x y) (g y x) 1.0e-13) (str x " " y)))))

(t/deftest extrapolation-continues-the-edge-patch
  ;; outside of the grid: the patch of the nearest border cell, no exception; the value at the
  ;; boundary line is continuous with the edge
  (let [[xs ys] (grids :nonuni7x6)
        f (sut/bicubic xs ys (values target xs ys))]
    (doseq [[x y] [[-1.0 2.0] [5.0 2.0] [2.0 -1.0] [2.0 5.0] [-1.0 -1.0] [5.0 5.0] [-1.0 5.0] [5.0 -1.0]]]
      (t/is (m/valid-double? (f x y)) (str x " " y)))
    (t/is (m/delta-eq (f -1.0e-9 2.0) (f 0.0 2.0) 1.0e-7))
    (t/is (m/delta-eq (f 4.0000000001 2.0) (f 4.0 2.0) 1.0e-7))
    ;; polynomial continuation: exact bilinear data is exact outside too
    ;; (far away the powers of the local coordinate amplify rounding: u up to 10 here)
    (let [b (sut/bicubic xs ys (values bilinear-fn xs ys))]
      (t/is (m/delta-eq (bilinear-fn -3.5 8.0) (b -3.5 8.0) 1.0e-8 1.0e-10)))))

(t/deftest huge-and-tiny-values
  ;; scaling of the values and of the coordinates is a symmetry
  (let [[xs ys] (grids :nonuni7x6)
        vss (values target xs ys)
        base (sut/bicubic xs ys vss)]
    (doseq [c [1.0e-9 1.0e9]
            :let [fv (sut/bicubic xs ys (mapv (fn [r] (mapv #(* c (double %)) r)) vss))
                  fc (sut/bicubic (mapv #(* c (double %)) xs) (mapv #(* c (double %)) ys) vss)]
            [x y] [[0.3 0.2] [1.7 2.5] [3.9 3.3] [-0.5 4.5]]]
      (t/is (m/delta-eq (* c (double (base x y))) (fv x y) 1.0e-12 1.0e-12) (str "values " c))
      (t/is (m/delta-eq (base x y) (fc (* c (double x)) (* c (double y))) 1.0e-9 1.0e-9) (str "coordinates " c)))))

(t/deftest nan-and-infinite
  (let [xs (linspace 0.0 5.0 6)
        ys (linspace 0.0 5.0 6)
        vss (assoc-in (values target xs ys) [0 0] ##NaN)
        f (sut/bicubic xs ys vss)]
    (t/is (m/nan? (f 0.5 0.5)) "cell with the NaN node")
    (t/is (m/valid-double? (f 4.5 4.5)) "far cell is not affected")
    (t/is (m/nan? ((sut/bicubic xs ys (assoc-in (values target xs ys) [3 3] ##Inf)) 3.5 3.5))))
  (let [f (sut/bicubic [0.0 1.0 2.0] [0.0 1.0 2.0] (values target [0.0 1.0 2.0] [0.0 1.0 2.0]))]
    (t/is (m/nan? (f ##NaN 1.0)))
    (t/is (m/nan? (f 1.0 ##NaN)))
    (t/is (m/nan? (f [##NaN ##NaN])))))

(t/deftest invalid-input
  (let [ok-v [[1.0 2.0 3.0] [2.0 3.0 4.0] [3.0 4.0 5.0]]
        ok-x [0.0 1.0 2.0]]
    (t/is (nil? (property #(sut/bicubic ok-x ok-x ok-v))))
    (t/testing "too few points on an axis"
      (t/is (= :too-few-points (property #(sut/bicubic [] ok-x ok-v))))
      (t/is (= :too-few-points (property #(sut/bicubic [0.0] ok-x [[1.0 2.0 3.0]]))))
      (t/is (= :too-few-points (property #(sut/bicubic ok-x [] [[] [] []]))))
      (t/is (= :too-few-points (property #(sut/bicubic ok-x [1.0] [[1.0] [2.0] [3.0]]))))
      (t/is (nil? (property #(sut/bicubic [0.0 1.0] [0.0 1.0] [[1.0 2.0] [3.0 4.0]])))))
    (t/testing "not sorted, duplicates, not finite"
      (t/is (= :not-increasing (property #(sut/bicubic [0.0 2.0 1.0] ok-x ok-v))))
      (t/is (= :not-increasing (property #(sut/bicubic ok-x [0.0 1.0 1.0] ok-v))))
      (t/is (= :not-increasing (property #(sut/bicubic [2.0 1.0 0.0] ok-x ok-v))))
      (t/is (= :not-increasing (property #(sut/bicubic [1.0 1.0] ok-x [[1.0 2.0 3.0] [2.0 3.0 4.0]]))))
      (t/is (= :not-finite (property #(sut/bicubic [0.0 ##NaN 2.0] ok-x ok-v))))
      (t/is (= :not-finite (property #(sut/bicubic ok-x [0.0 1.0 ##Inf] ok-v)))))
    (t/testing "shape of the values"
      (t/is (= :length-mismatch (property #(sut/bicubic ok-x ok-x (butlast ok-v)))))
      (t/is (= :length-mismatch (property #(sut/bicubic ok-x ok-x (conj ok-v [1.0 2.0 3.0])))))
      (t/is (= :length-mismatch (property #(sut/bicubic ok-x ok-x []))))
      (t/is (= :length-mismatch (property #(sut/bicubic ok-x ok-x [[1.0 2.0 3.0] [2.0 3.0] [3.0 4.0 5.0]]))))
      (t/is (= :length-mismatch (property #(sut/bicubic ok-x ok-x [[1.0 2.0 3.0] [2.0 3.0 4.0 9.0] [3.0 4.0 5.0]]))))
      (t/is (= :length-mismatch (property #(sut/bicubic ok-x ok-x [[] [] []]))))
      ;; transposed shape of a non square grid
      (t/is (= :length-mismatch (property #(sut/bicubic [0.0 1.0] [0.0 1.0 2.0] [[1.0 2.0] [3.0 4.0] [5.0 6.0]])))))
    (t/testing "ex-data"
      (let [d (try (sut/bicubic [0.0 1.0 1.0] ok-x ok-v) nil
                   (catch clojure.lang.ExceptionInfo e (ex-data e)))]
        (t/is (= {:type :invalid-input :method :bicubic :property :not-increasing :axis :xs :index 2 :value 1.0} d))))))

(t/deftest large-grid
  (let [xs (linspace 0.0 10.0 300)
        ys (linspace -3.0 3.0 250)
        f (sut/bicubic xs ys (values target xs ys))]
    (t/is (m/delta-eq (target 3.3333 1.1111) (f 3.3333 1.1111) 1.0e-5))
    (t/is (m/delta-eq (target (xs 150) (ys 100)) (f (xs 150) (ys 100)) 1.0e-12))))

(t/deftest keeps-the-input-sequences
  ;; lazy sequences and vectors of vectors, seqs of seqs
  (let [xs (map double (range 4))
        ys (map double (range 5))
        vss (map (fn [x] (map (fn [y] (target x y)) ys)) xs)
        f (sut/bicubic xs ys vss)]
    (t/is (m/delta-eq (target 2.0 3.0) (f 2.0 3.0) 1.0e-13))
    (t/is (m/delta-eq (f 1.3 2.2) (f [1.3 2.2]) 0.0))))
