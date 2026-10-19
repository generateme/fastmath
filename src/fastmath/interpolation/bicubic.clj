(ns fastmath.interpolation.bicubic
  "Bicubic Hermite interpolation on a rectilinear grid.

  Every grid cell is a C1 bicubic patch defined by the values, the first partial derivatives and the mixed derivative at its four corners. The derivatives are estimated from the grid values by finite differences."
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(deftype Grid [^doubles xs ^doubles ys
               ^doubles f ^doubles fx ^doubles fy ^doubles fxy])

;; Derivative stencils: node i uses three (index, weight) pairs stored at 3i, 3i+1, 3i+2.
;; Unused pairs have weight 0 and point to a node of the same cell.

(defn- stencils
  "Returns `[idx wts]`, the finite difference stencils of the first derivative at every node of the strictly increasing `xs`.

  Interior nodes: central difference of the neighbours `(f[i+1] - f[i-1]) / (x[i+1] - x[i-1])`, the rule of Apache Commons Math `BicubicInterpolator`. Boundary nodes: second order one-sided three point formula on nonuniform spacing. An axis with two nodes: the secant."
  [^doubles xs]
  (let [n (alength xs)
        idx (long-array (m/* 3 n))
        wts (double-array (m/* 3 n))
        put! (fn [^long i ^long s ^long j ^double w]
               (Array/aset idx (m/+ (m/* 3 i) s) j)
               (Array/aset wts (m/+ (m/* 3 i) s) w))]
    (if (m/== n 2)
      (let [h (m/- (Array/aget xs 1) (Array/aget xs 0))]
        (dotimes [i 2]
          (put! i 0 0 (m// -1.0 h))
          (put! i 1 1 (m// 1.0 h))
          (put! i 2 1 0.0)))
      (do
        (let [h1 (m/- (Array/aget xs 1) (Array/aget xs 0))
              h2 (m/- (Array/aget xs 2) (Array/aget xs 1))
              h12 (m/+ h1 h2)]
          (put! 0 0 0 (m/- (m// (m/+ (m/* 2.0 h1) h2) (m/* h1 h12))))
          (put! 0 1 1 (m// h12 (m/* h1 h2)))
          (put! 0 2 2 (m/- (m// h1 (m/* h2 h12)))))
        (dotimes [c (m/- n 2)]
          (let [i (m/inc c)
                dx (m/- (Array/aget xs (m/inc i)) (Array/aget xs (m/dec i)))]
            (put! i 0 (m/dec i) (m// -1.0 dx))
            (put! i 1 i 0.0)
            (put! i 2 (m/inc i) (m// 1.0 dx))))
        (let [l (m/dec n)
              h1 (m/- (Array/aget xs l) (Array/aget xs (m/dec l)))
              h2 (m/- (Array/aget xs (m/dec l)) (Array/aget xs (m/- l 2)))
              h12 (m/+ h1 h2)]
          (put! l 0 (m/- l 2) (m// h1 (m/* h2 h12)))
          (put! l 1 (m/dec l) (m/- (m// h12 (m/* h1 h2))))
          (put! l 2 l (m// (m/+ (m/* 2.0 h1) h2) (m/* h1 h12))))))
    [idx wts]))

(defn- derivatives
  "Returns `[fx fy fxy]`, flat arrays (row `i`, column `j` at `i*ny + j`) of the partial derivatives estimated at all nodes. The mixed derivative is the product of both stencils."
  [^doubles f ^long ny xstencil ystencil]
  (let [nx (m/long-quot (alength f) ny)
        [^longs xi ^doubles xw] xstencil
        [^longs yi ^doubles yw] ystencil
        fx (double-array (m/* nx ny))
        fy (double-array (m/* nx ny))
        fxy (double-array (m/* nx ny))]
    (dotimes [i nx]
      (dotimes [j ny]
        (let [pos (m/+ (m/* i ny) j)]
          (dotimes [s 3]
            (let [a (Array/aget xi (m/+ (m/* 3 i) s))
                  wa (Array/aget xw (m/+ (m/* 3 i) s))
                  b (Array/aget yi (m/+ (m/* 3 j) s))
                  wb (Array/aget yw (m/+ (m/* 3 j) s))]
              (Array/add fx pos (m/* wa (Array/aget f (m/+ (m/* a ny) j))))
              (Array/add fy pos (m/* wb (Array/aget f (m/+ (m/* i ny) b))))
              (dotimes [r 3]
                (let [b2 (Array/aget yi (m/+ (m/* 3 j) r))
                      wb2 (Array/aget yw (m/+ (m/* 3 j) r))]
                  (Array/add fxy pos (m/* wa wb2 (Array/aget f (m/+ (m/* a ny) b2)))))))))))
    [fx fy fxy]))

(defn- check-shape!
  "Checks that `vss` has a row for every x node and a value for every y node in every row."
  [vss ^long nx ^long ny]
  (when-not (m/== (count vss) nx)
    (ic/invalid-input! :bicubic :length-mismatch
                       (str "vss needs one row per x node: " nx " rows expected, got " (count vss))
                       {:axis :xs :expected nx :count (count vss)}))
  (loop [i 0
         rows (seq vss)]
    (when rows
      (let [c (count (first rows))]
        (when-not (m/== c ny)
          (ic/invalid-input! :bicubic :length-mismatch
                             (str "every row of vss needs one value per y node: " ny " values expected, got " c " in row " i)
                             {:axis :ys :index i :expected ny :count c})))
      (recur (m/inc i) (next rows)))))

(defn- flatten-values
  ^doubles [vss ^long nx ^long ny]
  (let [f (double-array (m/* nx ny))]
    (loop [i 0
           rows (seq vss)]
      (when rows
        (let [^doubles row (m/seq->double-array (first rows))]
          (System/arraycopy row 0 f (m/* i ny) ny))
        (recur (m/inc i) (next rows))))
    f))

(defn- grid-value
  "Value of the patch of the cell containing (or nearest to) `x`, `y`."
  ^double [^Grid g ^double x ^double y]
  (let [^doubles xs (.xs g)
        ^doubles ys (.ys g)
        ny (alength ys)
        i (ic/binary-search xs x)
        j (ic/binary-search ys y)
        hx (m/- (Array/aget xs (m/inc i)) (Array/aget xs i))
        hy (m/- (Array/aget ys (m/inc j)) (Array/aget ys j))
        u (m// (m/- x (Array/aget xs i)) hx)
        v (m// (m/- y (Array/aget ys j)) hy)
        u2 (m/* u u)
        u3 (m/* u2 u)
        v2 (m/* v v)
        v3 (m/* v2 v)
        ;; Hermite basis: values at the corners 0 and 1 and derivatives (scaled by the cell size)
        a0 (m/+ (m/- (m/* 2.0 u3) (m/* 3.0 u2)) 1.0)
        a1 (m/- (m/* 3.0 u2) (m/* 2.0 u3))
        da0 (m/* hx (m/+ (m/- u3 (m/* 2.0 u2)) u))
        da1 (m/* hx (m/- u3 u2))
        c0 (m/+ (m/- (m/* 2.0 v3) (m/* 3.0 v2)) 1.0)
        c1 (m/- (m/* 3.0 v2) (m/* 2.0 v3))
        dc0 (m/* hy (m/+ (m/- v3 (m/* 2.0 v2)) v))
        dc1 (m/* hy (m/- v3 v2))
        p00 (long (m/+ (m/* i ny) j))
        p01 (long (m/inc p00))
        p10 (long (m/+ p00 ny))
        p11 (long (m/inc p10))
        ^doubles f (.f g)
        ^doubles fx (.fx g)
        ^doubles fy (.fy g)
        ^doubles fxy (.fxy g)]
    (m/+ (m/+ (m/* a0 (m/+ (m/+ (m/* c0 (Array/aget f p00)) (m/* c1 (Array/aget f p01)))
                           (m/+ (m/* dc0 (Array/aget fy p00)) (m/* dc1 (Array/aget fy p01)))))
              (m/* a1 (m/+ (m/+ (m/* c0 (Array/aget f p10)) (m/* c1 (Array/aget f p11)))
                           (m/+ (m/* dc0 (Array/aget fy p10)) (m/* dc1 (Array/aget fy p11))))))
         (m/+ (m/* da0 (m/+ (m/+ (m/* c0 (Array/aget fx p00)) (m/* c1 (Array/aget fx p01)))
                            (m/+ (m/* dc0 (Array/aget fxy p00)) (m/* dc1 (Array/aget fxy p01)))))
              (m/* da1 (m/+ (m/+ (m/* c0 (Array/aget fx p10)) (m/* c1 (Array/aget fx p11)))
                            (m/+ (m/* dc0 (Array/aget fxy p10)) (m/* dc1 (Array/aget fxy p11)))))))))

(defn bicubic
  "Creates a C1 bicubic Hermite interpolator of values given on a rectilinear grid.

  The surface passes through all grid values. Every cell is a bicubic patch built from the values, the first partial derivatives `fx`, `fy` and the mixed derivative `fxy` at its corners. Derivatives are estimated by finite differences of the grid values (nonuniform spacing is supported):

  - interior nodes: central differences between the two neighbours, `fx = (f[i+1] - f[i-1]) / (x[i+1] - x[i-1])` and the same for `fy`; `fxy` is the product of both differences. This is the rule of Apache Commons Math `BicubicInterpolator`.
  - boundary nodes: second order one-sided three point differences on the actual spacing, instead of zero derivatives, so the cells at the border are as accurate as the inner ones.
  - an axis with two nodes: the secant of the two values.

  The result is continuous with continuous first derivatives (C1) across the cell borders, the second derivatives jump. It reproduces exactly every function `a + b x + c y + d x y` on all cells including the border cells and the extrapolation.

  Parameters:

  - `xs` (sequence of numbers): nodes of the first axis, at least `2`, strictly increasing and finite
  - `ys` (sequence of numbers): nodes of the second axis, at least `2`, strictly increasing and finite
  - `vss` (sequence of sequences of numbers): values, `vss[i][j]` is the value at `(xs[i], ys[j])`; `count xs` rows of `count ys` values each. `NaN` or infinite values make the cells around the node `NaN`

  Returns a function of two arguments `(f x y)` or one vector argument `(f [x y])`. Outside the grid the patch of the nearest border cell is continued (extrapolation, no exception, which differs from `OutOfRangeException` of the Apache implementation). `NaN` coordinates give `NaN`.

  A grid of `2 x 2` points gives the bilinear interpolation. Throws `ex-info` (`:type :invalid-input`) when an axis has less than two nodes, is not finite or not strictly increasing, or when `vss` does not have the shape `count xs` times `count ys`.

  See also [[fastmath.interpolation.cubic/cubic-2d]] (natural tensor spline, C2), [[fastmath.interpolation/bilinear]]."
  [xs ys vss]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)]
    (ic/check-axis! :bicubic :xs xs {:min-points 2 :sorted? true})
    (ic/check-axis! :bicubic :ys ys {:min-points 2 :sorted? true})
    (let [nx (alength xs)
          ny (alength ys)]
      (check-shape! vss nx ny)
      (let [f (flatten-values vss nx ny)
            [fx fy fxy] (derivatives f ny (stencils xs) (stencils ys))
            g (Grid. xs ys f fx fy fxy)]
        (fn bicubic-interpolator
          (^double [[^double x ^double y]] (grid-value g x y))
          (^double [^double x ^double y] (grid-value g x y)))))))
