(ns fastmath.interpolation.monotone
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- secants
  "Slopes of the segments between consecutive nodes."
  ^doubles [^doubles xs ^doubles ys]
  (let [k (m/dec (alength xs))
        d (double-array k)]
    (dotimes [i k]
      (Array/aset d i (m// (m/- (Array/aget ys (m/inc i)) (Array/aget ys i))
                           (m/- (Array/aget xs (m/inc i)) (Array/aget xs i)))))
    d))

(defn- opposite-or-zero?
  "True when two neighbouring secants have different signs or one of them is zero."
  [^double a ^double b]
  (or (m/zero? a) (m/zero? b) (not= (m/pos? a) (m/pos? b))))

(defn- fc-tangents
  "Fritsch-Carlson tangents for the secants `d`.

  Interior tangents are arithmetic means of the neighbouring secants and are set to zero when the secants change
  sign or one of them is zero. End tangents are the secants of the end segments. Then every segment is limited:
  the pair (alpha, beta) of tangent to secant ratios is scaled back to the circle of radius 3."
  ^doubles [^doubles d]
  (let [k (alength d)
        ms (double-array (m/inc k))]
    (Array/aset ms 0 (Array/aget d 0))
    (Array/aset ms k (Array/aget d (m/dec k)))
    (dotimes [i (m/dec k)]
      (let [a (Array/aget d i)
            b (Array/aget d (m/inc i))]
        (Array/aset ms (m/inc i) (if (opposite-or-zero? a b) 0.0 (m/* 0.5 (m/+ a b))))))
    (dotimes [i k]
      (let [di (Array/aget d i)]
        (if (m/zero? di)
          (do (Array/aset ms i 0.0)
              (Array/aset ms (m/inc i) 0.0))
          (let [alpha (m// (Array/aget ms i) di)
                beta (m// (Array/aget ms (m/inc i)) di)
                h (Math/hypot alpha beta)]
            (when (m/> h 3.0)
              (let [tau (m// 3.0 h)]
                (Array/aset ms i (m/* tau alpha di))
                (Array/aset ms (m/inc i) (m/* tau beta di))))))))
    ms))

(defn- pchip-edge
  "Non-centered, shape-preserving three point end tangent (h0 is the end segment width, m0 its secant)."
  ^double [^double h0 ^double h1 ^double m0 ^double m1]
  (let [d (m// (m/- (m/* (m/+ (m/* 2.0 h0) h1) m0) (m/* h0 m1)) (m/+ h0 h1))
        sd (Math/signum d)
        s0 (Math/signum m0)]
    (cond
      (m/not== sd s0) 0.0
      (and (m/not== s0 (Math/signum m1)) (m/> (m/abs d) (m/* 3.0 (m/abs m0)))) (m/* 3.0 m0)
      :else d)))

(defn- pchip-tangents
  "Fritsch-Butland tangents (weighted harmonic mean of the neighbouring secants, zero where the secants change
  sign or vanish) with the three point end formulas."
  ^doubles [^doubles xs ^doubles d]
  (let [k (alength d)
        ms (double-array (m/inc k))]
    (if (m/== k 1)
      (do (Array/aset ms 0 (Array/aget d 0))
          (Array/aset ms 1 (Array/aget d 0)))
      (do
        (dotimes [i (m/dec k)]
          (let [a (Array/aget d i)
                b (Array/aget d (m/inc i))]
            (when-not (opposite-or-zero? a b)
              (let [h0 (m/- (Array/aget xs (m/inc i)) (Array/aget xs i))
                    h1 (m/- (Array/aget xs (m/+ i 2)) (Array/aget xs (m/inc i)))
                    w1 (m/+ (m/* 2.0 h1) h0)
                    w2 (m/+ h1 (m/* 2.0 h0))
                    whmean (m// (m/+ (m// w1 a) (m// w2 b)) (m/+ w1 w2))]
                (Array/aset ms (m/inc i) (m// 1.0 whmean))))))
        (let [h0 (m/- (Array/aget xs 1) (Array/aget xs 0))
              h1 (m/- (Array/aget xs 2) (Array/aget xs 1))]
          (Array/aset ms 0 (pchip-edge h0 h1 (Array/aget d 0) (Array/aget d 1))))
        (let [h0 (m/- (Array/aget xs k) (Array/aget xs (m/dec k)))
              h1 (m/- (Array/aget xs (m/dec k)) (Array/aget xs (m/- k 2)))]
          (Array/aset ms k (pchip-edge h0 h1 (Array/aget d (m/dec k)) (Array/aget d (m/- k 2)))))))
    ms))

(defn- hermite-interpolator
  "Piecewise cubic Hermite function with node values `ys` and node tangents `ms`."
  [^doubles xs ^doubles ys ^doubles ms]
  (let [n (alength xs)]
    (fn ^double [^double v]
      (let [b (java.util.Arrays/binarySearch xs v)]
        (if (m/not-neg? b)
          (Array/aget ys b)
          (let [i (ic/binary-search-id b n)
                mxi (Array/aget xs i)
                h (m/- (Array/aget xs (m/inc i)) mxi)
                t (m// (m/- v mxi) h)]
            (m/+ (m/* (m/sq (m/- 1.0 t))
                      (m/+ (m/* h t (Array/aget ms i))
                           (m/* (m/inc (m/* 2.0 t)) (Array/aget ys i))))
                 (m/* t t (m/+ (m/* h (m/dec t) (Array/aget ms (m/inc i)))
                               (m/* (m/- 3.0 t t) (Array/aget ys (m/inc i))))))))))))

(defn monotone
  "Creates a 1D monotone piecewise cubic Hermite interpolator (Fritsch-Carlson).

  The interpolant passes through all data points and preserves the monotonicity of the data: on every interval it stays between the two neighbouring values, so it never overshoots and has a local extremum only at a data point. This avoids the oscillations of a natural cubic spline.

  The tangent at an interior node is the arithmetic mean of the two neighbouring secants, or zero when the secants have different signs or one of them is zero. The end tangents are the secants of the first and the last interval. Then every interval (including the last one) is limited: the ratios of the two tangents to the secant of the interval are scaled back into the circle of radius 3, which is sufficient for monotonicity (Fritsch and Carlson, 1980). The interpolant is C1.

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and strictly increasing; at least 2 points.
  - `ys` (sequence of numbers): values, the same length as `xs`.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`. `f` returns the data value exactly at a node. Outside the domain the cubic of the first or the last interval is continued (extrapolation, not limited). A `NaN` in `ys` produces `NaN` around it.

  Throws `ex-info` (with `:type :invalid-input`, `:method :monotone` and a `:property`) when the lengths differ, there are fewer than 2 points, `xs` contains `NaN` or an infinity, or `xs` is not strictly increasing (unsorted or duplicated).

  See also [[pchip]] (harmonic mean tangents), [[fastmath.interpolation.cubic/cubic]], [[fastmath.interpolation.linear/linear]]."
  [xs ys]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)]
    (ic/validate-xy! :monotone xs ys {:min-points 2})
    (hermite-interpolator xs ys (fc-tangents (secants xs ys)))))

(defn pchip
  "Creates a 1D shape-preserving piecewise cubic Hermite interpolator (PCHIP, Fritsch-Butland).

  The interpolant passes through all data points, is C1 and does not overshoot: on every interval it stays between the two neighbouring values. It equals `scipy.interpolate.PchipInterpolator`.

  The tangent at an interior node is the weighted harmonic mean of the two neighbouring secants (weights depend on the widths of the intervals, so a nonuniform grid is handled), or zero when the secants have different signs or one of them is zero. The end tangents use the non-centered three point formula, set to zero when its sign differs from the end secant and limited to three times the end secant when the secants of the end intervals have different signs. With two points the interpolant is the straight line.

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and strictly increasing; at least 2 points.
  - `ys` (sequence of numbers): values, the same length as `xs`.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`. `f` returns the data value exactly at a node. Outside the domain the cubic of the first or the last interval is continued (extrapolation, not limited).

  Throws `ex-info` (with `:type :invalid-input`, `:method :pchip` and a `:property`) when the lengths differ, there are fewer than 2 points, `xs` contains `NaN` or an infinity, or `xs` is not strictly increasing.

  See also [[monotone]] (arithmetic mean tangents, Fritsch-Carlson), [[fastmath.interpolation.cubic/cubic]]."
  [xs ys]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)]
    (ic/validate-xy! :pchip xs ys {:min-points 2})
    (hermite-interpolator xs ys (pchip-tangents xs (secants xs ys)))))
