(ns fastmath.interpolation.akima
  "Akima and modified Akima (makima) piecewise cubic interpolation of one dimensional data."
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^:private ^:const break-mult
  ;; the node slope is defined by the weights only when the sum of the weights is
  ;; larger than this fraction of the largest sum (the same cutoff as scipy)
  1.0e-9)

(defn- slope-weight
  "Weight of the extended secant slopes `sl[k]` and `sl[k+1]`: their absolute difference, plus for makima
  half of the absolute sum."
  ^double [^doubles sl makima? k]
  (let [k (long k)
        a (Array/aget sl k)
        b (Array/aget sl (m/inc k))
        d (m/abs (m/- b a))]
    (if makima?
      (m/+ d (m/* 0.5 (m/abs (m/+ b a))))
      d)))

(defn- max-weight-sum
  "Largest sum of the two weights of a node; NaN sums are skipped, so a NaN value does not affect the distant nodes."
  ^double [^doubles sl makima? ^long n]
  (loop [i (long 0)
         mx ##-Inf]
    (if (m/< i n)
      (let [f12 (m/+ (slope-weight sl makima? (m/+ i 2))
                     (slope-weight sl makima? i))]
        (recur (m/inc i) (if (m/nan? f12) mx (m/max mx f12))))
      mx)))

(defn- node-slopes
  "Slopes of the interpolant at the nodes.

  `sl[j + 2]` is the secant slope of interval `j`; two extended slopes are added at both ends by
  the linear extrapolation `sl[1] = 2 sl[2] - sl[3]`, `sl[0] = 2 sl[1] - sl[2]` (and symmetrically). The slope
  at node `i` is `sl[i+1] + f2 / (f1 + f2) (sl[i+2] - sl[i+1])` with weights `f1 = |sl[i+3] - sl[i+2]|`
  and `f2 = |sl[i+1] - sl[i]|` (makima adds half of the absolute sums). When `f1 + f2` is not larger than
  `break-mult` times the largest sum, the mean of `sl[i]` and `sl[i+3]` is used."
  ^doubles [^doubles xs ^doubles ys makima?]
  (let [n (alength xs)]
    (if (m/== n 2)
      (let [s (m// (m/- (Array/aget ys 1) (Array/aget ys 0))
                   (m/- (Array/aget xs 1) (Array/aget xs 0)))]
        (double-array [s s]))
      (let [sl (double-array (m/+ n 3))]
        (dotimes [j (m/dec n)]
          (Array/aset sl (m/+ j 2) (m// (m/- (Array/aget ys (m/inc j)) (Array/aget ys j))
                                        (m/- (Array/aget xs (m/inc j)) (Array/aget xs j)))))
        (Array/aset sl 1 (m/- (m/* 2.0 (Array/aget sl 2)) (Array/aget sl 3)))
        (Array/aset sl 0 (m/- (m/* 2.0 (Array/aget sl 1)) (Array/aget sl 2)))
        (Array/aset sl (m/inc n) (m/- (m/* 2.0 (Array/aget sl n)) (Array/aget sl (m/dec n))))
        (Array/aset sl (m/+ n 2) (m/- (m/* 2.0 (Array/aget sl (m/inc n))) (Array/aget sl n)))
        (let [cutoff (m/* ^double break-mult (max-weight-sum sl makima? n))
              t (double-array n)]
          (dotimes [i n]
            (let [f1 (slope-weight sl makima? (m/+ i 2))
                  f2 (slope-weight sl makima? i)
                  f12 (m/+ f1 f2)
                  m1 (Array/aget sl (m/inc i))]
              (Array/aset t i (if (m/> f12 cutoff)
                                (m/+ m1 (m/* (m// f2 f12) (m/- (Array/aget sl (m/+ i 2)) m1)))
                                (m/* 0.5 (m/+ (Array/aget sl (m/+ i 3)) (Array/aget sl i)))))))
          t)))))

(defn- hermite
  "Builds the piecewise cubic Hermite interpolant through the nodes with the slopes of [[node-slopes]].
  Outside of the domain the end cubics continue."
  [^doubles xs ^doubles ys makima?]
  (let [n (alength xs)
        ns (m/dec n)
        t (node-slopes xs ys makima?)
        c2 (double-array ns)
        c3 (double-array ns)]
    (dotimes [i ns]
      (let [h (m/- (Array/aget xs (m/inc i)) (Array/aget xs i))
            s (m// (m/- (Array/aget ys (m/inc i)) (Array/aget ys i)) h)
            t0 (Array/aget t i)
            t1 (Array/aget t (m/inc i))]
        (Array/aset c2 i (m// (m/- (m/- (m/* 3.0 s) (m/* 2.0 t0)) t1) h))
        (Array/aset c3 i (m// (m/- (m/+ t0 t1) (m/* 2.0 s)) (m/* h h)))))
    (fn ^double [^double x]
      (let [i (ic/binary-search xs x)
            d (m/- x (Array/aget xs i))]
        (m/+ (Array/aget ys i)
             (m/* d (m/+ (Array/aget t i)
                         (m/* d (m/+ (Array/aget c2 i)
                                     (m/* d (Array/aget c3 i)))))))))))

(def ^:private valid-methods #{:akima :makima})

(defn akima
  "Creates a one dimensional Akima (or modified Akima) piecewise cubic interpolator.

  The interpolant is a cubic Hermite spline. The slope at every node is a weighted average of the two neighbouring secant slopes, where the weights are the absolute differences of the more distant secant slopes. This makes the curve local and much less prone to overshoot and wiggles than the natural cubic spline, but it is only C1 continuous. The slopes at the ends are computed from two additional secant slopes obtained by linear extrapolation of the nearest ones, so the end intervals are treated the same way as the interior ones. The result equals `scipy.interpolate.Akima1DInterpolator`.

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and strictly increasing, at least 5 points
  - `ys` (sequence of numbers): values at the nodes, the same count as `xs`; `NaN` and infinite values are allowed and make the curve `NaN` only near them (up to three intervals on both sides), the rest of the curve is unaffected
  - `opts` (map, optional):
      - `:method` (keyword, default: `:akima`): `:akima` for the original weights or `:makima` for the modified Akima algorithm, in which each weight additionally contains half of the absolute sum of the two secant slopes. `:makima` has no jump of the slope when the weights nearly vanish and gives smaller overshoot near plateaus. It equals `Akima1DInterpolator(method='makima')`. `nil` means the default.

  Returns a function `f` of one argument `x` returning the interpolated `double` value. The curve passes through all nodes. Outside of the domain the first or the last cubic polynomial is continued (extrapolation like scipy), `NaN` input gives `NaN`. When the weights of a node vanish (for example four equal consecutive slopes) the slope is the mean of the two surrounding more distant secant slopes. Data lying exactly on a line is reproduced exactly.

  Throws `ex-info` (with `:type` `:invalid-input`, `:method` `:akima` and a `:property`) when the lengths differ (`:length-mismatch`), there are fewer than 5 points (`:too-few-points`), `xs` contains a non-finite value (`:not-finite`), `xs` is not strictly increasing (`:not-increasing`) or `:method` is unknown (`:unknown-method`).

  See also [[cubic]], [[monotone]], [[linear]], [[interpolation]]."
  ([xs ys] (akima xs ys nil))
  ([xs ys {:keys [method]}]
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         method (if (nil? method) :akima method)]
     (when-not (contains? valid-methods method)
       (ic/invalid-input! :akima :unknown-method
                          (str "unknown method " (pr-str method) ", expected :akima or :makima")
                          {:value method :expected valid-methods}))
     (ic/validate-xy! :akima xs ys {:min-points 5})
     (hermite xs ys (identical? method :makima)))))
