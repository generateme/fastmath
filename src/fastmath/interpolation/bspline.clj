(ns fastmath.interpolation.bspline
  "Interpolating B-spline of arbitrary degree, parametrised in x.

  The spline passes through all data points. Default knots are the ones of `scipy.interpolate.make_interp_spline`: not-a-knot for odd degree (de Boor) and midpoints between data points for even degree. The spline is evaluated with de Boor's algorithm, outside the data range the end polynomial is continued."
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^:private ^:const default-degree 3)

(defn- default-knots
  "Knots of the interpolating B-spline for strictly increasing `xs` and degree `k`: odd degree not-a-knot, even degree midpoints (as scipy), with `k+1` copies of both ends."
  ^doubles [^doubles xs ^long k]
  (let [n (alength xs)
        odd? (m/odd? k)
        k2 (if odd? (m/long-quot (m/long-inc k) 2) (m/long-quot k 2))
        ;; interior candidates: data points (odd) or midpoints (even)
        cnt (long (if odd? n (m/dec n)))
        ni (m/- cnt (m/* 2 k2))
        t (double-array (m/+ n k 1))
        xe0 (Array/aget xs 0)
        xe1 (Array/aget xs (m/dec n))]
    (dotimes [i (m/inc k)]
      (Array/aset t i xe0)
      (Array/aset t (m/- (m/+ n k) i) xe1))
    (dotimes [i ni]
      (let [j (m/+ i k2)]
        (Array/aset t (m/+ i k 1)
                    (if odd?
                      (Array/aget xs j)
                      (m// (m/+ (Array/aget xs j) (Array/aget xs (m/inc j))) 2.0)))))
    t))

(defn- check-knots!
  "Validates user supplied knots `t` for nodes `xs` and degree `k`."
  [^doubles xs ^doubles t ^long k]
  (let [n (alength xs)
        nt (alength t)]
    (when-not (m/== nt (m/+ n k 1))
      (ic/invalid-input! :b-spline :knots-count
                         (str "knots need " (m/+ n k 1) " values (points + degree + 1), got " nt)
                         {:expected (m/+ n k 1) :count nt}))
    (dotimes [i nt]
      (let [a (Array/aget t i)]
        (when-not (m/valid-double? a)
          (ic/invalid-input! :b-spline :not-finite
                             (str "knots have to be finite, got " a " at index " i)
                             {:axis :knots :index i :value a}))
        (when (and (m/pos? i) (m/< a (Array/aget t (m/dec i))))
          (ic/invalid-input! :b-spline :knots-not-sorted
                             (str "knots have to be non-decreasing, got " a " after "
                                  (Array/aget t (m/dec i)) " at index " i)
                             {:axis :knots :index i :value a}))))
    (when-not (m/< (Array/aget t k) (Array/aget t n))
      (ic/invalid-input! :b-spline :knots-domain
                         "the knot interval between the indices degree and points is empty"
                         {:axis :knots}))
    (when (or (m/< (Array/aget xs 0) (Array/aget t k))
              (m/> (Array/aget xs (m/dec n)) (Array/aget t n)))
      (ic/invalid-input! :b-spline :knots-domain
                         "the nodes have to lie in the knot interval between the indices degree and points"
                         {:axis :knots}))))

(defn- find-interval
  "Index `ell` in `[k, n-1]` of the non-empty knot interval used to evaluate at `x`: the last one with `t[ell] <= x`, the first or the last one outside of the domain."
  ^long [^doubles t ^long k ^double x]
  (let [n (m/- (alength t) k 1)
        ell (long (loop [lo k
                         hi (m/dec n)]
                    (if (m/< lo hi)
                      (let [mid (m/long-quot (m/long-inc (m/long-add lo hi)) 2)]
                        (if (m/<= (Array/aget t mid) x)
                          (recur mid hi)
                          (recur lo (m/dec mid))))
                      lo)))]
    (if (m/== (Array/aget t ell) (Array/aget t (m/inc ell)))
      (if (m/< x (Array/aget t ell))
        (loop [e ell]
          (if (m/== (Array/aget t e) (Array/aget t (m/inc e)))
            (recur (m/inc e))
            e))
        (loop [e ell]
          (if (m/== (Array/aget t e) (Array/aget t (m/inc e)))
            (recur (m/dec e))
            e)))
      ell)))

(defn- basis-values!
  "Fills `nb` with the `k+1` non-zero B-spline values at `x` in the interval `ell` (Cox - de Boor)."
  [^doubles nb ^doubles t ^long ell ^double x]
  (let [k (m/dec (alength nb))
        left (double-array (m/inc k))
        right (double-array (m/inc k))]
    (Array/aset nb 0 1.0)
    (loop [j 1]
      (when (m/<= j k)
        (Array/aset left j (m/- x (Array/aget t (m/- (m/inc ell) j))))
        (Array/aset right j (m/- (Array/aget t (m/+ ell j)) x))
        (let [saved (double (loop [r 0
                                   saved 0.0]
                              (if (m/< r j)
                                (let [temp (m// (Array/aget nb r)
                                                (m/+ (Array/aget right (m/inc r))
                                                     (Array/aget left (m/- j r))))]
                                  (Array/aset nb r (m/+ saved (m/* (Array/aget right (m/inc r)) temp)))
                                  (recur (m/inc r) (m/* (Array/aget left (m/- j r)) temp)))
                                saved)))]
          (Array/aset nb j saved))
        (recur (m/inc j))))
    nb))

(defn- collocation
  "Banded collocation matrix of the nodes `xs`: row `i` in a flat array with `2k+1` entries per row, entry `(i, j)` at `i*(2k+1) + j - i + k`. Throws when the Schoenberg-Whitney conditions are not satisfied."
  ^doubles [^doubles xs ^doubles t ^long k]
  (let [n (alength xs)
        w (m/inc (m/* 2 k))
        ab (double-array (m/* n w))
        nb (double-array (m/inc k))]
    (dotimes [i n]
      (let [x (Array/aget xs i)
            ell (find-interval t k x)
            first-col (m/- ell k)]
        (when (or (m/< i first-col) (m/> i ell))
          (ic/invalid-input! :b-spline :schoenberg-whitney
                             (str "knots and nodes do not satisfy the Schoenberg-Whitney conditions at index " i)
                             {:axis :knots :index i}))
        (basis-values! nb t ell x)
        (when-not (m/pos? (Array/aget nb (m/- i first-col)))
          (ic/invalid-input! :b-spline :schoenberg-whitney
                             (str "knots and nodes do not satisfy the Schoenberg-Whitney conditions at index " i)
                             {:axis :knots :index i}))
        (dotimes [j (m/inc k)]
          (Array/aset ab (m/+ (m/* i w) (m/- (m/+ first-col j) i) k) (Array/aget nb j)))))
    ab))

(defn- band-solve!
  "Solves the banded system (`k` sub and super diagonals) in place without pivoting; the collocation matrix is totally positive so the elimination is stable. Returns the solution (the modified `rhs`)."
  ^doubles [^doubles ab ^doubles rhs ^long k]
  (let [n (alength rhs)
        w (m/inc (m/* 2 k))]
    ;; elimination
    (dotimes [c n]
      (let [pivot (Array/aget ab (m/+ (m/* c w) k))]
        (when-not (m/pos? pivot)
          (ic/invalid-input! :b-spline :singular
                             "the collocation matrix is singular"
                             {:index c}))
        (let [rmax (m/long-min (m/long-dec n) (m/long-add c k))]
          (loop [r (m/inc c)]
            (when (m/<= r rmax)
              (let [f (m// (Array/aget ab (m/+ (m/* r w) (m/- c r) k)) pivot)]
                (when-not (m/zero? f)
                  (loop [j c]
                    (when (m/<= j rmax)
                      (let [ir (m/+ (m/* r w) (m/- j r) k)]
                        (Array/aset ab ir (m/- (Array/aget ab ir)
                                               (m/* f (Array/aget ab (m/+ (m/* c w) (m/- j c) k))))))
                      (recur (m/inc j))))
                  (Array/aset rhs r (m/- (Array/aget rhs r) (m/* f (Array/aget rhs c))))))
              (recur (m/inc r)))))))
    ;; back substitution
    (loop [i (m/dec n)]
      (when (m/>= i 0)
        (let [jmax (m/long-min (m/long-dec n) (m/long-add i k))
              s (double (loop [j (m/inc i)
                               s (Array/aget rhs i)]
                          (if (m/<= j jmax)
                            (recur (m/inc j)
                                   (m/- s (m/* (Array/aget ab (m/+ (m/* i w) (m/- j i) k))
                                               (Array/aget rhs j))))
                            s)))]
          (Array/aset rhs i (m// s (Array/aget ab (m/+ (m/* i w) k)))))
        (recur (m/dec i))))
    rhs))

(defn- de-boor
  "Value of the spline with coefficients `c`, knots `t` and degree `k` at `x`."
  ^double [^doubles c ^doubles t ^long k ^double x]
  (let [ell (find-interval t k x)
        base (m/- ell k)
        d (double-array (m/inc k))]
    (dotimes [j (m/inc k)]
      (Array/aset d j (Array/aget c (m/+ base j))))
    (loop [r 1]
      (when (m/<= r k)
        (loop [j k]
          (when (m/>= j r)
            (let [i (m/+ base j)
                  ti (Array/aget t i)
                  alpha (m// (m/- x ti) (m/- (Array/aget t (m/- (m/+ i k 1) r)) ti))]
              (Array/aset d j (m/+ (m/* (m/- 1.0 alpha) (Array/aget d (m/dec j)))
                                   (m/* alpha (Array/aget d j)))))
            (recur (m/dec j))))
        (recur (m/inc r))))
    (Array/aget d k)))

(defn b-spline
  "Creates an interpolating B-spline of a given degree, parametrised in x.

  The spline passes through all data points and equals `scipy.interpolate.make_interp_spline(x, y, k=degree)` (default knots) up to rounding errors. It replaces the former control polygon form: `xs` are the nodes of the interpolation and not the control points.

  Parameters:

  - `xs` (sequence of numbers): nodes, strictly increasing and finite
  - `ys` (sequence of numbers): values, same length as `xs`; `NaN` or infinite values give `NaN` results
  - `opts` (map, optional):
      - `:degree` (long, default: `3`): degree of the spline, an integer `>= 1`; `1` is the linear interpolation, `3` the cubic spline
      - `:knots` (sequence of numbers, default: `nil`): full knot vector of `count xs + degree + 1` non-decreasing values; `nil` selects the default knots

  Default knots: for an odd degree the not-a-knot condition (de Boor): the data points without the `(degree + 1) / 2` nearest to each end; for an even degree the midpoints between data points without the `degree / 2` nearest to each end (as scipy). The end points are repeated `degree + 1` times. With these knots the cubic spline is the not-a-knot spline, so for `count xs = degree + 1` it is the interpolating polynomial.

  Custom `:knots` have to be non-decreasing and finite, the interval between the knots with indices `degree` and `count xs` has to contain all nodes and the Schoenberg-Whitney conditions have to be satisfied, otherwise the spline would not exist. The spline is defined by the knots on that interval and it is continued by the end polynomials outside of it.

  Returns a function of one argument `x` which returns the value of the spline. Outside of the data range the end polynomial is continued (no clamping), `NaN` gives `NaN`. The result is a polynomial of degree `degree` on every knot interval and it is continuously differentiable `degree - 1` times where knots are simple.

  Throws `ex-info` (`:type :invalid-input`) when the number of points is smaller than `degree + 1`, when `xs` and `ys` have different lengths, when `xs` is not finite or not strictly increasing, when `:degree` is not an integer `>= 1` or when the knots are invalid.

  See also [[fastmath.interpolation/interpolation]], [[fastmath.interpolation.cubic/cubic]]."
  ([xs ys] (b-spline xs ys nil))
  ([xs ys {:keys [degree knots]}]
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         degree (or degree default-degree)]
     (when-not (and (number? degree)
                    (m/valid-double? (double degree))
                    (m/== (double degree) (Math/rint (double degree)))
                    (m/>= (double degree) 1.0))
       (ic/invalid-input! :b-spline :degree
                          (str ":degree has to be an integer >= 1, got " degree)
                          {:option :degree :value degree :expected "integer >= 1"}))
     (let [k (long (m/min (double degree) 1.0e9))]
       (ic/validate-xy! :b-spline xs ys {:min-points (m/inc k) :sorted? true})
       (let [t (if knots
                 (let [t (m/seq->double-array knots)]
                   (check-knots! xs t k)
                   t)
                 (default-knots xs k))
             c (band-solve! (collocation xs t k) (aclone ys) k)]
         (fn ^double [^double x] (de-boor c t k x)))))))
