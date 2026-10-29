;; Numerical recipes, p. 128
(ns fastmath.interpolation.barycentric
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- weights
  "Floater-Hormann barycentric weights of order `d` (all weights share the factor `scale^d`, which is
  irrelevant for the interpolant and keeps the products of node distances in the double range)."
  ^doubles [^doubles xs ^long d ^double scale]
  (let [n (alength xs)
        n-d (m/- n d)
        w (double-array n)]
    (dotimes [k n]
      (let [xk (Array/aget xs k)
            imin (m/max 0 (m/- k d))
            imax (m/min k (m/dec n-d))
            s (double (loop [i imin
                             s 0.0]
                        (if (m/> i imax)
                          s
                          (recur (m/inc i)
                                 (m/+ s (m// 1.0 (double (loop [j i
                                                                p 1.0]
                                                           (cond
                                                             (m/> j (m/+ i d)) p
                                                             (m/== j k) (recur (m/inc j) p)
                                                             :else (recur (m/inc j)
                                                                          (m/* p scale (m/abs (m/- xk (Array/aget xs j))))))))))))))]
        (Array/aset w k (if (m/even? (m/long-sub k d)) s (m/- s)))))
    w))

(defn- sort-nodes
  "Returns `[xs ys]` sorted by `xs` (the arrays are returned unchanged when `xs` is already sorted)."
  [^doubles xs ^doubles ys]
  (let [n (alength xs)]
    (if (loop [i 1]
          (cond
            (m/>= i n) true
            (m/< (Array/aget xs i) (Array/aget xs (m/dec i))) false
            :else (recur (m/inc i))))
      [xs ys]
      (let [order (vec (sort-by (fn [^long i] (Array/aget xs i)) (range n)))]
        [(double-array (map (fn [^long i] (Array/aget xs i)) order))
         (double-array (map (fn [^long i] (Array/aget ys i)) order))]))))

(defn- order-valid?
  [order ^long n]
  (and (number? order)
       (m/valid-double? (double order))
       (m/== (double order) (Math/floor (double order)))
       (m/>= (double order) 0.0)
       (m/<= (double order) (m/dec n))))

(defn- barycentric*
  "Barycentric rational interpolator over arrays; `order` is validated here."
  [method xs ys order]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)
        _ (ic/validate-xy! method xs ys {:min-points 1 :sorted? false})
        n (alength xs)]
    (when-not (order-valid? order n)
      (ic/invalid-input! method :invalid-order
                         (str "order has to be an integer between 0 and " (m/dec n) " (count of xs - 1), got " order)
                         {:order order :count n}))
    (let [d (long order)
          [^doubles xs ^doubles ys] (sort-nodes xs ys)]
      (loop [i (long 1)]
        (when (m/< i n)
          (when (m/== (Array/aget xs i) (Array/aget xs (m/dec i)))
            (ic/invalid-input! method :duplicate-nodes
                               (str "xs has to be distinct, got " (Array/aget xs i) " twice")
                               {:value (Array/aget xs i)}))
          (recur (m/inc i))))
      (let [span (m/- (Array/aget xs (m/dec n)) (Array/aget xs 0))
            scale (if (m/pos? span) (m// 4.0 span) 1.0)
            ^doubles w (weights xs d scale)]
        (fn ^double [^double x]
          (double
           (loop [i (long 0)
                  num 0.0
                  den 0.0]
             (if (m/== i n)
               (m// num den)
               (let [xi (Array/aget xs i)]
                 (if (m/== x xi)
                   (Array/aget ys i)
                   (let [t (m// (Array/aget w i) (m/- x xi))]
                     (recur (m/inc i)
                            (m/+ num (m/* t (Array/aget ys i)))
                            (m/+ den t)))))))))))))

(defn barycentric
  "Creates a 1D barycentric rational interpolator (Floater-Hormann).

  The interpolant is the barycentric rational function with the Floater-Hormann weights of the given order: a blend of local polynomial interpolants of degree `order` through `order + 1` neighbouring nodes. It has no poles on the real line for any distribution of nodes and passes through all data points. The approximation order is `order + 1` for smooth data, so a higher `order` gives a better fit for smooth data, up to the polynomial interpolant for `order` equal to `(count xs) - 1`, which is the least robust choice. Order 0 gives the Berrut rational interpolant with alternating weights. For order close to the number of points on equispaced nodes the usual instability of high degree polynomial interpolation returns.

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and distinct; they do not need to be sorted (the nodes are sorted internally, the weights are defined by the neighbourhood in the sorted order).
  - `ys` (sequence of numbers): values, the same length as `xs`.
  - `opts` (map, optional):
      - `:order` (integer, default: 1): blending degree, an integer-valued number with `0 <= order <= (count xs) - 1`. A `double` with an integer value is accepted.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`. A query exactly equal to a node returns the value of that node. Outside the domain the rational function is evaluated as is (extrapolation, accuracy decreases quickly).

  Throws `ex-info` (with `:type :invalid-input`, `:method :barycentric` and a `:property`) when the lengths differ, there are no points (`:too-few-points`), `xs` contains `NaN` or an infinity (`:not-finite`), `xs` has duplicates (`:duplicate-nodes`), or `:order` is not an integer-valued number in the range above (`:invalid-order`, that includes an order for which there are too few points).

  See also [[polynomial]], [[fastmath.interpolation.cubic/cubic]]."
  ([xs ys] (barycentric xs ys nil))
  ([xs ys {:keys [order]}]
   (barycentric* :barycentric xs ys (if (nil? order) 1 order))))

(defn polynomial
  "Creates a 1D polynomial interpolator in the barycentric form.

  The interpolant is the unique polynomial of degree `(count xs) - 1` through all data points. It is the [[barycentric]] interpolant with `:order` equal to `(count xs) - 1`, evaluated with the barycentric formula, which is numerically stable for well distributed nodes (for example Chebyshev nodes). On equispaced nodes with many points the polynomial interpolation itself is ill-conditioned (Runge's phenomenon), whatever the algorithm.

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and distinct; at least 2 points; they do not need to be sorted.
  - `ys` (sequence of numbers): values, the same length as `xs`.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`. A query exactly equal to a node returns the value of that node. Outside the domain the polynomial is evaluated as is (extrapolation).

  Throws `ex-info` (with `:type :invalid-input`, `:method :polynomial` and a `:property`) when the lengths differ, there are fewer than 2 points, `xs` contains `NaN` or an infinity, or `xs` has duplicates (`:duplicate-nodes`).

  See also [[barycentric]]."
  [xs ys]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)]
    (ic/validate-xy! :polynomial xs ys {:min-points 2 :sorted? false})
    (barycentric* :polynomial xs ys (m/dec (alength xs)))))
