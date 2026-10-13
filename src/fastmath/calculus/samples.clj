(ns fastmath.calculus.samples
  "Numerical integration of sampled (tabulated) data.

  The integrand is not a function but a sequence of samples `ys`, taken either on a uniform grid with a known spacing or at arbitrary abscissae. The integral over the whole sampled range (or its running value at every sample) is computed with the composite trapezoid rule or the composite Simpson rule. The behaviour follows `scipy.integrate.trapezoid`, `simpson`, `cumulative_trapezoid` and `cumulative_simpson`.

  The namespace contains a single function, [[integrate-samples]]."
  (:require [fastmath.core :as m])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^:private methods-set #{:trapezoid :simpson})

(defn- samples->array
  "Converts a sequence or an array of numbers to `double[]`, throws for anything else."
  ^doubles [v option]
  (cond
    (nil? v) (double-array 0)
    (or (sequential? v) (.isArray (class v))) (m/seq->double-array v)
    :else (throw (ex-info "Samples must be a sequence or an array of numbers."
                          {:option option :value v :expected "sequence of numbers"}))))

(defn- spacings
  "Returns the array of the `n-1` spacings between consecutive samples."
  ^doubles [^doubles xs ^double dx ^long n]
  (let [nh (if (pos? n) (m/long-dec n) 0)
        h (double-array nh)]
    (dotimes [i nh]
      (Array/aset h i (if xs
                        (m/- (Array/aget xs (m/long-inc i)) (Array/aget xs i))
                        dx)))
    h))

;; trapezoid

(defn- trapezoid-total
  ^double [^doubles ys ^doubles h]
  (loop [i 0
         s 0.0]
    (if (< i (alength h))
      (recur (m/long-inc i)
             (m/+ s (m// (m/* (Array/aget h i)
                              (m/+ (Array/aget ys (m/long-inc i)) (Array/aget ys i)))
                         2.0)))
      s)))

(defn- trapezoid-cumulative
  ^doubles [^doubles ys ^doubles h]
  (let [n (alength ys)
        res (double-array n)]
    (loop [i 0
           s 0.0]
      (when (< i (alength h))
        (let [s (m/+ s (m// (m/* (Array/aget h i)
                                 (m/+ (Array/aget ys (m/long-inc i)) (Array/aget ys i)))
                            2.0))]
          (Array/aset res (m/long-inc i) s)
          (recur (m/long-inc i) s))))
    res))

;; simpson

(defn- simpson-pair
  "Simpson integral over two consecutive intervals starting at `i`, uniform spacing."
  ^double [^doubles ys ^long i]
  (m/+ (Array/aget ys i)
       (m/* 4.0 (Array/aget ys (m/long-inc i)))
       (Array/aget ys (m/long-add i 2))))

(defn- simpson-pair-nonuniform
  "Simpson integral over two consecutive intervals starting at `i`, arbitrary spacings.
  Division by a zero spacing contributes zero."
  ^double [^doubles ys ^doubles h ^long i]
  (let [h0 (Array/aget h i)
        h1 (Array/aget h (m/long-inc i))
        hsum (m/+ h0 h1)
        hprod (m/* h0 h1)
        h0divh1 (if (m/zero? h1) 0.0 (m// h0 h1))
        c0 (m/- 2.0 (if (m/zero? h0divh1) 0.0 (m// 1.0 h0divh1)))
        c1 (m/* hsum (if (m/zero? hprod) 0.0 (m// hsum hprod)))
        c2 (m/- 2.0 h0divh1)]
    (m/* (m// hsum 6.0)
         (m/+ (m/* (Array/aget ys i) c0)
              (m/* (Array/aget ys (m/long-inc i)) c1)
              (m/* (Array/aget ys (m/long-add i 2)) c2)))))

(defn- simpson-pairs
  "Sum of the Simpson pairs, uniform spacing (to be multiplied by dx/3)."
  ^double [^doubles ys ^long npairs]
  (loop [j 0
         s 0.0]
    (if (< j npairs)
      (recur (m/long-inc j) (m/+ s (simpson-pair ys (m/long-mult 2 j))))
      s)))

(defn- simpson-pairs-nonuniform
  "Sum of the Simpson pairs, arbitrary spacings."
  ^double [^doubles ys ^doubles h ^long npairs]
  (loop [j 0
         s 0.0]
    (if (< j npairs)
      (recur (m/long-inc j) (m/+ s (simpson-pair-nonuniform ys h (m/long-mult 2 j))))
      s)))

(defn- safe-div
  ^double [^double num ^double den]
  (if (m/zero? den) 0.0 (m// num den)))

(defn- simpson-total
  ^double [^doubles ys ^doubles h uniform?]
  (let [n (alength ys)]
    (if (< n 3)
      (if (== n 2) (trapezoid-total ys h) 0.0)
      (let [npairs (m/long-quot (m/long-dec n) 2)
            main (if uniform?
                   (m/* (m// (Array/aget h 0) 3.0) (simpson-pairs ys npairs))
                   (simpson-pairs-nonuniform ys h npairs))]
        (if (m/odd? n)
          main
          ;; even number of samples: correction for the last interval (Cartwright)
          (let [h0 (Array/aget h (m/long-sub n 3))
                h1 (Array/aget h (m/long-sub n 2))
                alpha (safe-div (m/+ (m/* 2.0 h1 h1) (m/* 3.0 h0 h1))
                                (m/* 6.0 (m/+ h1 h0)))
                beta (safe-div (m/+ (m/* h1 h1) (m/* 3.0 h0 h1))
                               (m/* 6.0 h0))
                eta (safe-div (m/* h1 h1 h1)
                              (m/* 6.0 h0 (m/+ h0 h1)))]
            (m/- (m/+ main
                      (m/* alpha (Array/aget ys (m/long-sub n 1)))
                      (m/* beta (Array/aget ys (m/long-sub n 2))))
                 (m/* eta (Array/aget ys (m/long-sub n 3))))))))))

(defn- simpson-sub
  "Integral over the first interval (of width `x21`) of the parabola through three points
  with values `fs` = [f1 f2 f3]. `x32` is the width of the second interval."
  ^double [^double x21 ^double x32 ^doubles fs uniform?]
  (let [f1 (Array/aget fs 0)
        f2 (Array/aget fs 1)
        f3 (Array/aget fs 2)]
    (if uniform?
      (m/* (m// x21 3.0)
           (m/- (m/+ (m// (m/* 5.0 f1) 4.0) (m/* 2.0 f2)) (m// f3 4.0)))
      (let [x31 (m/+ x21 x32)
            x21-x31 (m// x21 x31)
            x21-x32 (m// x21 x32)
            x21x21-x31x32 (m/* x21-x31 x21-x32)
            coeff1 (m/- 3.0 x21-x31)
            coeff2 (m/+ 3.0 x21x21-x31x32 x21-x31)
            coeff3 (m/- x21x21-x31x32)]
        (m/* (m// x21 6.0)
             (m/+ (m/* coeff1 f1) (m/* coeff2 f2) (m/* coeff3 f3)))))))

(defn- simpson-cumulative
  ^doubles [^doubles ys ^doubles h uniform?]
  (let [n (alength ys)]
    (if (< n 3)
      (trapezoid-cumulative ys h)
      (let [res (double-array n)
            fs (double-array 3)
            last-interval (m/long-sub n 2)]
        (loop [k 0
               s 0.0]
          (when (<= k last-interval)
            (let [;; odd intervals and the last one are integrated over the second interval of the
                  ;; triple (k-1, k, k+1), the other ones over the first interval of (k, k+1, k+2);
                  ;; the former is the latter for the reversed triple
                  second? (or (m/odd? k) (== k last-interval))
                  i (if second? (m/long-dec k) k)
                  ya (Array/aget ys i)
                  yb (Array/aget ys (m/long-inc i))
                  yc (Array/aget ys (m/long-add i 2))
                  ha (Array/aget h i)
                  hb (Array/aget h (m/long-inc i))
                  part (if second?
                         (do (Array/aset fs 0 yc)
                             (Array/aset fs 1 yb)
                             (Array/aset fs 2 ya)
                             (simpson-sub hb ha fs uniform?))
                         (do (Array/aset fs 0 ya)
                             (Array/aset fs 1 yb)
                             (Array/aset fs 2 yc)
                             (simpson-sub ha hb fs uniform?)))
                  s (m/+ s part)]
              (Array/aset res (m/long-inc k) s)
              (recur (m/long-inc k) s))))
        res))))

(defn integrate-samples
  "Integrates sampled data with the composite trapezoid or Simpson rule.

  The samples `ys` are taken either at equally spaced points with a known spacing, or at arbitrary (not necessarily uniform) abscissae. The result is the integral over the whole sampled range, or, in the cumulative mode, the running integral at every sample. The results match `scipy.integrate.trapezoid`, `simpson`, `cumulative_trapezoid` and `cumulative_simpson` (the last two with `initial=0`).

  Parameters:

  - `ys` (sequence or `double[]` of numbers): the samples.
  - `x-or-dx` (optional): either a sequence or `double[]` of abscissae, of the same length as `ys`, or a number `dx` being the spacing between consecutive samples. Defaults to the spacing `1.0` (also when `nil`). A negative or zero spacing is allowed and is used as is, as a decreasing `x`.
  - `options` (optional map):
    - `:method` (`:trapezoid` or `:simpson`, default `:trapezoid`): the quadrature rule.
    - `:cumulative?` (boolean, default `false`): when true the running integral is returned instead of the total.

  Returns a `double` when `:cumulative?` is false. Otherwise returns a vector of doubles of the same length as `ys`, where the first element is `0.0` and the i-th element is the integral from the first to the i-th sample.

  The trapezoid rule is exact for linear data. The Simpson rule fits a parabola through each three consecutive samples; the total is exact for polynomials up to degree 3 on a uniform grid with an odd number of samples, and up to degree 2 otherwise. For an even number of samples the last interval is treated with the three-point correction of Cartwright, so any number of samples is handled. The cumulative Simpson integral is computed over every interval from the parabola through the nearest three samples. Simpson with fewer than 3 samples is the trapezoid rule.

  Corner cases:

  - no samples: returns `0.0`, or `[]` in the cumulative mode.
  - one sample: returns `0.0`, or `[0.0]` in the cumulative mode.
  - zero spacings in `x` are not rejected: the trapezoid rule gives zero for such an interval, the total Simpson rule treats the undefined terms as zero, the cumulative Simpson rule gives NaN or infinite values.
  - throws `ex-info` with `:ys-count` and `:x-count` in the ex-data when `x-or-dx` is a sequence of a different length than `ys`.
  - throws `ex-info` with `:option`, `:value` and `:expected` in the ex-data for an unknown `:method`, and for `ys` or `x-or-dx` which are neither numbers nor sequences as required.

  See also [[fastmath.calculus/integrate]] (integration of functions)."
  ([ys] (integrate-samples ys nil nil))
  ([ys x-or-dx] (integrate-samples ys x-or-dx nil))
  ([ys x-or-dx options]
   (let [{:keys [method cumulative?] :or {method :trapezoid}} options
         _ (when-not (methods-set method)
             (throw (ex-info "Unknown integration method."
                             {:option :method :value method :expected methods-set})))
         ys (samples->array ys :ys)
         n (alength ys)
         [xs dx] (cond
                   (nil? x-or-dx) [nil 1.0]
                   (number? x-or-dx) [nil (double x-or-dx)]
                   :else [(samples->array x-or-dx :x-or-dx) 1.0])
         xs (when xs
              (when-not (== (alength ^doubles xs) n)
                (throw (ex-info "Lengths of ys and x differ."
                                {:ys-count n :x-count (alength ^doubles xs)})))
              xs)
         uniform? (nil? xs)
         h (spacings xs dx n)]
     (if cumulative?
       (if (zero? n)
         []
         (vec (if (= method :simpson)
                (simpson-cumulative ys h uniform?)
                (trapezoid-cumulative ys h))))
       (if (= method :simpson)
         (simpson-total ys h uniform?)
         (trapezoid-total ys h))))))
