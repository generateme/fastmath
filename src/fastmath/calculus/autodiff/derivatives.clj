(ns ^:no-doc fastmath.calculus.autodiff.derivatives
  "Internal: seeding of the inputs, extraction of the partials and the passes of [[fastmath.calculus.autodiff/derivative]], [[fastmath.calculus.autodiff/gradient]], [[fastmath.calculus.autodiff/jacobian]], [[fastmath.calculus.autodiff/hessian]], [[fastmath.calculus.autodiff/jvp]] and [[fastmath.calculus.autodiff/hvp]]."
  (:require [fastmath.core :as m]
            [fastmath.calculus.common :as cc]
            [fastmath.calculus.autodiff.dual :as d])
  (:import [fastmath.calculus.autodiff.dual Dual]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; ## Extraction

(defn- result-partials
  "The `k` partials (an array) of the layer `t` of one element of the result of a differentiated function.

  A number or a dual with an older tag does not depend on the inputs of the layer: `k` zeros. A dual with a newer tag escaped from a finished derivative call."
  ^objects [^long t ^long k result]
  (cond
    (instance? Dual result) (let [rt (.tag ^Dual result)]
                              (cond
                                (m/== rt t) (.partials ^Dual result)
                                (m/< rt t) (d/zero-partials k)
                                :else (throw (ex-info "A dual escaped from its derivative call into the result."
                                                      {:tag t :result-tag rt}))))
    (number? result) (d/zero-partials k)
    (fn? result) (throw (ex-info "Derivative of a function returning a function is not supported."
                                 {:tag t :result-type (type result)}))
    :else (throw (ex-info "Derivative of a function returning neither a number nor a dual is not supported."
                          {:tag t :result-type (type result)}))))

;; ## Derivative

(defn- first-derivative
  [f]
  (fn [x]
    (let [t (d/fresh-tag)]
      (aget ^objects (result-partials t 1 (f (d/new-dual t (d/leaf x) (d/basis-partials 1 0)))) 0))))

(defn derivative
  "The `n`-th derivative of `f` as `n` nested first derivatives, see [[fastmath.calculus.autodiff/derivative]]."
  ([f] (derivative f 1))
  ([f ^long n]
   (when (m/neg? n) (cc/invalid-option! :n n "an integer >= 0"))
   (loop [g f
          k n]
     (if (m/zero? k)
       g
       (recur (first-derivative g) (m/long-dec k))))))

;; ## Gradient, Jacobian, Hessian

(def ^:private ^:const ^long chunk-threshold 8)

(defn- ceil-quot ^long [^long a ^long b] (m/long-quot (m/long-add a (m/long-dec b)) b))

(defn- default-chunk-size
  "Balanced chunks of at most `chunk-threshold` inputs: `n` for small `n`, otherwise `n` split evenly into the least number of chunks."
  ^long [^long n]
  (if (m/<= n chunk-threshold)
    n
    (ceil-quot n (ceil-quot n chunk-threshold))))

(defn- validated-chunk-size
  "`:chunk-size` of the options as a long (truncated), `nil` when missing. A value which is not a number throws
  `ClassCastException`; a value below 1, `NaN` or an infinity throws `ex-info`; a finite value above the `long` range
  is `Long/MAX_VALUE` (it is clamped to the number of inputs anyway)."
  [options]
  (let [k (:chunk-size options)]
    (when-not (nil? k)
      (let [x (double k)]
        (when-not (and (m/>= x 1.0) (not (Double/isInfinite x)))
          (cc/invalid-option! :chunk-size k "an integer >= 1"))
        (if (m/>= x 9.223372036854775807E18) Long/MAX_VALUE (long x))))))

(defn- chunk-length
  "Inputs per pass for `n` inputs: the given chunk size clamped to `n`, or the default."
  ^long [chunk-size ^long n]
  (if chunk-size (m/long-min (long chunk-size) n) (default-chunk-size n)))

(defn- point
  "A point (any sequence or a `double[]`) as a vector, numbers coerced to doubles, duals kept."
  [v]
  (mapv d/leaf v))

(defn- seed-chunk
  "The point `x` with the `k` inputs from `start` replaced by duals of tag `t` with basis partials; the other inputs stay constants."
  [x ^long t ^long start ^long k]
  (loop [i 0
         xs (transient x)]
    (if (m/< i k)
      (let [j (m/long-add start i)]
        (recur (m/long-inc i) (assoc! xs j (d/new-dual t (x j) (d/basis-partials k i)))))
      (persistent! xs))))

(defn- gradient-at
  [f x ^long chunk]
  (let [n (count x)]
    (loop [start 0
           acc (transient [])]
      (if (m/< start n)
        (let [k (m/long-min chunk (m/long-sub n start))
              t (d/fresh-tag)
              ps (result-partials t k (f (seed-chunk x t start k)))]
          (recur (m/long-add start k) (reduce conj! acc ps)))
        (persistent! acc)))))

(defn- result-values
  "The values of a vector valued function as a vector (a lazy result is realized here, inside the pass); throws when the result is not sequential or an element is neither a number nor a dual."
  [r]
  (when-not (sequential? r)
    (throw (ex-info "Jacobian of a function returning something else than a sequential of numbers or duals is not supported."
                    {:result-type (type r)})))
  (let [ys (vec r)]
    (doseq [y ys]
      (when-not (or (number? y) (instance? Dual y))
        (throw (ex-info "Jacobian of a function returning a sequential with an element which is neither a number nor a dual is not supported."
                        {:element-type (type y)}))))
    ys))

(defn- jacobian-at
  [f x ^long chunk]
  (let [n (count x)]
    ;; `rows` is nil before the first pass: the first pass always runs, also for the empty point (n = 0),
    ;; where it is a pass over 0 inputs which calls `f` once to get the number of rows m
    (loop [start 0
           rows nil]
      (if (or (m/< start n) (nil? rows))
        (let [k (m/long-min chunk (m/long-sub n start))
              t (d/fresh-tag)
              ys (result-values (f (seed-chunk x t start k)))
              block (mapv (fn [y] (vec (result-partials t k y))) ys)]
          (when (and rows (not (m/== (count rows) (count block))))
            (throw (ex-info "The function returned a different number of values."
                            {:expected (count rows) :actual (count block)})))
          (recur (m/long-add start k) (if rows (mapv into rows block) block)))
        rows))))

(defn gradient
  "The gradient function of `f`, see [[fastmath.calculus.autodiff/gradient]]."
  ([f] (gradient f nil))
  ([f options]
   (let [chunk-size (validated-chunk-size options)]
     (fn local-gradient [v]
       (let [x (point v)]
         (gradient-at f x (chunk-length chunk-size (count x))))))))

(defn jacobian
  "The Jacobian function of `f`, see [[fastmath.calculus.autodiff/jacobian]]."
  ([f] (jacobian f nil))
  ([f options]
   (let [chunk-size (validated-chunk-size options)]
     (fn local-jacobian [v]
       (let [x (point v)]
         (jacobian-at f x (chunk-length chunk-size (count x))))))))

(defn hessian
  "The Hessian function of `f` as the Jacobian of the gradient, see [[fastmath.calculus.autodiff/hessian]]."
  ([f] (hessian f nil))
  ([f options]
   (jacobian (gradient f options) options)))

;; ## Derivative products

(defn- seed-direction
  "The point `x` with every input replaced by a dual of tag `t` with the single partial `v_i` (the direction); both
  vectors have the same length."
  [x v ^long t]
  (let [n (count x)]
    (loop [i 0
           xs (transient x)]
      (if (m/< i n)
        (let [p (object-array 1)]
          (aset p 0 (v i))
          (recur (m/long-inc i) (assoc! xs i (d/new-dual t (x i) p))))
        (persistent! xs)))))

(defn- directional-derivative
  "The partial of the layer `t` of the result `r` of one pass: a number for a number or a dual, a vector for a sequential."
  [^long t r]
  (if (sequential? r)
    (mapv (fn [y] (aget ^objects (result-partials t 1 y) 0)) (result-values r))
    (aget ^objects (result-partials t 1 r) 0)))

(defn jvp
  "The Jacobian-vector product function of `f`, one pass, see [[fastmath.calculus.autodiff/jvp]]."
  [f]
  (fn local-jvp [x v]
    (let [x (point x)
          v (point v)]
      (when-not (m/== (count x) (count v))
        (throw (ex-info "The point and the direction have different lengths."
                        {:count-x (count x) :count-v (count v)})))
      (let [t (d/fresh-tag)]
        (directional-derivative t (f (seed-direction x v t)))))))

(defn hvp
  "The Hessian-vector product function of `f` as the Jacobian-vector product of the gradient, see [[fastmath.calculus.autodiff/hvp]]."
  ([f] (hvp f nil))
  ([f options]
   (jvp (gradient f options))))
