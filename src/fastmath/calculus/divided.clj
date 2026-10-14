(ns fastmath.calculus.divided
  "Divided differences of a function at given nodes, with repeated nodes (Hermite, confluent divided differences).

  The divided difference `f[x0, x1, ..., xn]` is the leading coefficient of the Newton interpolation polynomial through the nodes. When nodes are repeated, the polynomial interpolates also the derivatives and the divided difference of `k + 1` equal nodes is `f^(k)(x) / k!`. The derivatives are given by the user or approximated with finite differences.

  Functions: [[divided]] (the highest divided difference only) and [[divided-all]] (all levels of the table)."
  (:require [fastmath.core :as m]
            [fastmath.calculus.finite :as finite])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- invalid!
  [option value expected]
  (throw (ex-info (str "Invalid " option ": " (pr-str value) ", expected " expected)
                  {:option option :value value :expected expected})))

(defn- function-like?
  "A function (including vars, partials, etc.), but not a collection, keyword or symbol."
  [f]
  (and (ifn? f) (not (coll? f)) (not (keyword? f)) (not (symbol? f))))

(defn- run-lengths
  "Lengths of the runs of equal nodes, throws when a node value appears in two separate runs."
  [^doubles axs xs]
  (let [n (alength axs)]
    (loop [i 0
           seen #{}
           runs []]
      (if (m/< i n)
        (let [x (Array/aget axs i)
              end (long (loop [j (m/long-inc i)]
                          (if (and (m/< j n) (m/== (Array/aget axs j) x)) (recur (m/long-inc j)) j)))
              key (m/+ x 0.0)] ;; -0.0 is 0.0
          (when (contains? seen key)
            (invalid! :xs xs "equal nodes must be adjacent"))
          (recur end (conj seen key) (conj runs (m/long-sub end i))))
        runs))))

(defn- finite-derivative
  "Derivative of order `k` by finite differences with Richardson extrapolation. When the starting step 0.5 gives
  a non-finite value (e.g. `f` is not defined so far from the node), the starting step is divided by 10, up to 8 times."
  [f ^long k]
  (fn ^double [^double x]
    (loop [init-h 0.5
           attempt (long 1)]
      (let [v ((finite/derivative f k {:extrapolate? {:init-h init-h}}) x)]
        (if (or (m/valid-double? v) (m/>= attempt 8))
          v
          (recur (m// init-h 10.0) (m/long-inc attempt)))))))

(defn- node-values
  "For every node a vector of `f^(k)(x) / k!` for k = 0 .. (length of the run) - 1, shared by the nodes of a run."
  [f ^doubles axs xs derivatives runs]
  (let [n (alength axs)
        derivs (vec derivatives)
        needed (m/long-dec (long (reduce m/long-max 1 runs)))]
    (when (and (m/pos? needed) (nil? derivatives) (not (function-like? f)))
      (invalid! :f f "a function when nodes are repeated and :derivatives is not given"))
    (when (and (m/pos? needed) derivatives (m/< (count derivs) needed))
      (invalid! :derivatives derivatives (str "a sequence of at least " needed " functions (f' f'' ...)")))
    (let [dfns (when (m/pos? needed)
                 (if derivatives
                   derivs
                   (mapv (fn [^long k] (finite-derivative f (m/long-inc k)))
                         (range needed))))]
      (loop [i 0
             runs runs
             res (transient [])]
        (if (m/< i n)
          (let [len (long (first runs))
                x (nth xs i)
                dx (Array/aget axs i)
                vals (reduce (fn [acc ^long k]
                               (conj acc (m// (double ((dfns (m/long-dec k)) dx)) (m/factorial k))))
                             [(double (f x))]
                             (range 1 len))]
            (recur (m/long-add i len)
                   (rest runs)
                   (reduce (fn [r _] (conj! r vals)) res (range len))))
          (persistent! res))))))

(defn- prepare
  "Validates the input, returns `[ys axs values]`: function values, nodes and per-node values."
  [f xs derivatives]
  (when-not (ifn? f)
    (invalid! :f f "a function (or a vector or a map) of the node"))
  (let [xs (if (sequential? xs) xs (seq xs))]
    (when-not (seq xs)
      (invalid! :xs xs "a non-empty sequence of numbers"))
    (when-not (or (nil? derivatives) (sequential? derivatives))
      (invalid! :derivatives derivatives "a sequence of functions (f' f'' ...)"))
    (let [xs (vec xs)
          axs (m/seq->double-array xs)
          runs (run-lengths axs xs)
          values (node-values f axs xs derivatives runs)]
      [axs values])))

(defn divided-all
  "Divided differences of `f` at the nodes `xs`: all levels of the divided differences table.

  The level `k` consists of the divided differences `f[x_i, ..., x_(i+k)]` of `k + 1` consecutive nodes. Repeated nodes are allowed and give the Hermite (confluent) divided differences: for `k + 1` equal nodes the value is `f^(k)(x) / k!`.

  Parameters:

  - `f` (function): the function of a node. A vector or a map (lookup by the node) is also accepted as `f`, but then the nodes have to be distinct or `:derivatives` has to be given.
  - `xs` (sequence of numbers): the nodes, at least one. Equal nodes must be adjacent.
  - `options` (map, optional):
    - `:derivatives` - sequence of functions `[f' f'' ...]` giving the derivatives of `f`; the element `k - 1` is the `k`-th derivative. It has to contain the derivatives up to the order of the longest run of equal nodes minus one. When it is not given, the derivatives are approximated with [[fastmath.calculus.finite/derivative]] with Richardson extrapolation (`:extrapolate? true`) started with the step `0.5` (divided by 10, up to 8 times, when `f` is not finite so far from the node). It is accurate to about 1e-6 or better for the first and the second derivative of smooth functions, and gets worse for higher orders. `f` has to be a smooth function of a double.

  Returns a list of levels in the order from the highest to the lowest: the first element is a list with the single divided difference of all nodes, the last element is the list of the values of `f` at the nodes. Every level is a list of doubles. For a single node the result is `((f(x0)))`. The first (highest) element is the value returned by [[divided]].

  Throws `ex-info` with `{:option :value :expected}` when `xs` is empty, equal nodes are not adjacent, `f` is not callable, `f` is not a function while nodes are repeated and `:derivatives` is not given, or `:derivatives` is too short.

  See also [[divided]]."
  ([f xs] (divided-all f xs nil))
  ([f xs {:keys [derivatives]}]
   (let [[^doubles axs values] (prepare f xs derivatives)
         n (alength axs)
         ys (double-array (map first values))]
     (loop [level 1
            prev ys
            buff (list (seq ys))]
       (if (m/< level n)
         (let [^doubles prev prev
               cnt (m/long-sub n level)
               cur (double-array cnt)]
           (dotimes [i cnt]
             (let [x0 (Array/aget axs i)
                   x1 (Array/aget axs (m/long-add i level))]
               (Array/aset cur i (if (m/== x0 x1)
                                   (double ((values i) level))
                                   (m// (m/- (Array/aget prev (m/long-inc i))
                                             (Array/aget prev i))
                                        (m/- x1 x0))))))
           (recur (m/long-inc level) cur (conj buff (seq cur))))
         buff)))))

(defn divided
  "Divided difference `f[x0, x1, ..., xn]` of `f` at the nodes `xs`, the highest level of the table.

  It is the leading coefficient of the Newton interpolation polynomial through the nodes. Repeated nodes are allowed and give the Hermite (confluent) divided differences: for `k + 1` equal nodes the value is `f^(k)(x) / k!`.

  Parameters:

  - `f` (function): the function of a node. A vector or a map (lookup by the node) is also accepted as `f`, but then the nodes have to be distinct or `:derivatives` has to be given.
  - `xs` (sequence of numbers): the nodes, at least one. Equal nodes must be adjacent.
  - `options` (map, optional):
    - `:derivatives` - sequence of functions `[f' f'' ...]` giving the derivatives of `f`; the element `k - 1` is the `k`-th derivative. It has to contain the derivatives up to the order of the longest run of equal nodes minus one. When it is not given, the derivatives are approximated with [[fastmath.calculus.finite/derivative]] with Richardson extrapolation (`:extrapolate? true`) started with the step `0.5` (divided by 10, up to 8 times, when `f` is not finite so far from the node). It is accurate to about 1e-6 or better for the first and the second derivative of smooth functions, and gets worse for higher orders. `f` has to be a smooth function of a double.

  Returns the divided difference as a double. For a single node it is `f(x0)`.

  Throws `ex-info` with `{:option :value :expected}` when `xs` is empty, equal nodes are not adjacent, `f` is not callable, `f` is not a function while nodes are repeated and `:derivatives` is not given, or `:derivatives` is too short.

  See also [[divided-all]]."
  (^double [f xs] (divided f xs nil))
  (^double [f xs {:keys [derivatives]}]
   (let [[^doubles axs values] (prepare f xs derivatives)
         n (alength axs)
         ys (double-array (map first values))]
     (loop [level 1]
       (when (m/< level n)
         (dotimes [j (m/long-sub n level)]
           (let [x0 (Array/aget axs j)
                 x1 (Array/aget axs (m/long-add j level))]
             (Array/aset ys j (if (m/== x0 x1)
                                (double ((values j) level))
                                (m// (m/- (Array/aget ys (m/long-inc j))
                                          (Array/aget ys j))
                                     (m/- x1 x0))))))
         (recur (m/long-inc level))))
     (Array/aget ys 0))))
