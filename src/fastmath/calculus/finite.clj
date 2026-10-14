(ns fastmath.calculus.finite
  "Finite difference derivatives and Richardson extrapolation.

  Univariate derivatives of any order and any order of accuracy with central, forward and backward stencils ([[derivative]]), the stencil weights ([[fd-coeffs]], [[fd-coeffs-for-accuracy]]), multivariate derivatives ([[gradient]], [[hessian]], [[hessian-diagonal]], [[jacobian]], [[laplacian]]) and Richardson extrapolation ([[extrapolate]]).

  When the step `:h` is not given (or is `0`) it is chosen automatically from the machine epsilon and the magnitude of the point, which balances the truncation and the rounding errors. An explicit `:h` is always an absolute step.

  Options: a `nil` option is a missing option. Numeric options are coerced to `long` (an integer option) or `double` and only their ranges are checked, so a fractional `n`, `:acc` or `:max-evals` is truncated and a non-number gives a `ClassCastException`. The step given to a function returned by [[derivative]], [[gradient]], [[jacobian]], [[hessian]], [[hessian-diagonal]] or [[laplacian]] at call time is not checked: a step which is not positive means the automatic step."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.calculus.common :as cc])
  (:import [fastmath.vector Vec3]
           [fastmath.java Array]
           [clojure.lang IFn$DDD]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; validation

(defn- invalid-option!
  [option value expected]
  (throw (ex-info (str "Invalid option " option ": " (pr-str value) ", expected " expected)
                  {:option option :value value :expected expected})))

(def ^:private methods-set #{:central :forward :backward})

(defn- check-n!
  [^long n]
  (when (m/neg? n)
    (invalid-option! :n n "an integer >= 0")))

(defn- check-acc!
  [^long acc]
  (when-not (m/>= acc 1)
    (invalid-option! :acc acc "an integer >= 1")))

(defn- check-method!
  [method]
  (when-not (contains? methods-set method)
    (invalid-option! :method method "one of :central, :forward, :backward")))

(defn- check-h!
  [^double h]
  (when-not (and (m/valid-double? h)
                 (m/not-neg? h))
    (invalid-option! :h h "a finite number >= 0 (0 means automatic step)")))

;; extrapolation

;; https://github.com/JuliaMath/Richardson.jl/blob/master/src/Richardson.jl
(defn extrapolate
  "Richardson extrapolation of a function `g(x, h)` to `h -> 0`, returns the extrapolated function of `x`.

  The function `g` is evaluated for a geometric sequence of steps `h` (started at `:init-h` and multiplied by `:contract` each time). The values are combined with Neville's scheme, which eliminates consecutive powers of `h` from the error expansion. The result is the combination with the smallest error estimate. The method is a port of Richardson.jl.

  Parameters:

  - `g` (function): the function of `x` and `h`, both doubles, returning a double. It has to be defined (and finite) for the steps used.
  - `options` (map, optional):
    - `:contract` - contraction factor of the step, must be in the open interval `(0, 1)` (default: `0.5`).
    - `:power` - the power of `h` in the leading error term; the next terms are assumed to be in powers of `contract^power` steps, e.g. `2.0` for the central differences (default: `1.0`).
    - `:init-h` - the first step, a nonzero number (default: `0.5`).
    - `:rel` - relative tolerance of the error estimate (default: square root of the ulp of `:init-h`).
    - `:abs` - absolute tolerance of the error estimate (default: machine epsilon).
    - `:max-evals` - maximum number of evaluations of `g`, an integer >= 1 (default: `Integer/MAX_VALUE`).
    - `:tol` - the iteration stops when the best error estimate grows by more than this factor with respect to the previous one (default: `2.0`).

  The extrapolation stops when the error estimate is not larger than `max(rel * |result|, abs)`, when the error estimate starts to grow, or when `:max-evals` is reached. In the last case the best value found so far is returned.

  For an infinite `x` the function `g(0, 1/h)` is evaluated for `x = ##Inf` and `g(0, -1/h)` for `x = ##-Inf`, which allows to extrapolate limits at infinity. `NaN` gives `NaN`.

  The stopping rule differs from Richardson.jl: the error estimate of the previous iteration is tested against the tolerance (one iteration later than in the original), therefore `g` can be evaluated once more. The result is not less accurate.

  Throws `ex-info` with `{:option :value :expected}` when `:contract` is not in `(0, 1)`, `:max-evals` is less than 1 or `:init-h` is zero or not finite.

  See also [[derivative]] (the `:extrapolate?` option)."
  ([g] (extrapolate g nil))
  ([g options]
   (let [{:keys [^double contract ^double power ^double init-h rel ^double abs ^long max-evals ^double tol]
          :or {contract 0.5 power 1.0 init-h 0.5
               abs m/MACHINE-EPSILON
               max-evals Integer/MAX_VALUE tol 2.0}} (cc/drop-nil-options options)]
     (when-not (m/< 0.0 contract 1.0)
       (invalid-option! :contract contract "a number in the open interval (0, 1)"))
     (when-not (m/>= max-evals 1)
       (invalid-option! :max-evals max-evals "an integer >= 1"))
     (when-not (and (m/valid-double? init-h)
                    (m/not-zero? init-h))
       (invalid-option! :init-h init-h "a finite nonzero number"))
     (let [rel (double (or rel (m/sqrt (m/ulp init-h))))
           invcontract (m/pow (m// contract) power)
           nf (fn ^double [^double x0 ^double h] (cond
                                                  (m/valid-double? x0) (g x0 h)
                                                  (m/pos-inf? x0) (g 0.0 (m// 1.0 h))
                                                  (m/neg-inf? x0) (g 0.0 (m// -1.0 h))
                                                  :else ##NaN))]
       (fn ^double [^double x0]
         (loop [i (long 1)
                h init-h
                ;; keep in 3-tuple of doubles to avoid destructuring
                ^Vec3 f0+err+invc (Vec3. (nf x0 init-h) ##Inf invcontract) ;; curr val, err, contraction den
                neville (list f0+err+invc)]
           (if (m/>= i max-evals)
             (.x f0+err+invc)
             (let [nh (m/* h contract)
                   nf0 (nf x0 nh)
                   res (reductions (fn [^Vec3 curr ^Vec3 ni+err+invc]
                                     (let [ni (.x ni+err+invc)
                                           ni+1 (.x curr)
                                           c (.z curr)
                                           n (m/+ ni+1 (m// (m/- ni+1 ni) (m/dec c)))]
                                       (Vec3. n (m/abs (m/- n ni)) (m/* c invcontract))))
                                   (Vec3. nf0 ##Inf invcontract) neville)
                   ^Vec3 minimal (reduce (fn [^Vec3 a ^Vec3 b]
                                           (if (m/< (.y a) (.y b)) a b)) res)
                   minimal' (.y minimal)
                   err (.y f0+err+invc)
                   ^Vec3 minimal (if (or (m/invalid-double? minimal')
                                         (m/< err minimal')) f0+err+invc minimal)]
               (if (or (m/invalid-double? minimal')
                       (m/> minimal' (m/* tol err))
                       (m/<= err (m/max (m/* rel (m/abs (.x minimal))) abs)))
                 (.x minimal)
                 (recur (m/long-inc i)
                        nh
                        minimal
                        res))))))))))

;; coefficients

;; generic (boxed) arithmetic on exact ratios, the boxed math warning is switched off on purpose
(set! *unchecked-math* false)

;; the promoting operators: the numerators and denominators of high orders do not fit in a long
(defn- rat- [a b] (-' a b))
(defn- rat* [a b] (*' a b))
(defn- rat-div [a b] (/ a b))

(set! *unchecked-math* :warn-on-boxed)

(defn- fornberg-update-row
  "Updates the row of the weights of the node `j` after the node `i` was added, `mn` is `min(i, n)`."
  [row ^long mn c4 c3]
  (reduce (fn [r ^long k]
            (assoc r k (rat-div (rat- (rat* c4 (row k))
                                      (rat* k (row (m/long-dec k))))
                                c3)))
          (assoc row 0 (rat-div (rat* c4 (row 0)) c3))
          (range 1 (m/long-inc mn))))

(defn- fornberg-new-row
  "Row of the weights of the added node `i`."
  [prev-row mn c1 c5 c2]
  (reduce (fn [r ^long k]
            (assoc r k (rat-div (rat* c1 (rat- (rat* k (prev-row (m/long-dec k)))
                                               (rat* c5 (prev-row k))))
                                c2)))
          (assoc (vec (repeat (count prev-row) 0))
                 0 (rat-div (rat- 0 (rat* c1 (rat* c5 (prev-row 0)))) c2))
          (range 1 (m/long-inc (long mn)))))

(defn- fd-weights-exact
  "Weights of the n-th derivative at 0 for the (rational, distinct) nodes `xs` as exact ratios (Fornberg's algorithm)."
  [^long n xs]
  (let [xs (vec xs)
        cnt (count xs)
        row0 (vec (repeat (m/long-inc n) 0))
        c-init (assoc (vec (repeat cnt row0)) 0 (assoc row0 0 1))]
    (loop [i 1
           c1 (rat* 1 1)
           c4 (xs 0)
           c c-init]
      (if (m/== i cnt)
        (mapv (fn [row] (row n)) c)
        (let [mn (m/long-min i n)
              xi (xs i)
              [c2 c] (reduce (fn [[c2 c] ^long j]
                               (let [c3 (rat- xi (xs j))
                                     c2 (rat* c2 c3)
                                     c (if (m/== j (m/long-dec i))
                                         (assoc c i (fornberg-new-row (c j) mn c1 c4 c2))
                                         c)]
                                 [c2 (assoc c j (fornberg-update-row (c j) mn xi c3))]))
                             [1 c]
                             (range i))]
          (recur (m/long-inc i) c2 xi c))))))

(defn- ratio->double
  ^double [x]
  (if (ratio? x)
    (let [^java.math.BigInteger bn (biginteger (numerator x))
          ^java.math.BigInteger bd (biginteger (denominator x))]
      (if (and (m/<= (.bitLength bn) 53) (m/<= (.bitLength bd) 53))
        (m// (.doubleValue bn) (.doubleValue bd)) ;; correctly rounded
        (double x)))
    (double x)))

(defn fd-coeffs
  "Weights of the finite difference approximation of the `n`-th derivative at a point for given offsets of the stencil.

  The weights are computed with Fornberg's algorithm in exact rational arithmetic and converted to doubles, so they are accurate to the last digit for any order and stencil. The derivative is approximated by `sum(coeff_i * f(x + offset_i * h)) / h^n`.

  Parameters:

  - `n` (integer >= 0): the order of the derivative.
  - `offsets` (sequence of numbers): distinct offsets of the nodes in units of the step `h`, e.g. `[-1 0 1]`. Doubles are converted to exact ratios.

  Returns `[offsets coeffs]`, where `coeffs` is a vector of doubles of the same length as `offsets`. When there are not more than `n` offsets, the derivative is not determined and all weights are `0.0`.

  Throws `ex-info` with `{:option :value :expected}` when `n` is negative or when `offsets` is empty, contains duplicates or non-finite numbers.

  See also [[fd-coeffs-for-accuracy]], [[derivative]]."
  [n offsets]
  (check-n! n)
  (when-not (and (sequential? offsets)
                 (seq offsets)
                 (every? #(and (number? %) (m/valid-double? (double %))) offsets)
                 (apply distinct? (map rationalize offsets)))
    (invalid-option! :offsets offsets "a non-empty sequence of distinct finite numbers"))
  [offsets (mapv ratio->double (fd-weights-exact (long n) (map rationalize offsets)))])

(defn fd-coeffs-for-accuracy
  "Offsets and weights of the finite difference stencil for the `n`-th derivative with a given order of accuracy.

  Parameters:

  - `n` (integer >= 0): the order of the derivative (default: `1`).
  - `acc` (integer >= 1): the order of accuracy, the error is `O(h^acc)` (default: `2`). The central stencil gives an even order, an odd `acc` results in the next even accuracy.
  - `kind` (keyword): `:central` (default), `:forward` or `:backward`.

  Returns `[offsets coeffs]` as [[fd-coeffs]] does: the integer offsets in units of the step `h` and the weights of the stencil (divide the weighted sum by `h^n`).

  Throws `ex-info` with `{:option :value :expected}` when `n` is negative, `acc` is less than 1 or `kind` is unknown.

  See also [[fd-coeffs]], [[derivative]]."
  ([^long n] (fd-coeffs-for-accuracy n 2))
  ([^long n ^long acc] (fd-coeffs-for-accuracy n acc :central))
  ([^long n ^long acc kind]
   (check-n! n)
   (check-acc! acc)
   (check-method! kind)
   (let [num-central (m/long-add (m/long-sub (m/long-mult 2 (m/long-quot (m/long-inc n) 2)) 1) acc)
         num-side (m/long-quot num-central 2)
         num-coeffs (if (m/even? n) (m/long-inc num-central) num-central)]
     (fd-coeffs n (case kind
                    :forward (range num-coeffs)
                    :backward (range (m/long-inc (m/long-sub 0 num-coeffs)) 1)
                    (range (m/long-sub 0 num-side) (m/long-inc num-side)))))))

;; derivative

(defn- effective-acc
  "Order of accuracy really obtained: the next even number for the central stencil."
  ^long [^long acc method]
  (if (and (identical? method :central) (m/odd? acc))
    (m/long-inc acc)
    acc))

(defn- auto-step-base
  "Automatic step for `|x| <= 10`: `eps^(1/(n+acc))`."
  ^double [^long n ^long acc method]
  (m/pow m/MACHINE-EPSILON (m// 1.0 (m/long-add n (effective-acc acc method)))))

(def ^{:private true :const true :tag 'double} step-scale-start
  "The automatic step is the same up to this magnitude of the point. Functions of order-one variability (exp, sin) need no
  larger step there, while a step proportional to the point (the right one for powers and logarithms) is chosen above it.
  Measured on exp and sin for `|x|` up to 5, orders 1-4 and accuracies 2, 4, 6 the step never has a larger error than the
  constant step `eps^(1/(n+2))`; scaling from 1 instead of 10 made it up to 27 times worse for `acc` 2."
  10.0)

(defn- scaled-step
  "Automatic step for the point `x`: `base` for `|x| <= 10`, proportional to `|x|` above it."
  ^double [^double base ^double x]
  (m/* base (m/max 1.0 (m// (m/abs x) step-scale-start))))

(defmacro ^:private produce-symbols
  "Stencil `(fn [x h])` returning the weighted sum of the values (not divided by `h^n`)."
  [n acc kind]
  (let [[offsets coeffs] (fd-coeffs-for-accuracy n acc kind)
        x (gensym "x")
        h (gensym "h")
        args (with-meta [(with-meta x {:tag 'double}) (with-meta h {:tag 'double})] {:tag 'double})
        fname (gensym (str "local-" n "-" acc "-" (name kind)))
        point (fn [id] (case (long id)
                         0 x
                         1 `(m/+ ~x ~h)
                         -1 `(m/- ~x ~h)
                         `(m/+ ~x (m/* ~(double id) ~h))))
        terms (for [[id coeff] (map vector offsets coeffs)
                    :when (not (m/zero? (double coeff)))]
                `(m/* ~(double coeff) (double (~'f ~(point id)))))]
    `(fn ~fname ~args
       (m/+ ~@terms))))

(defmacro ^:private produce-case-for-methods
  [n acc]
  `(case ~'method
     :forward (produce-symbols ~n ~acc :forward)
     :backward (produce-symbols ~n ~acc :backward)
     (produce-symbols ~n ~acc :central)))

(defn- generic-stencil
  "Stencil of any order and accuracy."
  [f ^long n ^long acc method]
  (let [[offsets coeffs] (fd-coeffs-for-accuracy n acc method)
        ^doubles offs (m/seq->double-array offsets)
        ^doubles cfs (m/seq->double-array coeffs)
        cnt (alength offs)]
    (fn ^double [^double x ^double h]
      (loop [i 0
             s 0.0]
        (if (m/< i cnt)
          (recur (m/long-inc i)
                 (m/+ s (m/* (Array/aget cfs i) (double (f (m/+ x (m/* (Array/aget offs i) h)))))))
          s)))))

(defmacro ^:private produce-cases-for-accuracy
  [n]
  `(case (unchecked-int ~'acc)
     1 (produce-case-for-methods ~n 1)
     2 (produce-case-for-methods ~n 2)
     3 (produce-case-for-methods ~n 3)
     4 (produce-case-for-methods ~n 4)
     5 (produce-case-for-methods ~n 5)
     6 (produce-case-for-methods ~n 6)
     (generic-stencil ~'f ~n ~'acc ~'method)))

(defn derivative
  "Creates the `n`-th derivative of a univariate function `f` approximated with finite differences.

  Parameters:

  - `f` (function): double -> double.
  - `n` (integer >= 0): the order of the derivative (default: `1`). For `0` the function `f` itself is returned.
  - `options` (map, optional):
    - `:acc` - order of accuracy, an integer >= 1 (default: `2`). The error of the approximation is `O(h^acc)`; an odd `acc` for the central stencil gives the next even accuracy.
    - `:h` - step, a number >= 0 (default: `0`). `0` means the automatic step `eps^(1/(n+acc)) * max(1, |x|/10)`, where `eps` is the machine epsilon and `acc` is the accuracy (for the central method rounded up to even); this balances the truncation and rounding errors and grows with the magnitude of the point above `10`. An explicit positive `:h` is an absolute step.
    - `:method` - `:central` (default), `:forward` or `:backward`.
    - `:extrapolate?` - when true or a map with options of [[extrapolate]] the derivative is Richardson-extrapolated to `h -> 0` (default: `false`); the extrapolation starts with the step given in its options (`:init-h`, default `0.5`) and `:power` is `2` for the central and `1` otherwise, unless given.

  Returns a function of `x` returning the derivative as a double. The function has also the arities `(d x h)` with the absolute step `h` (a step which is not positive means automatic) and `(d x h hn)` with the step and its precomputed power `h^n`. It is the function `(d x h)` that [[extrapolate]] uses.

  Throws `ex-info` with `{:option :value :expected}` when `n` is negative, `:acc` is less than 1, `:method` is unknown or `:h` is negative or not finite.

  Higher orders and high accuracies need smaller relative steps and are more sensitive to rounding errors: the attainable relative error grows quickly with `n`.

  See also [[gradient]], [[hessian]], [[extrapolate]], [[fd-coeffs-for-accuracy]]."
  ([f] (derivative f 1))
  ([f ^long n] (derivative f n nil))
  ([f ^long n options]
   (let [{:keys [^long acc ^double h method extrapolate?]
          :or {acc 2 h 0.0 method :central}} (cc/drop-nil-options options)]
     (check-n! n)
     (check-acc! acc)
     (check-method! method)
     (check-h! h)
     (if (m/zero? n)
       f
       (let [hn (m/fpow h n)
             fixed? (m/pos? h)
             base (auto-step-base n acc method)
             ^IFn$DDD stencil (case n
                                1 (produce-cases-for-accuracy 1)
                                2 (produce-cases-for-accuracy 2)
                                3 (produce-cases-for-accuracy 3)
                                4 (produce-cases-for-accuracy 4)
                                5 (produce-cases-for-accuracy 5)
                                6 (produce-cases-for-accuracy 6)
                                (generic-stencil f n acc method))
             result (fn local-derivative
                      (^double [^double x]
                       (if fixed?
                         (m// (.invokePrim stencil x h) hn)
                         (let [h (scaled-step base x)]
                           (m// (.invokePrim stencil x h) (m/fpow h n)))))
                      (^double [^double x ^double h2]
                       (let [h2 (if (m/pos? h2) h2 (if fixed? h (scaled-step base x)))]
                         (m// (.invokePrim stencil x h2) (m/fpow h2 n))))
                      (^double [^double x ^double h2 ^double hn2]
                       (m// (.invokePrim stencil x h2) hn2)))]
         (if extrapolate?
           (let [eoptions (if (map? extrapolate?) extrapolate? {})]
             (extrapolate result (if (:power eoptions)
                                   eoptions
                                   (assoc eoptions :power (if (= method :central) 2 1)))))
           result))))))

;; multivariate

(def ^{:private true :const true :tag 'double} hessian-base (m/pow m/MACHINE-EPSILON 0.25))

(defn- central-stencil
  "Nonzero offsets and weights of the central stencil as `[offsets coeffs]` double arrays."
  [^long n ^long acc]
  (let [[offsets coeffs] (fd-coeffs-for-accuracy n acc :central)
        pairs (filter (fn [[_ ^double c]] (m/not-zero? c)) (map vector offsets coeffs))]
    [(m/seq->double-array (map first pairs))
     (m/seq->double-array (map second pairs))]))

(defn- partial-derivative
  "Partial derivative of scalar valued `f` for the i-th component of the vector `v` (vector of doubles)."
  [f v i hi ^doubles offs ^doubles cfs]
  (let [vi (double (v i))
        hi (double hi)
        cnt (alength offs)]
    (loop [k (long 0)
           s 0.0]
      (if (m/< k cnt)
        (recur (m/long-inc k)
               (m/+ s (m/* (Array/aget cfs k)
                           (double (f (assoc v i (m/+ vi (m/* (Array/aget offs k) hi))))))))
        (m// s hi)))))

(defn- values-vector
  "Values of the vector valued `f` at `v` as a vector, throws `ex-info` when `f` does not return a sequence or an array."
  [f v]
  (let [r (f v)]
    (when-not (and (some? r) (or (sequential? r) (.isArray (class r))))
      (invalid-option! :f r "a function returning a sequence of numbers"))
    (vec r)))

(defn- partial-derivative-vector
  "Partial derivatives of the vector valued `f` for the i-th component of `v`, throws `ex-info` when the number of values changes."
  [f v i hi ^doubles offs ^doubles cfs]
  (let [vi (double (v i))
        hi (double hi)
        cnt (alength offs)]
    (loop [k 0
           acc nil]
      (if (m/< k cnt)
        (let [c (Array/aget cfs k)
              fv (values-vector f (assoc v i (m/+ vi (m/* (Array/aget offs k) hi))))
              scaled (v/mult fv c)]
          (when (and acc (not= (count acc) (count scaled)))
            (invalid-option! :f (count scaled) (str "a function returning " (count acc) " values")))
          (recur (m/long-inc k)
                 (if acc
                   (v/add acc scaled)
                   scaled)))
        (v/div acc hi)))))

(defn gradient
  "Creates the gradient (vector of the first partial derivatives) of a multivariate scalar function `f`, approximated with central finite differences.

  Parameters:

  - `f` (function): a vector of doubles -> double.
  - `options` (map, optional):
    - `:acc` - order of accuracy, an integer >= 1 (default: `2`); an odd value gives the next even accuracy.
    - `:h` - step, a number >= 0 (default: `0`). `0` means the automatic step for each component `i`: `eps^(1/(1+acc)) * max(1, |x_i|/10)` (with `acc` rounded up to even). An explicit positive `:h` is an absolute step used for every component.

  Returns a function of a point `v` (any sequence of numbers or a `double[]`) returning the gradient as a vector of doubles. The function has also the arity `(g v h)` with an explicit absolute step (a step which is not positive means automatic). The function `f` receives a vector and is called when the gradient is called, the result is not lazy. The gradient of the empty point is `[]`.

  Throws `ex-info` with `{:option :value :expected}` when `:acc` is less than 1 or `:h` is negative or not finite.

  See also [[hessian]], [[jacobian]], [[derivative]]."
  ([f] (gradient f nil))
  ([f options]
   (let [{:keys [^double h ^long acc]
          :or {h 0.0 acc 2}} (cc/drop-nil-options options)]
     (check-acc! acc)
     (check-h! h)
     (let [[^doubles offs ^doubles cfs] (central-stencil 1 acc)
           base (auto-step-base 1 acc :central)]
       (fn local-gradient
         ([v] (local-gradient v h))
         ([v ^double h2]
          (let [v (vec v)]
            (mapv (fn [^long i]
                    (partial-derivative f v i (if (m/pos? h2) h2 (scaled-step base (v i))) offs cfs))
                  (range (count v))))))))))

(defn jacobian
  "Creates the Jacobian matrix of a vector valued function `f`, approximated with central finite differences.

  Parameters:

  - `f` (function): a vector of `n` doubles -> a sequence (or array) of `m` numbers.
  - `options` (map, optional): as in [[gradient]]
    - `:acc` - order of accuracy, an integer >= 1 (default: `2`).
    - `:h` - step, a number >= 0 (default: `0` = automatic step for each component, `eps^(1/(1+acc)) * max(1, |x_i|/10)`); an explicit positive value is an absolute step.

  Returns a function of a point `v` (any sequence of numbers or a `double[]`) returning the matrix as a vector of `m` rows, each a vector of `n` partial derivatives `d f_j / d x_i`. The function has also the arity `(j v h)` with an explicit absolute step (a step which is not positive means automatic). The function `f` is called when the Jacobian is called. For the empty point the result is a vector of `m` empty rows.

  Throws `ex-info` with `{:option :value :expected}` when `:acc` is less than 1, `:h` is negative or not finite, `f` returns a number instead of a sequence or the number of returned values differs between the calls.

  See also [[gradient]], [[hessian]]."
  ([f] (jacobian f nil))
  ([f options]
   (let [{:keys [^double h ^long acc]
          :or {h 0.0 acc 2}} (cc/drop-nil-options options)]
     (check-acc! acc)
     (check-h! h)
     (let [[^doubles offs ^doubles cfs] (central-stencil 1 acc)
           base (auto-step-base 1 acc :central)]
       (fn local-jacobian
         ([v] (local-jacobian v h))
         ([v ^double h2]
          (let [v (vec v)
                columns (mapv (fn [^long i]
                                (partial-derivative-vector f v i (if (m/pos? h2) h2 (scaled-step base (v i))) offs cfs))
                              (range (count v)))]
            (if (seq columns)
              (do (when-not (apply = (map count columns))
                    (invalid-option! :f (mapv count columns) "a function returning the same number of values at every point"))
                  (apply mapv vector columns))
              (mapv (constantly []) (values-vector f v))))))))))

(defn hessian
  "Creates the Hessian matrix of a multivariate scalar function `f`, approximated with central finite differences of the second order.

  Parameters:

  - `f` (function): a vector of doubles -> double.
  - `options` (map, optional):
    - `:h` - step, a number >= 0 (default: `0`). `0` means the automatic step for each component `i`: `eps^(1/4) * max(1, |x_i|/10)`. An explicit positive `:h` is an absolute step used for every component.

  Returns a function of a point `v` (any sequence of numbers or a `double[]`) returning the symmetric matrix as a vector of vectors of doubles. The function has also the arity `(hs v h)` with an explicit absolute step (a step which is not positive means automatic). The diagonal term is `(f(x+h_i) - 2 f(x) + f(x-h_i)) / h_i^2` and the mixed term is `(f(++) - f(+-) - f(-+) + f(--)) / (4 h_i h_j)`. The function `f` is called when the Hessian is called. The Hessian of the empty point is `[]`.

  Throws `ex-info` with `{:option :value :expected}` when `:h` is negative or not finite.

  See also [[hessian-diagonal]], [[laplacian]], [[gradient]]."
  ([f] (hessian f nil))
  ([f options]
   (let [{:keys [^double h] :or {h 0.0}} (cc/drop-nil-options options)]
     (check-h! h)
     (fn local-hessian
       ([v] (local-hessian v h))
       ([v ^double h2]
        (let [v (vec v)
              cv (count v)
              hs (mapv (fn [x] (if (m/pos? h2) h2 (scaled-step hessian-base x))) v)
              buff (double-array (m/long-mult cv cv))]
          (when (m/pos? cv)
            (let [fv-2 (m/* -2.0 (double (f v)))]
              (dotimes [i cv]
                (dotimes [j cv]
                  (when (m/<= i j)
                    (let [^double hi (hs i)
                          ^double hj (hs j)
                          ^double vi (v i)
                          ^double vj (v j)]
                      (if (m/== i j)
                        (let [^double x1 (f (assoc v i (m/+ vi hi)))
                              ^double x2 (f (assoc v i (m/- vi hi)))]
                          (Array/set buff (m/long-add (m/long-mult i cv) i)
                                     (m/* (m/+ x1 fv-2 x2) (m// (m/* hi hi)))))
                        (let [^double x1 (f (assoc v i (m/+ vi hi) j (m/+ vj hj)))
                              ^double x2 (f (assoc v i (m/+ vi hi) j (m/- vj hj)))
                              ^double x3 (f (assoc v i (m/- vi hi) j (m/+ vj hj)))
                              ^double x4 (f (assoc v i (m/- vi hi) j (m/- vj hj)))
                              val (m/* (m/- (m/+ x1 x4) x2 x3) (m// (m/* 4.0 hi hj)))]
                          (Array/set buff (m/long-add (m/long-mult i cv) j) val)
                          (Array/set buff (m/long-add (m/long-mult j cv) i) val)))))))))
          (mapv (fn [^long i]
                  (mapv (fn [^long j] (Array/aget buff (m/long-add (m/long-mult i cv) j)))
                        (range cv)))
                (range cv))))))))

(defn hessian-diagonal
  "Creates the function returning the diagonal of the Hessian matrix of a multivariate scalar function `f`, approximated with central finite differences of the second order.

  Parameters:

  - `f` (function): a vector of doubles -> double.
  - `options` (map, optional):
    - `:h` - step, a number >= 0 (default: `0`). `0` means the automatic step for each component `i`: `eps^(1/4) * max(1, |x_i|/10)`. An explicit positive `:h` is an absolute step used for every component.

  Returns a function of a point `v` (any sequence of numbers or a `double[]`) returning the diagonal `(f(x+h_i) - 2 f(x) + f(x-h_i)) / h_i^2` as a `double[]`. The function has also the arity `(hd v h)` with an explicit absolute step (a step which is not positive means automatic). The function `f` is called when the diagonal is computed.

  Throws `ex-info` with `{:option :value :expected}` when `:h` is negative or not finite.

  See also [[hessian]], [[laplacian]]."
  ([f] (hessian-diagonal f nil))
  ([f options]
   (let [{:keys [^double h] :or {h 0.0}} (cc/drop-nil-options options)]
     (check-h! h)
     (fn local-hessian-diagonal
       ([v] (local-hessian-diagonal v h))
       ([v ^double h2]
        (let [v (vec v)
              cv (count v)
              buff (double-array cv)]
          (when (m/pos? cv)
            (let [fv-2 (m/* -2.0 (double (f v)))]
              (dotimes [i cv]
                (let [^double vi (v i)
                      hi (if (m/pos? h2) h2 (scaled-step hessian-base vi))
                      ^double x1 (f (assoc v i (m/+ vi hi)))
                      ^double x2 (f (assoc v i (m/- vi hi)))]
                  (Array/set buff i (m/* (m/+ x1 fv-2 x2) (m// (m/* hi hi))))))))
          buff))))))

(defn laplacian
  "Creates the Laplacian (the sum of the second partial derivatives) of a multivariate scalar function `f`, approximated with central finite differences.

  Parameters:

  - `f` (function): a vector of doubles -> double.
  - `options` (map, optional): as in [[hessian-diagonal]]
    - `:h` - step, a number >= 0 (default: `0` = automatic step for each component, `eps^(1/4) * max(1, |x_i|/10)`); an explicit positive value is an absolute step.

  Returns a function of a point `v` (any sequence of numbers or a `double[]`) returning the Laplacian as a double. The function has also the arity `(l v h)` with an explicit absolute step (a step which is not positive means automatic). The Laplacian of the empty point is `0.0`.

  Throws `ex-info` with `{:option :value :expected}` when `:h` is negative or not finite.

  See also [[hessian-diagonal]], [[hessian]]."
  ([f] (laplacian f nil))
  ([f options]
   (let [hd (hessian-diagonal f options)]
     (fn local-laplacian
       (^double [v] (v/sum (hd v)))
       (^double [v ^double h] (v/sum (hd v h)))))))
