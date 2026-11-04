(ns ^:no-doc fastmath.calculus.autodiff.rules
  "Internal: the chain rule of the derivative rules of [[fastmath.calculus.autodiff]] (the `chain` macro) and the functions made with [[fastmath.calculus.autodiff/rule]]."
  (:require [fastmath.core :as m]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :as ops])
  (:import [fastmath.calculus.autodiff.dual Dual]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; ## Dispatch

(defmacro dispatch
  "A dual in any argument (symbols) selects `(dual-rule args...)`, otherwise the arguments are coerced to doubles and
  passed to `double-fn` (a function of `fastmath.core` or `fastmath.special`)."
  [dual-rule double-fn & args]
  `(if (or ~@(map (fn [a] `(instance? Dual ~a)) args))
     (~dual-rule ~@args)
     (~double-fn ~@(map (fn [a] `(.doubleValue ~(with-meta a {:tag 'Number}))) args))))

;; ## Chain rule helpers

(defn outside-domain?
  "Is the value `y` `NaN` while no input is `NaN`? Then the function is not defined at the inputs and neither are its
  derivatives, which are made `NaN` (a finite derivative formula would give a number outside the domain)."
  ([y a] (and (m/nan? (d/scalar y)) (not (m/nan? (d/scalar a)))))
  ([y a b] (and (m/nan? (d/scalar y)) (not (m/nan? (d/scalar a))) (not (m/nan? (d/scalar b)))))
  ([y a b c] (and (m/nan? (d/scalar y)) (not (m/nan? (d/scalar a))) (not (m/nan? (d/scalar b))) (not (m/nan? (d/scalar c)))))
  ([y a b c & more] (and (outside-domain? y a b c) (not-any? (fn [x] (m/nan? (d/scalar x))) more))))

(defn outside-domain-seq?
  "[[outside-domain?]] of the value `y` and the inputs `xs` (a sequence)."
  [y xs]
  (and (m/nan? (d/scalar y)) (not-any? (fn [x] (m/nan? (d/scalar x))) xs)))

(defn guard
  "The derivative factor `df`, times `NaN` outside the domain. The factor is scaled, not the partials, so the
  perturbations of the outer layers carried by the factor are `NaN` too and all higher derivatives are `NaN`."
  [df outside?]
  (if outside? (ops/mul df ##NaN) df))

(defn no-derivative!
  "An argument without a derivative (a shape parameter) accepts a constant or a dual with zero partials in the layer,
  otherwise throws."
  [function argument p]
  (when-not (d/p-zero? p)
    (throw (ex-info (str "The derivative of " (name function) " in the argument " (name argument) " is not available.")
                    {:function function :argument argument}))))

(defmacro chain
  "The dual result of a function of a fixed number of arguments, at least one of them a dual, in the layer of the newest tag.

  `(chain [a x, b y, value (f a b)] :f dfa dfb)`: each argument pair binds the value of the argument in the layer
  (a constant is itself), the last pair binds the value of the function. The body is the name of the function (a
  keyword, for the error data) and the derivatives in the order of the arguments, expressions of the bound values
  written with the operations on duals, so the rules nest. A derivative is evaluated only when its argument has
  partials in the layer; `nil` marks an argument without a derivative ([[no-derivative!]]). Outside the domain the
  derivatives are `NaN` ([[outside-domain?]]).

  Shaped as `let` for the linters: the bindings are in scope of the body."
  [bindings function & derivatives]
  (let [pairs (partition 2 bindings)
        args (butlast pairs)
        [y value-expr] (last pairs)
        n (count args)
        gs (repeatedly n #(gensym "arg"))
        ps (repeatedly n #(gensym "partials"))
        t (gensym "tag")
        out (gensym "outside")
        acc (gensym "acc")]
    (assert (and (<= 1 n) (= n (count derivatives))) "at least one argument, one derivative for each")
    `(let [~@(interleave gs (map second args))
           ~t ~(reduce (fn [a b] `(m/long-max ~a ~b)) (map (fn [g] `(d/tag-of ~g)) gs))
           ~@(mapcat (fn [[b] g] [b `(d/value-at ~g ~t)]) args gs)
           ~@(mapcat (fn [p g] [p `(d/partials-at ~g ~t)]) ps gs)
           ~y ~value-expr
           ~out (outside-domain? ~y ~@(map first args))
           ~acc nil
           ~@(mapcat (fn [[b] p df]
                       [acc (if (nil? df)
                              `(do (no-derivative! ~function ~(keyword (name b)) ~p) ~acc)
                              `(if ~p (ops/p-add ~acc (ops/p-scale ~p (guard ~df ~out))) ~acc))])
                     args ps derivatives)]
       (d/new-dual ~t ~y (or ~acc (d/p-zeros (or ~@ps)))))))

;; ## Chain rule of a function of a sequence

(defn- tag-of-seq
  ^long [xs]
  (reduce (fn [^long t x] (m/long-max t (d/tag-of x))) 0 xs))

(defn chain-seq
  "The dual result of a function of a sequence `xs` (the arguments with a derivative) and of `params` (the arguments
  without a derivative), at least one element of them a dual, in the layer of the newest tag.

  - `function` (keyword): the name of the function, for the error data.
  - `xs` (non-empty vector of numbers and duals).
  - `params` (vector of `[name argument]` pairs): each argument accepts a constant or a dual with zero partials in the
    layer ([[no-derivative!]]).
  - `value-fn`: `(value-fn vs pvs)`, the value of the function at the values `vs` of `xs` and `pvs` of the parameters
    in the layer (numbers or duals of the older layers).
  - `derivative-fn`: `(derivative-fn vs pvs y)` returns the derivative of the function in an element as a function of
    the value of that element; it is called once, only when an element has partials in the layer, and its result only
    for such elements.

  Outside the domain the derivatives are `NaN` ([[outside-domain-seq?]])."
  [function xs params value-fn derivative-fn]
  (let [t (m/long-max (tag-of-seq xs) (tag-of-seq (map second params)))
        vs (mapv (fn [x] (d/value-at x t)) xs)
        ps (mapv (fn [x] (d/partials-at x t)) xs)
        pvs (mapv (fn [[argument x]]
                    (no-derivative! function argument (d/partials-at x t))
                    (d/value-at x t)) params)
        y (value-fn vs pvs)
        out (outside-domain-seq? y vs)
        acc (when (some some? ps)
              (let [df (derivative-fn vs pvs y)]
                (reduce (fn [acc i]
                          (if-let [p (ps i)]
                            (ops/p-add acc (ops/p-scale p (guard (df (vs i)) out)))
                            acc))
                        nil (range (count vs)))))]
    (d/new-dual t y (or acc (d/p-zeros (or (some identity ps)
                                           (some (fn [[_ x]] (d/partials-at x t)) params)))))))

;; ## Rules defined by the user (first order)

(defn- first-order-value
  "The value of `x` in the layer `t` as a double; a dual there is a nested derivative, which needs the derivative of
  the derivative, and throws."
  ^double [x ^long t]
  (let [v (d/value-at x t)]
    (if (instance? Dual v)
      (throw (ex-info "A function made with rule supports first-order derivatives only." {:tag t}))
      (.doubleValue ^Number v))))

(defn- user-term
  "Partials `p` times the derivative `(df vx vy)` added to `acc`; `nil` derivative is an argument without a derivative."
  [acc p df vx vy outside? argument]
  (cond
    (nil? p) acc
    (nil? df) (do (no-derivative! :rule argument p) acc)
    :else (ops/p-add acc (ops/p-scale p (guard (if (nil? vy) (df vx) (df vx vy)) outside?)))))

(defn rule1
  "A function of one argument with the value `f` and the derivative `df` (functions of doubles), see [[fastmath.calculus.autodiff/rule]]."
  [f df]
  (fn rule-function [x]
    (if (instance? Dual x)
      (let [^Dual x x
            t (.tag x)
            v (first-order-value x t)
            p (.partials x)
            y (.doubleValue ^Number (f v))]
        (d/new-dual t y (or (user-term nil p df v nil (outside-domain? y v) :x) (d/p-zeros p))))
      (f (d/leaf x)))))

(defn rule2
  "A function of two arguments with the value `f` and the partial derivatives `dfx`, `dfy` (functions of doubles), see [[fastmath.calculus.autodiff/rule]]."
  [f dfx dfy]
  (fn rule-function [x y]
    (if (or (instance? Dual x) (instance? Dual y))
      (let [t (d/top-tag x y)
            vx (first-order-value x t)
            vy (first-order-value y t)
            px (d/partials-at x t)
            py (d/partials-at y t)
            z (.doubleValue ^Number (f vx vy))
            outside? (outside-domain? z vx vy)
            acc (-> nil
                    (user-term px dfx vx vy outside? :x)
                    (user-term py dfy vx vy outside? :y))]
        (d/new-dual t z (or acc (d/p-zeros (or px py)))))
      (f (d/leaf x) (d/leaf y)))))
