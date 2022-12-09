;; # Namespace scope
;;
;; Collection of math function:
;;
;; * Several constants from Java, C, Processing, etc.
;; * Functions based on FastMath exposed as macros or functions (trigonometry, powers/logarithms/roots, rounding)
;; * Primitive operators (as in primitive-math package)
;; * Additional math functions (signum, constrain, interpolation)
;; * Statistics

(ns fastmath.core
  "Collection of fast math functions and plethora of constants known from other math libraries.

  ### Primitive math operators

  Based on [Primitive Math by Zach Tellman](https://github.com/ztellman/primitive-math) several operators are introduced and replace `clojure.core` functions. All operators are macros and can't be used as functions. List includes:

  Known from Clojure: `*` `+` `-` `/` `>` `<` `>=` `<=` `==` `rem` `quot` `mod` `bit-or` `bit-and` `bit-xor` `bit-not` `bit-shift-left` `bit-shift-right` `unsigned-bit-shift-right` `inc` `dec` `zero?` `neg?` `pos?` `min` `max` `even?` `odd?` `abs`

  And additionally:

  * `bool-and` - `and` working on booleans
  * `bool-or` - boolean `or`
  * `bool-xor` - boolean `xor`
  * `bool-not` - boolean `not`
  * `<<` - bit shift left
  * `>>` - signed bit shift right
  * `>>>` - unsigned bit shift right
  * `not==` - not equal

  Warning: All `bool-` evaluate all parameters.
  
  To turn on primitive math on your namespace call [[use-primitive-operators]].
  To turn off and revert original versions call [[unuse-primitive-operators]]

  ### Fast Math

  Almost all math functions are backed by [FastMath](https://github.com/jeffhain/jafama) library. Most of them are macros. Some of them are wrapped in Clojure functions. Almost all operates on primitive `double` and returns `double` (with an exception [[round]] or [[qround]] which returns `long`).

  ### Other functions

  Additionally namespace contains functions which are common in frameworks like OpenFrameworks and Processing."
  {:metadoc/categories {:trig "Trigonometry"
                        :pow "Powers / logarithms"
                        :conv "Conversions"
                        :err "Error"
                        :dist "Distance"
                        :round "Rounding"
                        :sign "Sign"
                        :bitwise "Bitwise"
                        :mod "Mod"
                        :compare "Comparison"
                        :prim "Primitive"
                        :special "Special functions"
                        :stat "Statistical"
                        :bool "Boolean"
                        :rank "Rank"
                        :seq "Primitive <-> Seq converters"
                        :sample "Sampling"}}
  (:refer-clojure
   :exclude [* + - / > < >= <= == rem quot mod bit-or bit-and bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd? abs])
  (:import [net.jafama FastMath]
           [fastmath.java PrimitiveMath]
           [org.apache.commons.math3.util Precision]
           [org.apache.commons.math3.special Gamma])
  (:require [fastmath.core :as m]))

(set! *unchecked-math* :warn-on-boxed)

;; metadoc sets

(def ^:const ^:private -prim-set- #{:prim})
(def ^:const ^:private -trig-set- #{:trig})
(def ^:const ^:private -conv-set- #{:conv})
(def ^:const ^:private -err-set- #{:err})
(def ^:const ^:private -pow-set- #{:pow})
(def ^:const ^:private -bitwise-set- #{:bitwise})
(def ^:const ^:private -compare-set- #{:compare})
(def ^:const ^:private -round-set- #{:round})
(def ^:const ^:private -special-set- #{:special})
(def ^:const ^:private -mod-set- #{:mod})

;; which java?

(def ^:private jvm-version-type (if (= "1." (subs (System/getProperty "java.version") 0 2)) :old :new))

;; ## Macros

(defmacro ^:private javaclass-proxy
  "Wrapps operation into macro"
  ([class arity alt-name name]
   (let [cf (str class "/" name)
         f (symbol cf)
         x (symbol "x")
         y (symbol "y")
         z (symbol "z")
         doc (or (:doc (meta alt-name)) (str cf " function wrapped in macro."))
         arity-1 `([~x] (list '~f ~x))
         arity-2 `([~x ~y] (list '~f ~x ~y))
         arity-3 `([~x ~y ~z] (list '~f ~x ~y ~z))]
     (condp = arity
       :two `(defmacro ~alt-name ~doc ~arity-2)
       :onetwo `(defmacro ~alt-name ~doc ~arity-1 ~arity-2)
       :three `(defmacro ~alt-name ~doc ~arity-3)
       `(defmacro ~alt-name ~doc ~arity-1))))
  ([class arity name]
   `(javaclass-proxy ~class ~arity ~name ~name)))

(defmacro ^:private fastmath-proxy [& rest] `(javaclass-proxy "net.jafama.FastMath" ~@rest))
(defmacro ^:private primitivemath-proxy [& rest] `(javaclass-proxy "fastmath.java.PrimitiveMath" ~@rest))
(defmacro ^:private erf-proxy [& rest] `(javaclass-proxy "org.apache.commons.math3.special.Erf" ~@rest))
(defmacro ^:private gamma-proxy [& rest] `(javaclass-proxy "org.apache.commons.math3.special.Gamma" ~@rest))
(defmacro ^:private beta-proxy [& rest] `(javaclass-proxy "org.apache.commons.math3.special.Beta" ~@rest))
(defmacro ^:private besselj-proxy [& rest] `(javaclass-proxy "org.apache.commons.math3.special.BesselJ" ~@rest))

(defmacro ^:private variadic-proxy
  "Creates left-associative variadic forms for any operator.
  https://github.com/ztellman/primitive-math/blob/master/src/primitive_math.clj#L10"
  ([name]
   `(variadic-proxy ~name ~name))
  ([name fn]
   `(variadic-proxy ~name ~fn identity))
  ([name fn single-arg-form]
   (let [x (symbol "x")
         y (symbol "y")
         rest (symbol "rest")
         fname (symbol (str "fastmath.java.PrimitiveMath/" fn))
         qname (symbol (str "fastmath.core/" name))
         doc (or (:doc (meta name)) (str "A primitive math version of `" name "`"))]
     `(defmacro ~name
        ~doc
        ([~x]
         ~((eval single-arg-form) x))
        ([~x ~y]
         (list '~fname ~x ~y))
        ([~x ~y ~'& ~rest]
         (list* '~qname (list '~fname ~x ~y) ~rest))))))


(defmacro ^:private variadic-predicate-proxy
  "Turns variadic predicates into multiple pair-wise comparisons.
  https://github.com/ztellman/primitive-math/blob/master/src/primitive_math.clj#L27"
  ([name]
   `(variadic-predicate-proxy ~name ~name))
  ([name fn]
   `(variadic-predicate-proxy ~name ~fn (constantly true)))
  ([name fn single-arg-form]
   (let [x (symbol "x")
         y (symbol "y")
         rest (symbol "rest")
         fname (symbol (str "fastmath.java.PrimitiveMath/" fn))
         qname (symbol (str "fastmath.core/" name))
         doc (or (:doc (meta name)) (str "A primitive math version of `" name "`"))]
     `(defmacro ~name
        ~doc
        ([~x]
         ~((eval single-arg-form) x))
        ([~x ~y]
         (list '~fname ~x ~y))
        ([~x ~y ~'& ~rest]
         (list 'fastmath.java.PrimitiveMath/and (list '~fname ~x ~y) (list* '~qname ~y ~rest)))))))

;; ## Basic operations

(variadic-proxy ^{:metadoc/categories -prim-set-} + add)
(variadic-proxy ^{:metadoc/categories -prim-set-} - subtract (fn [x] `(list 'fastmath.java.PrimitiveMath/negate ~x)))
(variadic-proxy ^{:metadoc/categories -prim-set-} * multiply)
(variadic-proxy ^{:metadoc/categories -prim-set-} / divide (fn [x] `(list 'fastmath.java.PrimitiveMath/reciprocal ~x)))
(primitivemath-proxy :one ^{:metadoc/categories -prim-set-} inc)
(primitivemath-proxy :one ^{:metadoc/categories -prim-set-} dec)
(primitivemath-proxy :two ^{:metadoc/categories -mod-set-} rem remainder)
(primitivemath-proxy :two ^{:metadoc/categories -mod-set-} quot quotient)
(primitivemath-proxy :two ^{:metadoc/categories -mod-set-} mod modulus)
(variadic-proxy ^{:metadoc/categories -bitwise-set-} bit-and bitAnd)
(variadic-proxy ^{:metadoc/categories -bitwise-set-} bit-or bitOr)
(variadic-proxy ^{:metadoc/categories -bitwise-set-} bit-xor bitXor)
(primitivemath-proxy :one ^{:metadoc/categories -bitwise-set-} bit-not bitNot)
(variadic-proxy ^{:metadoc/categories #{:bool}
                  :deprecated true} bool-and and)
(variadic-proxy ^{:metadoc/categories #{:bool}
                  :deprecated true} bool-or or)
(variadic-proxy ^{:metadoc/categories #{:bool}
                  :deprecated true} bool-xor xor)
(primitivemath-proxy :one ^{:metadoc/categories #{:bool}} bool-not not)
(variadic-proxy ^{:metadoc/categories -prim-set-} min)
(variadic-proxy ^{:metadoc/categories -prim-set-} max)
(primitivemath-proxy :one ^{:metadoc/categories -compare-set-} zero? isZero)
(primitivemath-proxy :one ^{:metadoc/categories -compare-set-} one? isOne)
(primitivemath-proxy :one ^{:metadoc/categories -compare-set-} neg? isNeg)
(primitivemath-proxy :one ^{:metadoc/categories -compare-set-} pos? isPos)
(primitivemath-proxy :one ^{:metadoc/categories -compare-set-} even? isEven)
(primitivemath-proxy :one ^{:metadoc/categories -compare-set-} odd? isOdd)
(primitivemath-proxy :two ^{:metadoc/categories -bitwise-set-} << shiftLeft)
(primitivemath-proxy :two ^{:metadoc/categories -bitwise-set-} >> shiftRight)
(primitivemath-proxy :two ^{:metadoc/categories -bitwise-set-} >>> unsignedShiftRight)
(primitivemath-proxy :two ^{:metadoc/categories -bitwise-set-} bit-shift-left shiftLeft)
(primitivemath-proxy :two ^{:metadoc/categories -bitwise-set-} bit-shift-right shiftRight)
(primitivemath-proxy :two ^{:metadoc/categories -bitwise-set-} unsigned-bit-shift-right unsignedShiftRight)

(variadic-predicate-proxy ^{:metadoc/categories -compare-set-} < lt)
(variadic-predicate-proxy ^{:metadoc/categories -compare-set-} > gt)
(variadic-predicate-proxy ^{:metadoc/categories -compare-set-} <= lte)
(variadic-predicate-proxy ^{:metadoc/categories -compare-set-} >= gte)
(variadic-predicate-proxy ^{:doc "Equality. See also [[eq]] for function version." :metadoc/categories -compare-set-} == eq)
(variadic-predicate-proxy ^{:metadoc/categories -compare-set-} not== neq)

(defn fast+ {:inline (fn [x y] `(+ ~x ~y)) :inline-arities #{2}
             :doc "Primitive `+` for two arguments as function (not macro)."
             :metadoc/categories -prim-set-} ^double [^double a ^double b] (+ a b))
(defn fast- {:inline (fn [x y] `(- ~x ~y)) :inline-arities #{2}
             :doc "Primitive `-` for two arguments as function (not macro)."
             :metadoc/categories -prim-set-} ^double [^double a ^double b] (- a b))
(defn fast* {:inline (fn [x y] `(* ~x ~y)) :inline-arities #{2}
             :doc "Primitive `*` for two arguments as function (not macro)."
             :metadoc/categories -prim-set-} ^double [^double a ^double b] (* a b))
(defn fast-max {:inline (fn [x y] `(+ ~x ~y)) :inline-arities #{2}
                :doc "Primitive `max` for two arguments as function (not macro)."
                :metadoc/categories -prim-set-} ^double [^double a ^double b] (max a b))
(defn fast-min {:inline (fn [x y] `(+ ~x ~y)) :inline-arities #{2}
                :doc "Primitive `min` for two arguments as function (not macro)."
                :metadoc/categories -prim-set-} ^double [^double a ^double b] (min a b))
(defn fast-identity {:inline (fn [x] `~x) :inline-arities #{1}
                     :doc "Identity on double."
                     :metadoc/categories -prim-set-} ^double [^double a] a)

;; Primitive math eq
(defn eq 
  "Primitive math equality function for doubles. See [[==]]."
  {:metadoc/categories -compare-set-}
  ([_] true)
  ([^double a ^double b]
   (== a b))
  ([^double a ^double b ^double c]
   (and (== a b) (== b c)))
  ([^double a ^double b ^double c ^double d]
   (and (== a b) (== b c) (== c d))))

;; macros for polynomials

(defmacro muladd
  "`[x y z]` -> `(+ z (* x y))` or `Math/fma` for java 9+"
  [x y z]
  (if (= :new jvm-version-type)
    `(Math/fma ~x ~y ~z)
    `(+ ~z (* ~x ~y))))

(defmacro fma
  "`[x y z]` -> `(+ z (* x y))` or `Math/fma` for java 9+"
  [x y z]
  `(muladd ~x ~y ~z))

(defmacro negmuladd
  "`[x y z]` -> `(+ z (* -1.0 x y)`"
  [x y z]
  `(muladd (- ~x) ~y ~z))

(defmacro mevalpoly
  "Evaluate polynomial macro version"
  [x & coeffs]
  (let [cnt (count coeffs)]
    (condp clojure.core/= cnt
      0 `0.0
      1 `~(first coeffs)
      2 (let [[z y] coeffs]
          `(muladd ~x ~y ~z))
      `(muladd ~x (mevalpoly ~x ~@(rest coeffs)) ~(first coeffs)))))

(defn evalpoly
  "Evaluate polynomial"
  [x & coeffs]
  (if-not (seq coeffs)
    0.0
    (let [rc (reverse coeffs)]
      (loop [rcoeffs (rest rc)
             ^double ex (first rc)]
        (if-not (seq rcoeffs)
          ex
          (recur (rest rcoeffs)
                 (muladd ^double x ex ^double (first rcoeffs))))))))

(defn makepoly
  "Create polynomial function for given coefficients"
  [coeffs]
  (cond
    (not (seq coeffs)) (constantly 0.0)
    (= 1 (count coeffs)) (constantly (first coeffs))
    :else (let [rc (reverse coeffs)]
            (fn [^double x]
              (loop [rcoeffs (rest rc)
                     ^double ex (first rc)]
                (if-not (seq rcoeffs)
                  ex
                  (recur (rest rcoeffs)
                         (muladd x ex ^double (first rcoeffs)))))))))
;; some stuff from pbrt
(defn difference-of-products
  ^double [^double a ^double b ^double c ^double d]
  (let [cd (* c d)]
    (+ (fma a b (- cd)) (fma (- c) d cd))))

(defn sum-of-products
  ^double [^double a ^double b ^double c ^double d]
  (let [cd (* c d)]
    (+ (fma a b cd) (fma c d (- cd)))))


;; Processing math constants
(def ^:const ^double ^{:doc "Value of \\\\(\\pi\\\\)"} PI Math/PI)
(def ^:const ^double ^{:doc "Value of \\\\(\\frac{\\pi}{2}\\\\)"} HALF_PI (* PI 0.5))
(def ^:const ^double ^{:doc "Value of \\\\(\\frac{\\pi}{3}\\\\)"} THIRD_PI (/ PI 3.0))
(def ^:const ^double ^{:doc "Value of \\\\(\\frac{\\pi}{4}\\\\)"} QUARTER_PI (* PI 0.25))
(def ^:const ^double ^{:doc "Value of \\\\(2 {\\pi}\\\\)"} TWO_PI (+ PI PI))
(def ^:const ^double ^{:doc "Alias for [[TWO_PI]]"} TAU TWO_PI)
(def ^:const ^double ^{:doc "Value of \\\\(e\\\\)"} E Math/E)

(def ^:const ^double ^{:doc "Value of \\\\(\\-pi\\\\)"} -PI (- Math/PI))
(def ^:const ^double ^{:doc "Value of \\\\(-\\frac{\\pi}{2}\\\\)"} -HALF_PI (* PI -0.5))
(def ^:const ^double ^{:doc "Value of \\\\(-\\frac{\\pi}{3}\\\\)"} -THIRD_PI (/ -PI -3.0))
(def ^:const ^double ^{:doc "Value of \\\\(-\\frac{\\pi}{4}\\\\)"} -QUARTER_PI (* PI -0.25))
(def ^:const ^double ^{:doc "Value of \\\\(-2 {\\pi}\\\\)"} -TWO_PI (- TWO_PI))
(def ^:const ^double ^{:doc "Alias for [[TWO_PI-]]"} -TAU -TWO_PI)
(def ^:const ^double ^{:doc "Value of \\\\(e\\\\)"} -E (- Math/E))

(def ^:const ^double ^{:doc "Value of \\\\(\\frac{1}{\\pi}\\\\)"} INV_PI (/ PI))
(def ^:const ^double ^{:doc "Value of \\\\(\\frac{2}{\\pi}\\\\)"} TWO_INV_PI (/ 2.0 PI))
(def ^:const ^double ^{:doc "Value of \\\\(\\frac{4}{\\pi}\\\\)"} FOUR_INV_PI (/ 4.0 PI))
(def ^:const ^double ^{:doc "Value of \\\\(\\frac{1}{2 \\pi}\\\\)"} INV_TWO_PI (/ TWO_PI))

(def ^:const ^double ^{:doc "Very small number \\\\(\\varepsilon\\\\)"} EPSILON 1.0e-10)

(def ^:const ^double ^{:doc "Euler-Mascheroni constant"} GAMMA Gamma/GAMMA)
(def ^:const ^double ^{:doc "Lanchos approximation `g` constant"} LANCZOS_G Gamma/LANCZOS_G)

(def ^:const ^double ^{:doc "Smallest machine number. Value is calculated during evaluation and may differ on different processors."}
  MACHINE-EPSILON (* 0.5 (double (loop [d (double 1.0)]
                                   (if (not== 1.0 (+ 1.0 (* d 0.5)))
                                     (recur (* d 0.5))
                                     d)))))

(def ^:const ^double ^{:doc "Value of \\\\(\\frac{1}{3}\\\\)"} THIRD (/ 3.0))
(def ^:const ^double ^{:doc "Value of \\\\(\\frac{2}{3}\\\\)"} TWO_THIRD (/ 2.0 3.0))
(def ^:const ^double ^{:doc "Value of \\\\(\\frac{1}{6}\\\\)"} SIXTH (/ 6.0))

;; Trigonometry
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} sin)
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} cos)
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} tan)
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} asin)
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} acos)
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} atan)
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} sinh)
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} cosh)
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} tanh)
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} asinh)
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} acosh)
(fastmath-proxy :one ^{:metadoc/categories -trig-set-} atanh)

(fastmath-proxy :one ^{:doc "Fast and less accurate [[sin]]." :metadoc/categories -trig-set-} qsin sinQuick)
(fastmath-proxy :one ^{:doc "Fast and less accurate [[cos]]." :metadoc/categories -trig-set-} qcos cosQuick)

;; Additional trigonometry functions
(defn ^{:doc "Cotangent" :metadoc/categories -trig-set-} cot ^double [^double v] (FastMath/tan (- HALF_PI v)))
(defn ^{:doc "Secant" :metadoc/categories -trig-set-} sec ^double [^double v] (/ (FastMath/cos v)))
(defn ^{:doc "Cosecant" :metadoc/categories -trig-set-} csc ^double [^double v] (/ (FastMath/sin v)))

;; Additional cyclometric functions
(defn ^{:doc "Arccotangent" :metadoc/categories -trig-set-} acot ^double [^double v] (- HALF_PI (FastMath/atan v)))
(defn ^{:doc "Arcsecant" :metadoc/categories -trig-set-} asec ^double [^double v] (FastMath/acos (/ 1.0 v)))
(defn ^{:doc "Arccosecant" :metadoc/categories -trig-set-} acsc ^double [^double v] (FastMath/asin (/ 1.0 v)))
(fastmath-proxy :two ^{:metadoc/categories -trig-set-} atan2)

;; Additional hyperbolic functions
(defn ^{:doc "Hyperbolic cotangent" :metadoc/categories -trig-set-} coth ^double [^double v] (/ (FastMath/tanh v)))
(defn ^{:doc "Hyperbolic secant" :metadoc/categories -trig-set-} sech ^double [^double v] (/ (FastMath/cosh v)))
(defn ^{:doc "Hyperbolic cosecant" :metadoc/categories -trig-set-} csch ^double [^double v] (/ (FastMath/sinh v)))

;; Additional inverse hyperbolic functions
(defn ^{:doc "Area hyperbolic cotangent" :metadoc/categories -trig-set-} acoth ^double [^double v] (FastMath/atanh (/ v)))
(defn ^{:doc "Area hyperbolic secant" :metadoc/categories -trig-set-} asech ^double [^double v] (FastMath/acosh (/ v)))
(defn ^{:doc "Area hyperbolic cosecant" :metadoc/categories -trig-set-} acsch ^double [^double v] (FastMath/asinh (/ v)))

;; haversine

(defn ^{:doc "Haversine formula" :metadoc/categories -trig-set-} haversine
  (^double [^double v] (* 0.5 (- 1.0 (FastMath/cos v))))
  (^double [[^double lat1 ^double lon1] [^double lat2 ^double lon2]]
   (haversine lat1 lon1 lat2 lon2))
  (^double [^double lat1 ^double lon1 ^double lat2 ^double lon2]
   (+ (haversine (- lat2 lat1))
      (* (FastMath/cos lat1)
         (FastMath/cos lat2)
         (haversine (- lon2 lon1))))))

(defn ^{:doc "Haversine distance `d` for `r=1`" :metadoc/categories -trig-set-} haversine-dist
  (^double [[^double lat1 ^double lon1] [^double lat2 ^double lon2]]
   (haversine-dist lat1 lon1 lat2 lon2))
  (^double [^double lat1 ^double lon1 ^double lat2 ^double lon2]
   (* 2 (FastMath/asin (FastMath/sqrt (haversine lat1 lon1 lat2 lon2))))))

;; exp and log
(fastmath-proxy :one ^{:metadoc/categories -pow-set-} exp)
(fastmath-proxy :one ^{:metadoc/categories -pow-set-} log)
(fastmath-proxy :one ^{:doc "\\\\(\\ln_{10}{x}\\\\)" :metadoc/categories -pow-set-} log10)
;; Alias for natural logarithm
(fastmath-proxy :one ^{:metadoc/categories -pow-set-} ln log)

(fastmath-proxy :one ^{:metadoc/categories -pow-set-} log1p)
(fastmath-proxy :one ^{:metadoc/categories -pow-set-} expm1)

(defn ^{:doc "log(1+exp(x))"
        :metadoc/categories -pow-set-}
  log1pexp ^double [^double x] (FastMath/log1p (FastMath/exp x)))

(defn ^{:doc "log(1-exp(x))"
        :metadoc/categories -pow-set-}
  log1mexp ^double [^double x] (FastMath/log1p (- (FastMath/exp x))))

(defn ^{:doc "log(1+x^2))"
        :metadoc/categories -pow-set-}
  log1psq ^double [^double x] (FastMath/log1p (* x x)))

(defn ^{:doc "log(exp(x)-1))"
        :metadoc/categories -pow-set-}
  logexpm1 ^double [^double x] (FastMath/log (FastMath/expm1 x)))

;; from julia
(defn- log1pmx-ker
  ^double [^double x]
  (let [r (/ x (+ 2.0 x))
        t (* r r)
        w (mevalpoly t 6.66666666666666667e-1 4.00000000000000000e-1 2.85714285714285714e-1 2.22222222222222222e-1
                     1.81818181818181818e-1 1.53846153846153846e-1 1.33333333333333333e-1 1.17647058823529412e-1)
        hxsq (* 0.5 x x)]
    (- (* r (+ hxsq (* w t))) hxsq)))

(defn ^{:doc "log(1+x)-x"
        :metadoc/categories -pow-set-}
  log1pmx ^double [^double x]
  (cond
    (not (< -0.7 x 0.9)) (- (FastMath/log1p x) x)
    (> x 0.315) (let [u (/ (- x 0.5) 1.5)]
                  (- (log1pmx-ker u) 9.45348918918356180e-2 (* 0.5 u)))
    (> x -0.227) (log1pmx-ker x)
    (> x -0.4) (let [u (/ (+ x 0.25) 0.75)]
                 (+ (log1pmx-ker u) -3.76820724517809274e-2 (* 0.25 u)))
    (> x -0.6) (let [u (* (+ x 0.5) 2.0)]
                 (+ (log1pmx-ker u) -1.93147180559945309e-1 (* 0.5 u)))
    :else (let [u (/ (+ x 0.625) 0.375)]
            (+ (log1pmx-ker u) -3.55829253011726237e-1 (* 0.625 u)))))

(defn ^{:doc "log(x)-x+1"
        :metadoc/categories -pow-set-}
  logmxp1 ^double [^double x]
  (cond
    (<= x 0.3) (- (inc (FastMath/log x)) x)
    (<= x 0.4) (let [u (/ (- x 0.375) 0.375)]
                 (+ (log1pmx-ker u) -3.55829253011726237e-1 (* 0.625 u)))
    (<= x 0.6) (let [u (* (- x 0.5) 2.0)]
                 (+ (log1pmx-ker u) -1.93147180559945309e-1 (* 0.5 u)))
    :else (log1pmx (dec x))))


;; Roots (square and cubic)
(fastmath-proxy :one ^{:doc "\\\\(\\sqrt{x}\\\\)" :metadoc/categories -pow-set-} sqrt)
(fastmath-proxy :one ^{:doc "\\\\(\\sqrt[3]{x}\\\\)" :metadoc/categories -pow-set-} cbrt)

;; Quick version of exponential \\(e^x\\)
(fastmath-proxy :one ^{:doc "Quick and less accurate version of [[exp]]." :metadoc/categories -pow-set-} qexp expQuick)

;; Radians to degrees (and opposite) conversions
(def ^:const ^double ^{:doc "\\\\(\\frac{180}{\\pi}\\\\)"} rad-in-deg (/ 180.0 PI))
(def ^:const ^double ^{:doc "\\\\(\\frac{\\pi}{180}\\\\)"} deg-in-rad (/ PI 180.0))
(defn ^{:doc "Convert degrees into radians."
        :metadoc/categories -conv-set-}
  radians ^double [^double deg] (* deg-in-rad deg))
(defn ^{:doc "Convert radians into degrees."
        :metadoc/categories -conv-set-}
  degrees ^double [^double rad] (* rad-in-deg rad))

;; Erf
(erf-proxy :onetwo ^{:doc "Error function. For two arguments return difference between `(erf x)` and `(erf y)`." :metadoc/categories -err-set-} erf)
(erf-proxy :one ^{:doc "Complementary error function." :metadoc/categories -err-set-} erfc)
(erf-proxy :one ^{:doc "Inverse [[erf]]." :metadoc/categories -err-set-} inv-erf erfInv)
(erf-proxy :one ^{:doc "Inverse [[erfc]]." :metadoc/categories -err-set-} inv-erfc erfcInv)

;; Gamma

(gamma-proxy :one ^{:doc "Gamma function \\\\(\\Gamma(x)\\\\)" :metadoc/categories -special-set-} gamma gamma)
(gamma-proxy :one ^{:doc "Gamma function \\\\(\\ln\\Gamma(x)\\\\)" :metadoc/categories -special-set-} log-gamma logGamma)
(gamma-proxy :one ^{:doc "Gamma function \\\\(\\ln\\Gamma(1+x)\\\\)" :metadoc/categories -special-set-} log-gamma-1p logGamma1p)
(gamma-proxy :one ^{:doc "Logarithmic derivative of \\\\(\\Gamma\\\\)." :metadoc/categories -special-set-} digamma)
(gamma-proxy :one ^{:doc "Derivative of [[digamma]]." :metadoc/categories -special-set-} trigamma)
(gamma-proxy :one ^{:doc "\\\\(\\frac{1}{\\Gamma(1+x)}\\\\)." :metadoc/categories -special-set-} inv-gamma-1pm1 invGamma1pm1)
(gamma-proxy :two ^{:doc "Regularized `gamma` P" :metadoc/categories -special-set-} regularized-gamma-p regularizedGammaP)
(gamma-proxy :two ^{:doc "Regularized `gamma` Q" :metadoc/categories -special-set-} regularized-gamma-q regularizedGammaQ)

(defn minkowski
  "Minkowski's question mark function ?(x)"
  {:metadoc/categories -special-set-}
  ^double [^double x]
  (loop [it (long 0) p 0.0 q 1.0 r 1.0 s 1.0 d 1.0 y 0.0]
    (if (< it 20)
      (let [d (* d 0.5)
            m (+ p r)
            n (+ q s)
            p? (< x (/ m n))]
        (recur (inc it)
               (if p? p m)
               (if p? q n)
               (if p? m r)
               (if p? n s)
               d
               (if p? y (+ y d))))
      (+ y d))))

;; Beta

(beta-proxy :two ^{:doc "Logarithm of Beta function." :metadoc/categories -special-set-} log-beta logBeta)
(beta-proxy :three ^{:doc "Regularized `Beta`." :metadoc/categories -special-set-} regularized-beta regularizedBeta)

;; BesselJ
(besselj-proxy :two ^{:doc "Bessel J function value for given order and argument." :metadoc/categories -special-set-} bessel-j value)

(def ^:private ^:const ^double jinc-c4 (/ (* PI PI PI PI) 192.0))
(def ^:private ^:const ^double jinc-c2 (/ (* PI PI) -8.0))

(defn jinc
  "Besselj1 devided by `x`"
  {:metadoc/categories -special-set-}
  ^double [^double x]
  (if (< (FastMath/abs x) 0.002)
    (let [x2 (* x x)]
      (mevalpoly x2 1.0 jinc-c2 jinc-c4))
    (let [pix (* PI x)]
      (* 2.0 (/ (bessel-j 1 pix) pix)))))

;; I0

(defn I0
  "Modified Bessel function of the first kind, order 0."
  {:metadoc/categories -special-set-}
  ^double [^double x]
  (let [x2 (* x x)
        ;; i=1
        val (inc (/ x2 4))
        x2i (* x2 x2)
        ;; i=2
        val (+ val (/ x2i 64))
        x2i (* x2i x2)
        ;; i=3
        val (+ val (/ x2i 2304))
        x2i (* x2i x2)
        ;; i=4
        val (+ val (/ x2i 147456))
        x2i (* x2i x2)
        ;; i=5
        val (+ val (/ x2i 14745600))
        x2i (* x2i x2)
        ;; i=6
        val (+ val (/ x2i 2123366400))
        x2i (* x2i x2)
        ;; i=7
        val (+ val (/ x2i 416179814400))
        x2i (* x2i x2)
        ;; i=8
        val (+ val (/ x2i 106542032486400))
        x2i (* x2i x2)]
    (+ val (/ x2i 34519618525593600))))

(defn logI0
  "Log of [[I0]]."
  {:metadoc/categories -special-set-}
  ^double [^double x]
  (m/log (I0 x)))

;; Sinc
(defn sinc
  "Sinc function."
  {:metadoc/categories -trig-set-}
  ^double [^double v]
  (let [x (* PI (FastMath/abs v))]
    (if (< x 1.0e-5) 1.0
        (/ (FastMath/sin x) x))))

;;
(defn sigmoid
  "Sigmoid function"
  {:metadoc/categories -pow-set-}
  ^double [^double x]
  (/ (inc (FastMath/exp (- x)))))

(def ^{:metadoc/categories -pow-set-
       :doc "Alias for [[sigmoid]]"}
  logistic sigmoid)

(defn logit
  "Logit function"
  {:metadoc/categories -pow-set-}
  ^double [^double x]
  (FastMath/log (/ x (- 1.0 x))))

(def ^:const ^double ^{:doc "\\\\(\\ln{2}\\\\)"} LN2 (log 2.0))
(def ^:const ^double ^{:doc "\\\\(\\frac{1}{\\ln{2}}\\\\)"} INV_LN2 (/ LN2))
(def ^:const ^double ^{:doc "\\\\(\\frac{\\ln{2}}{2}\\\\)"} LN2_2 (* 0.5 LN2))
(def ^:const ^double ^{:doc "\\\\(\\ln{10}\\\\)"} LN10 (log 10.0))
(def ^:const ^double ^{:doc "\\\\(\\frac{1}{\\ln{0.5}}\\\\)"} INV_LOG_HALF (/ (log 0.5)))
(def ^:const ^double ^{:doc "\\\\(\\ln{0.5}\\\\)"} LOG_HALF (log 0.5))
(def ^:const ^double ^{:doc "\\\\(\\ln{\\pi}\\\\)"} LOG_PI (log PI))
(def ^:const ^double ^{:doc "\\\\(\\ln{2 \\pi}\\\\)"} LOG_TWO_PI (log TWO_PI))

(defn xlogx
  "x * log(x)"
  {:metadoc/categories -pow-set-}
  ^double [^double x]
  (if (zero? x) 0.0 (* x (log x))))

(defn xlogy
  "x * log(y)"
  {:metadoc/categories -pow-set-}
  ^double [^double x ^double y]
  (if (and (zero? x)
           (not (Double/isNaN y))) 0.0 (* x (log y))))


(defn log2
  "Logarithm with base 2.

  \\\\(\\ln_2{x}\\\\)"
  {:metadoc/categories -pow-set-}
  ^double [^double x]
  (* (FastMath/log x) INV_LN2))

;; \\(\log_b x\\)
(defn logb
  "Logarithm with base `b`.

  \\\\(\\ln_b{x}\\\\)"
  {:metadoc/categories -pow-set-}
  ^double [^double b ^double x]
  (/ (FastMath/log x) (FastMath/log b)))

;; Quick logarithm
(fastmath-proxy :one ^{:doc "Fast and less accurate version of [[log]]." :metadoc/categories -pow-set-} qlog logQuick)

;; \\(\log_2 e\\)
(def ^:const ^double ^{:doc "\\\\(\\log_{2}{e}\\\\)"} LOG2E (log2 E))

;; \\(\log_{10} e\\)
(def ^:const ^double ^{:doc "\\\\(\\log_{10}{e}\\\\)"} LOG10E (log10 E))

;; Powers (normal, quick)
(fastmath-proxy :two ^{:metadoc/categories -pow-set-} pow)
(fastmath-proxy :two ^{:doc "Fast and less accurate version of [[pow]]." :metadoc/categories -pow-set-} qpow powQuick)

;; Fast version of power, second parameter should be integer
(fastmath-proxy :two ^{:doc "Fast version of pow where exponent is integer." :metadoc/categories -pow-set-} fpow powFast)

(def ^:private factorial20-table [1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600
                                6227020800 87178291200 1307674368000 20922789888000
                                355687428096000 6402373705728000 121645100408832000
                                2432902008176640000])

(defn factorial20
  "Factorial table up to 20!"
  {:metadoc/categories -pow-set-}
  ^long [^long n]
  (factorial20-table n))

;; Square and cubic
(defn sq "Same as [[pow2]]. \\\\(x^2\\\\)" {:metadoc/categories -pow-set-} ^double [^double x] (* x x))
(defn pow2 "Same as [[sq]]. \\\\(x^2\\\\)" {:metadoc/categories -pow-set-} ^double [^double x] (* x x))
(defn pow3 "\\\\(x^3\\\\)" {:metadoc/categories -pow-set-} ^double [^double x] (* x (* x x)))
(defn cb "\\\\(x^3\\\\)" {:metadoc/categories -pow-set-} ^double [^double x] (* x (* x x)))

(defn safe-sqrt
  "Safe sqrt, for value <= 0 result is 0.

  \\\\(
  \\left\\\\{
  \\begin{array}{lr}
  0 & : x \\leq 0\\\\\\\\
  \\sqrt{x} & : x > 0
  \\end{array}
  \\\\right.
  \\\\)"
  {:metadoc/categories -pow-set-}
  ^double [^double value]
  (if (neg? value) 0.0 (sqrt value)))

;; Approximated sqrt via binary operations (error 1.0E-2)
(fastmath-proxy :one ^{:doc "Approximated [[sqrt]] using binary operations with error `1.0E-2`." :metadoc/categories -pow-set-} qsqrt sqrtQuick)
(fastmath-proxy :one ^{:doc "Inversed version of [[qsqrt]]. Quick and less accurate." :metadoc/categories -pow-set-} rqsqrt invSqrtQuick)

(defn hypot
  "Hypot.
  See also [[hypot-sqrt]]."
  {:metadoc/categories #{:dist}}
  (^double [^double x ^double y]
   (FastMath/hypot x y))
  (^double [^double x ^double y ^double z]
   (FastMath/hypot x y z)))

(defn hypot-sqrt
  "Hypot, sqrt version: \\\\(\\sqrt{x^2+y^2}\\\\) or \\\\(\\sqrt{x^2+y^2+z^2}\\\\).
  Should be faster than [[hypot]]."
  {:metadoc/categories #{:dist}}
  (^double [^double x ^double y]
   (sqrt (+ (* x x) (* y y))))
  (^double [^double x ^double y ^double z]
   (sqrt (+ (* x x) (* y y) (* z z)))))

;; distance
(defn dist
  "Euclidean distance between points `(x1,y1)` and `(x2,y2)`. See [[fastmath.vector]] namespace to see other metrics which work on vectors."
  {:metadoc/categories #{:dist}}
  ^double [^double x1 ^double y1 ^double x2 ^double y2]
  (sqrt (+ (sq (- x2 x1)) (sq (- y2 y1)))))

(defn qdist
  "Quick version of Euclidean distance between points. [[qsqrt]] is used instead of [[sqrt]]."
  {:metadoc/categories #{:dist}}
  ^double [^double x1 ^double y1 ^double x2 ^double y2]
  (qsqrt (+ (sq (- x2 x1)) (sq (- y2 y1)))))

;; Rounding functions
(defn floor
  "\\\\(\\lfloor x \\rfloor\\\\). See: [[qfloor]].

  Rounding is done to a multiply of scale value (when provided)."
  {:metadoc/categories -round-set-}
  (^double [^double x] (FastMath/floor x))
  (^double [^double x ^double scale] (* (FastMath/floor (/ x scale)) scale)))

(defn ceil
  "\\\\(\\lceil x \\rceil\\\\). See: [[qceil]].

  Rounding is done to a multiply of scale value (when provided)."
  {:metadoc/categories -round-set-}
  (^double [^double x] (FastMath/ceil x))
  (^double [^double x ^double scale] (* (FastMath/ceil (/ x scale)) scale)))

(defn ^{:doc "Round to `long`. See: [[rint]], [[qround]]."
     :metadoc/categories -round-set-} round ^long [^double x] (FastMath/round x))
(defn ^{:doc "Round to `double`. See [[round]], [[qround]]."
        :metadoc/categories -round-set-} rint ^double [^double x] (FastMath/rint x))
(defn round-even
  "Round evenly (like in round in R), IEEE / IEC rounding"
  {:metadoc/categories -round-set-}
  ^long [^double x] (FastMath/roundEven x))

(primitivemath-proxy :one ^{:doc "Fast version of [[floor]]. Returns `long`. See: [[floor]]." :metadoc/categories -round-set-} qfloor fastFloor)
(primitivemath-proxy :one ^{:doc "Fast version of [[ceil]]. Returns `long`. See: [[ceil]]." :metadoc/categories -round-set-} qceil fastCeil)
(primitivemath-proxy :one ^{:doc "Fast version of [[round]]. Returns `long`. See: [[rint]], [[round]]." :metadoc/categories -round-set-} qround fastRound)


(fastmath-proxy :two ^{:doc "From `FastMath` doc: returns dividend - divisor * n,
where n is the mathematical integer closest to dividend/divisor. Returned value in `[-|divisor|/2,|divisor|/2]`"
                       :metadoc/categories -mod-set-} remainder)

(defn abs "\\\\(|x|\\\\) - `double` version. See [[iabs]]." {:metadoc/categories -round-set-} ^double [^double x] (FastMath/abs x))
(defn iabs "\\\\(|x|\\\\) - `long` version. See [[abs]]." {:metadoc/categories -round-set-} ^long [^long x] (if (neg? x) (- x) x))

(defn trunc
  "Truncate fractional part, keep sign. Returns `double`."
  {:metadoc/categories -round-set-}
  ^double [^double v] (if (neg? v) (ceil v) (floor v)))

(defn itrunc
  "Truncate fractional part, keep sign. Returns `long`."
  {:metadoc/categories -round-set-}
  ^long [^double v] (if (neg? v) (qceil v) (qfloor v)))

;; return approximate value
(defn approx
  "Round `v` to specified (default: 2) decimal places. Be aware of `double` number accuracy."
  {:metadoc/categories -round-set-}
  (^double [^double v] (Precision/round v (int 2)))
  (^double [^double v ^long digits] (Precision/round v (int digits))))

(defn approx-eq
  "Checks equality approximately. See [[approx]]."
  {:metadoc/categories -round-set-}
  ([^double a ^double b] (== (approx a) (approx b)))
  ([^double a ^double b ^long digits] (== (approx a digits)
                                          (approx b digits))))

(defn delta-eq
  "Checks equality for given accuracy (default `10.0e-6`)."
  ([^double a ^double b] (delta-eq a b 10.0e-6))
  ([^double a ^double b ^double accuracy]
   (< (m/abs (- a b)) accuracy)))

(def ^{:metadoc/categories -round-set-
     :doc "Alias for [[approx-eq]]"} approx= approx-eq)

(def ^{:metadoc/categories -round-set-
     :doc "Alias for [[delta-eq]]"} delta= delta-eq)

(defn frac
  "Fractional part, always returns values from 0.0 to 1.0 (exclusive). See [[sfrac]] for signed version."
  {:metadoc/categories -round-set-}
  ^double [^double v] (abs (- v (unchecked-long v))))

(defn sfrac
  "Fractional part, always returns values from -1.0 to 1.0 (exclusive). See [[frac]] for unsigned version."
  {:metadoc/categories -round-set-}
  ^double [^double v] (- v (trunc v)))

;; Find power of 2 exponent for double number where  
;; \\(2^(n-1)\leq x\leq 2^n\\)  
;; where n-1 is result of `low-2-exp` and n is result of `high-2-exp`
;; `(low-2-exp TWO_PI) => 2` \\(2^2\eq 4\leq 6.28\\)  
;; `(high-2-exp TWO_PI) => 3` \\(6.28\leq 2^3\eq 8\\)
(defn low-2-exp
  "Find greatest exponent (power of 2) which is lower or equal `x`. See [[high-2-exp]]."
  {:metadoc/categories -pow-set-}
  ^long [^double x] (-> x log2 floor unchecked-long))

(defn high-2-exp
  "Find lowest exponent (power of 2) which is greater or equal `x`. See [[low-2-exp]]."
  {:metadoc/categories -pow-set-}
  ^long [^double v] (-> v log2 ceil unchecked-long))

(defn low-exp
  "Find greatest exponent for base `b` which is lower or equal `x`. See also [[high-exp]]."
  {:metadoc/categories -pow-set-}
  ^long [^double b ^double x] (->> x (logb b) floor unchecked-long))

(defn high-exp
  "Find lowest exponent for base `b` which is higher or equal`x`. See also [[low-exp]]."
  {:metadoc/categories -pow-set-}
  ^long [^double b ^double x] (->> x (logb b) ceil unchecked-long))

(defn round-up-pow2
  "Round long to the next power of 2"
  {:metadoc/categories -round-set-}
  ^long [^long v]
  (as-> (dec v) v
    (bit-or v (>> v 1))
    (bit-or v (>> v 2))
    (bit-or v (>> v 4))
    (bit-or v (>> v 8))
    (bit-or v (>> v 16))
    (bit-or v (>> v 32))
    (inc v)))

(defn next-double
  "Next double value. Optional value `delta` sets step amount."
  (^double [^double v]
   (FastMath/nextUp v))
  (^double [^double v ^long delta]
   (nth (iterate next-double v) delta)))

(defn prev-double
  "Next double value. Optional value `delta` sets step amount."
  (^double [^double v]
   (FastMath/nextDown v))
  (^double [^double v ^long delta]
   (nth (iterate prev-double v) delta)))

(defn double-high-bits
  "Returns high word from double as bits"
  ^long [^double v]
  (bit-and (>>> (Double/doubleToRawLongBits v) 32) 0xffffffff))

(defn double-low-bits
  "Returns low word from double as bits"
  ^long [^double v]
  (bit-and (Double/doubleToRawLongBits v) 0xffffffff))

;; More constants

(def ^:const ^double ^{:doc "Value of 0x1.fffffffffffffp-1d = 0.(9)"}
  double-one-minus-epsilon (Double/parseDouble "0x1.fffffffffffffp-1d"))

;; \\(\sqrt{2}\\)
(def ^:const ^double ^{:doc "\\\\(\\sqrt{2}\\\\)"} SQRT2 (sqrt 2.0))
(def ^:const ^double ^{:doc "\\\\(\\frac{\\sqrt{2}}{2}\\\\)"} SQRT2_2 (* 0.5 SQRT2))

;; \\(\sqrt{3}\\)
(def ^:const ^double ^{:doc "\\\\(\\sqrt{3}\\\\)"} SQRT3 (sqrt 3.0))
(def ^:const ^double ^{:doc "\\\\(\\frac{\\sqrt{3}}{2}\\\\)"} SQRT3_2 (* 0.5 (sqrt 3.0)))
(def ^:const ^double ^{:doc "\\\\(\\frac{\\sqrt{3}}{3}\\\\)"} SQRT3_3 (/ (sqrt 3.0) 3.0))
(def ^:const ^double ^{:doc "\\\\(\\frac{\\sqrt{3}}{4}\\\\)"} SQRT3_4 (/ (sqrt 3.0) 4.0))

;; \\(\sqrt{5}\\)
(def ^:const ^double ^{:doc "\\\\(\\sqrt{5}\\\\)"} SQRT5 (sqrt 5.0))

;; \\(\sqrt{\pi}\\)
(def ^:const ^double ^{:doc "\\\\(\\sqrt{\\pi}\\\\)"} SQRTPI (sqrt PI))
(def ^:const ^double ^{:doc "\\\\(\\sqrt{2\\pi}\\\\)"} SQRT2PI (sqrt TWO_PI))
(def ^:const ^double ^{:doc "\\\\(\\sqrt{\\frac{1}{2}\\pi}\\\\)"} SQRT_HALFPI (sqrt HALF_PI))

;; 
(def ^:const ^double ^{:doc "Golden ratio \\\\(\\varphi\\\\)"} PHI (* (inc SQRT5) 0.5))

;; math.h predefined constants names
(def ^:const ^double ^{:doc "\\\\(e\\\\)"} M_E E)
(def ^:const ^double ^{:doc "\\\\(\\log_{2}{e}\\\\)"} M_LOG2E LOG2E)
(def ^:const ^double ^{:doc "\\\\(\\log_{10}{e}\\\\)"} M_LOG10E LOG10E)
(def ^:const ^double ^{:doc "\\\\(\\ln{2}\\\\)"} M_LN2 LN2)
(def ^:const ^double ^{:doc "\\\\(\\ln{10}\\\\)"} M_LN10 LN10)
(def ^:const ^double ^{:doc "\\\\(\\pi\\\\)"} M_PI PI)
(def ^:const ^double ^{:doc "\\\\(\\frac{\\pi}{2}\\\\)"} M_PI_2 HALF_PI)
(def ^:const ^double ^{:doc "\\\\(\\frac{\\pi}{4}\\\\)"} M_PI_4 QUARTER_PI)
(def ^:const ^double ^{:doc "\\\\(\\frac{1}{\\pi}\\\\)"} M_1_PI (/ PI))
(def ^:const ^double ^{:doc "\\\\(\\frac{2}{\\pi}\\\\)"} M_2_PI (/ 2.0 PI))
(def ^:const ^double ^{:doc "\\\\(\\frac{2}{\\sqrt\\pi}\\\\)"} M_2_SQRTPI (/ 2.0 SQRTPI))
(def ^:const ^double ^{:doc "\\\\(\\frac{1}{\\sqrt\\pi}\\\\)"} INV_SQRT2PI (/ 1.0 SQRT2PI))
(def ^:const ^double ^{:doc "\\\\(\\sqrt{2}\\\\)"} M_SQRT2 SQRT2)
(def ^:const ^double ^{:doc "\\\\(\\frac{1}{\\sqrt{2}}\\\\)"} M_SQRT1_2 (/ SQRT2))
(def ^:const ^double ^{:doc "\\\\(\\frac{1}{\\sqrt{2}}\\\\)"} INV_SQRT_2 M_SQRT1_2)


(def ^:const ^double ^{:doc "\\\\(2\\pi\\\\)"} M_TWOPI TWO_PI)
(def ^:const ^double ^{:doc "\\\\(\\frac{3\\pi}{4}\\\\)"} M_3PI_4 (* PI 0.75))
(def ^:const ^double ^{:doc "\\\\(\\sqrt\\pi\\\\)"} M_SQRT_PI SQRTPI)
(def ^:const ^double ^{:doc "\\\\(\\sqrt{3}\\\\)"} M_SQRT3 SQRT3)
(def ^:const ^double ^{:doc "\\\\(\\frac{1}{\\ln{10}}\\\\)"} M_IVLN10 (/ LN10))
(def ^:const ^double ^{:doc "\\\\(\\ln{2}\\\\)"} M_LOG2_E LN2)
(def ^:const ^double ^{:doc "\\\\(\\frac{1}{\\ln{2}}\\\\)"} M_INVLN2 (/ LN2))

(defn signum
  "Return 1 if `value` is > 0, 0 if it is 0, -1 otherwise. See also [[sgn]].

  \\\\(
  \\left\\\\{
  \\begin{array}{lr}
  1.0 & : x > 0\\\\\\\\
  -1.0 & : x < 0\\\\\\\\
  0.0 & : x = 0
  \\end{array}
  \\\\right.
  \\\\)"
  {:metadoc/categories #{:sign}}
  ^double [^double value]
  (cond (pos? value) 1.0
        (neg? value) -1.0
        :else 0.0))

(defn sgn
  "Return -1 when `value` is negative, 1 otherwise. See also [[signum]].

  \\\\(
  \\left\\\\{
  \\begin{array}{lr}
  1.0 & : x \\geq 0\\\\\\\\
  -1.0 & : x < 0\\\\\\\\
  \\end{array}
  \\\\right.
  \\\\)"
  {:metadoc/categories #{:sign}}
  ^double [^double value]
  (if (neg? value) -1.0 1.0))


;; copy-sign

(fastmath-proxy :two ^{:doc "Returns a value with a magnitude of first arguemnt and sign of second."
                       :metadoc/categories #{:sign}} copy-sign copySign)

(defmacro constrain
  "Clamp `value` to the range `[mn,mx]`."
  {:metadoc/categories -conv-set-}
  [value mn mx]
  `(max (min ~value ~mx) ~mn))

(defn norm
  "Normalize `v` from the range `[start,stop]` to the range `[0,1]` or map `v` from the range `[start1,stop1]` to the range `[start2,stop2]`. See also [[make-norm]]."
  {:inline (fn
             ([v start stop] `(PrimitiveMath/norm ~v ~start ~stop))
             ([v start1 stop1 start2 stop2] `(PrimitiveMath/norm ~v ~start1 ~stop1 ~start2 ~stop2)))
   :inline-arities #{3 5}
   :metadoc/categories -conv-set-}
  (^double [^double v ^double start ^double stop] ;; norm
   (PrimitiveMath/norm v start stop))
  ([v start1 stop1 start2 stop2] ;; map
   (PrimitiveMath/norm v start1 stop1 start2 stop2)))

(defmacro mnorm
  "Macro version of [[norm]]."
  {:metadoc/categoriewes -conv-set-}
  ([v start stop]
   `(PrimitiveMath/norm ~v ~start ~stop))
  ([v start1 stop1 start2 stop2]
   `(PrimitiveMath/norm ~v ~start1 ~stop1 ~start2 ~stop2)))

(defn make-norm
  "Make [[norm]] function for given range. Resulting function accepts `double` value (with optional target `[dstart,dstop]` range) and returns `double`."
  {:metadoc/categories -conv-set-}
  ([^double start ^double stop]
   (fn ^double [^double v ^double dstart ^double dstop]
     (PrimitiveMath/norm v start stop dstart dstop)))
  ([^double start ^double stop ^double dstart ^double dstop]
   (fn ^double [^double v]
     (PrimitiveMath/norm v start stop dstart dstop))))

(defn cnorm
  "Constrained version of norm. Result of [[norm]] is applied to [[constrain]] to `[0,1]` or `[start2,stop2]` ranges."
  {:metadoc/categories -conv-set-}
  ([v start1 stop1 start2 stop2]
   (constrain ^double (PrimitiveMath/norm v start1 stop1 start2 stop2) ^double start2 ^double stop2))
  (^double [v ^double start ^double stop]
   (constrain ^double (PrimitiveMath/norm v start stop) 0.0 1.0)))

;;; Interpolation functions

;; Linear interpolation between `start` and `stop`.
(defn lerp
  "Linear interpolation between `start` and `stop` for amount `t`. See also [[mlerp]], [[cos-interpolation]], [[quad-interpolation]] or [[smooth-interpolation]]."
  {:metadoc/categories -conv-set-}
  ^double [^double start ^double stop ^double t]
  (+ start (* t (- stop start))))

(defmacro mlerp
  "[[lerp]] as macro. For inline code. See also [[lerp]], [[cos-interpolation]], [[quad-interpolation]] or [[smooth-interpolation]]."
  {:metadoc/categories -conv-set-}
  [start stop t]
  `(+ ~start (* ~t (- ~stop ~start))))

;; Cosine interpolation between `start` and `stop`
(defn cos-interpolation
  "oF interpolateCosine interpolation. See also [[lerp]]/[[mlerp]], [[quad-interpolation]] or [[smooth-interpolation]]."
  {:metadoc/categories -conv-set-}
  ^double [^double start ^double stop ^double t]
  (mlerp start stop (* 0.5 (- 1.0 (cos (* t PI))))))

(defn smooth-interpolation
  "Smoothstep based interpolation. See also [[lerp]]/[[mlerp]], [[quad-interpolation]] or [[cos-interpolation]]."
  {:metadoc/categories -conv-set-}
  ^double [^double start ^double stop ^double t]
  (mlerp start stop (* t t (- 3.0 (* 2.0 t)))))

(defn quad-interpolation
  "Quad interpolation. See also [[lerp]]/[[mlerp]], [[cos-interpolation]] or [[smooth-interpolation]]."
  {:metadoc/categories -conv-set-}
  ^double [^double start ^double stop ^double t]
  (mlerp start stop (let [t' (* 2.0 t)]
                      (if (< t' 1.0)
                        (* 0.5 (* t' t'))
                        (* -0.5 (dec (* (dec t') (- t' 3.0))))))))

(defn smoothstep
  "GL [smoothstep](https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/smoothstep.xhtml)."
  {:metadoc/categories -conv-set-}
  ^double [^double edge0 ^double edge1 ^double x]
  (let [t (cnorm x edge0 edge1)]
    (* t t (- 3.0 (* 2.0 t)))))

;;`(wrap 0 -1 1) => 0.0`  
;;`(wrap -1.1 -1 1) => 0.8999999999999999`  
;;`(wrap 1.1 -1 1) => -0.8999999999999999`
(defn wrap
  "Wrap overflowed value into the range, similar to [ofWrap](http://openframeworks.cc/documentation/math/ofMath/#!show_ofWrap)."
  {:metadoc/categories -conv-set-}
  ^double [^double start ^double stop ^double value]
  (let [p (> start stop)
        from (if p stop start)
        to (if p start stop)
        cycle (- to from)]
    (if (zero? cycle)
      to
      (->> cycle
           (/ (- value from))
           (floor)
           (* cycle)
           (- value)))))

;;

(defn nan?
  "Check if number is NaN"
  {:inline (fn [v] `(Double/isNaN ~v)) :inline-arities #{1}
   :metadoc/categories #{:bool :compare}}
  [^double v]
  (Double/isNaN v))

(defn inf?
  "Check if number is infinite"
  {:inline (fn [v] `(Double/isInfinite ~v)) :inline-arities #{1}
   :metadoc/categories #{:bool :compare}}
  [^double v]
  (Double/isInfinite v))

(defn pos-inf?
  "Check if number is positively infinite"
  {:inline (fn [v] `(== ~v ##Inf)) :inline-arities #{1}
   :metadoc/categories #{:bool :compare}}
  [^double v]
  (== v ##Inf))

(defn neg-inf?
  "Check if number is negatively infinite"
  {:inline (fn [v] `(== ~v ##-Inf)) :inline-arities #{1}
   :metadoc/categories #{:bool :compare}}
  [^double v]
  (== v ##-Inf))

(defn invalid-double?
  "Check if number is invalid"
  {:inline (fn [v] `(bool-not (Double/isFinite ~v))) :inline-arities #{1}
   :metadoc/categories #{:bool :compare}}
  [^double v]
  (bool-not (Double/isFinite v)))

(defn valid-double?
  "Check if number is invalid"
  {:inline (fn [v] `(Double/isFinite ~v)) :inline-arities #{1}
   :metadoc/categories #{:bool :compare}}
  [^double v]
  (Double/isFinite v))

(defn between?
  "Check if given number is within the range [x,y]."
  {:metadoc/categories #{:bool :compare}}
  ([[^double x ^double y] ^double v] (<= x v y))
  ([^double x ^double y ^double v] (<= x v y)))

(defn between-?
  "Check if given number is within the range (x,y]."
  {:metadoc/categories #{:bool :compare}}
  ([[^double x ^double y] ^double v] (and (< x v) (<= v y)))
  ([^double x ^double y ^double v] (and (< x v) (<= v y))))

;; intervals

(defn slice-range 
  "Slice range to get `cnt` number of points evenly distanced."
  ([^double start ^double end ^long cnt] (if (= cnt 1)
                                           (list (+ start (* 0.5 (- end start))))
                                           (map (make-norm 0.0 (dec cnt) start end) (range cnt))))
  ([^long cnt] (slice-range 0.0 1.0 cnt)))

(defn cut
  "Cut range or sequence into intervals"
  ([data ^long breaks]
   (let [d (sort (remove invalid-double? data))]
     (cut (first d) (last d) breaks)))
  ([^double x1 ^double x2 ^long breaks]
   (let [[[^double start end] & r] (->> (slice-range x1 x2 (inc breaks))
                                        (partition 2 1))]
     (conj r (list (- start EPSILON) end)))))

(defn co-intervals
  "Divide sequence to overlaping intervals containing similar number of values. Same as R's `co.intervals()`"
  ([data] (co-intervals data 6))
  ([data ^long number] (co-intervals data number 0.5))
  ([data ^long number ^double overlap]
   (let [o- (- 1.0 overlap)
         x (vec (sort (remove invalid-double? data)))
         n (count x)
         r (/ n (+ (* number o-) overlap))
         ii-mult (*  o- r)
         ii (mapv #(* ^long % ii-mult) (range number))
         x1 (map #(x (round-even %)) ii)
         xr (map #(x (dec (round-even (+ r ^double %)))) ii)
         diffs (filter #(pos? ^double %) (mapv (fn [[^double x ^double y]] (- y x)) (partition 2 1 x)))
         eps (* 0.5 (double (if (seq diffs) (reduce fast-min diffs) 0.0)))]
     (for [[[^double px ^double cx] [^double py ^double cy]] (map vector
                                                                  (partition 2 1 (conj x1 (dec ^double (first x1))))
                                                                  (partition 2 1 (conj xr (dec ^double (first xr)))))
           :when (or (pos? (- cx px))
                     (pos? (- cy py)))]
       [(- cx eps) (+ cy eps)]))))

(defn group-by-intervals
  "Group sequence of values into given intervals.

  If `intervals` are missing, use [[co-intervals]] to find some."
  ([coll] (group-by-intervals (co-intervals coll) coll))
  ([intervals coll]
   (into {} (map (fn [[^double x1 ^double x2 :as i]]
                   [i (filter #(between-? x1 x2 %) coll)]) intervals))))

;; gcd

(defn- gcd-
  "Input is unsigned!"
  ^long [^long a ^long b]
  (cond
    (== a b) a
    (zero? a) b
    (zero? b) a
    (and (even? a) (even? b)) (<< (gcd- (>> a 1) (>> b 1)) 1)
    (and (even? a) (odd? b)) (recur (>> a 1) b)
    (and (odd? a) (even? b)) (recur a (>> b 1))
    (and (odd? a) (odd? a)) (if (> a b)
                              (recur (>> (- a b) 1) b)
                              (recur (>> (- b a) 1) a))))

(defn gcd
  "Fast binary greatest common divisor (Stein's algorithm)"
  {:metadoc/categories -mod-set-}
  ^long [^long a ^long b]
  (gcd- (iabs a) (iabs b)))

(defn lcm
  "Fast binary least common multiplier."
  {:metadoc/categories -mod-set-}
  ^long [^long a ^long b]
  (/ (* a b) (gcd- (iabs a) (iabs b))))

;;

(defn sample
  "Sample function `f` and return sequence of values.

  `range-min` defaults to 0.0, `range-max` to 1.0.

  Range is inclusive.

  When optional `domain?` is set to true (default: false) function returns pairs `[x,(f x)]`."
  {:metadoc/categories #{:sample}}
  ([f number-of-values]
   (sample f 0.0 1.0 number-of-values false))
  ([f ^long number-of-values domain?]
   (sample f 0.0 1.0 number-of-values domain?))
  ([f range-min range-max number-of-values]
   (sample f range-min range-max number-of-values false))
  ([f range-min range-max number-of-values domain?]
   (let [n- (dec ^long number-of-values)
         f (if domain? #(vector % (f %)) f)]
     (->> (range number-of-values)
          (map #(norm % 0.0 n- range-min range-max))
          (map f)))))

;; rank/order

(defn rank
  "Sample ranks. See [R docs](https://www.rdocumentation.org/packages/base/versions/3.6.1/topics/rank).

  Rank uses 0 based indexing.
  
  Possible tie strategies: `:average`, `:first`, `:last`, `:random`, `:min`, `:max`, `:dense`.

  `:dense` is the same as in `data.table::frank` from R"
  {:metadoc/categories #{:rank}}
  ([vs] (rank vs :average))
  ([vs ties] (rank vs ties false))
  ([vs ties desc?]
   (let [cmp (if desc? #(compare %2 %1) compare)
         indexed-sorted-map (group-by second (map-indexed vector (sort cmp vs)))]
     (if (#{:first :last :random} ties)
       (let [tie-sort (case ties
                        :first (partial sort-by first clojure.core/<)
                        :last (partial sort-by first clojure.core/>)
                        :random shuffle)
             sorted2-map (into {} (map (fn [[k v]] [k (tie-sort v)]) indexed-sorted-map))]
         (first (reduce (fn [[res curr] v]
                          (let [lst (curr v)]
                            [(conj res (ffirst lst))
                             (assoc curr v (rest lst))])) [[] sorted2-map] vs)))
       (let [tie-fn (case ties
                      :min ffirst
                      :dense ffirst
                      :max (comp first last)
                      (fn ^double [v] (/ ^double (reduce #(+ ^double %1 ^double %2) (map first v)) (count v))))
             m (map (fn [[k v]] [k (tie-fn v)]) indexed-sorted-map)
             m (if (= ties :dense)
                 (map-indexed (fn [id [k _]]
                                [k id]) (sort-by second m))
                 m)]
         (map (into {} m) vs))))))

(defn order
  "Ordering permutation. See [R docs](https://www.rdocumentation.org/packages/base/versions/3.6.1/topics/order)

  Order uses 0 based indexing."
  {:metadoc/categories #{:rank}}
  ([vs] (order vs false))
  ([vs decreasing?]
   (->> (map-indexed vector vs)
        (sort-by second (if decreasing?
                          clojure.core/>
                          clojure.core/<))
        (map first))))

;; sinint / cosint

(defn Si
  "Integral of sin(t)/t from 0 to x"
  ^double [^double x]
  (if (nan? x)
    ##NaN
    (let [t (* x x)]
      (cond
        (<= t 36.0) (* x (/ (mevalpoly t 1.00000000000000000000E0 -0.44663998931312457298E-1 0.11209146443112369449E-2
                                       -0.13276124407928422367E-4 0.85118014179823463879E-7 -0.29989314303147656479E-9
                                       0.55401971660186204711E-12 -0.42406353433133212926E-15)
                            (mevalpoly t 1.00000000000000000000E0 0.10891556624243098264E-1 0.59334456769186835896E-4
                                       0.21231112954641805908E-6 0.54747121846510390750E-9 0.10378561511331814674E-11
                                       0.13754880327250272679E-14 0.10223981202236205703E-17)))
        (<= t 144.0) (let [invt (/ t)
                           p (if (neg? x) -HALF_PI HALF_PI)]
                       (- p
                          (* (cos x) (/ (mevalpoly invt 0.99999999962173909991E0 0.36451060338631902917E3
                                                   0.44218548041288440874E5 0.22467569405961151887E7
                                                   0.49315316723035561922E8 0.43186795279670283193E9
                                                   0.11847992519956804350E10 0.45573267593795103181E9)
                                        (* x (mevalpoly invt 1.00000000000000000000E0 0.36651060273229347594E3
                                                        0.44927569814970692777E5 0.23285354882204041700E7
                                                        0.53117852017228262911E8 0.50335310667241870372E9
                                                        0.16575285015623175410E10 0.11746532837038341076E10))))
                          (* (sin x) invt (/ (mevalpoly invt 0.99999999920484901956E0 0.51385504875307321394E3
                                                        0.92293483452013810811E5 0.74071341863359841727E7
                                                        0.28142356162841356551E9 0.49280890357734623984E10
                                                        0.35524762685554302472E11 0.79194271662085049376E11
                                                        0.17942522624413898907E11)
                                             (mevalpoly invt 1.00000000000000000000E0 0.51985504708814870209E3
                                                        0.95292615508125947321E5 0.79215459679762667578E7
                                                        0.31977567790733781460E9 0.62273134702439012114E10
                                                        0.54570971054996441467E11 0.18241750166645704670E12
                                                        0.15407148148861454434E12)))))
        (< t ##Inf) (let [invt (/ t)
                          p (if (neg? x) -HALF_PI HALF_PI)]
                      (- p
                         (* (/ (cos x) x) (- 1.0 (* invt
                                                    (/ (mevalpoly invt 0.19999999999999978257E1 0.22206119380434958727E4
                                                                  0.84749007623988236808E6 0.13959267954823943232E9
                                                                  0.10197205463267975592E11 0.30229865264524075951E12
                                                                  0.27504053804288471142E13 0.21818989704686874983E13)
                                                       (mevalpoly invt 1.00000000000000000000E0 0.11223059690217167788E4
                                                                  0.43685270974851313242E6 0.74654702140658116258E8
                                                                  0.58580034751805687471E10 0.20157980379272098841E12
                                                                  0.26229141857684496445E13 0.87852907334918467516E13)))))
                         (* (sin x) invt (- 1.0 (* invt
                                                   (/ (mevalpoly invt 0.59999999999999993089E1 0.96527746044997139158E4
                                                                 0.56077626996568834185E7 0.15022667718927317198E10
                                                                 0.19644271064733088465E12 0.12191368281163225043E14
                                                                 0.31924389898645609533E15 0.25876053010027485934E16
                                                                 0.12754978896268878403E16)
                                                      (mevalpoly invt 1.00000000000000000000E0 0.16287957674166143196E4
                                                                 0.96636303195787870963E6 0.26839734750950667021E9
                                                                 0.37388510548029219241E11 0.26028585666152144496E13
                                                                 0.85134283716950697226E14 0.11304079361627952930E16
                                                                 0.42519841479489798424E16)))))))
        
        :else (if (neg? x) -HALF_PI HALF_PI)))))

(def ^:private ^:const ^double ci-r0 0.616505485620716233797110404100)
(def ^:private ^:const ^double ci-r1 3.384180422551186426397851146402)
(def ^:private ^:const ^double ci-r01 0.6162109375)
(def ^:private ^:const ^double ci-r02 0.29454812071623379711E-3)
(def ^:private ^:const ^double ci-r11 3.3837890625)
(def ^:private ^:const ^double ci-r12 0.39136005118642639785E-3)

(defn Ci
  "Integral of cos(t)/t from -x to inf"
  ^double [^double x]
  (assert (not (neg? x)) "x must be non-negative")
  (if (nan? x)
    ##NaN
    (let [t (* x x)]
      (cond
        (<= x 3.0) (+ (log (/ x ci-r0)) (* (- (- x ci-r01) ci-r02)
                                           (+ x ci-r0)
                                           (/ (mevalpoly t -0.24607411378767540707E0 0.72113492241301534559E-2
                                                         -0.11867127836204767056E-3 0.90542655466969866243E-6
                                                         -0.34322242412444409037E-8 0.51950683460656886834E-11)
                                              (mevalpoly t 1.00000000000000000000E0 0.12670095552700637845E-1
                                                         0.78168450570724148921E-4 0.29959200177005821677E-6
                                                         0.73191677761328838216E-9 0.94351174530907529061E-12))))
        (<= x 6.0) (+ (log (/ x ci-r1)) (* (- (- x ci-r11) ci-r12)
                                           (+ x ci-r1)
                                           (/ (mevalpoly t -0.15684781827145408780E0 0.66253165609605468916E-2
                                                         -0.12822297297864512864E-3 0.12360964097729408891E-5
                                                         -0.66450975112876224532E-8 0.20326936466803159446E-10
                                                         -0.33590883135343844613E-13 0.23686934961435015119E-16)
                                              (mevalpoly t 1.00000000000000000000E0 0.96166044388828741188E-2
                                                         0.45257514591257035006E-4 0.13544922659627723233E-6
                                                         0.27715365686570002081E-9 0.37718676301688932926E-12
                                                         0.27706844497155995398E-15))))
        (<= x 12.0) (let [invt (/ t)]
                      (- (* (sin x) (/ (mevalpoly invt 0.99999999962173909991E0 0.36451060338631902917E3
                                                  0.44218548041288440874E5 0.22467569405961151887E7
                                                  0.49315316723035561922E8 0.43186795279670283193E9
                                                  0.11847992519956804350E10 0.45573267593795103181E9)
                                       (* x (mevalpoly invt 1.00000000000000000000E0 0.36651060273229347594E3
                                                       0.44927569814970692777E5 0.23285354882204041700E7
                                                       0.53117852017228262911E8 0.50335310667241870372E9
                                                       0.16575285015623175410E10 0.11746532837038341076E10))))
                         (* (cos x) invt (/ (mevalpoly invt 0.99999999920484901956E0 0.51385504875307321394E3
                                                       0.92293483452013810811E5 0.74071341863359841727E7
                                                       0.28142356162841356551E9 0.49280890357734623984E10
                                                       0.35524762685554302472E11 0.79194271662085049376E11
                                                       0.17942522624413898907E11)
                                            (mevalpoly invt 1.00000000000000000000E0 0.51985504708814870209E3
                                                       0.95292615508125947321E5 0.79215459679762667578E7
                                                       0.31977567790733781460E9 0.62273134702439012114E10
                                                       0.54570971054996441467E11 0.18241750166645704670E12
                                                       0.15407148148861454434E12)))))
        (< x ##Inf) (let [invt (/ t)]
                      (- (* (/ (sin x) x) (- 1.0 (* invt
                                                    (/ (mevalpoly invt 0.19999999999999978257E1 0.22206119380434958727E4
                                                                  0.84749007623988236808E6 0.13959267954823943232E9
                                                                  0.10197205463267975592E11 0.30229865264524075951E12
                                                                  0.27504053804288471142E13 0.21818989704686874983E13)
                                                       (mevalpoly invt 1.00000000000000000000E0 0.11223059690217167788E4
                                                                  0.43685270974851313242E6 0.74654702140658116258E8
                                                                  0.58580034751805687471E10 0.20157980379272098841E12
                                                                  0.26229141857684496445E13 0.87852907334918467516E13)))))
                         (* (cos x) invt (- 1.0 (* invt
                                                   (/ (mevalpoly invt 0.59999999999999993089E1 0.96527746044997139158E4
                                                                 0.56077626996568834185E7 0.15022667718927317198E10
                                                                 0.19644271064733088465E12 0.12191368281163225043E14
                                                                 0.31924389898645609533E15 0.25876053010027485934E16
                                                                 0.12754978896268878403E16)
                                                      (mevalpoly invt 1.00000000000000000000E0 0.16287957674166143196E4
                                                                 0.96636303195787870963E6 0.26839734750950667021E9
                                                                 0.37388510548029219241E11 0.26028585666152144496E13
                                                                 0.85134283716950697226E14 0.11304079361627952930E16
                                                                 0.42519841479489798424E16)))))))        
        :else 0.0))))



;;

(def ^:const double-array-type (Class/forName "[D"))
(def ^:const double-double-array-type (Class/forName "[[D"))

(def ^{:doc "Convert double array into sequence.

  Alias for `seq`."
       :metadoc/categories #{:seq}} double-array->seq seq)

#_(defmacro ^:private seq->any-array
    [primitive-type]
    (let [atype (symbol (str primitive-type "-array"))
          clazz (symbol (str atype "-type"))
          fname (symbol (str "seq->" atype))
          docs (str "Convert sequence to " atype ".")]
      `(defn ~fname ~docs [vs#]
         (cond
           (nil? vs#) nil
           (= (type vs#) ~clazz) vs#
           (seqable? vs#) (~atype vs#)
           :else (let [arr# (~atype 1)] 
                   (aset arr# 0 (~primitive-type vs#))
                   arr#)))))

(defn seq->double-array
  "Convert sequence to double array. Returns input if `vs` is double array already."
  {:metadoc/categories #{:seq}}
  ^doubles [vs]
  (cond
    (= (type vs) double-array-type) vs
    (nil? vs) nil
    (seqable? vs) (double-array vs)
    :else (let [arr (double-array 1)] 
            (aset arr 0 (double vs))
            arr)))

(defn double-double-array->seq
  "Convert double array of double arrays into sequence of sequences. "
  {:metadoc/categories #{:seq}}
  [res]
  (seq (map seq res)))

(defn seq->double-double-array
  "Convert sequence to double-array of double-arrays.
  
  If sequence is double-array of double-arrays returns `vss`"
  {:metadoc/categories #{:seq}}
  #^"[[D" [vss]
  (cond 
    (= (type vss) double-double-array-type) vss
    (nil? vss) nil
    :else (into-array (mapv seq->double-array vss))))

;; ## Copy of primitive math machinery
;;
;; Simplified to be used after `ns` is defined.

(def ^:private vars-to-exclude
  '[* + - / > < >= <= == abs rem quot mod bit-or bit-and bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd? bool-and bool-or bool-xor bool-not << >> >>> not==])

(defn- using-primitive-operators? []
  (= #'fastmath.core/+ (resolve '+)))

(defn use-primitive-operators
  "Replaces Clojure's arithmetic and number coercion functions with primitive equivalents.  These are
   defined as macros, so they cannot be used as higher-order functions. This is an idempotent operation. Undo with [[unuse-primitive-operators]]."
  {:metadoc/categories -prim-set-}
  ([] (use-primitive-operators #{}))
  ([skip-set]
   (when-not (using-primitive-operators?)
     (let [v2e (remove skip-set vars-to-exclude)]
       (doseq [v v2e]
         (ns-unmap *ns* v))
       (require ['fastmath.core :refer v2e])))))

(defn unuse-primitive-operators
  "Undoes the work of [[use-primitive-operators]]. This is idempotent."
  {:metadoc/categories -prim-set-}
  []
  (when (using-primitive-operators?)
    (doseq [v vars-to-exclude]
      (ns-unmap *ns* v))
    (refer 'clojure.core)))

