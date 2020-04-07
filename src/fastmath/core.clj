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

  Known from Clojure: `*` `+` `-` `/` `>` `<` `>=` `<=` `==` `rem` `quot` `mod` `bit-or` `bit-and` `bit-xor` `bit-not` `bit-shift-left` `bit-shift-right` `unsigned-bit-shift-right` `inc` `dec` `zero?` `neg?` `pos?` `min` `max` `even?` `odd?`

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
   :exclude [* + - / > < >= <= == rem quot mod bit-or bit-and bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd?])
  (:import [net.jafama FastMath]
           [fastmath.java PrimitiveMath]
           [clojure.lang Numbers]
           [org.apache.commons.math3.util Precision]
           [org.apache.commons.math3.special Erf Gamma Beta BesselJ])
  (:require [fastmath.core :as m]))

(set! *warn-on-reflection* true)
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
         doc (or (:doc (meta name)) (str "A primitive math version of `" name "`"))]
     `(defmacro ~name
        ~doc
        ([~x]
         ~((eval single-arg-form) x))
        ([~x ~y]
         (list '~fname ~x ~y))
        ([~x ~y ~'& ~rest]
         (list* '~name (list '~name ~x ~y) ~rest))))))

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
         doc (or (:doc (meta name)) (str "A primitive math version of `" name "`"))]
     `(defmacro ~name
        ~doc
        ([~x]
         ~((eval single-arg-form) x))
        ([~x ~y]
         (list '~fname ~x ~y))
        ([~x ~y ~'& ~rest]
         (list 'fastmath.java.PrimitiveMath/and (list '~name ~x ~y) (list* '~name ~y ~rest)))))))

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
(variadic-proxy ^{:metadoc/categories #{:bool}} bool-and and)
(variadic-proxy ^{:metadoc/categories #{:bool}} bool-or or)
(variadic-proxy ^{:metadoc/categories #{:bool}} bool-xor xor)
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
   (bool-and (== a b) (== b c)))
  ([^double a ^double b ^double c ^double d]
   (bool-and (== a b) (== b c) (== c d))))

;; Processing math constants
(def ^:const ^double ^{:doc "Value of \\\\(\\pi\\\\)"} PI Math/PI)
(def ^:const ^double ^{:doc "Value of \\\\(\\frac{\\pi}{2}\\\\)"} HALF_PI (/ PI 2.0))
(def ^:const ^double ^{:doc "Value of \\\\(\\frac{\\pi}{3}\\\\)"} THIRD_PI (/ PI 3.0))
(def ^:const ^double ^{:doc "Value of \\\\(\\frac{\\pi}{4}\\\\)"} QUARTER_PI (/ PI 4.0))
(def ^:const ^double ^{:doc "Value of \\\\(2 {\\pi}\\\\)"} TWO_PI (+ PI PI))
(def ^:const ^double ^{:doc "Alias for [[TWO_PI]]"} TAU TWO_PI)
(def ^:const ^double ^{:doc "Value of \\\\(e\\\\)"} E Math/E)

(def ^:const ^double ^{:doc "Value of \\\\(\\-pi\\\\)"} -PI (- Math/PI))
(def ^:const ^double ^{:doc "Value of \\\\(-\\frac{\\pi}{2}\\\\)"} -HALF_PI (/ -PI 2.0))
(def ^:const ^double ^{:doc "Value of \\\\(-\\frac{\\pi}{3}\\\\)"} -THIRD_PI (/ -PI 3.0))
(def ^:const ^double ^{:doc "Value of \\\\(-\\frac{\\pi}{4}\\\\)"} -QUARTER_PI (/ -PI 4.0))
(def ^:const ^double ^{:doc "Value of \\\\(-2 {\\pi}\\\\)"} -TWO_PI (+ -PI -PI))
(def ^:const ^double ^{:doc "Alias for [[TWO_PI-]]"} -TAU -TWO_PI)
(def ^:const ^double ^{:doc "Value of \\\\(e\\\\)"} -E (- Math/E))

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

(defn
  ^{:doc "log(1+exp(x))"
    :metadoc/categories -pow-set-}
  log1pexp ^double [^double x] (log1p (exp x)))

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

;; Beta

(beta-proxy :two ^{:doc "Logarithm of Beta function." :metadoc/categories -special-set-} log-beta logBeta)
(beta-proxy :three ^{:doc "Regularized `Beta`." :metadoc/categories -special-set-} regularized-beta regularizedBeta)

;; BesselJ
(besselj-proxy :two ^{:doc "Bessel J function value for given order and argument." :metadoc/categories -special-set-} bessel-j value)

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
(defn floor "\\\\(\\lfloor x \\rfloor\\\\). See: [[qfloor]]." {:metadoc/categories -round-set-} ^double [^double x] (FastMath/floor x))
(defn ceil "\\\\(\\lceil x \\rceil\\\\). See: [[qceil]]." {:metadoc/categories -round-set-} ^double [^double x] (FastMath/ceil x))
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
   (let [ui (Double/doubleToRawLongBits (if (zero? v) 0.0 v))]
     (Double/longBitsToDouble (if (neg? v) (dec ui) (inc ui)))))
  (^double [^double v ^long delta]
   (let [ui (Double/doubleToRawLongBits (if (zero? v) 0.0 v))]
     (Double/longBitsToDouble (if (neg? v) (- ui delta) (+ ui delta))))))

(defn prev-double
  "Previous double value. Optional value `delta` sets step amount."
  (^double [^double v]
   (let [ui (Double/doubleToRawLongBits (if (zero? v) 0.0 v))]
     (Double/longBitsToDouble (if (pos? v) (dec ui) (inc ui)))))
  (^double [^double v ^long delta]
   (let [ui (Double/doubleToRawLongBits (if (zero? v) 0.0 v))]
     (Double/longBitsToDouble (if (pos? v) (- ui delta) (+ ui delta))))))

;; More constants

;; \\(\sqrt{2}\\)
(def ^:const ^double ^{:doc "\\\\(\\sqrt{2}\\\\)"} SQRT2 (sqrt 2.0))
(def ^:const ^double ^{:doc "\\\\(\\frac{\\sqrt{2}}{2}\\\\)"} SQRT2_2 (* 0.5 SQRT2))

;; \\(\sqrt{3}\\)
(def ^:const ^double ^{:doc "\\\\(\\sqrt{3}\\\\)"} SQRT3 (sqrt 3.0))
(def ^:const ^double ^{:doc "\\\\(\\frac{\\sqrt{3}}{2}\\\\)"} SQRT3_2 (* 0.5 (sqrt 3.0)))
(def ^:const ^double ^{:doc "\\\\(\\frac{\\sqrt{3}}{3}\\\\)"} SQRT3_3 (/ (sqrt 3.0) 3.0))

;; \\(\sqrt{5}\\)
(def ^:const ^double ^{:doc "\\\\(\\sqrt{5}\\\\)" }SQRT5 (sqrt 5.0))

;; \\(\sqrt{\pi}\\)
(def ^:const ^double ^{:doc "\\\\(\\sqrt{\\pi}\\\\)"} SQRTPI (sqrt PI))
(def ^:const ^double ^{:doc "\\\\(\\sqrt{2\\pi}\\\\)"} SQRT2PI (sqrt TWO_PI))

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
(def ^:const ^double ^{:doc "\\\\(\\sqrt{2}\\\\)"} M_SQRT2 SQRT2)
(def ^:const ^double ^{:doc "\\\\(\\frac{1}{\\sqrt{2}}\\\\)"} M_SQRT1_2 (/ SQRT2))

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
                                           (map #(m/mnorm % 0.0 (dec cnt) start end) (range cnt))))
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
         n- (dec n)
         r (/ n (+ (* number o-) overlap))
         ii-mult (*  o- r)
         ii (mapv #(* ^long % ii-mult) (range number))
         x1 (map #(x (round-even %)) ii)
         xr (map #(x (dec (round-even (+ r ^double %)))) ii)
         _ (println x n r ii x1 xr)
         diffs (filter #(pos? ^double %) (mapv (fn [[^double x ^double y]] (- y x)) (partition 2 1 x)))
         eps (* 0.5 (double (if (seq diffs) (reduce fast-min diffs) 0.0)))]
     (println diffs eps)
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

  Possible tie strategies: `:average`, `:first`, `:last`, `:random`, `:min`, `:max`"
  {:metadoc/categories #{:rank}}
  ([vs] (rank vs :average))
  ([vs ties]
   (let [indexed-sorted-map (group-by second (map-indexed (fn [^long idx v] [(inc idx) v]) (sort vs)))]
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
                      :max (comp first last)
                      (fn ^double [v] (/ ^double (reduce #(+ ^double %1 ^double %2) (map first v)) (count v))))
             m (into {} (map (fn [[k v]] [k (tie-fn v)]) indexed-sorted-map))]
         (map m vs))))))

(defn order
  "Ordering permutation. See [R docs](https://www.rdocumentation.org/packages/base/versions/3.6.1/topics/order)"
  {:metadoc/categories #{:rank}}
  ([vs] (order vs false))
  ([vs decreasing?]
   (->> (map-indexed vector vs)
        (sort-by second (if decreasing?
                          clojure.core/>
                          clojure.core/<))
        (map (comp #(inc ^long %) first)))))

;;

(def ^:const double-array-type (Class/forName "[D"))
(def ^:const double-double-array-type (Class/forName "[[D"))

(def ^{:doc "Convert double array into sequence.

  Alias for `seq`."
       :metadoc/categories #{:seq}} double-array->seq seq)

(defmacro ^:private seq->any-array
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
  '[* + - / > < >= <= == rem quot mod bit-or bit-and bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd? bool-and bool-or bool-xor bool-not << >> >>> not==])

(defn- using-primitive-operators? []
  (= #'fastmath.core/+ (resolve '+)))

(defn use-primitive-operators
  "Replaces Clojure's arithmetic and number coercion functions with primitive equivalents.  These are
   defined as macros, so they cannot be used as higher-order functions. This is an idempotent operation. Undo with [[unuse-primitive-operators]]."
  {:metadoc/categories -prim-set-}
  []
  (when-not (using-primitive-operators?)
    (doseq [v vars-to-exclude]
      (ns-unmap *ns* v))
    (require (vector 'fastmath.core :refer vars-to-exclude))))

(defn unuse-primitive-operators
  "Undoes the work of [[use-primitive-operators]]. This is idempotent."
  {:metadoc/categories -prim-set-}
  []
  (when (using-primitive-operators?)
    (doseq [v vars-to-exclude]
      (ns-unmap *ns* v))
    (refer 'clojure.core)))

