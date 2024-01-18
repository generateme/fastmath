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

  Known from Clojure: `*` `+` `-` `/` `>` `<` `>=` `<=` `==` `rem` `quot` `mod` `bit-or` `bit-and` `bit-xor` `bit-and-not` `bit-set` `bit-clear` `bit-test` `bit-flip` `bit-not` `bit-shift-left` `bit-shift-right` `unsigned-bit-shift-right` `inc` `dec` `zero?` `neg?` `pos?` `min` `max` `even?` `odd?` `abs`

  And additionally:

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
  (:refer-clojure
   :exclude [* + - / > < >= <= == rem quot mod bit-or bit-and bit-and-not bit-set bit-clear bit-test bit-flip bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd? abs])
  (:import [net.jafama FastMath]
           [fastmath.java PrimitiveMath]
           [org.apache.commons.math3.util Precision]
           [org.apache.commons.math3.special Gamma]))

(set! *unchecked-math* :warn-on-boxed)

;; which java?

(def ^{:const true :tag 'long} jvm-version
  (->> (System/getProperty "java.version")
       (re-seq #"\d+")
       (first)
       (Long/parseLong)))

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

(variadic-proxy + add)
(variadic-proxy - subtract (fn [x] `(list 'fastmath.java.PrimitiveMath/negate ~x)))
(variadic-proxy * multiply)
(variadic-proxy / divide (fn [x] `(list 'fastmath.java.PrimitiveMath/reciprocal ~x)))
(primitivemath-proxy :one inc)
(primitivemath-proxy :one dec)
(primitivemath-proxy :two rem remainder)
(primitivemath-proxy :two quot quotient)
(primitivemath-proxy :two mod modulus)
(variadic-proxy bit-and bitAnd)
(variadic-proxy bit-nand bitNand)
(variadic-proxy bit-and-not bitAndNot)
(variadic-proxy bit-or bitOr)
(variadic-proxy bit-nor bitNor)
(variadic-proxy bit-xor bitXor)
(primitivemath-proxy :one bit-not bitNot)
(primitivemath-proxy :two bit-set bitSet)
(primitivemath-proxy :two bit-clear bitClear)
(primitivemath-proxy :two bit-flip bitFlip)
(primitivemath-proxy :two bit-test bitTest)
(variadic-proxy ^{:deprecated true} bool-and and)
(variadic-proxy ^{:deprecated true} bool-or or)
(variadic-proxy ^{:deprecated true} bool-xor xor)
(primitivemath-proxy :one bool-not not)
(variadic-proxy min)
(variadic-proxy max)
(primitivemath-proxy :one zero? isZero)
(primitivemath-proxy :one one? isOne)
(primitivemath-proxy :one neg? isNeg)
(primitivemath-proxy :one pos? isPos)
(primitivemath-proxy :one not-neg? isNNeg)
(primitivemath-proxy :one not-pos? isNPos)
(primitivemath-proxy :one even? isEven)
(primitivemath-proxy :one odd? isOdd)
(primitivemath-proxy :two << shiftLeft)
(primitivemath-proxy :two >> shiftRight)
(primitivemath-proxy :two >>> unsignedShiftRight)
(primitivemath-proxy :two bit-shift-left shiftLeft)
(primitivemath-proxy :two bit-shift-right shiftRight)
(primitivemath-proxy :two unsigned-bit-shift-right unsignedShiftRight)

(variadic-predicate-proxy < lt)
(variadic-predicate-proxy > gt)
(variadic-predicate-proxy <= lte)
(variadic-predicate-proxy >= gte)
(variadic-predicate-proxy ^{:doc "Equality. See also [[eq]] for function version."} == eq)
(variadic-predicate-proxy not== neq)

(defn negative-zero?
  "Check if zero is negative, ie. -0.0"
  [^double x]
  (== (Double/doubleToLongBits x) -9223372036854775808)) ;; -0.0

(defn fast+
  {:inline (fn [x y] `(+ ~x ~y)) :inline-arities #{2}
   :doc "Primitive `+` for two doubles as function."}
  ^double [^double a ^double b] (+ a b))
(defn fast-
  {:inline (fn [x y] `(- ~x ~y)) :inline-arities #{2}
   :doc "Primitive `-` for two doubles as function."}
  ^double [^double a ^double b] (- a b))
(defn fast*
  {:inline (fn [x y] `(* ~x ~y)) :inline-arities #{2}
   :doc "Primitive `*` for two doubles as function."}
  ^double [^double a ^double b] (* a b))
(defn fast-max
  {:inline (fn [x y] `(max ~x ~y)) :inline-arities #{2}
   :doc "Primitive `max` for two doubles as function."}
  ^double [^double a ^double b] (max a b))
(defn fast-min
  {:inline (fn [x y] `(min ~x ~y)) :inline-arities #{2}
   :doc "Primitive `min` for two doubles as function"}
  ^double [^double a ^double b] (min a b))
(defn fast-identity
  {:inline (fn [x] `~x) :inline-arities #{1}
   :doc "Identity on double."}
  ^double [^double a] a)

;; Primitive math eq
(defn eq 
  "Primitive math equality function for doubles. See [[==]]."
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
  (if (< jvm-version 9)
    `(+ ~z (* ~x ~y))
    `(Math/fma ~x ~y ~z)))

(defmacro fma
  "`[x y z]` -> `(+ z (* x y))` or `Math/fma` for java 9+"
  [x y z]
  `(muladd ~x ~y ~z))

(defmacro negmuladd
  "`[x y z]` -> `(+ z (* -1.0 x y)`"
  [x y z]
  `(muladd (- ~x) ~y ~z))

(defmacro mevalpoly
  "Evaluate polynomial macro version in the form coeffs[0]+coeffs[1]*x+coeffs[2]*x^2+...."
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
  "Kahan's algorithm for (a*b)-(c*d) to avoid catastrophic cancellation."
  ^double [^double a ^double b ^double c ^double d]
  (let [cd (* c d)]
    (+ (fma a b (- cd)) (fma (- c) d cd))))

(defn sum-of-products
  "Kahan's algorithm for (a*b)+(c*d) to avoid catastrophic cancellation."
  ^double [^double a ^double b ^double c ^double d]
  (let [cd (* c d)]
    (+ (fma a b cd) (fma c d (- cd)))))

;; Processing math constants
(def ^{:const true :tag 'double :doc "Value of \\\\(\\pi\\\\)"} PI Math/PI)
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{\\pi}{2}\\\\)"} HALF_PI (* PI 0.5))
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{\\pi}{3}\\\\)"} THIRD_PI (/ PI 3.0))
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{\\pi}{4}\\\\)"} QUARTER_PI (* PI 0.25))
(def ^{:const true :tag 'double :doc "Value of \\\\(2 {\\pi}\\\\)"} TWO_PI (+ PI PI))
(def ^{:const true :tag 'double :doc "Alias for [[TWO_PI]]"} TAU TWO_PI)
(def ^{:const true :tag 'double :doc "Value of \\\\(e\\\\)"} E Math/E)
(def ^{:const true :tag 'double :doc "Value of \\\\(-\\pi\\\\)"} -PI (- Math/PI))
(def ^{:const true :tag 'double :doc "Value of \\\\(-\\frac{\\pi}{2}\\\\)"} -HALF_PI (* PI -0.5))
(def ^{:const true :tag 'double :doc "Value of \\\\(-\\frac{\\pi}{3}\\\\)"} -THIRD_PI (/ -PI -3.0))
(def ^{:const true :tag 'double :doc "Value of \\\\(-\\frac{\\pi}{4}\\\\)"} -QUARTER_PI (* PI -0.25))
(def ^{:const true :tag 'double :doc "Value of \\\\(-2 {\\pi}\\\\)"} -TWO_PI (- TWO_PI))
(def ^{:const true :tag 'double :doc "Alias for [[TWO_PI-]]"} -TAU -TWO_PI)
(def ^{:const true :tag 'double :doc "Value of \\\\(-e\\\\)"} -E (- Math/E))
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{1}{\\pi}\\\\)"} INV_PI (/ PI))
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{2}{\\pi}\\\\)"} TWO_INV_PI (/ 2.0 PI))
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{4}{\\pi}\\\\)"} FOUR_INV_PI (/ 4.0 PI))
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{1}{2 \\pi}\\\\)"} INV_TWO_PI (/ TWO_PI))
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{1}{4 \\pi}\\\\)"} INV_FOUR_PI (/ (* 2.0 TWO_PI)))
(def ^{:const true :tag 'double :doc "Very small number \\\\(\\varepsilon\\\\)"} EPSILON 1.0e-10)
(def ^{:const true :tag 'double :doc "Euler-Mascheroni constant"} GAMMA Gamma/GAMMA)
(def ^{:const true :tag 'double :doc "Lanchos approximation `g` constant"} LANCZOS_G Gamma/LANCZOS_G)
(def ^{:const true :tag 'double :doc "Catalan G"} CATALAN_G 0.915965594177219015054603514932384110774)

(defonce ^{:const true :tag 'double :doc "Smallest machine number. Value is calculated during evaluation and may differ on different processors."}
  MACHINE-EPSILON (* 0.5 (double (loop [d (double 1.0)]
                                   (if (not== 1.0 (+ 1.0 (* d 0.5)))
                                     (recur (* d 0.5))
                                     d)))))

(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{1}{3}\\\\)"} THIRD (/ 3.0))
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{2}{3}\\\\)"} TWO_THIRD (/ 2.0 3.0))
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{1}{6}\\\\)"} SIXTH (/ 6.0))

;; Trigonometry
(fastmath-proxy :one sin)
(fastmath-proxy :one cos)
(fastmath-proxy :one tan)
(fastmath-proxy :one asin)
(fastmath-proxy :one acos)
(fastmath-proxy :one atan)
(fastmath-proxy :one sinh)
(fastmath-proxy :one cosh)
(fastmath-proxy :one tanh)
(fastmath-proxy :one asinh)
(fastmath-proxy :one acosh)
(fastmath-proxy :one atanh)

(fastmath-proxy :one ^{:doc "Fast and less accurate [[sin]]."} qsin sinQuick)
(fastmath-proxy :one ^{:doc "Fast and less accurate [[cos]]."} qcos cosQuick)

;; Additional trigonometry functions
(defn cot "Cotangent" ^double [^double v] (/ (FastMath/tan v)
                                          #_(FastMath/tan (- HALF_PI v))))
(defn sec "Secant" ^double [^double v] (/ (FastMath/cos v)))
(defn csc "Cosecant" ^double [^double v] (/ (FastMath/sin v)))

;; Additional cyclometric functions
(defn acot "Arccotangent" ^double [^double v] (- HALF_PI (FastMath/atan v)))
(defn asec "Arcsecant" ^double [^double v] (FastMath/acos (/ 1.0 v)))
(defn acsc "Arcosecant" ^double [^double v] (FastMath/asin (/ 1.0 v)))

(fastmath-proxy :two atan2)

;; Additional hyperbolic functions
(defn coth "Hyperbolic cotangent"^double [^double v] (/ (FastMath/tanh v)))
(defn sech "Hyperbolic secant" ^double [^double v] (/ (FastMath/cosh v)))
(defn csch "Hyperbilic cosecant" ^double [^double v] (/ (FastMath/sinh v)))

;; Additional inverse hyperbolic functions
(defn acoth "Area hyperbolic cotangent" ^double [^double v] (FastMath/atanh (/ v)))
(defn asech "Area hyperbolic secant" ^double [^double v] (FastMath/acosh (/ v)))
(defn acsch "Area hyperbolic cosecant" ^double [^double v] (FastMath/asinh (/ v)))

;; historical

(defn crd "Chord" ^double [^double v] (* 2.0 (FastMath/sin (* 0.5 v))))
(defn acrd "Inverse chord" ^double [^double v] (* 2.0 (FastMath/asin (* 0.5 v))))

(defn versin "Versine" ^double [^double v] (- 1.0 (FastMath/cos v)))
(defn coversin "Coversine" ^double [^double v] (- 1.0 (FastMath/sin v)))
(defn vercos "Vercosine" ^double [^double v] (inc (FastMath/cos v)))
(defn covercos "Covercosine" ^double [^double v] (inc (FastMath/sin v)))

(defn aversin "Arc versine" ^double [^double v] (FastMath/acos (- 1.0 v)))
(defn acoversin "Arc coversine" ^double [^double v] (FastMath/asin (- 1.0 v)))
(defn avercos "Arc vecosine" ^double [^double v] (FastMath/acos (dec v)))
(defn acovercos "Arc covercosine" ^double [^double v] (FastMath/asin (dec v)))

(defn haversin
  "Haversine formula for value or lattitude and longitude pairs."
  (^double [^double v] (* 0.5 (- 1.0 (FastMath/cos v))))
  (^double [[^double lat1 ^double lon1] [^double lat2 ^double lon2]]
   (haversin lat1 lon1 lat2 lon2))
  (^double [^double lat1 ^double lon1 ^double lat2 ^double lon2]
   (+ (haversin (- lat2 lat1))
      (* (FastMath/cos lat1)
         (FastMath/cos lat2)
         (haversin (- lon2 lon1))))))

(def ^{:doc "Haversine ([[haversin]] alias)"} haversine haversin)
(defn hacoversin "Hacoversine" ^double [^double v] (* 0.5 (- 1.0 (FastMath/sin v))))
(defn havercos "Havercosine" ^double [^double v] (* 0.5 (inc (FastMath/cos v))))
(defn hacovercos "Hacovercosine" ^double [^double v] (* 0.5 (inc (FastMath/sin v))))

(defn ahaversin "Arc haversine" ^double [^double v] (FastMath/acos (- 1.0 (* 2.0 v))))
(defn ahacoversin "Arc hacoversine" ^double [^double v] (FastMath/asin (- 1.0 (* 2.0 v))))
(defn ahavercos "Arc havecosine" ^double [^double v] (FastMath/acos (dec (* 2.0 v))))
(defn ahacovercos "Arc hacovercosine" ^double [^double v] (FastMath/asin (dec (* 2.0 v))))

(defn exsec "Exsecant" ^double [^double v] (dec (sec v)))
(defn excsc "Excosecant" ^double [^double v] (dec (csc v)))
(defn aexsec "Arc exsecant" ^double [^double v] (asec (inc v)))
(defn aexcsc "Arc excosecant" ^double [^double v] (acsc (inc v)))

(defn haversine-dist
  "Haversine distance `d` for `r=1`"
  (^double [[^double lat1 ^double lon1] [^double lat2 ^double lon2]]
   (haversine-dist lat1 lon1 lat2 lon2))
  (^double [^double lat1 ^double lon1 ^double lat2 ^double lon2]
   (* 2.0 (FastMath/asin (FastMath/sqrt (haversin lat1 lon1 lat2 lon2))))))

;; exp and log
(fastmath-proxy :one exp)
(fastmath-proxy :one log)
(fastmath-proxy :one ^{:doc "\\\\(\\ln_{10}{x}\\\\)"} log10)
;; Alias for natural logarithm
(fastmath-proxy :one ln log)

(fastmath-proxy :one log1p)
(fastmath-proxy :one expm1)

(def ^{:const true :tag 'double :doc "\\\\(\\ln{2}\\\\)"} LN2 (log 2.0))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{1}{\\ln{2}}\\\\)"} INV_LN2 (/ LN2))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{\\ln{2}}{2}\\\\)"} LN2_2 (* 0.5 LN2))
(def ^{:const true :tag 'double :doc "\\\\(\\ln{10}\\\\)"} LN10 (log 10.0))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{1}{\\ln{0.5}}\\\\)"} INV_LOG_HALF (/ (log 0.5)))
(def ^{:const true :tag 'double :doc "\\\\(\\ln{0.5}\\\\)"} LOG_HALF (log 0.5))
(def ^{:const true :tag 'double :doc "\\\\(\\ln{\\pi}\\\\)"} LOG_PI (log PI))
(def ^{:const true :tag 'double :doc "\\\\(\\ln{2 \\pi}\\\\)"} LOG_TWO_PI (log TWO_PI))

(defn log1pexp
  "log(1+exp(x))"
  ^double [^double x]
  (cond
    (< x -745.1332191019412) 0.0
    (< x -36.7368005696771) (FastMath/exp x)
    (< x 18.021826694558577) (FastMath/log1p (FastMath/exp x))
    (< x 33.23111882352963) (+ x (FastMath/exp (- x)))
    :else x))

(defn log1mexp
  "log(1-exp(x))"
  ^double [^double x]
  (if (< x LOG_HALF)
    (FastMath/log1p (- (FastMath/exp x)))
    (FastMath/log (- (FastMath/expm1 x)))))

(defn log2mexp
  "log(2-exp(x))"
  ^double [^double x]
  (FastMath/log1p (- (FastMath/expm1 x))))

(defn log1psq
  "log(1+x^2))"
  ^double [^double x]
  (if (< x 9007199254740992)
    (FastMath/log1p (* x x))
    (* 2.0 (log x))))

(defn logexpm1
  "log(exp(x)-1))"
  ^double [^double x] (FastMath/log (FastMath/expm1 x)))

;; from julia
(defn- log1pmx-ker
  ^double [^double x]
  (let [r (/ x (+ 2.0 x))
        t (* r r)
        w (mevalpoly t 6.66666666666666667e-1 4.00000000000000000e-1 2.85714285714285714e-1 2.22222222222222222e-1
                     1.81818181818181818e-1 1.53846153846153846e-1 1.33333333333333333e-1 1.17647058823529412e-1)
        hxsq (* 0.5 x x)]
    (- (* r (+ hxsq (* w t))) hxsq)))

(defn log1pmx
  "log(1+x)-x"
  ^double [^double x]
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

(defn logmxp1
  "log(x)-x+1"
  ^double [^double x]
  (cond
    (<= x 0.3) (- (inc (FastMath/log x)) x)
    (<= x 0.4) (let [u (/ (- x 0.375) 0.375)]
                 (+ (log1pmx-ker u) -3.55829253011726237e-1 (* 0.625 u)))
    (<= x 0.6) (let [u (* (- x 0.5) 2.0)]
                 (+ (log1pmx-ker u) -1.93147180559945309e-1 (* 0.5 u)))
    :else (log1pmx (dec x))))

(defn logaddexp
  "log(exp(x)+exp(y))"
  ^double [^double x ^double y]
  (if (< x y)
    (+ y (log1pexp (- x y)))
    (+ (if-not (Double/isNaN y) x y)
       (log1pexp (- y x)))))

(defn logsubexp
  "log(abs(exp(x)-exp(y)))"
  ^double [^double x ^double y]
  (+ (max x y)
     (log1mexp (- (if (and (== x y)
                           (or (Double/isFinite x) (neg? x))) 0.0 (FastMath/abs (- x y)))))))


(defn logsumexp
  "log(exp(x1)+...+exp(xn))"
  ^double [xs]
  (loop [[^double x & rst] xs
         r 0.0
         alpha ##-Inf]
    (if (<= x alpha)
      (let [nr (+ r (FastMath/exp (- x alpha)))]
        (if-not (seq rst)
          (+ (FastMath/log nr) alpha)
          (recur rst nr alpha)))
      (let [nr (inc (* r (FastMath/exp (- alpha x))))]
        (if-not (seq rst)
          (+ (FastMath/log nr) x)
          (recur rst nr (double x)))))))

(defn xlogx
  "x * log(x)"
  ^double [^double x]
  (if (zero? x) 0.0 (* x (FastMath/log x))))

(defn xlogy
  "x * log(y)"
  ^double [^double x ^double y]
  (if (and (zero? x)
           (not (Double/isNaN y))) 0.0 (* x (log y))))

(defn xlog1py
  "x * log(1+y)"
  ^double [^double x ^double y]
  (if (and (zero? x)
           (not (Double/isNaN y))) 0.0 (* x (log1p y))))

(defn cloglog
  "log(-log(1-x))"
  ^double [^double x]
  (FastMath/log (- (FastMath/log1p (- x)))))

(defn xexpx
  "x * exp(x)"
  ^double [^double x]
  (let [expx (exp x)]
    (if (zero? expx) 0.0 (* x expx))))

(defn xexpy
  "x * exp(x)"
  ^double [^double x ^double y]
  (let [expy (exp y)]
    (if (and (zero? expy)
             (not (Double/isNaN x))) 0.0 (* x expy))))

(defn cexpexp
  "1-exp(-exp(x))"
  ^double [^double x]
  (- (FastMath/expm1 (- (FastMath/exp x)))))

;; Quick logarithm
(fastmath-proxy :one ^{:doc "Fast and less accurate version of [[log]]."} qlog logQuick)

;; Roots (square and cubic)
(fastmath-proxy :one ^{:doc "\\\\(\\sqrt{x}\\\\)"} sqrt)
(fastmath-proxy :one ^{:doc "\\\\(\\sqrt[3]{x}\\\\)"} cbrt)

;; Quick version of exponential \\(e^x\\)
(fastmath-proxy :one ^{:doc "Quick and less accurate version of [[exp]]."} qexp expQuick)

;; Radians to degrees (and opposite) conversions
(def ^{:const true :tag 'double :doc "\\\\(\\frac{180}{\\pi}\\\\)"} rad-in-deg (/ 180.0 PI))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{\\pi}{180}\\\\)"} deg-in-rad (/ PI 180.0))
(defn radians "Convert degrees into radians." ^double [^double deg] (* deg-in-rad deg))
(defn degrees "Convert radians into degrees." ^double [^double rad] (* rad-in-deg rad))

;; Erf
(erf-proxy :onetwo ^{:doc "Error function. For two arguments return difference between `(erf x)` and `(erf y)`."} erf)
(erf-proxy :one ^{:doc "Complementary error function."} erfc)
(erf-proxy :one ^{:doc "Inverse [[erf]]."} inv-erf erfInv)
(erf-proxy :one ^{:doc "Inverse [[erfc]]."} inv-erfc erfcInv)

;; Gamma

(gamma-proxy :one ^{:doc "Gamma function \\\\(\\Gamma(x)\\\\)"} gamma gamma)
(gamma-proxy :one ^{:doc "Log of Gamma function \\\\(\\ln\\Gamma(x)\\\\)"} log-gamma logGamma)
(gamma-proxy :one ^{:doc "Log of Gamma function \\\\(\\ln\\Gamma(1+x)\\\\)"} log-gamma-1p logGamma1p)
(gamma-proxy :one ^{:doc "Logarithmic derivative of \\\\(\\Gamma\\\\)."} digamma)
(gamma-proxy :one ^{:doc "Derivative of [[digamma]]."} trigamma)
(gamma-proxy :one ^{:doc "\\\\(\\frac{1}{\\Gamma(1+x)}-1\\\\)."} inv-gamma-1pm1 invGamma1pm1)
(gamma-proxy :two ^{:doc "Regularized `gamma` P"} regularized-gamma-p regularizedGammaP)
(gamma-proxy :two ^{:doc "Regularized `gamma` Q"} regularized-gamma-q regularizedGammaQ)

(defn minkowski
  "Minkowski's question mark function ?(x)"
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

(beta-proxy :two ^{:doc "Logarithm of Beta function."} log-beta logBeta)
(beta-proxy :three ^{:doc "Regularized `Beta`."} regularized-beta regularizedBeta)

;; BesselJ
(besselj-proxy :two ^{:doc "Bessel J function value for given order and argument."} bessel-j value)

;; jinc-c4 (/ (* PI PI PI PI) 192.0)
;; jinc-c2 (/ (* PI PI) -8.0)

(defn jinc
  "Besselj1 devided by `x`"
  ^double [^double x]
  (if (< (FastMath/abs x) 0.002)
    (let [x2 (* x x)]
      (mevalpoly x2 1.0 -1.2337005501361697 0.5073390158020964))
    (let [pix (* PI x)]
      (* 2.0 (/ (bessel-j 1 pix) pix)))))

;; I0

(defn I0
  "Modified Bessel function of the first kind, order 0."
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

(defn log-I0
  "Log of [[I0]]."
  ^double [^double x]
  (FastMath/log (I0 x)))

;; Sinc
(defn sinc
  "Sinc function."
  ^double [^double v]
  (let [x (* PI (FastMath/abs v))]
    (if (< x 1.0e-5) 1.0
        (/ (FastMath/sin x) x))))

;;
(defn sigmoid
  "Sigmoid function"
  ^double [^double x]
  (/ (inc (FastMath/exp (- x)))))

(def ^{:doc "Alias for [[sigmoid]]"} logistic sigmoid)

(defn logit
  "Logit function"
  ^double [^double x]
  (FastMath/log (/ x (- 1.0 x))))

(defn log2
  "Logarithm with base 2.

  \\\\(\\ln_2{x}\\\\)"
  ^double [^double x]
  (* (FastMath/log x) INV_LN2))

;; \\(\log_b x\\)
(defn logb
  "Logarithm with base `b`.

  \\\\(\\ln_b{x}\\\\)"
  ^double [^double b ^double x]
  (/ (FastMath/log x) (FastMath/log b)))

(defn logcosh
  "log(cosh(x))"
  ^double [^double x]
  (let [absx (FastMath/abs x)]
    (- (+ absx (log1pexp (* -2.0 absx))) LN2)))

;; \\(\log_2 e\\)
(def ^{:const true :tag 'double :doc "\\\\(\\log_{2}{e}\\\\)"} LOG2E (log2 E))

;; \\(\log_{10} e\\)
(def ^{:const true :tag 'double :doc "\\\\(\\log_{10}{e}\\\\)"} LOG10E (log10 E))

;; Powers (normal, quick)
(fastmath-proxy :two  pow)
(fastmath-proxy :two ^{:doc "Fast and less accurate version of [[pow]]."} qpow powQuick)

;; Fast version of power, second parameter should be integer
(fastmath-proxy :two ^{:doc "Fast version of pow where exponent is integer."} fpow powFast)

(def ^:private factorial20-table [1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600
                                6227020800 87178291200 1307674368000 20922789888000
                                355687428096000 6402373705728000 121645100408832000
                                2432902008176640000])

(defn factorial20
  "Factorial table up to 20!"
  ^long [^long n]
  (factorial20-table n))

(defn factorial "Factorial"
  ^double [^long n]
  (if (< n 21)
    (factorial20-table n)
    (exp (log-gamma (double (inc n))))))

(defn ^{:doc "Log factorial, alias to log-gamma"} log-factorial
  ^double [^long x] (log-gamma (double (inc x))))

(defn combinations
  "Binomial coefficient (n choose k)"
  ^double [^long n ^long k]
  (let [k (min k (- n k))]
    (cond
      (neg? k) 0.0
      (zero? k) 1.0
      (< k 30) (loop [j (long 2)
                      r (double n)]
                 (if (> j k)
                   r
                   (recur (inc j) (* r (/ (inc (- n j)) (double j))))))
      :else (exp (- (- (ln (inc n)))
                    (log-beta (inc (- n k)) (inc k)))))))

(defn log-combinations
  "Log of binomial coefficient (n choose k)"
  ^double [^long n ^long k]
  (let [k (min k (- n k))]
    (cond
      (neg? k) ##-Inf
      (zero? k) 0.0
      (one? k) (ln n)
      (< n k) ##-Inf
      (== n k) 0.0
      :else (- (- (ln (inc n)))
               (log-beta (inc (- n k)) (inc k))))))

;; Square and cubic
(defn sq "Same as [[pow2]]. \\\\(x^2\\\\)" ^double [^double x] (* x x))
(defn pow2 "Same as [[sq]]. \\\\(x^2\\\\)" ^double [^double x] (* x x))
(defn pow3 "\\\\(x^3\\\\)" ^double [^double x] (* x (* x x)))
(defn cb "\\\\(x^3\\\\)" ^double [^double x] (* x (* x x)))

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
  ^double [^double value]
  (if (neg? value) 0.0 (FastMath/sqrt value)))

;; Approximated sqrt via binary operations (error 1.0E-2)
(fastmath-proxy :one ^{:doc "Approximated [[sqrt]] using binary operations with error `1.0E-2`."} qsqrt sqrtQuick)
(fastmath-proxy :one ^{:doc "Inversed version of [[qsqrt]]. Quick and less accurate."} rqsqrt invSqrtQuick)

(defn hypot
  "Hypot.
  See also [[hypot-sqrt]]."
  (^double [^double x ^double y]
   (FastMath/hypot x y))
  (^double [^double x ^double y ^double z]
   (FastMath/hypot x y z)))

(defn hypot-sqrt
  "Hypot, sqrt version: \\\\(\\sqrt{x^2+y^2}\\\\) or \\\\(\\sqrt{x^2+y^2+z^2}\\\\).
  Should be faster than [[hypot]]."
  (^double [^double x ^double y]
   (FastMath/sqrt (+ (* x x) (* y y))))
  (^double [^double x ^double y ^double z]
   (FastMath/sqrt (+ (* x x) (* y y) (* z z)))))

;; distance
(defn dist
  "Euclidean distance between points `(x1,y1)` and `(x2,y2)`. See [[fastmath.vector]] namespace to see other metrics which work on vectors."
  (^double [[^double x1 ^double y1] [^double x2 ^double y2]] (dist x1 y1 x2 y2))
  (^double [^double x1 ^double y1 ^double x2 ^double y2]
   (FastMath/sqrt (+ (sq (- x2 x1)) (sq (- y2 y1))))))

(defn qdist
  "Quick version of Euclidean distance between points. [[qsqrt]] is used instead of [[sqrt]]."
  (^double [[^double x1 ^double y1] [^double x2 ^double y2]] (qdist x1 y1 x2 y2))
  (^double [^double x1 ^double y1 ^double x2 ^double y2]
   (FastMath/sqrtQuick (+ (sq (- x2 x1)) (sq (- y2 y1))))))

;; Rounding functions
(defn floor
  "\\\\(\\lfloor x \\rfloor\\\\). See: [[qfloor]].

  Rounding is done to a multiply of scale value (when provided)."
  (^double [^double x] (FastMath/floor x))
  (^double [^double x ^double scale] (* (FastMath/floor (/ x scale)) scale)))

(defn ceil
  "\\\\(\\lceil x \\rceil\\\\). See: [[qceil]].

  Rounding is done to a multiply of scale value (when provided)."
  (^double [^double x] (FastMath/ceil x))
  (^double [^double x ^double scale] (* (FastMath/ceil (/ x scale)) scale)))

(defn round "Round to `long`. See: [[rint]], [[qround]]." ^long [^double x] (FastMath/round x))

(defn rint
  "Round to `double`. See [[round]], [[qround]].

  Rounding is done to a multiply of scale value (when provided)."
  (^double [^double x] (FastMath/rint x))
  (^double [^double x ^double scale] (* (FastMath/rint (/ x scale)) scale)))

(defn round-even
  "Round evenly (like in round in R), IEEE / IEC rounding"
  ^long [^double x] (FastMath/roundEven x))

(primitivemath-proxy :one ^{:doc "Fast version of [[floor]]. Returns `long`. See: [[floor]]."} qfloor fastFloor)
(primitivemath-proxy :one ^{:doc "Fast version of [[ceil]]. Returns `long`. See: [[ceil]]."} qceil fastCeil)
(primitivemath-proxy :one ^{:doc "Fast version of [[round]]. Returns `long`. See: [[rint]], [[round]]."} qround fastRound)


(fastmath-proxy :two ^{:doc "From `FastMath` doc: returns dividend - divisor * n,
where n is the mathematical integer closest to dividend/divisor. Returned value in `[-|divisor|/2,|divisor|/2]`"} remainder)

(defn abs "\\\\(|x|\\\\) - `double` version. See [[iabs]]." ^double [^double x] (FastMath/abs x))
(defn iabs "\\\\(|x|\\\\) - `long` version. See [[abs]]." ^long [^long x] (if (neg? x) (- x) x))

(defn trunc
  "Truncate fractional part, keep sign. Returns `double`."
  ^double [^double v] (if (neg? v) (ceil v) (floor v)))

(defn itrunc
  "Truncate fractional part, keep sign. Returns `long`."
  ^long [^double v] (if (neg? v) (qceil v) (qfloor v)))

;; return approximate value
(defn approx
  "Round `v` to specified (default: 2) decimal places. Be aware of `double` number accuracy."
  (^double [^double v] (Precision/round v (int 2)))
  (^double [^double v ^long digits] (Precision/round v (int digits))))

(defn approx-eq
  "Checks equality approximately. See [[approx]]."
  ([^double a ^double b] (== (approx a) (approx b)))
  ([^double a ^double b ^long digits] (== (approx a digits)
                                          (approx b digits))))

(defn delta-eq
  "Checks equality for given absolute accuracy (default `1.0e-6`).

  Version with 4-arity accepts absolute and relative accuracy."
  ([^double a ^double b] (delta-eq a b 1.0e-6))
  ([^double a ^double b ^double accuracy]
   (< (abs (- a b)) accuracy))
  ([^double a ^double b ^double abs-tol ^double rel-tol]
   (< (abs (- a b)) (max abs-tol (* rel-tol (max (abs a) (abs b)))))))

(def ^{:doc "Alias for [[approx-eq]]"} approx= approx-eq)
(def ^{:doc "Alias for [[delta-eq]]"} delta= delta-eq)

(defn near-zero?
  "Checks if given value is near zero with absolute (default: `1.0e-6`) and/or relative (default `0.0`) tolerance."
  ([^double x] (near-zero? x 1.0e-6))
  ([^double x ^double abs-tol] (< (abs x) abs-tol))
  ([^double x ^double abs-tol ^double rel-tol]
   (let [ax (abs x)] (< ax (max abs-tol (* rel-tol ax)) ))))

(defn frac
  "Fractional part, always returns values from 0.0 to 1.0 (exclusive). See [[sfrac]] for signed version."
  ^double [^double v] (abs (- v (unchecked-long v))))

(defn sfrac
  "Fractional part, always returns values from -1.0 to 1.0 (exclusive). See [[frac]] for unsigned version."
  ^double [^double v] (- v (trunc v)))

;; Find power of 2 exponent for double number where  
;; \\(2^(n-1)\leq x\leq 2^n\\)  
;; where n-1 is result of `low-2-exp` and n is result of `high-2-exp`
;; `(low-2-exp TWO_PI) => 2` \\(2^2\eq 4\leq 6.28\\)  
;; `(high-2-exp TWO_PI) => 3` \\(6.28\leq 2^3\eq 8\\)
(defn low-2-exp
  "Find greatest exponent (power of 2) which is lower or equal `x`. See [[high-2-exp]]."
  ^long [^double x] (-> x log2 floor unchecked-long))

(defn high-2-exp
  "Find lowest exponent (power of 2) which is greater or equal `x`. See [[low-2-exp]]."
  ^long [^double v] (-> v log2 ceil unchecked-long))

(defn low-exp
  "Find greatest exponent for base `b` which is lower or equal `x`. See also [[high-exp]]."
  ^long [^double b ^double x] (->> x (logb b) floor unchecked-long))

(defn high-exp
  "Find lowest exponent for base `b` which is higher or equal`x`. See also [[low-exp]]."
  ^long [^double b ^double x] (->> x (logb b) ceil unchecked-long))

(defn round-up-pow2
  "Round long to the next power of 2"
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

(defn double-bits
  "Returns double as 64-bits (long)"
  ^long [^double v]
  (Double/doubleToRawLongBits v))

(defn bits->double
  "Convert 64 bits to double"
  ^double [^long v]
  (Double/longBitsToDouble v))

(defn double-exponent
  "Extract exponent information from double"
  ^long [^double v]
  (FastMath/getExponent v))

(defn double-significand
  "Extract significand from double"
  ^long [^double v]
  (bit-and (Double/doubleToRawLongBits v) 4503599627370495))

(defn log2int
  "Fast and integer version of log2, returns long"
  ^long [^double v]
  (if (< v 1.0)
    (- (log2int (/ v)))
    (let [s (double-significand v)]
      (+ (FastMath/getExponent v)
         (if (or (> s 1865452045155277) ;; (double-significand (pow 2 1.5))
                 (neg? s)) 1 0)))))

(defn ulp
  "Unit in the Last Place, distance between next value larger than `v` and `v`"
  ^double [^double v]
  (FastMath/ulp v))

;; More constants

(def ^{:const true :tag 'double :doc "Value of 0x1.fffffffffffffp-1d = 0.(9)"}
  double-one-minus-epsilon (Double/parseDouble "0x1.fffffffffffffp-1d"))

;; \\(\sqrt{2}\\)
(def ^{:const true :tag 'double :doc "\\\\(\\sqrt{2}\\\\)"} SQRT2 (sqrt 2.0))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{\\sqrt{2}}{2}\\\\)"} SQRT2_2 (* 0.5 SQRT2))

;; \\(\sqrt{3}\\)
(def ^{:const true :tag 'double :doc "\\\\(\\sqrt{3}\\\\)"} SQRT3 (sqrt 3.0))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{\\sqrt{3}}{2}\\\\)"} SQRT3_2 (* 0.5 (sqrt 3.0)))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{\\sqrt{3}}{3}\\\\)"} SQRT3_3 (/ (sqrt 3.0) 3.0))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{\\sqrt{3}}{4}\\\\)"} SQRT3_4 (/ (sqrt 3.0) 4.0))

;; \\(\sqrt{5}\\)
(def ^{:const true :tag 'double :doc "\\\\(\\sqrt{5}\\\\)"} SQRT5 (sqrt 5.0))

;; \\(\sqrt{\pi}\\)
(def ^{:const true :tag 'double :doc "\\\\(\\sqrt{\\pi}\\\\)"} SQRTPI (sqrt PI))
(def ^{:const true :tag 'double :doc "\\\\(\\sqrt{2\\pi}\\\\)"} SQRT2PI (sqrt TWO_PI))
(def ^{:const true :tag 'double :doc "\\\\(\\sqrt{\\frac{1}{2}\\pi}\\\\)"} SQRT_HALFPI (sqrt HALF_PI))

;; 
(def ^{:const true :tag 'double :doc "Golden ratio \\\\(\\phi\\\\)"} PHI (* (inc SQRT5) 0.5))
(def ^{:const true :tag 'double :doc "Silver ratio \\\\(\\delta_S\\\\)"} SILVER (inc SQRT2))

;; math.h predefined constants names
(def ^{:const true :tag 'double :doc "\\\\(e\\\\)"} M_E E)
(def ^{:const true :tag 'double :doc "\\\\(\\log_{2}{e}\\\\)"} M_LOG2E LOG2E)
(def ^{:const true :tag 'double :doc "\\\\(\\log_{10}{e}\\\\)"} M_LOG10E LOG10E)
(def ^{:const true :tag 'double :doc "\\\\(\\ln{2}\\\\)"} M_LN2 LN2)
(def ^{:const true :tag 'double :doc "\\\\(\\ln{10}\\\\)"} M_LN10 LN10)
(def ^{:const true :tag 'double :doc "\\\\(\\pi\\\\)"} M_PI PI)
(def ^{:const true :tag 'double :doc "\\\\(\\frac{\\pi}{2}\\\\)"} M_PI_2 HALF_PI)
(def ^{:const true :tag 'double :doc "\\\\(\\frac{\\pi}{4}\\\\)"} M_PI_4 QUARTER_PI)
(def ^{:const true :tag 'double :doc "\\\\(\\frac{1}{\\pi}\\\\)"} M_1_PI (/ PI))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{2}{\\pi}\\\\)"} M_2_PI (/ 2.0 PI))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{2}{\\sqrt\\pi}\\\\)"} M_2_SQRTPI (/ 2.0 SQRTPI))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{1}{\\sqrt\\pi}\\\\)"} INV_SQRT2PI (/ 1.0 SQRT2PI))
(def ^{:const true :tag 'double :doc "\\\\(\\sqrt{2}\\\\)"} M_SQRT2 SQRT2)
(def ^{:const true :tag 'double :doc "\\\\(\\frac{1}{\\sqrt{2}}\\\\)"} M_SQRT1_2 (/ SQRT2))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{1}{\\sqrt{2}}\\\\)"} INV_SQRT_2 M_SQRT1_2)


(def ^{:const true :tag 'double :doc "\\\\(2\\pi\\\\)"} M_TWOPI TWO_PI)
(def ^{:const true :tag 'double :doc "\\\\(\\frac{3\\pi}{4}\\\\)"} M_3PI_4 (* PI 0.75))
(def ^{:const true :tag 'double :doc "\\\\(\\sqrt\\pi\\\\)"} M_SQRT_PI SQRTPI)
(def ^{:const true :tag 'double :doc "\\\\(\\sqrt{3}\\\\)"} M_SQRT3 SQRT3)
(def ^{:const true :tag 'double :doc "\\\\(\\frac{1}{\\ln{10}}\\\\)"} M_IVLN10 (/ LN10))
(def ^{:const true :tag 'double :doc "\\\\(\\ln{2}\\\\)"} M_LOG2_E LN2)
(def ^{:const true :tag 'double :doc "\\\\(\\frac{1}{\\ln{2}}\\\\)"} M_INVLN2 (/ LN2))

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
  ^double [^double value]
  (if (neg? value) -1.0 1.0))


;; copy-sign

(fastmath-proxy :two ^{:doc "Returns a value with a magnitude of first argument and sign of second."}
                copy-sign copySign)

(defmacro constrain
  "Clamp `value` to the range `[mn,mx]`."
  [value mn mx]
  `(max (min ~value ~mx) ~mn))

(defn norm
  "Normalize `v` from the range `[start,stop]` to the range `[0,1]` or map `v` from the range `[start1,stop1]` to the range `[start2,stop2]`. See also [[make-norm]]."
  {:inline (fn
             ([v start stop] `(PrimitiveMath/norm ~v ~start ~stop))
             ([v start1 stop1 start2 stop2] `(PrimitiveMath/norm ~v ~start1 ~stop1 ~start2 ~stop2)))
   :inline-arities #{3 5}}
  (^double [^double v ^double start ^double stop] ;; norm
   (PrimitiveMath/norm v start stop))
  ([v start1 stop1 start2 stop2] ;; map
   (PrimitiveMath/norm v start1 stop1 start2 stop2)))

(defmacro mnorm
  "Macro version of [[norm]]."
  ([v start stop]
   `(PrimitiveMath/norm ~v ~start ~stop))
  ([v start1 stop1 start2 stop2]
   `(PrimitiveMath/norm ~v ~start1 ~stop1 ~start2 ~stop2)))

(defn make-norm
  "Make [[norm]] function for given range. Resulting function accepts `double` value (with optional target `[dstart,dstop]` range) and returns `double`."
  ([^double start ^double stop]
   (fn ^double [^double v ^double dstart ^double dstop]
     (PrimitiveMath/norm v start stop dstart dstop)))
  ([^double start ^double stop ^double dstart ^double dstop]
   (fn ^double [^double v]
     (PrimitiveMath/norm v start stop dstart dstop))))

(defn cnorm
  "Constrained version of norm. Result of [[norm]] is applied to [[constrain]] to `[0,1]` or `[start2,stop2]` ranges."
  ([v start1 stop1 start2 stop2]
   (constrain ^double (PrimitiveMath/norm v start1 stop1 start2 stop2) ^double start2 ^double stop2))
  (^double [v ^double start ^double stop]
   (constrain ^double (PrimitiveMath/norm v start stop) 0.0 1.0)))

;;; Interpolation functions

;; Linear interpolation between `start` and `stop`.
(defn lerp
  "Linear interpolation between `start` and `stop` for amount `t`. See also [[mlerp]], [[cos-interpolation]], [[quad-interpolation]] or [[smooth-interpolation]]."
  ^double [^double start ^double stop ^double t]
  (+ start (* t (- stop start))))

(defmacro mlerp
  "[[lerp]] as macro. For inline code. See also [[lerp]], [[cos-interpolation]], [[quad-interpolation]] or [[smooth-interpolation]]."
  [start stop t]
  `(+ ~start (* ~t (- ~stop ~start))))

;; Cosine interpolation between `start` and `stop`
(defn cos-interpolation
  "oF interpolateCosine interpolation. See also [[lerp]]/[[mlerp]], [[quad-interpolation]] or [[smooth-interpolation]]."
  ^double [^double start ^double stop ^double t]
  (mlerp start stop (* 0.5 (- 1.0 (cos (* t PI))))))

(defn smooth-interpolation
  "Smoothstep based interpolation. See also [[lerp]]/[[mlerp]], [[quad-interpolation]] or [[cos-interpolation]]."
  ^double [^double start ^double stop ^double t]
  (mlerp start stop (* t t (- 3.0 (* 2.0 t)))))

(defn quad-interpolation
  "Quad interpolation. See also [[lerp]]/[[mlerp]], [[cos-interpolation]] or [[smooth-interpolation]]."
  ^double [^double start ^double stop ^double t]
  (mlerp start stop (let [t' (* 2.0 t)]
                      (if (< t' 1.0)
                        (* 0.5 (* t' t'))
                        (* -0.5 (dec (* (dec t') (- t' 3.0))))))))

(defn smoothstep
  "GL [smoothstep](https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/smoothstep.xhtml)."
  ^double [^double edge0 ^double edge1 ^double x]
  (let [t (cnorm x edge0 edge1)]
    (* t t (- 3.0 (* 2.0 t)))))

;;`(wrap 0 -1 1) => 0.0`  
;;`(wrap -1.1 -1 1) => 0.8999999999999999`  
;;`(wrap 1.1 -1 1) => -0.8999999999999999`
(defn wrap
  "Wrap overflowed value into the range, similar to [ofWrap](http://openframeworks.cc/documentation/math/ofMath/#!show_ofWrap)."
  (^double [[^double start ^double stop] ^double value] (wrap start stop value))
  (^double [^double start ^double stop ^double value]
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
            (- value))))))

;;

(defn- scale-xs
  [xs ^double alpha]
  (map (fn [^double x] (* x alpha)) xs))

(defn- smooth-max-boltzmann
  ^double [xs ^double alpha]
  (let [eaxs (map (fn [^double x] (FastMath/exp (* alpha x))) xs)
        ^double den (reduce fast+ eaxs)]
    (reduce fast+ (map (fn [^double x ^double eax]
                         (/ (* x eax) den)) xs eaxs))))

(defn smooth-max
  "Smooth maximum function.

  A smooth function with `alpha` argument. When `alpha` goes to infinity, function returns maximum value of `xs`.

  Family:

  * `:lse` - LogSumExp (default)
  * `:boltzmann` - Boltzmann operator, works for small alpha values
  * `:mellowmax`
  * `:p-norm`
  * `:smu` - smooth maximum unit, epsilon = 1/alpha > 0"
  (^double [xs] (smooth-max xs 1.0))
  (^double [xs ^double alpha] (smooth-max xs alpha :lse))
  (^double [xs ^double alpha family]
   (case family
     :boltzmann (smooth-max-boltzmann xs alpha)
     :lse (/ (logsumexp (scale-xs xs alpha)) alpha)
     :mellowmax (/ (- (logsumexp (scale-xs xs alpha))
                      (log (count xs))) alpha)
     :p-norm (pow (reduce fast+ (map (fn [^double x]
                                       (pow (abs x) alpha)) xs)) (/ alpha))
     :smu (let [epsilon (/ alpha)]
            (reduce (fn [^double a ^double b]
                      (* 0.5 (+ a b (sqrt (+ (sq (- a b))
                                             epsilon))))) xs)))))

;;

(defn nan?
  "Check if number is NaN"
  {:inline (fn [v] `(Double/isNaN ~v)) :inline-arities #{1}}
  [^double v]
  (Double/isNaN v))

(defn inf?
  "Check if number is infinite"
  {:inline (fn [v] `(Double/isInfinite ~v)) :inline-arities #{1}}
  [^double v]
  (Double/isInfinite v))

(defn pos-inf?
  "Check if number is positively infinite"
  {:inline (fn [v] `(== ~v ##Inf)) :inline-arities #{1}}
  [^double v]
  (== v ##Inf))

(defn neg-inf?
  "Check if number is negatively infinite"
  {:inline (fn [v] `(== ~v ##-Inf)) :inline-arities #{1}}
  [^double v]
  (== v ##-Inf))

(defn invalid-double?
  "Check if number is invalid"
  {:inline (fn [v] `(bool-not (Double/isFinite ~v))) :inline-arities #{1}}
  [^double v]
  (bool-not (Double/isFinite v)))

(defn valid-double?
  "Check if number is invalid"
  {:inline (fn [v] `(Double/isFinite ~v)) :inline-arities #{1}}
  [^double v]
  (Double/isFinite v))

(defn between?
  "Check if given number is within the range [x,y]."
  ([[^double x ^double y] ^double v] (<= x v y))
  ([^double x ^double y ^double v] (<= x v y)))

(defn between-?
  "Check if given number is within the range (x,y]."
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
     (conj r (list (prev-double start) end)))))

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
  ^long [^long a ^long b]
  (gcd- (iabs a) (iabs b)))

(defn lcm
  "Fast binary least common multiplier."
  ^long [^long a ^long b]
  (if (> a b)
    (* b (/ a (gcd- (iabs a) (iabs b))))
    (* a (/ b (gcd- (iabs a) (iabs b))))))

;;

(defn sample
  "Sample function `f` and return sequence of values.

  `range-min` defaults to 0.0, `range-max` to 1.0.

  Range is inclusive.

  When optional `domain?` is set to true (default: false) function returns pairs `[x,(f x)]`."
  ([f number-of-values]
   (sample f 0.0 1.0 number-of-values false))
  ([f ^long number-of-values domain?]
   (sample f 0.0 1.0 number-of-values domain?))
  ([f domain-min domain-max number-of-values]
   (sample f domain-min domain-max number-of-values false))
  ([f domain-min domain-max number-of-values domain?]
   (let [n- (dec ^long number-of-values)
         f (if domain? #(vector % (f %)) f)]
     (->> (range number-of-values)
          (map #(norm % 0.0 n- domain-min domain-max))
          (map f)))))

;; rank/order

(defn rank
  "Sample ranks. See [R docs](https://www.rdocumentation.org/packages/base/versions/3.6.1/topics/rank).

  Rank uses 0 based indexing.
  
  Possible tie strategies: `:average`, `:first`, `:last`, `:random`, `:min`, `:max`, `:dense`.

  `:dense` is the same as in `data.table::frank` from R"
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
                      (fn ^double [v] (/ ^double (reduce fast+ (map first v)) (count v))))
             m (map (fn [[k v]] [k (tie-fn v)]) indexed-sorted-map)
             m (if (= ties :dense)
                 (map-indexed (fn [id [k _]]
                                [k id]) (sort-by second m))
                 m)]
         (map (into {} m) vs))))))

(def rank1 ^{:doc "[[rank]] with indexing statring from 1"}
  (comp (partial map clojure.core/inc) rank))

(defn order
  "Ordering permutation. See [R docs](https://www.rdocumentation.org/packages/base/versions/3.6.1/topics/order)

  Order uses 0 based indexing."
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

(def ^:private ^:const ^{:tag 'double} ci-r0 0.616505485620716233797110404100)
(def ^:private ^:const ^{:tag 'double} ci-r1 3.384180422551186426397851146402)
(def ^:private ^:const ^{:tag 'double} ci-r01 0.6162109375)
(def ^:private ^:const ^{:tag 'double} ci-r02 0.29454812071623379711E-3)
(def ^:private ^:const ^{:tag 'double} ci-r11 3.3837890625)
(def ^:private ^:const ^{:tag 'double} ci-r12 0.39136005118642639785E-3)

(defn Ci
  "Negative of integral of cos(t)/t from x to inf"
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

  Alias for `seq`."} double-array->seq seq)

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
  [res]
  (seq (map seq res)))

(defn seq->double-double-array
  "Convert sequence to double-array of double-arrays.
  
  If sequence is double-array of double-arrays returns `vss`"
  #^"[[D" [vss]
  (cond 
    (= (type vss) double-double-array-type) vss
    (nil? vss) nil
    :else (into-array (map seq->double-array vss))))

;; ## Copy of primitive math machinery
;;
;; Simplified to be used after `ns` is defined.

(def ^:private vars-to-exclude
  '[* + - / > < >= <= == abs rem quot mod bit-and-not bit-set bit-clear bit-test bit-flip bit-or bit-and bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd? bool-and bool-or bool-xor bool-not << >> >>> not==])

(defn- using-primitive-operators? []
  (= #'fastmath.core/+ (resolve '+)))

(defn use-primitive-operators
  "Replaces Clojure's arithmetic and number coercion functions with primitive equivalents.  These are
   defined as macros, so they cannot be used as higher-order functions. This is an idempotent operation. Undo with [[unuse-primitive-operators]]."
  ([] (use-primitive-operators #{}))
  ([skip-set]
   (when-not (using-primitive-operators?)
     (let [v2e (remove skip-set vars-to-exclude)]
       (doseq [v v2e]
         (ns-unmap *ns* v))
       (require ['fastmath.core :refer v2e])))))

(defn unuse-primitive-operators
  "Undoes the work of [[use-primitive-operators]]. This is idempotent."
  []
  (when (using-primitive-operators?)
    (doseq [v vars-to-exclude]
      (ns-unmap *ns* v))
    (refer 'clojure.core)))

