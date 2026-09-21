(ns fastmath.core
  "Core, high-performance mathematical functions and constants, specialized for primitive `double` and `long` types.

  Most functions here are inlined and operate directly on primitives to avoid boxing overhead, and are primarily backed by the FastMath (jafama) library, Apache Commons Math, and custom primitive implementations. Many mathematical constants (`PI`, `E`, roots, reciprocals, logarithms of common values, etc.) are also provided.

  Functions and macros defined in this namespace cover:

  - Primitive-typed arithmetic, comparison, predicate and bitwise operators (`+`, `-`, `*`, `/`, `==`, `<`, `zero?`, `pos?`, `bit-and`, `bit-shift-left`, etc.), including `long`-specific variants and fused multiply-add operations (`fma`, `muladd`).
  - Trigonometric, hyperbolic, and their inverse and reciprocal functions, plus less common variants (versine, haversine, exsecant family).
  - Exponentials, logarithms (including numerically stable variants such as `log1p`, `logsumexp`, `log1pexp`), and power functions (`pow`, `sq`, `cb`, `fpow`, `mpow`, `tpow`).
  - Combinatorics: factorials, falling/rising factorials, binomial coefficients.
  - Rounding, truncation, and floating-point precision helpers: `floor`/`ceil`/`round` variants, approximate equality, ulp stepping (`next-double`, `prev-double`), and raw bit manipulation of doubles.
  - Distance and hypotenuse calculations.
  - Interpolation and range mapping: `lerp`, `norm`, `make-norm`, `wrap`, `smoothstep`, `smooth-max`.
  - Range and interval utilities: `slice-range`, `cut`, `co-intervals`, `group-by-intervals`.
  - Predicates for special double values (`nan?`, `inf?`, `valid-double?`) and range checks (`between?`).
  - Other utilities: `gcd`, `lcm`, `agm` (arithmetic-geometric mean), `sample`, `rank`, `order`, error calculation, and `double`-array conversions.

  Primitive Math Operators:

  A set of inlined macros is provided to replace selected `clojure.core` arithmetic, comparison, and bitwise operators for potential performance gains with primitive arguments. These macros operate on `double` and `long` primitives and generally return primitive values.

  Replaced operators:

  - `* + - / > < >= <= == rem quot mod`
  - `bit-or bit-and bit-xor bit-not bit-and-not`
  - `bit-shift-left bit-shift-right unsigned-bit-shift-right`
  - `bit-set bit-clear bit-flip bit-test`
  - `inc dec`
  - `zero? neg? pos? even? odd?`
  - `min max`
  - `abs`
  - Additionally: `<< >> >>> not==`

  To enable these primitive operators in your namespace, call [[use-primitive-operators]].
  To revert to the original `clojure.core` functions, call [[unuse-primitive-operators]].
  Note that the `fastmath.core` versions are not a complete drop-in replacement due to their primitive-specific behavior (e.g., return types), and calling `unuse-primitive-operators` at the end of the namespace is recommended, especially in Clojure 1.12+."
  (:refer-clojure
   :exclude [* + - / > < >= <= == rem quot mod bit-or bit-and bit-and-not bit-set bit-clear bit-test bit-flip bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd? abs integer?])
  (:import [fastmath.java PrimitiveMath]
           [net.jafama FastMath]
           [org.apache.commons.math3.util Precision]
           [org.apache.commons.math3.special Gamma Beta]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; which java?

(def ^{:const true :tag 'long} jvm-version
  (->> (System/getProperty "java.version")
       (re-seq #"\d+")
       (first)
       (Long/parseLong)))

;;

(defn ^:private primitivemath-nary-inline
  ([op] (primitivemath-nary-inline op nil nil))
  ([op one] (primitivemath-nary-inline op one nil))
  ([op one v]
   (fn
     ([] `~v)
     ([x] (if one `(. PrimitiveMath (~one ~x)) `~x))
     ([x y] `(. PrimitiveMath (~op ~x ~y)))
     ([x y & more]
      (reduce
       (fn [a b] `(. PrimitiveMath (~op ~a ~b)))
       `(. PrimitiveMath (~op ~x ~y)) more)))))

(defn ^:private primitivemath-nary-inline-long
  ([op] (primitivemath-nary-inline-long op nil nil))
  ([op one] (primitivemath-nary-inline-long op one nil))
  ([op one v]
   (fn
     ([] `~v)
     ([x] (if one `(. PrimitiveMath (~one (long ~x))) `(long ~x)))
     ([x y] `(. PrimitiveMath (~op (long ~x) (long ~y))))
     ([x y & more]
      (reduce
       (fn [a b] `(. PrimitiveMath (~op (long ~a) (long ~b))))
       `(. PrimitiveMath (~op (long ~x) (long ~y))) more)))))

(defn ^:private >=2? [^long n] (fastmath.java.PrimitiveMath/gte n 2))
(defn ^:private >=1? [^long n] (fastmath.java.PrimitiveMath/gte n 1))
(defn ^:private >=0? [^long n] (fastmath.java.PrimitiveMath/gte n 0))

;; ## Basic operations

(defn +
  {:inline (primitivemath-nary-inline 'add nil 0.0)
   :inline-arities >=0?
   :doc "Adds numbers together. Primitive and inlined replacement for `clojure.core/+`.

  Parameters:

  - zero or more `double` values to sum.

  Returns the sum as a double. With no arguments returns `0.0` (the additive identity). With one argument returns it unchanged.

  See also [[use-primitive-operators]] (enables this replacement), [[long-add]] (long-coerced version), [[-]], [[*]], [[/]]."}
  (^double [] 0.0)
  (^double [^double a] a)
  (^double [^double a ^double b] (. PrimitiveMath (add a b)))
  (^double [^double a ^double b ^double c] (. PrimitiveMath (add (. PrimitiveMath (add a b)) c)))
  (^double [^double a ^double b ^double c ^double d]
   (. PrimitiveMath (add (. PrimitiveMath (add (. PrimitiveMath (add a b)) c)) d)))
  ([a b c d & r] (reduce + (+ (double a) (double b) (double c) (double d)) r)))

(defn long-add
  {:inline (primitivemath-nary-inline-long 'add nil 0)
   :inline-arities >=0?
   :doc "Adds numbers together, coercing arguments and the result to `long`. Primitive and inlined replacement for `clojure.core/+`.

  Parameters:

  - zero or more values coercible to `long`.

  Returns the sum as a long, following standard JVM two's complement overflow wraparound (no overflow checking). With no arguments returns `0`. With one argument returns it unchanged.

  See also [[+]] (double version), [[long-sub]], [[long-mult]], [[long-div]]."}
  (^long [] 0)
  (^long [^long a] a)
  (^long [^long a ^long b] (. PrimitiveMath (add a b)))
  (^long [^long a ^long b ^long c] (. PrimitiveMath (add (. PrimitiveMath (add a b)) c)))
  (^long [^long a ^long b ^long c ^long d]
   (. PrimitiveMath (add (. PrimitiveMath (add (. PrimitiveMath (add a b)) c)) d)))
  ([a b c d & r] (reduce long-add (long-add a b c d) r)))

(defn -
  {:inline (primitivemath-nary-inline 'subtract 'negate)
   :inline-arities >=1?
   :doc "Subtracts numbers, or negates a single value. Primitive and inlined replacement for `clojure.core/-`.

  Parameters:

  - one or more `double` values.

  Returns the negation of the sole argument when called with one value, otherwise the left-to-right cumulative subtraction of all arguments, as a double.

  See also [[long-sub]] (long-coerced version), [[+]], [[*]], [[/]]."}
  (^double [^double a] (. PrimitiveMath (negate a)))
  (^double [^double a ^double b] (. PrimitiveMath (subtract a b)))
  (^double [^double a ^double b ^double c] (. PrimitiveMath (subtract (. PrimitiveMath (subtract a b)) c)))
  (^double [^double a ^double b ^double c ^double d]
   (. PrimitiveMath (subtract (. PrimitiveMath (subtract (. PrimitiveMath (subtract a b)) c)) d)))
  ([a b c d & r] (reduce - (- (double a) (double b) (double c) (double d)) r)))

(defn long-sub
  {:inline (primitivemath-nary-inline-long 'subtract 'negate)
   :inline-arities >=1?
   :doc "Subtracts numbers, or negates a single value, coercing arguments and the result to `long`. Primitive and inlined replacement for `clojure.core/-`.

  Parameters:

  - one or more values coercible to `long`.

  Returns the negation of the sole argument when called with one value, otherwise the left-to-right cumulative subtraction of all arguments, as a long, following standard JVM two's complement overflow wraparound.

  See also [[-]] (double version), [[long-add]], [[long-mult]], [[long-div]]."}
  (^long [^long a] (. PrimitiveMath (negate a)))
  (^long [^long a ^long b] (. PrimitiveMath (subtract a b)))
  (^long [^long a ^long b ^long c] (. PrimitiveMath (subtract (. PrimitiveMath (subtract a b)) c)))
  (^long [^long a ^long b ^long c ^long d]
   (. PrimitiveMath (subtract (. PrimitiveMath (subtract (. PrimitiveMath (subtract a b)) c)) d)))
  ([a b c d & r] (reduce long-sub (long-sub a b c d) r)))

(defn *
  {:inline (primitivemath-nary-inline 'multiply nil 1.0)
   :inline-arities >=0?
   :doc "Multiplies numbers together. Primitive and inlined replacement for `clojure.core/*`.

  Parameters:

  - zero or more `double` values to multiply.

  Returns the product as a double. With no arguments returns `1.0` (the multiplicative identity). With one argument returns it unchanged.

  See also [[long-mult]] (long-coerced version), [[+]], [[-]], [[/]]."}
  (^double [] 1.0)
  (^double [^double a] a)
  (^double [^double a ^double b] (. PrimitiveMath (multiply a b)))
  (^double [^double a ^double b ^double c] (. PrimitiveMath (multiply (. PrimitiveMath (multiply a b)) c)))
  (^double [^double a ^double b ^double c ^double d]
   (. PrimitiveMath (multiply (. PrimitiveMath (multiply (. PrimitiveMath (multiply a b)) c)) d)))
  ([a b c d & r] (reduce * (* (double a) (double b) (double c) (double d)) r)))

(defn long-mult
  {:inline (primitivemath-nary-inline-long 'multiply nil 1)
   :inline-arities >=0?
   :doc "Multiplies numbers together, coercing arguments and the result to `long`. Primitive and inlined replacement for `clojure.core/*`.

  Parameters:

  - zero or more values coercible to `long`.

  Returns the product as a long, following standard JVM two's complement overflow wraparound. With no arguments returns `1`. With one argument returns it unchanged.

  See also [[*]] (double version), [[long-add]], [[long-sub]], [[long-div]]."}
  (^long [] 1)
  (^long [^long a] a)
  (^long [^long a ^long b] (. PrimitiveMath (multiply a b)))
  (^long [^long a ^long b ^long c] (. PrimitiveMath (multiply (. PrimitiveMath (multiply a b)) c)))
  (^long [^long a ^long b ^long c ^long d]
   (. PrimitiveMath (multiply (. PrimitiveMath (multiply (. PrimitiveMath (multiply a b)) c)) d)))
  ([a b c d & r] (reduce long-mult (long-mult a b c d) r)))

(defn /
  {:inline (primitivemath-nary-inline 'divide 'reciprocal)
   :inline-arities >=1?
   :doc "Divides numbers, or returns the reciprocal of a single value. Primitive and inlined replacement for `clojure.core//`.

  Parameters:

  - one or more `double` values.

  Returns the reciprocal of the sole argument when called with one value, otherwise the left-to-right cumulative division of all arguments, as a double. Division by zero follows IEEE 754 semantics (returns `##Inf`, `##-Inf` or `##NaN`, no exception thrown).

  See also [[long-div]] (long-coerced version), [[+]], [[-]], [[*]]."}
  (^double [^double a] (. PrimitiveMath (reciprocal a)))
  (^double [^double a ^double b] (. PrimitiveMath (divide a b)))
  (^double [^double a ^double b ^double c] (. PrimitiveMath (divide (. PrimitiveMath (divide a b)) c)))
  (^double [^double a ^double b ^double c ^double d]
   (. PrimitiveMath (divide (. PrimitiveMath (divide (. PrimitiveMath (divide a b)) c)) d)))
  ([a b c d & r] (reduce / (/ (double a) (double b) (double c) (double d)) r)))

(defn long-div
  {:inline (primitivemath-nary-inline-long 'divide 'reciprocal)
   :inline-arities >=1?
   :doc "Divides numbers, or returns the reciprocal of a single value, coercing arguments to `long`. Primitive and inlined replacement for `clojure.core//`.

  Parameters:

  - one or more values coercible to `long`.

  Returns the reciprocal of the sole argument as a double when called with one value (a long's reciprocal is generally not representable as a long, so this arity does not coerce its result). With two or more arguments, returns the left-to-right cumulative integer division of all arguments as a long. Throws an arithmetic exception on division by zero when two or more arguments are given.

  See also [[/]] (double version), [[long-add]], [[long-sub]], [[long-mult]]."}
  (^double [^long a] (. PrimitiveMath (reciprocal a)))
  (^long [^long a ^long b] (. PrimitiveMath (divide a b)))
  (^long [^long a ^long b ^long c] (. PrimitiveMath (divide (. PrimitiveMath (divide a b)) c)))
  (^long [^long a ^long b ^long c ^long d]
   (. PrimitiveMath (divide (. PrimitiveMath (divide (. PrimitiveMath (divide a b)) c)) d)))
  ([a b c d & r] (reduce long-div (long-div a b c d) r)))

(defn inc
  {:inline (fn [x] `(. PrimitiveMath (inc ~x)))
   :inline-arities #{1}
   :doc "Primitive and inlined `inc`"}
  ^double [^double x] (. PrimitiveMath (inc x)))

(defn long-inc
  {:inline (fn [x] `(. PrimitiveMath (inc (long ~x))))
   :inline-arities #{1}
   :doc "Primitive and inlined `inc` coerced to a long"}
  ^long [^long x] (. PrimitiveMath (inc x)))

(defn dec
  {:inline (fn [x] `(. PrimitiveMath (dec ~x)))
   :inline-arities #{1}
   :doc "Primitive and inlined `dec`"}
  ^double [^double x] (. PrimitiveMath (dec x)))

(defn long-dec
  {:inline (fn [x] `(. PrimitiveMath (dec (long ~x))))
   :inline-arities #{1}
   :doc "Primitive and inlined `dec` coerced to a long"}
  ^long [^long x] (. PrimitiveMath (dec x)))

(defn rem
  {:inline (fn [x y] `(. PrimitiveMath (remainder ~x ~y)))
   :inline-arities #{2}
   :doc "Primitive and inlined `rem`"}
  ^double [^double x ^double y] (. PrimitiveMath (remainder x y)))

(defn long-rem
  {:inline (fn [x y] `(. PrimitiveMath (remainder (long ~x) (long ~y))))
   :inline-arities #{2}
   :doc "Primitive and inlined `rem` coerced to longs"}
  ^long [^long x ^long y] (. PrimitiveMath (remainder x y)))

(defn quot
  {:inline (fn [x y] `(. PrimitiveMath (quotient ~x ~y)))
   :inline-arities #{2}
   :doc "Primitive and inlined `quot`"}
  ^double [^double x ^double y] (. PrimitiveMath (quotient x y)))

(defn long-quot
  {:inline (fn [x y] `(. PrimitiveMath (quotient (long ~x) (long ~y))))
   :inline-arities #{2}
   :doc "Primitive and inlined `quot` coerced to longs"}
  ^long [^long x ^long y] (. PrimitiveMath (quotient x y)))

(defn mod
  {:inline (fn [x y] `(. PrimitiveMath (modulus ~x ~y)))
   :inline-arities #{2}
   :doc "Primitive and inlined `mod`"}
  ^double [^double x ^double y] (. PrimitiveMath (modulus x y)))

(defn long-mod
  {:inline (fn [x y] `(. PrimitiveMath (modulus (long ~x) (long ~y))))
   :inline-arities #{2}
   :doc "Primitive and inlined `mod` coerced to longs"}
  ^long [^long x ^long y] (. PrimitiveMath (modulus x y)))

(defn min
  {:inline (primitivemath-nary-inline 'min)
   :inline-arities >=1?
   :doc "Returns the smallest of one or more values. Primitive and inlined replacement for `clojure.core/min`.

  Parameters:

  - one or more `double` values.

  Returns the minimum value as a double.

  See also [[long-min]] (long-coerced version), [[max]]."}
  (^double [^double a] a)
  (^double [^double a ^double b] (. PrimitiveMath (min a b)))
  (^double [^double a ^double b ^double c] (. PrimitiveMath (min (. PrimitiveMath (min a b)) c)))
  (^double [^double a ^double b ^double c ^double d]
   (. PrimitiveMath (min (. PrimitiveMath (min (. PrimitiveMath (min a b)) c)) d)))
  ([a b c d & r] (reduce min (min (double a) (double b) (double c) (double d)) r)))

(defn long-min
  {:inline (primitivemath-nary-inline-long 'min)
   :inline-arities >=1?
   :doc "Returns the smallest of one or more values, coercing arguments and the result to `long`. Primitive and inlined replacement for `clojure.core/min`.

  Parameters:

  - one or more values coercible to `long`.

  Returns the minimum value as a long.

  See also [[min]] (double version), [[long-max]]."}
  (^long [^long a] a)
  (^long [^long a ^long b] (. PrimitiveMath (min a b)))
  (^long [^long a ^long b ^long c] (. PrimitiveMath (min (. PrimitiveMath (min a b)) c)))
  (^long [^long a ^long b ^long c ^long d]
   (. PrimitiveMath (min (. PrimitiveMath (min (. PrimitiveMath (min a b)) c)) d)))
  ([a b c d & r] (reduce long-min (long-min a b c d) r)))

(defn max
  {:inline (primitivemath-nary-inline 'max)
   :inline-arities >=1?
   :doc "Returns the largest of one or more values. Primitive and inlined replacement for `clojure.core/max`.

  Parameters:

  - one or more `double` values.

  Returns the maximum value as a double.

  See also [[long-max]] (long-coerced version), [[min]]."}
  (^double [^double a] a)
  (^double [^double a ^double b] (. PrimitiveMath (max a b)))
  (^double [^double a ^double b ^double c] (. PrimitiveMath (max (. PrimitiveMath (max a b)) c)))
  (^double [^double a ^double b ^double c ^double d]
   (. PrimitiveMath (max (. PrimitiveMath (max (. PrimitiveMath (max a b)) c)) d)))
  ([a b c d & r] (reduce max (max (double a) (double b) (double c) (double d)) r)))

(defn long-max
  {:inline (primitivemath-nary-inline-long 'max)
   :inline-arities >=1?
   :doc "Returns the largest of one or more values, coercing arguments and the result to `long`. Primitive and inlined replacement for `clojure.core/max`.

  Parameters:

  - one or more values coercible to `long`.

  Returns the maximum value as a long.

  See also [[max]] (double version), [[long-min]]."}
  (^long [^long a] a)
  (^long [^long a ^long b] (. PrimitiveMath (max a b)))
  (^long [^long a ^long b ^long c] (. PrimitiveMath (max (. PrimitiveMath (max a b)) c)))
  (^long [^long a ^long b ^long c ^long d]
   (. PrimitiveMath (max (. PrimitiveMath (max (. PrimitiveMath (max a b)) c)) d)))
  ([a b c d & r] (reduce long-max (long-max a b c d) r)))

(defn zero?
  {:inline (fn [x] `(. PrimitiveMath (isZero ~x)))
   :inline-arities #{1}
   :doc "Primitive and inlined `zero?`"}
  [^double x] (. PrimitiveMath (isZero x)))

(defn not-zero?
  {:inline (fn [x] `(not (. PrimitiveMath (isZero ~x))))
   :inline-arities #{1}
   :doc "Primitive and inlined x<>0.0"}
  [^double x] (not (. PrimitiveMath (isZero x))))

(defn one?
  {:inline (fn [x] `(. PrimitiveMath (isOne ~x)))
   :inline-arities #{1}
   :doc "Primitive and inlined `one?` (x==1.0)"}
  [^double x] (. PrimitiveMath (isOne x)))

(defn neg?
  {:inline (fn [x] `(. PrimitiveMath (isNeg ~x)))
   :inline-arities #{1}
   :doc "Primitive and inlined `neg?`"}
  [^double x] (. PrimitiveMath (isNeg x)))

(defn pos?
  {:inline (fn [x] `(. PrimitiveMath (isPos ~x)))
   :inline-arities #{1}
   :doc "Primitive and inlined `pos?`"}
  [^double x] (. PrimitiveMath (isPos x)))

(defn not-neg?
  {:inline (fn [x] `(. PrimitiveMath (isNNeg ~x)))
   :inline-arities #{1}
   :doc "Primitive and inlined `not-neg?` (x>=0.0)"}
  [^double x] (. PrimitiveMath (isNNeg x)))

(defn not-pos?
  {:inline (fn [x] `(. PrimitiveMath (isNPos ~x)))
   :inline-arities #{1}
   :doc "Primitive and inlined `not-pos?` (x<=0.0)"}
  [^double x] (. PrimitiveMath (isNPos x)))

(defn even?
  {:inline (fn [x] `(. PrimitiveMath (isEven ~x)))
   :inline-arities #{1}
   :doc "Primitive and inlined `even?`"}
  [^long x] (. PrimitiveMath (isEven x)))

(defn odd?
  {:inline (fn [x] `(. PrimitiveMath (isOdd ~x)))
   :inline-arities #{1}
   :doc "Primitive and inlined `odd?`"}
  [^long x] (. PrimitiveMath (isOdd x)))

;;

(defn- primitivemath-nary-inline-predicate
  [op]
  (fn ([_] true)
    ([a b] `(. PrimitiveMath (~op ~a ~b)))
    ([a b & r] `(and (. PrimitiveMath (~op ~a ~b))
                     ~@(map (fn [[x y]] `(. PrimitiveMath (~op ~x ~y)))
                            (partition 2 1 (conj r b)))))))

(defn ==
  "Primitive math equality function.

  Parameters:

  - one or more `double` values.

  Returns `true` if all arguments are equal (chained pairwise: `a=b`, `b=c`, `c=d`, ...), `false` otherwise, as a Boolean. With one argument returns `true`.

  See also [[eq]] (alias), [[not==]], [[approx-eq]], [[delta-eq]]."
  {:inline (primitivemath-nary-inline-predicate 'eq)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (eq a b)))
  ([a b & r]
   (boolean (and (== (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (eq x y)) (reduced false) y)) b r)))))
(defn eq
  "Primitive math equality function. Alias for [[==]].

  Parameters:

  - one or more `double` values.

  Returns `true` if all arguments are equal (chained pairwise: `a=b`, `b=c`, `c=d`, ...), `false` otherwise, as a Boolean. With one argument returns `true`.

  See also [[==]], [[not==]], [[approx-eq]], [[delta-eq]]."
  {:inline (primitivemath-nary-inline-predicate 'eq)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (eq  a b)))
  ([a b & r]
   (boolean (and (eq (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (eq x y)) (reduced false) y)) b r)))))

(defn <
  "Primitive math less-than function.

  Parameters:

  - one or more `double` values.

  Returns `true` if arguments are in strictly increasing order (chained pairwise: `a<b`, `b<c`, `c<d`, ...), `false` otherwise, as a Boolean. With one argument returns `true`.

  See also [[>]], [[<=]], [[>=]]."
  {:inline (primitivemath-nary-inline-predicate 'lt)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (lt a b)))
  ([a b & r]
   (boolean (and (< (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (lt x y)) (reduced false) y)) b r)))))

(defn >
  "Primitive math greater-than function.

  Parameters:

  - one or more `double` values.

  Returns `true` if arguments are in strictly decreasing order (chained pairwise: `a>b`, `b>c`, `c>d`, ...), `false` otherwise, as a Boolean. With one argument returns `true`.

  See also [[<]], [[<=]], [[>=]]."
  {:inline (primitivemath-nary-inline-predicate 'gt)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (gt a b)))
  ([a b & r]
   (boolean (and (> (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (gt x y)) (reduced false) y)) b r)))))

(defn <=
  "Primitive math less-than-or-equal function.

  Parameters:

  - one or more `double` values.

  Returns `true` if arguments are in non-decreasing order (chained pairwise: `a<=b`, `b<=c`, `c<=d`, ...), `false` otherwise, as a Boolean. With one argument returns `true`.

  See also [[<]], [[>]], [[>=]]."
  {:inline (primitivemath-nary-inline-predicate 'lte)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (lte a b)))
  ([a b & r]
   (boolean (and (<= (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (lte x y)) (reduced false) y)) b r)))))

(defn >=
  "Primitive math greater-than-or-equal function.

  Parameters:

  - one or more `double` values.

  Returns `true` if arguments are in non-increasing order (chained pairwise: `a>=b`, `b>=c`, `c>=d`, ...), `false` otherwise, as a Boolean. With one argument returns `true`.

  See also [[<]], [[>]], [[<=]]."
  {:inline (primitivemath-nary-inline-predicate 'gte)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (gte a b)))
  ([a b & r]
   (boolean (and (>= (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (gte x y)) (reduced false) y)) b r)))))

(defn not==
  "Not equality. For more than two arguments, returns `true` when all values are unique.

  `(not== 1 2 1)` === `(and (not= 1 1) (not= 1 2))`

  Parameters:

  - one or more `double` values.

  Returns `true` when all arguments are pairwise distinct, `false` when any two arguments are equal, as a Boolean. With one argument returns `false`.

  See also [[==]], [[eq]]."
  {:inline (fn ([_] false)
             ([a b] `(. PrimitiveMath (neq ~a ~b)))
             ([a b & r] `(and ~@(map (fn [[x y]] `(. PrimitiveMath (neq ~x ~y)))
                                     (partition 2 1 (sort (conj r a b)))))))
   :inline-arities >=1?}
  ([_] false)
  ([^double a ^double b] (. PrimitiveMath (neq a b)))
  ([a b & r]
   (boolean (reduce (fn [^double x ^double y]
                      (if-not (. PrimitiveMath (neq x y)) (reduced false) y)) (sort (conj r a b))))))

;;;;;;;;;;;;;;

(defn bit-and
  "Bitwise AND (`x ∧ y`) of two or more integer values.

  Parameters:

  - one or more values coercible to `long`.

  Returns the bitwise AND of all arguments as a long, computed by folding pairwise from left to right. With one argument returns it unchanged.

  See also [[bit-or]], [[bit-xor]], [[bit-and-not]], [[bit-nand]]."
  {:inline (primitivemath-nary-inline-long 'bitAnd)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitAnd x y)))
  ([x y & r] (reduce bit-and (. PrimitiveMath (bitAnd x y)) r)))

(defn bit-nand
  "Bitwise NAND (`~(x ∧ y)`) of two values, or a left-to-right pairwise fold of NAND over more than two.

  Parameters:

  - one or more values coercible to `long`.

  Returns the bitwise NAND as a long. For three or more arguments, computes `(x nand y) nand z ...` -- a repeated pairwise fold, not the hardware-style `~(x ∧ y ∧ z ...)` multi-input NAND. With one argument returns it unchanged.

  See also [[bit-and]], [[bit-nor]], [[bit-xnor]]."
  {:inline (primitivemath-nary-inline-long 'bitNand)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitNand x y)))
  ([x y & r] (reduce bit-nand (. PrimitiveMath (bitNand x y)) r)))

(defn bit-and-not
  "Bitwise AND with complemented subsequent arguments (`x ∧ ~y`), folded left to right for more than two arguments.

  Parameters:

  - one or more values coercible to `long`.

  Returns `x` with all bits set in any of `y, z, ...` cleared, as a long, computed by folding pairwise from left to right (`(x and-not y) and-not z ...`). With one argument returns it unchanged.

  See also [[bit-and]], [[bit-clear]]."
  {:inline (primitivemath-nary-inline-long 'bitAndNot)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitAndNot x y)))
  ([x y & r] (reduce bit-and-not (. PrimitiveMath (bitAndNot x y)) r)))

(defn bit-or
  "Bitwise OR (`x ∨ y`) of two or more integer values.

  Parameters:

  - one or more values coercible to `long`.

  Returns the bitwise OR of all arguments as a long, computed by folding pairwise from left to right. With one argument returns it unchanged.

  See also [[bit-and]], [[bit-xor]], [[bit-nor]]."
  {:inline (primitivemath-nary-inline-long 'bitOr)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitOr x y)))
  ([x y & r] (reduce bit-or (. PrimitiveMath (bitOr x y)) r)))

(defn bit-nor
  "Bitwise NOR (`~(x ∨ y)`) of two values, or a left-to-right pairwise fold of NOR over more than two.

  Parameters:

  - one or more values coercible to `long`.

  Returns the bitwise NOR as a long. For three or more arguments, computes `(x nor y) nor z ...` -- a repeated pairwise fold, not the hardware-style `~(x ∨ y ∨ z ...)` multi-input NOR. With one argument returns it unchanged.

  See also [[bit-or]], [[bit-nand]], [[bit-xnor]]."
  {:inline (primitivemath-nary-inline-long 'bitNor)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitNor x y)))
  ([x y & r] (reduce bit-nor (. PrimitiveMath (bitNor x y)) r)))

(defn bit-xor
  "Bitwise XOR (`x⊕y`) of two or more integer values.

  Parameters:

  - one or more values coercible to `long`.

  Returns the bitwise XOR of all arguments as a long, computed by folding pairwise from left to right (associative, so fold order does not affect the result). With one argument returns it unchanged.

  See also [[bit-xnor]], [[bit-and]], [[bit-or]]."
  {:inline (primitivemath-nary-inline-long 'bitXor)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitXor x y)))
  ([x y & r] (reduce bit-xor (. PrimitiveMath (bitXor x y)) r)))

(defn bit-xnor
  "Bitwise XNOR (`~(x⊕y)`) of two values, or a left-to-right pairwise fold of XNOR over more than two.

  Parameters:

  - one or more values coercible to `long`.

  Returns the bitwise XNOR as a long. For three or more arguments, computes `(x xnor y) xnor z ...` -- a repeated pairwise fold, not the hardware-style `~(x⊕y⊕z...)` multi-input XNOR. With one argument returns it unchanged.

  See also [[bit-xor]], [[bit-nand]], [[bit-nor]]."
  {:inline (primitivemath-nary-inline-long 'bitXNor)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitXNor x y)))
  ([x y & r] (reduce bit-xnor (. PrimitiveMath (bitXNor x y)) r)))

(defn bit-not
  "~x - bitwise NOT"
  {:inline (fn [x] `(. PrimitiveMath (bitNot (long ~x))))
   :inline-arities #{1}}
  ^long [^long x] (. PrimitiveMath (bitNot x)))

(defn bit-set
  "Set bit (set to `1`)."
  {:inline (fn [x bit] `(. PrimitiveMath (bitSet (long ~x) (long ~bit))))
   :inline-arities #{2}}
  ^long [^long x ^long bit] (. PrimitiveMath (bitSet x bit)))

(defn bit-clear
  "Clear bit (set to `0`)."
  {:inline (fn [x bit] `(. PrimitiveMath (bitClear (long ~x) (long ~bit))))
   :inline-arities #{2}}
  ^long [^long x ^long bit] (. PrimitiveMath (bitClear x bit)))

(defn bit-flip
  "Flip bit (set to `0` when `1` or to `1` when `0`)."
  {:inline (fn [x bit] `(. PrimitiveMath (bitFlip (long ~x) (long ~bit))))
   :inline-arities #{2}}
  ^long [^long x ^long bit] (. PrimitiveMath (bitFlip x bit)))

(defn bit-test
  "Test bit (return to `true` when `1` or `false` when `0`)."
  {:inline (fn [x bit] `(. PrimitiveMath (bitTest (long ~x) (long ~bit))))
   :inline-arities #{2}}
  [^long x ^long bit] (. PrimitiveMath (bitTest x bit)))

(defn bit-count
  "Count set bits"
  {:inline (fn [x] `(Long/bitCount (long ~x)))
   :inline-arities #{1}}
  ^long [^long x] (Long/bitCount x))

(defn bit-shift-left
  "Shift bits left"
  {:inline (fn [x shift] `(. PrimitiveMath (shiftLeft (long ~x) (long ~shift))))
   :inline-arities #{2}}
  [^long x ^long shift] (. PrimitiveMath (shiftLeft x shift)))

(defn <<
  "Shift bits left"
  {:inline (fn [x shift] `(. PrimitiveMath (shiftLeft (long ~x) (long ~shift))))
   :inline-arities #{2}}
  [^long x ^long shift] (. PrimitiveMath (shiftLeft x shift)))

(defn bit-shift-right
  "Shift bits right and keep most significant bit unchanged"
  {:inline (fn [x shift] `(. PrimitiveMath (shiftRight (long ~x) (long ~shift))))
   :inline-arities #{2}}
  [^long x ^long shift] (. PrimitiveMath (shiftRight x shift)))

(defn >>
  "Shift bits right and keep most significant bit unchanged"
  {:inline (fn [x shift] `(. PrimitiveMath (shiftRight (long ~x) (long ~shift))))
   :inline-arities #{2}}
  [^long x ^long shift] (. PrimitiveMath (shiftRight x shift)))

(defn unsigned-bit-shift-right
  "Shift bits right and set most significant bit to `0`"
  {:inline (fn [x shift] `(. PrimitiveMath (unsignedShiftRight (long ~x) (long ~shift))))
   :inline-arities #{2}}
  [^long x ^long shift] (. PrimitiveMath (unsignedShiftRight x shift)))

(defn >>>
  "Shift bits right and set most significant bit to `0`"
  {:inline (fn [x shift] `(. PrimitiveMath (unsignedShiftRight (long ~x) (long ~shift))))
   :inline-arities #{2}}
  [^long x ^long shift] (. PrimitiveMath (unsignedShiftRight x shift)))

(defn bool-not
  "Primitive boolean not"
  {:inline (fn [x] `(. PrimitiveMath (not ~x)))
   :inline-arities #{1}}
  [x] (. PrimitiveMath (not (boolean x))))

(defn bool-xor
  "Primitive boolean XOR of two or more values.

  Parameters:

  - two or more values (coerced to boolean).

  Returns the exclusive-or of all arguments as a Boolean, computed by folding pairwise from left to right (associative, so fold order does not affect the result).

  See also [[xor]] (identical implementation), [[bool-not]]."
  {:inline (primitivemath-nary-inline 'xor)
   :inline-arities >=2?}
  ([x y] (. PrimitiveMath (xor (boolean x) (boolean y))))
  ([x y & r] (reduce bool-xor (. PrimitiveMath (xor (boolean x) (boolean y))) r)))

(defn xor
  "Primitive boolean XOR of two or more values. Identical implementation to [[bool-xor]].

  Parameters:

  - two or more values (coerced to boolean).

  Returns the exclusive-or of all arguments as a Boolean, computed by folding pairwise from left to right (associative, so fold order does not affect the result).

  See also [[bool-xor]], [[bool-not]]."
  {:inline (primitivemath-nary-inline 'xor)
   :inline-arities >=2?}
  ([x y] (. PrimitiveMath (xor (boolean x) (boolean y))))
  ([x y & r] (reduce xor (. PrimitiveMath (xor (boolean x) (boolean y))) r)))

;;;;

(defn negative-zero?
  "Checks whether a double is negative zero (`-0.0`), as distinct from positive zero (`0.0`).

  Parameters:

  - `x` (double): value to check.

  Returns `true` when `x` is bitwise equal to `-0.0`, `false` otherwise (including for `0.0`, despite `(== -0.0 0.0)` being `true` under normal floating-point equality).

  See also [[zero?]]."
  {:inline (fn [x] `(. PrimitiveMath (eq -9223372036854775808
                                        (Double/doubleToLongBits (double ~x)))))
   :inline-arities #{1}}
  [^double x]
  (== (Double/doubleToLongBits x) -9223372036854775808)) ;; -0.0

(defn identity-double
  {:inline (fn [x] `~x) :inline-arities #{1}
   :doc "Identity on double."}
  ^double [^double a] a)

(defn identity-long
  {:inline (fn [x] `(long ~x)) :inline-arities #{1}
   :doc "Identity on long."}
  ^long [^long a] a)

(defn integer?
  "Checks if a given real number is a mathematical integer (has zero fractional part).

  Parameters:

  - `v` (double): value to check.

  Returns `true` when `v` equals its own rounded value (`v == rint(v)`), `false` otherwise. Shadows `clojure.core/integer?`, which instead checks the value's Java type; this predicate checks the numeric value, so `(integer? 5.0)` returns `true`.

  See also [[frac]], [[sfrac]]."
  {:inline (fn [v] `(== ~v (FastMath/rint (double ~v)))) :inline-arities #{1}}
  [^double v]
  (== v (FastMath/rint v)))

;; macros for polynomials

(defn ^:private ->fma
  ([] (->fma false))
  ([n?]
   (if-not n?
     (if (< jvm-version 9)
       (fn [x y z] `(+ ~z (* ~x ~y)))
       (fn [x y z] `(Math/fma (double ~x) (double ~y) (double ~z))))
     (if (< jvm-version 9)
       (fn [x y z] `(+ ~z (* (- ~x) ~y)))
       (fn [x y z] `(Math/fma (- (double ~x)) (double ~y) (double ~z)))))))

(defmacro ^:private fma-macro
  [x y z]
  (if (< jvm-version 9)
    `(+ ~z (* ~x ~y))
    `(Math/fma ~x ~y ~z)))

(defn muladd
  "Computes `x*y + z` using fused multiply-add when available.

  Parameters:

  - `x`, `y`, `z` (doubles): multiplicands and addend.

  Returns `x*y + z` as a double. On Java 9+, uses `Math/fma` for a single correctly-rounded fused multiply-add (more accurate than a separate multiply and add). On earlier JVMs, falls back to a plain `(+ z (* x y))` with intermediate rounding.

  See also [[fma]] (identical implementation), [[negmuladd]], [[difference-of-products]], [[sum-of-products]]."
  {:inline (->fma)
   :inline-arities #{3}}
  ^double [^double x ^double y ^double z]
  (fma-macro x y z))

(defn fma
  "Computes `x*y + z` using fused multiply-add when available. Identical implementation to [[muladd]].

  Parameters:

  - `x`, `y`, `z` (doubles): multiplicands and addend.

  Returns `x*y + z` as a double. On Java 9+, uses `Math/fma` for a single correctly-rounded fused multiply-add (more accurate than a separate multiply and add). On earlier JVMs, falls back to a plain `(+ z (* x y))` with intermediate rounding.

  See also [[muladd]], [[negmuladd]], [[difference-of-products]], [[sum-of-products]]."
  {:inline (->fma)
   :inline-arities #{3}}
  ^double [^double x ^double y ^double z]
  (fma-macro x y z))

(defn negmuladd
  "Computes `-x*y + z` (negated fused multiply-add) using fused multiply-add when available.

  Parameters:

  - `x`, `y`, `z` (doubles): multiplicands (negated) and addend.

  Returns `-x*y + z` as a double. On Java 9+, uses `Math/fma` for a single correctly-rounded fused multiply-add. On earlier JVMs, falls back to a plain `(+ z (* (- x) y))` with intermediate rounding.

  See also [[muladd]], [[fma]]."
  {:inline (->fma true)
   :inline-arities #{3}}
  ^double [^double x ^double y ^double z]
  (fma-macro (- x) y z))

;; some stuff from pbrt
(defn difference-of-products
  "Computes `a*b - c*d` using Kahan's two-product algorithm to avoid catastrophic cancellation when `a*b` and `c*d` are close in magnitude.

  Parameters:

  - `a`, `b`, `c`, `d` (doubles): factors of the two products.

  Returns `a*b - c*d` as a double, computed via two fused multiply-adds plus a correction term rather than a direct subtraction of two separately-rounded products. This achieves near full double precision even when the naive `(- (* a b) (* c d))` loses many significant digits to cancellation. Requires a true hardware/JVM fused multiply-add (Java 9+) to realize the accuracy benefit; on earlier JVMs it is numerically equivalent to the naive computation.

  See also [[sum-of-products]], [[fma]]."
  ^double [^double a ^double b ^double c ^double d]
  (let [cd (* c d)]
    (+ (fma a b (- cd)) (fma (- c) d cd))))

(defn sum-of-products
  "Computes `a*b + c*d` using Kahan's two-product algorithm to avoid catastrophic cancellation when `a*b` and `c*d` have opposite signs and are close in magnitude.

  Parameters:

  - `a`, `b`, `c`, `d` (doubles): factors of the two products.

  Returns `a*b + c*d` as a double, computed via two fused multiply-adds rather than a direct sum of two separately-rounded products. This achieves near full double precision even when the naive `(+ (* a b) (* c d))` loses many significant digits to cancellation. Requires a true hardware/JVM fused multiply-add (Java 9+) to realize the accuracy benefit; on earlier JVMs it is numerically equivalent to the naive computation.

  See also [[difference-of-products]], [[fma]]."
  ^double [^double a ^double b ^double c ^double d]
  (let [cd (* c d)]
    (+ (fma a b cd) (fma c d (- cd)))))

;; Processing math constants
(def ^{:const true :tag 'double :doc "Value of $\\pi$"} PI Math/PI)
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\pi}{2}$"} HALF_PI (* PI 0.5))
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\pi}{3}$"} THIRD_PI (/ PI 3.0))
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\pi}{4}$"} QUARTER_PI (* PI 0.25))
(def ^{:const true :tag 'double :doc "Value of $2\\pi$"} TWO_PI (+ PI PI))
(def ^{:const true :tag 'double :doc "Value of $2\\pi$"} TAU TWO_PI)
(def ^{:const true :tag 'double :doc "Value of $\\mathrm{e}$"} E Math/E)
(def ^{:const true :tag 'double :doc "Value of $-\\pi$"} -PI (- Math/PI))
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\pi}{2}$"} -HALF_PI (* PI -0.5))
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\pi}{3}$"} -THIRD_PI (/ -PI 3.0))
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\pi}{4}$"} -QUARTER_PI (* PI -0.25))
(def ^{:const true :tag 'double :doc "Value of $-2\\pi$"} -TWO_PI (- TWO_PI))
(def ^{:const true :tag 'double :doc "Value of $-2\\pi$"} -TAU -TWO_PI)
(def ^{:const true :tag 'double :doc "Value of $-\\mathrm{e}$"} -E (- Math/E))
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{\\pi}$"} INV_PI (/ PI))
(def ^{:const true :tag 'double :doc "Value of $\\frac{2}{\\pi}$"} TWO_INV_PI (/ 2.0 PI))
(def ^{:const true :tag 'double :doc "Value of $\\frac{4}{\\pi}$"} FOUR_INV_PI (/ 4.0 PI))
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{2\\pi}$"} INV_TWO_PI (/ TWO_PI))
(def ^{:const true :tag 'double :doc "Value of $\\frac{2}{2\\pi}$"} INV_FOUR_PI (/ (* 2.0 TWO_PI)))
(def ^{:const true :tag 'double :doc "$\\varepsilon$, a small number"} EPSILON 1.0e-10)
(def ^{:const true :tag 'double :doc "$\\gamma$, Euler-Mascheroni constant"} GAMMA Gamma/GAMMA)
(def ^{:const true :tag 'double :doc "Lanchos approximation of `g` constant"} LANCZOS_G Gamma/LANCZOS_G)
(def ^{:const true :tag 'double :doc "Catalan G"} CATALAN_G 0.91596559417721901505)
(def ^{:const true :tag 'double :doc "Value of $\\pi^2$"} PI2 (* Math/PI Math/PI))

(defonce ^{:const true :tag 'double :doc "ulp(1)/2"}
  MACHINE-EPSILON (* 0.5 (FastMath/ulp 1.0)))

(def ^{:const true :tag 'double :doc "5ulp(1)"}
  MACHINE-EPSILON10 (* 10.0 MACHINE-EPSILON))

(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{3}$"} THIRD      0.333333333333333333333333)
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{3}$"} ONE_THIRD  0.333333333333333333333333)
(def ^{:const true :tag 'double :doc "Value of $\\frac{2}{3}$"} TWO_THIRD  0.666666666666666666666666)
(def ^{:const true :tag 'double :doc "Value of $\\frac{2}{3}$"} TWO_THIRDS 0.666666666666666666666666)
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{6}$"} SIXTH      0.166666666666666666666666)
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{6}$"} ONE_SIXTH  0.166666666666666666666666)

(defn signum
  "Returns 1 if `value` is > 0, 0 if it is 0, -1 otherwise. See also [[sgn]]."
  {:inline (fn [v] `(if (pos? (double ~v)) 1.0 (if (neg? (double ~v)) -1.0 0.0)))
   :inline-arities #{1}}
  ^double [^double value]
  (cond (pos? value) 1.0
        (neg? value) -1.0
        :else 0.0))

(defn sgn
  "Returns -1 when `value` is negative, 1 otherwise. See also [[signum]]."
  {:inline (fn [v] `(if (neg? (double ~v)) -1.0 1.0))
   :inline-arities #{1}}
  ^double [^double value]
  (if (neg? value) -1.0 1.0))

;; copy-sign

(defn copy-sign
  "Returns a value with a magnitude of first argument and sign of second."
  {:inline (fn [magnitude sign] `(FastMath/copySign (double ~magnitude) (double ~sign)))
   :inline-arities #{2}}
  ^double [^double magnitude ^double sign]
  (FastMath/copySign magnitude sign))

;; trigonometry

(defn sin
  "sin(x)"
  {:inline (fn [x] `(. FastMath (sin (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (sin x)))

(defn sinpi
  "Computes `sin(π·x)` -- the sine of `x` expressed in half-turns (units of π) rather than radians.

  Parameters:

  - `x` (double): value in half-turns to take the sine of.

  Returns `sin(π·x)` as a double, equivalent to `(sin (* PI x))`. Near-exact identity values (e.g. `x=0.5` giving `1.0`) are not guaranteed to be bit-exact, since `π` is only finitely represented as a double.

  See also [[sin]], [[cospi]], [[tanpi]]."
  {:inline (fn [x] `(. FastMath (sin (* PI (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (sin (* PI x))))

(defn cos
  "cos(x)"
  {:inline (fn [x] `(. FastMath (cos (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (cos x)))

(defn cospi
  "Computes `cos(π·x)` -- the cosine of `x` expressed in half-turns (units of π) rather than radians.

  Parameters:

  - `x` (double): value in half-turns to take the cosine of.

  Returns `cos(π·x)` as a double, equivalent to `(cos (* PI x))`. Near-exact identity values (e.g. `x=1.0` giving `-1.0`) are not guaranteed to be bit-exact, since `π` is only finitely represented as a double.

  See also [[cos]], [[sinpi]], [[tanpi]]."
  {:inline (fn [x] `(. FastMath (cos (* PI (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (cos (* PI x))))

(defn tan
  "tan(x)"
  {:inline (fn [x] `(. FastMath (tan (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (tan x)))

(defn tanpi
  "Computes `tan(π·x)` -- the tangent of `x` expressed in half-turns (units of π) rather than radians.

  Parameters:

  - `x` (double): value in half-turns to take the tangent of.

  Returns `tan(π·x)` as a double, equivalent to `(tan (* PI x))`. Diverges to very large magnitudes near `x = k+0.5` for integer `k` (where `cos(π·x)` is near zero), matching plain [[tan]]'s behavior at its own poles.

  See also [[tan]], [[sinpi]], [[cospi]]."
  {:inline (fn [x] `(. FastMath (tan (* PI (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (tan (* PI x))))

(defn asin
  "asin(x)"
  {:inline (fn [x] `(. FastMath (asin (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (asin x)))

(defn acos
  "acos(x)"
  {:inline (fn [x] `(. FastMath (acos (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (acos x)))

(defn atan
  "atan(x)"
  {:inline (fn [x] `(. FastMath (atan (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (atan x)))

(defn sinh
  "sinh(x)"
  {:inline (fn [x] `(. FastMath (sinh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (sinh x)))

(defn cosh
  "cosh(x)"
  {:inline (fn [x] `(. FastMath (cosh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (cosh x)))

(defn tanh
  "tanh(x)"
  {:inline (fn [x] `(. FastMath (tanh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (tanh x)))

(defn asinh
  "asinh(x)"
  {:inline (fn [x] `(. FastMath (asinh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (asinh x)))

(defn acosh
  "acosh(x)"
  {:inline (fn [x] `(. FastMath (acosh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (acosh x)))

(defn atanh
  "atanh(x)"
  {:inline (fn [x] `(. FastMath (atanh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (atanh x)))

(defn qsin
  "Fast and less accurate sin(x)"
  {:inline (fn [x] `(. FastMath (sinQuick (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (sinQuick x)))

(defn qcos
  "Fast and less accurate cos(x)"
  {:inline (fn [x] `(. FastMath (cosQuick (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (cosQuick x)))

;; Additional trigonometry functions

(defn cot
  "Computes the cotangent of `x`, `cot(x) = 1/tan(x)`.

  Parameters:

  - `x` (double): angle in radians.

  Returns `cot(x)` as a double. Diverges where `tan(x)` is zero, i.e. at integer multiples of π.

  See also [[tan]], [[cotpi]], [[sec]], [[csc]]."
  {:inline (fn [x] `(/ (tan (double ~x))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/tan x))))

(defn cotpi
  "Computes the cotangent of `x` expressed in half-turns (units of π), `cot(π·x) = 1/tan(π·x)`.

  Parameters:

  - `x` (double): value in half-turns.

  Returns `cot(π·x)` as a double. Diverges near integer `x` (where `tan(π·x)` is near zero).

  See also [[cot]], [[tanpi]], [[secpi]], [[cscpi]]."
  {:inline (fn [x] `(/ (tan (* PI (double ~x)))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/tan (* PI x)))))

(defn sec
  "Computes the secant of `x`, `sec(x) = 1/cos(x)`.

  Parameters:

  - `x` (double): angle in radians.

  Returns `sec(x)` as a double. Diverges where `cos(x)` is zero, i.e. at odd multiples of π/2.

  See also [[cos]], [[secpi]], [[cot]], [[csc]]."
  {:inline (fn [x] `(/ (cos (double ~x))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/cos x))))

(defn secpi
  "Computes the secant of `x` expressed in half-turns (units of π), `sec(π·x) = 1/cos(π·x)`.

  Parameters:

  - `x` (double): value in half-turns.

  Returns `sec(π·x)` as a double. Diverges near half-integer `x` (where `cos(π·x)` is near zero).

  See also [[sec]], [[cospi]], [[cotpi]], [[cscpi]]."
  {:inline (fn [x] `(/ (cos (* PI (double ~x)))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/cos (* PI x)))))

(defn csc
  "Computes the cosecant of `x`, `csc(x) = 1/sin(x)`.

  Parameters:

  - `x` (double): angle in radians.

  Returns `csc(x)` as a double. Diverges where `sin(x)` is zero, i.e. at integer multiples of π.

  See also [[sin]], [[cscpi]], [[cot]], [[sec]]."
  {:inline (fn [x] `(/ (sin (double ~x))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/sin x))))

(defn cscpi
  "Computes the cosecant of `x` expressed in half-turns (units of π), `csc(π·x) = 1/sin(π·x)`.

  Parameters:

  - `x` (double): value in half-turns.

  Returns `csc(π·x)` as a double. Diverges near integer `x` (where `sin(π·x)` is near zero).

  See also [[csc]], [[sinpi]], [[cotpi]], [[secpi]]."
  {:inline (fn [x] `(/ (sin (* PI (double ~x)))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/sin (* PI x)))))

;; Additional cyclometric functions

(defn acot
  "Computes the inverse cotangent of `x`, `acot(x) = π/2 - atan(x)`.

  Parameters:

  - `x` (double): value to take the inverse cotangent of.

  Returns `acot(x)` as a double, in the range `(0, π)`.

  See also [[atan]], [[asec]], [[acsc]]."
  {:inline (fn [x] `(- HALF_PI (atan (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (- HALF_PI (FastMath/atan x)))

(defn asec
  "Computes the inverse secant of `x`, `asec(x) = acos(1/x)`.

  Parameters:

  - `x` (double): value to take the inverse secant of; must satisfy `|x| >= 1`.

  Returns `asec(x)` as a double, in the range `[0, π]`. Returns `##NaN` for `|x| < 1` (outside the domain).

  See also [[acos]], [[acot]], [[acsc]]."
  {:inline (fn [x] `(acos (/ ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acos (/ 1.0 x)))

(defn acsc
  "Computes the inverse cosecant of `x`, `acsc(x) = asin(1/x)`.

  Parameters:

  - `x` (double): value to take the inverse cosecant of; must satisfy `|x| >= 1`.

  Returns `acsc(x)` as a double, in the range `[-π/2, π/2]`. Returns `##NaN` for `|x| < 1` (outside the domain).

  See also [[asin]], [[acot]], [[asec]]."
  {:inline (fn [x] `(asin (/ ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/asin (/ 1.0 x)))

(defn atan2
  "atan2(x,y)"
  {:inline (fn [x y] `(. FastMath (atan2 (double ~x) (double ~y))))
   :inline-arities #{2}}
  ^double [^double x ^double y] (FastMath/atan2 x y))

;; Additional hyperbolic functions
(defn coth
  "Computes the hyperbolic cotangent of `x`, `coth(x) = 1/tanh(x)`.

  Parameters:

  - `x` (double): value to take the hyperbolic cotangent of.

  Returns `coth(x)` as a double. Diverges at `x=0` (where `tanh(x)` is zero).

  See also [[tanh]], [[sech]], [[csch]]."
  {:inline (fn [x] `(/ (tanh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (/ (FastMath/tanh x)))

(defn sech
  "Computes the hyperbolic secant of `x`, `sech(x) = 1/cosh(x)`.

  Parameters:

  - `x` (double): value to take the hyperbolic secant of.

  Returns `sech(x)` as a double, always in `(0, 1]`.

  See also [[cosh]], [[coth]], [[csch]]."
  {:inline (fn [x] `(/ (cosh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (/ (FastMath/cosh x)))

(defn csch
  "Computes the hyperbolic cosecant of `x`, `csch(x) = 1/sinh(x)`.

  Parameters:

  - `x` (double): value to take the hyperbolic cosecant of.

  Returns `csch(x)` as a double. Diverges at `x=0` (where `sinh(x)` is zero).

  See also [[sinh]], [[coth]], [[sech]]."
  {:inline (fn [x] `(/ (sinh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (/ (FastMath/sinh x)))

;; Additional inverse hyperbolic functions
(defn acoth
  "Computes the area (inverse) hyperbolic cotangent of `x`, `acoth(x) = atanh(1/x)`.

  Parameters:

  - `x` (double): value to take the area hyperbolic cotangent of; must satisfy `|x| > 1`.

  Returns `acoth(x)` as a double. Returns `##NaN` for `|x| <= 1` (outside the domain).

  See also [[atanh]], [[asech]], [[acsch]]."
  {:inline (fn [x] `(atanh (/ ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/atanh (/ x)))

(defn asech
  "Computes the area (inverse) hyperbolic secant of `x`, `asech(x) = acosh(1/x)`.

  Parameters:

  - `x` (double): value to take the area hyperbolic secant of; must satisfy `0 < x <= 1`.

  Returns `asech(x)` as a double, always non-negative. Returns `##NaN` outside the domain.

  See also [[acosh]], [[acoth]], [[acsch]]."
  {:inline (fn [x] `(acosh (/ ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acosh (/ x)))

(defn acsch
  "Computes the area (inverse) hyperbolic cosecant of `x`, `acsch(x) = asinh(1/x)`.

  Parameters:

  - `x` (double): value to take the area hyperbolic cosecant of; must be nonzero.

  Returns `acsch(x)` as a double. Returns `##Inf`/`##-Inf` at `x=0` per IEEE 754 division semantics (`1/0 = ##Inf`).

  See also [[asinh]], [[acoth]], [[asech]]."
  {:inline (fn [x] `(asinh (/ ~x)))
   :inline-arities #{1}}
  ^double [^double v] (FastMath/asinh (/ v)))

;; historical

(defn crd
  "Computes the chord length of an arc, `crd(x) = 2*sin(x/2)`, for a unit circle.

  Parameters:

  - `x` (double): central angle in radians.

  Returns the chord length as a double.

  See also [[acrd]], [[sin]]."
  {:inline (fn [x] `(* 2.0 (sin (* 0.5 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (* 2.0 (FastMath/sin (* 0.5 x))))

(defn acrd
  "Computes the inverse chord function, `acrd(x) = 2*asin(x/2)`.

  Parameters:

  - `x` (double): chord length; must satisfy `|x| <= 2` for a real result.

  Returns the central angle in radians as a double. Returns `##NaN` outside the domain.

  See also [[crd]], [[asin]]."
  {:inline (fn [x] `(* 2.0 (asin (* 0.5 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (* 2.0 (FastMath/asin (* 0.5 x))))

(defn versin
  "Computes the versine (versed sine) of `x`, `versin(x) = 1 - cos(x)`.

  Parameters:

  - `x` (double): angle in radians.

  Returns the versine as a double, always in `[0, 2]`.

  See also [[aversin]], [[coversin]], [[vercos]], [[haversin]]."
  {:inline (fn [x] `(- 1.0 (cos (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (- 1.0 (FastMath/cos x)))

(defn coversin
  "Computes the coversine (coversed sine) of `x`, `coversin(x) = 1 - sin(x)`.

  Parameters:

  - `x` (double): angle in radians.

  Returns the coversine as a double, always in `[0, 2]`.

  See also [[acoversin]], [[versin]], [[covercos]]."
  {:inline (fn [x] `(- 1.0 (sin (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (- 1.0 (FastMath/sin x)))

(defn vercos
  "Computes the vercosine (versed cosine) of `x`, `vercos(x) = 1 + cos(x)`.

  Parameters:

  - `x` (double): angle in radians.

  Returns the vercosine as a double, always in `[0, 2]`.

  See also [[avercos]], [[versin]], [[covercos]]."
  {:inline (fn [x] `(inc (cos (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (inc (FastMath/cos x)))

(defn covercos
  "Computes the covercosine (coversed cosine) of `x`, `covercos(x) = 1 + sin(x)`.

  Parameters:

  - `x` (double): angle in radians.

  Returns the covercosine as a double, always in `[0, 2]`.

  See also [[acovercos]], [[coversin]], [[vercos]]."
  {:inline (fn [x] `(inc (sin (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (inc (FastMath/sin x)))

(defn aversin
  "Computes the arc (inverse) versine of `x`, `aversin(x) = acos(1 - x)`.

  Parameters:

  - `x` (double): versine value; must satisfy `0 <= x <= 2` for a real result.

  Returns the angle in radians as a double, in `[0, π]`. Returns `##NaN` outside the domain.

  See also [[versin]], [[acoversin]], [[avercos]]."
  {:inline (fn [x] `(acos (- 1.0 ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acos (- 1.0 x)))

(defn acoversin
  "Computes the arc (inverse) coversine of `x`, `acoversin(x) = asin(1 - x)`.

  Parameters:

  - `x` (double): coversine value; must satisfy `0 <= x <= 2` for a real result.

  Returns the angle in radians as a double, in `[-π/2, π/2]`. Returns `##NaN` outside the domain.

  See also [[coversin]], [[aversin]], [[acovercos]]."
  {:inline (fn [x] `(asin (- 1.0 ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/asin (- 1.0 x)))

(defn avercos
  "Computes the arc (inverse) vercosine of `x`, `avercos(x) = acos(x - 1)`.

  Parameters:

  - `x` (double): vercosine value; must satisfy `0 <= x <= 2` for a real result.

  Returns the angle in radians as a double, in `[0, π]`. Returns `##NaN` outside the domain.

  See also [[vercos]], [[aversin]], [[acovercos]]."
  {:inline (fn [x] `(acos (dec ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acos (dec x)))

(defn acovercos
  "Computes the arc (inverse) covercosine of `x`, `acovercos(x) = asin(x - 1)`.

  Parameters:

  - `x` (double): covercosine value; must satisfy `0 <= x <= 2` for a real result.

  Returns the angle in radians as a double, in `[-π/2, π/2]`. Returns `##NaN` outside the domain.

  See also [[covercos]], [[acoversin]], [[avercos]]."
  {:inline (fn [x] `(asin (dec ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/asin (dec x)))

(defn haversin
  "Computes the haversine of `x`, `haversin(x) = (1 - cos(x))/2`, or the haversine-formula central angle between two latitude/longitude points (in radians).

  Parameters:

  - `x` (double): angle in radians, for the 1-arity form.
  - `[lat1 lon1]`, `[lat2 lon2]` (pairs of doubles): coordinates in radians, for the 2-arity form.
  - `lat1`, `lon1`, `lat2`, `lon2` (doubles): coordinates in radians, for the 4-arity form.

  Returns the haversine value as a double, always in `[0, 1]`. For the 2- and 4-arity forms, returns the haversine of the central angle between the two points, via `hav(dlat) + cos(lat1)*cos(lat2)*hav(dlon)`.

  See also [[haversine]] (alias), [[haversine-dist]], [[versin]]."
  {:inline (fn [x] `(* 0.5 (- 1.0 (cos (double ~x)))))
   :inline-arities #{1}}
  (^double [^double x] (* 0.5 (- 1.0 (FastMath/cos x))))
  (^double [[^double lat1 ^double lon1] [^double lat2 ^double lon2]]
   (haversin lat1 lon1 lat2 lon2))
  (^double [^double lat1 ^double lon1 ^double lat2 ^double lon2]
   (+ (haversin (- lat2 lat1))
      (* (FastMath/cos lat1)
         (FastMath/cos lat2)
         (haversin (- lon2 lon1))))))

(def ^{:doc "Haversine ([[haversin]] alias)"} haversine haversin)

(defn hacoversin
  "Computes the hacoversine of `x`, `hacoversin(x) = (1 - sin(x))/2`.

  Parameters:

  - `x` (double): angle in radians.

  Returns the hacoversine as a double, always in `[0, 1]`.

  See also [[ahacoversin]], [[coversin]], [[haversin]]."
  {:inline (fn [x] `(* 0.5 (- 1.0 (sin (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x] (* 0.5 (- 1.0 (FastMath/sin x))))

(defn havercos
  "Computes the havercosine of `x`, `havercos(x) = (1 + cos(x))/2`.

  Parameters:

  - `x` (double): angle in radians.

  Returns the havercosine as a double, always in `[0, 1]`.

  See also [[ahavercos]], [[vercos]], [[haversin]]."
  {:inline (fn [x] `(* 0.5 (inc (cos (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x] (* 0.5 (inc (FastMath/cos x))))

(defn hacovercos
  "Computes the hacovercosine of `x`, `hacovercos(x) = (1 + sin(x))/2`.

  Parameters:

  - `x` (double): angle in radians.

  Returns the hacovercosine as a double, always in `[0, 1]`.

  See also [[ahacovercos]], [[covercos]], [[haversin]]."
  {:inline (fn [x] `(* 0.5 (inc (sin (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x] (* 0.5 (inc (FastMath/sin x))))

(defn ahaversin
  "Computes the arc (inverse) haversine of `x`, `ahaversin(x) = acos(1 - 2x)`.

  Parameters:

  - `x` (double): haversine value; must satisfy `0 <= x <= 1` for a real result.

  Returns the angle in radians as a double, in `[0, π]`. Returns `##NaN` outside the domain.

  See also [[haversin]], [[ahacoversin]], [[ahavercos]]."
  {:inline (fn [x] `(acos (- 1.0 (* 2.0 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acos (- 1.0 (* 2.0 x))))

(defn ahacoversin
  "Computes the arc (inverse) hacoversine of `x`, `ahacoversin(x) = asin(1 - 2x)`.

  Parameters:

  - `x` (double): hacoversine value; must satisfy `0 <= x <= 1` for a real result.

  Returns the angle in radians as a double, in `[-π/2, π/2]`. Returns `##NaN` outside the domain.

  See also [[hacoversin]], [[ahaversin]], [[ahacovercos]]."
  {:inline (fn [x] `(asin (- 1.0 (* 2.0 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/asin (- 1.0 (* 2.0 x))))

(defn ahavercos
  "Computes the arc (inverse) havercosine of `x`, `ahavercos(x) = acos(2x - 1)`.

  Parameters:

  - `x` (double): havercosine value; must satisfy `0 <= x <= 1` for a real result.

  Returns the angle in radians as a double, in `[0, π]`. Returns `##NaN` outside the domain.

  See also [[havercos]], [[ahaversin]], [[ahacovercos]]."
  {:inline (fn [x] `(acos (dec (* 2.0 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acos (dec (* 2.0 x))))

(defn ahacovercos
  "Computes the arc (inverse) hacovercosine of `x`, `ahacovercos(x) = asin(2x - 1)`.

  Parameters:

  - `x` (double): hacovercosine value; must satisfy `0 <= x <= 1` for a real result.

  Returns the angle in radians as a double, in `[-π/2, π/2]`. Returns `##NaN` outside the domain.

  See also [[hacovercos]], [[ahacoversin]], [[ahavercos]]."
  {:inline (fn [x] `(asin (dec (* 2.0 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/asin (dec (* 2.0 x))))

(defn exsec
  "Computes the exsecant of `x`, `exsec(x) = sec(x) - 1`.

  Parameters:

  - `x` (double): angle in radians.

  Returns the exsecant as a double.

  See also [[aexsec]], [[sec]], [[excsc]]."
  {:inline (fn [x] `(dec (sec ~x)))
   :inline-arities #{1}}
  ^double [^double x] (dec (sec x)))

(defn excsc
  "Computes the excosecant of `x`, `excsc(x) = csc(x) - 1`.

  Parameters:

  - `x` (double): angle in radians.

  Returns the excosecant as a double.

  See also [[aexcsc]], [[csc]], [[exsec]]."
  {:inline (fn [x] `(dec (csc ~x)))
   :inline-arities #{1}}
  ^double [^double x] (dec (csc x)))

(defn aexsec
  "Computes the arc (inverse) exsecant of `x`, `aexsec(x) = asec(x + 1)`.

  Parameters:

  - `x` (double): exsecant value; must satisfy `x >= 0` or `x <= -2` for a real result.

  Returns the angle in radians as a double, in `[0, π]`. Returns `##NaN` outside the domain.

  See also [[exsec]], [[asec]], [[aexcsc]]."
  {:inline (fn [x] `(asec (inc ~x)))
   :inline-arities #{1}}
  ^double [^double x] (asec (inc x)))

(defn aexcsc
  "Computes the arc (inverse) excosecant of `x`, `aexcsc(x) = acsc(x + 1)`.

  Parameters:

  - `x` (double): excosecant value; must satisfy `x >= 0` or `x <= -2` for a real result.

  Returns the angle in radians as a double, in `[-π/2, π/2]`. Returns `##NaN` outside the domain.

  See also [[excsc]], [[acsc]], [[aexsec]]."
  {:inline (fn [x] `(acsc (inc ~x)))
   :inline-arities #{1}}
  ^double [^double x] (acsc (inc x)))

(defn haversine-dist
  "Computes the great-circle distance between two latitude/longitude points on a unit sphere (`r=1`), using the haversine formula.

  Parameters:

  - `[lat1 lon1]`, `[lat2 lon2]` (pairs of doubles): coordinates in radians, for the 2-arity form.
  - `lat1`, `lon1`, `lat2`, `lon2` (doubles): coordinates in radians, for the 4-arity form.

  Returns the great-circle distance in radians as a double (the central angle between the two points). Multiply by a sphere's radius to get the distance in that radius's units (e.g. multiply by Earth's mean radius, ~6371 km, for a distance in kilometers).

  See also [[haversin]]."
  (^double [[^double lat1 ^double lon1] [^double lat2 ^double lon2]]
   (haversine-dist lat1 lon1 lat2 lon2))
  (^double [^double lat1 ^double lon1 ^double lat2 ^double lon2]
   (* 2.0 (FastMath/asin (FastMath/sqrt (haversin lat1 lon1 lat2 lon2))))))

;; exp and log

(defn exp
  "exp(x) = e^x"
  {:inline (fn [x] `(. FastMath (exp (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (exp x)))

(defn exp2
  "exp2(x) = 2^x"
  {:inline (fn [x] `(. Math (pow 2.0 (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. Math (pow 2.0 x)))

(defn exp10
  "exp10(x) = 10^x"
  {:inline (fn [x] `(. Math (pow 10.0 (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. Math (pow 10.0 x)))

(defn log
  "Computes the natural logarithm of `x`, or the logarithm of `x` with a given `base`.

  Parameters:

  - `x` (double): value to take the logarithm of; must be positive for a real result.
  - `base`, `x` (doubles), 2-arity: computes `log(x)/log(base)` (change of base).

  Returns the logarithm as a double. Returns `##NaN` for `x < 0`, and `##-Inf` for `x = 0`, matching `FastMath/log`'s IEEE 754 behavior.

  See also [[ln]] (alias, 1-arity only), [[log10]], [[log2]], [[logb]], [[log1p]]."
  {:inline (fn ([x] `(. FastMath (log (double ~x))))
             ([base x] `(/ (. FastMath (log (double ~x))) (. FastMath (log (double ~base))))))
   :inline-arities #{1 2}}
  (^double [^double x] (. FastMath (log x)))
  (^double [^double base ^double x] (/ (. FastMath (log x)) (. FastMath (log base)))))

(defn ln
  "log(x)=ln(x)"
  {:inline (fn [x] `(. FastMath (log (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (log x)))

(defn log10
  "log_10(x)"
  {:inline (fn [x] `(. FastMath (log10 (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (log10 x)))

(defn log1p
  "log(1+x) for small x"
  {:inline (fn [x] `(. FastMath (log1p (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (log1p x)))

(defn expm1
  "exp(x)-1 for small x"
  {:inline (fn [x] `(. FastMath (expm1 (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (expm1 x)))

(defn exprel
  "Computes `(exp(x) - 1) / x`, the relative rate of change of `exp` -- numerically stable near `x=0`, where the naive formula suffers catastrophic cancellation.

  Parameters:

  - `x` (double): value to evaluate.

  Returns the result as a double. Returns exactly `1.0` for `|x|` below `10` times machine epsilon (the correct limiting value as `x -> 0`), and `##Inf` for `x > 717.0` (avoiding `exp` overflow before it would occur), otherwise computes `expm1(x)/x` directly.

  See also [[expm1]], [[exp]]."
  {:inline (fn [x] `(let [x# (double ~x)]
                     (cond
                       (< (. FastMath abs x#) MACHINE-EPSILON10) 1.0
                       (> x# 717.0) ##Inf
                       :else (/ (. FastMath (expm1 x#)) x#))))
   :inline-arities #{1}}
  ^double [^double x]
  (cond
    (< (. FastMath abs x) MACHINE-EPSILON10) 1.0
    (> x 717.0) ##Inf
    :else (/ (. FastMath (expm1 x)) x)))

(def ^{:const true :tag 'double :doc "Value of $\\ln{2}$"} LN2 (log 2.0))
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{\\ln{2}}$"} INV_LN2 (/ LN2))
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\ln{2}}{2}$"} LN2_2 (* 0.5 LN2))
(def ^{:const true :tag 'double :doc "Value of $\\ln{10}$"} LN10 (log 10.0))
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{\\ln{\\frac{1}{2}}}$"} INV_LOG_HALF (/ (log 0.5)))
(def ^{:const true :tag 'double :doc "Value of $\\ln{\\frac{1}{2}}$"} LOG_HALF (log 0.5))
(def ^{:const true :tag 'double :doc "Value of $\\ln{\\pi}$"} LOG_PI (log PI))
(def ^{:const true :tag 'double :doc "Value of $\\ln{2\\pi}$"} LOG_TWO_PI (log TWO_PI))

(defn log1pexp
  "Computes `log(1+exp(x))` (the softplus function), using a numerically stable piecewise approximation to avoid overflow for large `x` and underflow for very negative `x`.

  Parameters:

  - `x` (double): value to evaluate.

  Returns the result as a double. Returns `0.0` for very negative `x`, and `x` itself for very large `x` (where `exp(x)` would overflow but the correction term becomes negligible), with a stable computation in between.

  See also [[log1mexp]], [[logaddexp]], [[expm1]]."
  ^double [^double x]
  (cond
    (< x -745.1332191019412) 0.0
    (< x -36.7368005696771) (FastMath/exp x)
    (< x 18.021826694558577) (FastMath/log1p (FastMath/exp x))
    (< x 33.23111882352963) (+ x (FastMath/exp (- x)))
    :else x))

(defn log1mexp
  "Computes `log(1-exp(x))` for `x < 0`, using a numerically stable form to avoid catastrophic cancellation near `x=0`.

  Parameters:

  - `x` (double): value to evaluate; must be negative.

  Returns the result as a double, always non-positive.

  See also [[log1pexp]], [[log2mexp]]."
  ^double [^double x]
  (if (< x LOG_HALF)
    (FastMath/log1p (- (FastMath/exp x)))
    (FastMath/log (- (FastMath/expm1 x)))))

(defn log2mexp
  "Computes `log(2-exp(x))`, equivalent to `log1p(-expm1(x))`.

  Parameters:

  - `x` (double): value to evaluate.

  Returns the result as a double.

  See also [[log1mexp]], [[expm1]]."
  {:inline (fn [x] `(FastMath/log1p (- (FastMath/expm1 (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x]
  (FastMath/log1p (- (FastMath/expm1 x))))

(defn log1psq
  "Computes `log(1+x^2)`, switching to a direct formula for very large `x` to avoid unnecessary precision loss in `x*x`.

  Parameters:

  - `x` (double): value to evaluate.

  Returns the result as a double, always non-negative.

  See also [[log1p]]."
  ^double [^double x]
  (if (< x 9007199254740992)
    (FastMath/log1p (* x x))
    (* 2.0 (log x))))

(defn logexpm1
  "Computes `log(exp(x)-1)`, the inverse of [[log1pexp]] (softplus).

  Parameters:

  - `x` (double): value to evaluate; must be positive for a real result (`exp(x)-1` must be positive).

  Returns the result as a double. Returns `##NaN` (via `log` of a non-positive value) for `x <= 0`.

  See also [[log1pexp]], [[expm1]]."
  {:inline (fn [x] `(log (expm1 (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/log (FastMath/expm1 x)))

;; from julia
(defn- log1pmx-ker
  ^double [^double x]
  (let [r (/ x (+ 2.0 x))
        t (* r r)
        w (muladd t (muladd t (muladd t (muladd t (muladd t (muladd t (muladd t 0.11764705882352941 0.13333333333333333) 0.15384615384615385) 0.18181818181818182) 0.2222222222222222) 0.2857142857142857) 0.4) 0.6666666666666666)
        hxsq (* 0.5 x x)]
    (- (* r (+ hxsq (* w t))) hxsq)))

(defn log1pmx
  "Computes `log(1+x) - x`, using a stable kernel-based polynomial approximation near `x=0` to avoid catastrophic cancellation (ported from Julia's `Base.Math`).

  Parameters:

  - `x` (double): value to evaluate; must satisfy `x > -1` for a real result.

  Returns the result as a double, always non-positive.

  See also [[logmxp1]], [[log1p]]."
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
  "Computes `log(x) - x + 1`, using a stable formulation for `x` near `1` (reusing [[log1pmx]]'s kernel).

  Parameters:

  - `x` (double): value to evaluate; must be positive for a real result.

  Returns the result as a double, always non-positive.

  See also [[log1pmx]]."
  ^double [^double x]
  (cond
    (<= x 0.3) (- (inc (FastMath/log x)) x)
    (<= x 0.4) (let [u (/ (- x 0.375) 0.375)]
                 (+ (log1pmx-ker u) -3.55829253011726237e-1 (* 0.625 u)))
    (<= x 0.6) (let [u (* (- x 0.5) 2.0)]
                 (+ (log1pmx-ker u) -1.93147180559945309e-1 (* 0.5 u)))
    :else (log1pmx (dec x))))

(defn logaddexp
  "Computes `log(exp(x)+exp(y))`, the numerically stable 2-argument log-sum-exp.

  Parameters:

  - `x`, `y` (doubles): values to combine.

  Returns the result as a double, computed without overflowing `exp(x)`/`exp(y)` for large arguments.

  See also [[logsumexp]] (n-ary version), [[logsubexp]], [[log1pexp]]."
  ^double [^double x ^double y]
  (if (< x y)
    (+ y (log1pexp (- x y)))
    (+ (if-not (Double/isNaN y) x y)
       (log1pexp (- y x)))))

(defn logsubexp
  "Computes `log(abs(exp(x)-exp(y)))`, the numerically stable log-difference-of-exponentials.

  Parameters:

  - `x`, `y` (doubles): values to combine.

  Returns the result as a double. Returns `##-Inf` when `x` equals `y` (since the difference is zero).

  See also [[logaddexp]], [[log1mexp]]."
  ^double [^double x ^double y]
  (+ (PrimitiveMath/max x y)
     (log1mexp (- (if (and (== x y)
                           (or (Double/isFinite x) (neg? x))) 0.0 (Math/abs (- x y)))))))

(defn logsumexp
  "Computes `log(exp(x1)+...+exp(xn))`, the numerically stable n-ary log-sum-exp, using an online (single-pass) running-maximum algorithm.

  Parameters:

  - `xs` (sequence of doubles): values to combine.

  Returns the result as a double, computed without overflowing any individual `exp(xi)` for large arguments.

  See also [[logaddexp]] (2-arity version)."
  ^double [xs]
  (loop [xs xs
         r 0.0
         alpha ##-Inf]
    (let [x (double (first xs))
          rst (rest xs)]
      (if (<= x alpha)
        (let [nr (+ r (FastMath/exp (- x alpha)))]
          (if-not (seq rst)
            (+ (FastMath/log nr) alpha)
            (recur rst nr alpha)))
        (let [nr (inc (* r (FastMath/exp (- alpha x))))]
          (if-not (seq rst)
            (+ (FastMath/log nr) x)
            (recur rst nr x)))))))

(defn xlogx
  "Computes `x * log(x)`, with the convention `0 * log(0) = 0` (the limiting value, rather than `##NaN`).

  Parameters:

  - `x` (double): value; must be non-negative for a real result.

  Returns the result as a double.

  See also [[xlogy]]."
  ^double [^double x]
  (if (zero? x) 0.0 (* x (FastMath/log x))))

(defn xlogy
  "Computes `x * log(y)`, with the convention `0 * log(y) = 0` (rather than `##NaN`) whenever `x` is zero and `y` is not `##NaN`.

  Parameters:

  - `x`, `y` (doubles): values; `y` must be non-negative for a real result.

  Returns the result as a double.

  See also [[xlogx]], [[xlog1py]]."
  ^double [^double x ^double y]
  (if (and (zero? x)
           (not (Double/isNaN y))) 0.0 (* x (log y))))

(defn xlog1py
  "Computes `x * log(1+y)`, with the convention `0 * log1p(y) = 0` (rather than `##NaN`) whenever `x` is zero and `y` is not `##NaN`.

  Parameters:

  - `x`, `y` (doubles): values; `y` must be greater than `-1` for a real result.

  Returns the result as a double.

  See also [[xlogy]], [[log1p]]."
  ^double [^double x ^double y]
  (if (and (zero? x)
           (not (Double/isNaN y))) 0.0 (* x (log1p y))))

(defn cloglog
  "Computes the complementary log-log function, `cloglog(x) = log(-log(1-x))`, used as a link function for binary/count models.

  Parameters:

  - `x` (double): probability-like value; must satisfy `0 < x < 1` for a real result.

  Returns the result as a double.

  See also [[loglog]]."
  {:inline (fn [x] `(FastMath/log (- (FastMath/log1p (- (double ~x))))))
   :inline-arities #{1}}
  ^double [^double x]
  (FastMath/log (- (FastMath/log1p (- x)))))

(defn loglog
  "Computes the log-log function, `loglog(x) = -log(-log(x))`, used as a link function for binary/count models.

  Parameters:

  - `x` (double): probability-like value; must satisfy `0 < x < 1` for a real result.

  Returns the result as a double.

  See also [[cloglog]]."
  {:inline (fn [x] `(- (FastMath/log (- (FastMath/log ~x)))))
   :inline-arities #{1}}
  ^double [^double x]
  (- (FastMath/log (- (FastMath/log x)))))

(defn xexpx
  "Computes `x * exp(x)`, with the convention that the result is `0.0` whenever `exp(x)` underflows to zero (rather than propagating a spurious `0.0 * x`).

  Parameters:

  - `x` (double): value.

  Returns the result as a double.

  See also [[xexpy]]."
  ^double [^double x]
  (let [expx (exp x)]
    (if (zero? expx) 0.0 (* x expx))))

(defn xexpy
  "Computes `x * exp(y)`, with the convention that the result is `0.0` whenever `exp(y)` underflows to zero and `x` is not `##NaN`.

  Parameters:

  - `x`, `y` (doubles): values.

  Returns the result as a double.

  See also [[xexpx]]."
  ^double [^double x ^double y]
  (let [expy (exp y)]
    (if (and (zero? expy)
             (not (Double/isNaN x))) 0.0 (* x expy))))

(defn cexpexp
  "Computes `1 - exp(-exp(x))`, the complementary Gumbel CDF form.

  Parameters:

  - `x` (double): value.

  Returns the result as a double, always in `[0, 1]`.

  See also [[expexp]]."
  {:inline (fn [x] `(- (FastMath/expm1 (- (FastMath/exp (double ~x))))))
   :inline-arities #{1}}
  ^double [^double x]
  (- (FastMath/expm1 (- (FastMath/exp x)))))

(defn expexp
  "Computes `exp(-exp(-x))`, the Gumbel CDF form.

  Parameters:

  - `x` (double): value.

  Returns the result as a double, always in `[0, 1]`.

  See also [[cexpexp]]."
  {:inline (fn [x] `(FastMath/exp (- (FastMath/exp (- (double ~x))))))
   :inline-arities #{1}}
  ^double [^double x]
  (FastMath/exp (- (FastMath/exp (- x)))))

;; Quick logarithm
(defn qlog
  "Fast and less accurate version of [[log]]."
  {:inline (fn [x] `(. FastMath (logQuick (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (logQuick x)))

;; Roots (square and cubic)
(defn sqrt
  "square root, sqrt(x)"
  {:inline (fn [x] `(. FastMath (sqrt (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (sqrt x)))

(defn cbrt
  "cubic root, cbrt(x)"
  {:inline (fn [x] `(. FastMath (cbrt (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (cbrt x)))

;; Quick version of exponential \\(e^x\\)
(defn qexp
  "Quick and less accurate version of [[exp]]."
  {:inline (fn [x] `(. FastMath (expQuick (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (expQuick x)))

;; Radians to degrees (and opposite) conversions
(def ^{:const true :tag 'double :doc "$\\frac{180}{\\pi}$"} rad-in-deg (/ 180.0 PI))
(def ^{:const true :tag 'double :doc "$\\frac{\\pi}{180}$"} deg-in-rad (/ PI 180.0))

(defn radians
  "Converts an angle from degrees to radians.

  Parameters:

  - `deg` (double): angle in degrees.

  Returns the angle in radians as a double.

  See also [[degrees]]."
  {:inline (fn [deg] `(* deg-in-rad ~deg))
   :inline-arities #{1}}
  ^double [^double deg] (* deg-in-rad deg))

(defn degrees
  "Converts an angle from radians to degrees.

  Parameters:

  - `rad` (double): angle in radians.

  Returns the angle in degrees as a double.

  See also [[radians]]."
  {:inline (fn [rad] `(* rad-in-deg ~rad))
   :inline-arities #{1}}
  ^double [^double rad] (* rad-in-deg rad))

;; Sinc
(defn sinc
  "Computes the normalized sinc function, `sinc(x) = sin(pi*x)/(pi*x)`, with the removable singularity at `x=0` handled explicitly.

  Parameters:

  - `v` (double): value to evaluate.

  Returns the result as a double, in `[-1, 1]`. Returns exactly `1.0` for `|pi*v|` below `1.0e-8` (avoiding division by a near-zero denominator).

  See also [[sin]]."
  ^double [^double v]
  (let [x (* PI (Math/abs v))]
    (if (< x 1.0e-8) 1.0
        (/ (FastMath/sin x) x))))

;;
(defn sigmoid
  "Computes the sigmoid (standard logistic) function, `sigmoid(x) = 1/(1+exp(-x))`.

  Parameters:

  - `x` (double): value to evaluate.

  Returns the result as a double, always in `(0, 1)`.

  See also [[logistic]] (alias), [[logit]] (inverse)."
  {:inline (fn [x] `(/ (inc (FastMath/exp (- (double ~x))))))
   :inline-arities #{1}}
  ^double [^double x]
  (/ (inc (FastMath/exp (- x)))))

(def ^{:doc "Alias for [[sigmoid]]"} logistic sigmoid)

(defn logit
  "Computes the logit function (log-odds), `logit(x) = log(x/(1-x))`, the inverse of [[sigmoid]]. Uses a numerically stable, `log1p`-based reformulation near `x=0.5` to avoid catastrophic cancellation.

  Parameters:

  - `x` (double): probability-like value; must satisfy `0 < x < 1` for a finite result.

  Returns the result as a double.

  See also [[sigmoid]] (inverse)."
  {:inline (fn [x] `(let [x# (double ~x)]
                     (if (< 0.3 x# 0.65)
                       (let [s# (* 2.0 (- x# 0.5))]
                         (- (FastMath/log1p s#)
                            (FastMath/log1p (- s#))))
                       (FastMath/log (/ x# (- 1.0 x#))))))
   :inline-arities #{1}}
  ^double [^double x]
  (if (< 0.3 x 0.65)
    (let [s (* 2.0 (- x 0.5))]
      (- (FastMath/log1p s) (FastMath/log1p (- s))))
    (FastMath/log (/ x (- 1.0 x)))))

(defn log2
  "Computes the base-2 logarithm of `x`.

  Parameters:

  - `x` (double): value to take the logarithm of; must be positive for a real result.

  Returns the result as a double.

  See also [[log]], [[logb]], [[LOG2E]]."
  {:inline (fn [x] `(* (FastMath/log (double ~x)) INV_LN2))
   :inline-arities #{1}}
  ^double [^double x]
  (* (FastMath/log x) INV_LN2))

;; \\(\log_b x\\)
(defn logb
  "Computes the logarithm of `x` with an explicit base `b`, `logb(b,x) = log(x)/log(b)`.

  Parameters:

  - `b` (double): logarithm base; must be positive and not equal to `1.0`.
  - `x` (double): value to take the logarithm of; must be positive for a real result.

  Returns the result as a double.

  See also [[log]] (2-arity form), [[log2]]."
  {:inline (fn [b x] `(/ (FastMath/log (double ~x)) (FastMath/log (double ~b))))
   :inline-arities #{2}}
  ^double [^double b ^double x]
  (/ (FastMath/log x) (FastMath/log b)))

(defn logcosh
  "Computes `log(cosh(x))`, using a numerically stable form (`|x| + log1pexp(-2|x|) - ln(2)`) to avoid overflow in `cosh(x)` for large `|x|`.

  Parameters:

  - `x` (double): value to evaluate.

  Returns the result as a double, always non-negative.

  See also [[cosh]], [[log1pexp]]."
  {:inline (fn [x] `(let [absx# (Math/abs (double ~x))]
                     (- (+ absx# (log1pexp (* -2.0 absx#))) LN2)))
   :inline-arities #{1}}
  ^double [^double x]
  (let [absx (Math/abs x)]
    (- (+ absx (log1pexp (* -2.0 absx))) LN2)))

;; \\(\log_2 e\\)
(def ^{:const true :tag 'double :doc "$\\log_{2}{\\mathrm{e}}$"} LOG2E (log2 E))

;; \\(\log_{10} e\\)
(def ^{:const true :tag 'double :doc "$\\log_{10}{\\mathrm{e}}$"} LOG10E (log10 E))

;; Powers (normal, quick)

;; using Math here due to some fastmath innacuracies

(defn pow
  "Power of a number"
  {:inline (fn [x exponent] `(. Math (pow (double ~x) (double ~exponent))))
   :inline-arities #{2}}
  ^double [^double x ^double exponent] (. Math (pow x exponent)))

(defn spow
  "Computes the symmetric power of `x`, `spow(x,e) = sign(x) * |x|^e`, preserving the sign of `x` (e.g. allows fractional or even exponents on negative bases without producing `##NaN`).

  Parameters:

  - `x` (double): base.
  - `exponent` (double): power to raise `|x|` to.

  Returns the result as a double, with the sign of `x` (or `0.0` when `x` is `0.0`).

  See also [[pow]], [[qpow]]."
  {:inline (fn [x exponent] `(let [v# (double ~x)]
                              (* (sgn v#) (. Math (pow (abs v#) (double ~exponent))))))
   :inline-arities #{2}}
  ^double [^double x ^double exponent] (* (sgn x) (. Math (pow (Math/abs x) exponent))))

(defn qpow
  "Fast and less accurate version of [[pow]]."
  {:inline (fn [x exponent] `(. FastMath (powQuick (double ~x) (double ~exponent))))
   :inline-arities #{2}}
  ^double [^double x ^double exponent] (. FastMath (powQuick x exponent)))

(defn fpow
  "Fast version of pow where exponent is integer."
  {:inline (fn [x exponent] `(. FastMath (powFast (double ~x) (long ~exponent))))
   :inline-arities #{2}}
  ^double [^double x ^long exponent] (. FastMath (powFast x exponent)))

(defn mpow
  "Calculates modular exponentiation, that is `x` raised to the power `e`, reduced modulo `m`.

  Uses binary (square-and-multiply) exponentiation, so it runs efficiently even for large exponents.

  Parameters:

  - `x` (long): Base.
  - `e` (long): Exponent, must be non-negative.
  - `m` (long): Modulus, must be positive.

  Returns the result as a long in the range `[0, m-1]`. Returns `0` when `m` is `1`.

  See also [[fpow]], [[pow]]."
  ^long [^long x ^long e ^long m]
  (if (one? m)
    0
    (loop [r (long 1)
           b (mod x m)
           e e]
      (if (zero? e)
        r
        (recur (if (odd? e) (mod (* r b) m) r)
               (mod (* b b) m)
               (>> e 1))))))

(defn tpow
  "Calculates the truncated power function of `x`.

  This is the truncated power basis function commonly used to construct splines, defined as `(x-shift)^exponent` when `x` is greater than `shift`, and `0.0` otherwise.

  Parameters:

  - `x` (double): Input value.
  - `exponent` (double): Power to raise `(x-shift)` to.
  - `shift` (double): Truncation point. Defaults to `0.0`.

  Returns `(x-shift)^exponent` as a double when `x > shift`, `0.0` when `x <= shift`."
  (^double [^double x ^double exponent] (tpow x exponent 0.0))
  (^double [^double x ^double exponent ^double shift]
   (let [diff (- x shift)]
     (if (pos? diff) (Math/pow diff exponent) 0.0))))

;;

(set! *unchecked-math* true)

(defn bernoulli
  "Calculates the `n`-th Bernoulli number, `B_n`.

  This implementation uses the `B_1 = +1/2` sign convention, giving the sequence `B_0 = 1`, `B_1 = 1/2`, `B_2 = 1/6`, `B_3 = 0`, `B_4 = -1/30`, and so on.

  Parameters:

  - `n` (long): Index of the Bernoulli number to compute, must be non-negative.

  Returns `B_n` as a double. All odd-indexed Bernoulli numbers above `B_1` (i.e. `B_3`, `B_5`, `B_7`, ...) are `0.0`.

  See also [[factorial]]."
  ^double [^long n]
  (loop [m 0
         j m
         buff (mapv (fn [^long m] (clojure.core// 1 m)) (range 1 (+ n 2)))]
    (cond
      (pos? j) (let [j- (dec j)]
                 (recur m j- (assoc buff j- (clojure.core/* j (clojure.core/- (buff j-) (buff j))))))
      (< m n) (let [m+ (inc m)]
                (recur m+ m+ buff))
      :else (buff 0))))

(set! *unchecked-math* :warn-on-boxed)

;;
(def ^:private factorial20-table [1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600
                                  6227020800 87178291200 1307674368000 20922789888000
                                  355687428096000 6402373705728000 121645100408832000
                                  2432902008176640000])

(defn factorial20
  "Looks up `n!` from a precomputed table for `n` in `[0, 20]`.

  Parameters:

  - `n` (long): index into the factorial table; must be in `[0, 20]` (`20!` is the largest factorial exactly representable as a `long`).

  Returns `n!` as a long.

  See also [[factorial]]."
  ^long [^long n]
  (factorial20-table n))

(defn factorial
  "Computes `x!`, using an exact table lookup for non-negative integers below `21`, and the gamma function (`exp(logGamma(x+1))`) otherwise.

  Parameters:

  - `x` (double): value to compute the factorial of.

  Returns the result as a double.

  See also [[factorial20]], [[inv-factorial]], [[log-factorial]], [[falling-factorial]], [[rising-factorial]]."
  ^double [^double x]
  (if (and (integer? x) (< x 21))
    (factorial20-table (long x))
    (exp (Gamma/logGamma (inc x)))))

(defn inv-factorial
  "Computes `1/x!`, the reciprocal of [[factorial]].

  Parameters:

  - `x` (double): value to compute the inverse factorial of.

  Returns the result as a double.

  See also [[factorial]]."
  ^double [^double x]
  (if (and (integer? x) (< x 21))
    (/ 1.0 (long (factorial20-table (long x))))
    (exp (- (Gamma/logGamma (inc x))))))

(defn stirling-factorial
  "Approximates `x!` using Stirling's asymptotic series with a 6-term correction, without relying on the gamma function.

  Parameters:

  - `x` (double): value to approximate the factorial of; accuracy improves for larger `x` (already near double-precision-level accuracy for `x >= 5`).

  Returns the approximate result as a double.

  See also [[factorial]], [[log-stirling-factorial]]."
  ^double [^double x]
  (let [x2 (* x x)
        x3 (* x x2)
        x5 (* x2 x3)
        x7 (* x2 x5)
        x9 (* x2 x7)]
    (* (sqrt (* TWO_PI x))
       (pow (/ x E) x)
       (exp (+ (/ 0.08333333333333333 x)
               (/ -0.002777777777777778 x3)
               (/ 7.936507936507937E-4 x5)
               (/ -5.952380952380953E-4 x7)
               (/ 8.417508417508417E-4 x9)
               (/ -0.0019175269175269176 (* x2 x9)))))))

(defn log-stirling-factorial
  "Approximates `log(x!)` using Stirling's asymptotic series with a 6-term correction, without relying on the gamma function.

  Parameters:

  - `x` (double): value to approximate the log-factorial of; accuracy improves for larger `x`.

  Returns the approximate result as a double.

  See also [[log-factorial]], [[stirling-factorial]]."
  ^double [^double x]
  (let [lx (log x)
        x2 (* x x)
        x3 (* x x2)
        x5 (* x2 x3)
        x7 (* x2 x5)
        x9 (* x2 x7)]
    (+ (* 0.5 (+ LOG_TWO_PI lx))
       (* x (dec lx))
       (+ (/ 0.08333333333333333 x)
          (/ -0.002777777777777778 x3)
          (/ 7.936507936507937E-4 x5)
          (/ -5.952380952380953E-4 x7)
          (/ 8.417508417508417E-4 x9)
          (/ -0.0019175269175269176 (* x2 x9))))))

(defn log-factorial
  "Computes `log(x!)`, equivalent to `log-gamma(x+1)`.

  Parameters:

  - `x` (double): value to compute the log-factorial of.

  Returns the result as a double.

  See also [[factorial]], [[log-stirling-factorial]], [[log-combinations]]."
  {:inline (fn [x] `(Gamma/logGamma (inc (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (Gamma/logGamma (inc x)))

(defn falling-factorial-int
  "Computes the falling (descending) factorial of `x` to the power of an integer `n`, `x*(x-1)*...*(x-n+1)`.

  Parameters:

  - `n` (long): number of descending terms; may be negative, in which case the result is the reciprocal of the falling factorial with `n` negated and `x` shifted by `n`.
  - `x` (double): starting value.

  Returns the falling factorial as a double.

  See also [[falling-factorial]], [[rising-factorial-int]]."
  ^double [^long n ^double x]
  (if (not-neg? n)
    (loop [i (long 0)
           v 1.0]
      (if (== i n) v
          (recur (inc i) (* v (- x i)))))
    (/ (falling-factorial-int (long-sub n) (- x n)))))

(defn falling-factorial
  "Calculates the falling (descending) factorial of `x` to the power `n`.

  Defined as the product `x*(x-1)*(x-2)*...*(x-n+1)` for a non-negative integer `n`. When `n` is not a non-negative integer, the definition is extended using the gamma function.

  Parameters:

  - `n` (double): Number of descending terms.
  - `x` (double): Starting value.

  Returns the falling factorial as a double.

  See also [[rising-factorial]], [[falling-factorial-int]], [[factorial]], [[combinations]]."
  ^double [^double n ^double x]
  (if (integer? n)
    (falling-factorial-int (long n) x)
    (let [x+ (inc x)]
      (/ (Gamma/gamma x+)
         (Gamma/gamma (- x+ n))))))

(defn rising-factorial-int
  "Computes the rising (ascending, Pochhammer) factorial of `x` to the power of an integer `n`, `x*(x+1)*...*(x+n-1)`.

  Parameters:

  - `n` (long): number of ascending terms; may be negative, in which case the result is the reciprocal of the rising factorial with `n` negated and `x` shifted by `n`.
  - `x` (double): starting value.

  Returns the rising factorial as a double.

  See also [[rising-factorial]], [[falling-factorial-int]]."
  ^double [^long n ^double x]
  (if (not-neg? n)
    (loop [i (long 0)
           v 1.0]
      (if (== i n) v
          (recur (inc i) (* v (+ x i)))))
    (/ (rising-factorial-int (long-sub n) (+ x n)))))

(defn rising-factorial
  "Calculates the rising (ascending) factorial of `x` to the power `n`, also known as the Pochhammer symbol.

  Defined as the product `x*(x+1)*(x+2)*...*(x+n-1)` for a non-negative integer `n`. When `n` is not a non-negative integer, the definition is extended using the gamma function.

  Parameters:

  - `n` (double): Number of ascending terms.
  - `x` (double): Starting value.

  Returns the rising factorial as a double.

  See also [[falling-factorial]], [[rising-factorial-int]], [[factorial]], [[combinations]]."
  ^double [^double n ^double x]
  (if (integer? n)
    (rising-factorial-int (long n) x)
    (/ (Gamma/gamma (+ x n))
       (Gamma/gamma x))))

(defn combinations
  "Computes the binomial coefficient, `n choose k`.

  Parameters:

  - `n` (long): total number of items.
  - `k` (long): number of items to choose.

  Returns the binomial coefficient as a double, `0.0` when `k` is negative or greater than `n`. Uses a direct iterative product for `k < 30`, and a log-beta-based formula for larger `k` to avoid overflow.

  See also [[log-combinations]], [[falling-factorial]]."
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
                    (Beta/logBeta (inc (- n k)) (inc k)))))))

(defn log-combinations
  "Computes the logarithm of the binomial coefficient, `log(n choose k)`.

  Parameters:

  - `n` (long): total number of items.
  - `k` (long): number of items to choose.

  Returns the result as a double, `##-Inf` when `k` is negative or greater than `n`, `0.0` when `k` is `0` or `k` equals `n`.

  See also [[combinations]], [[log-factorial]]."
  ^double [^long n ^long k]
  (let [k (min k (- n k))]
    (cond
      (neg? k) ##-Inf
      (zero? k) 0.0
      (one? k) (ln n)
      (< n k) ##-Inf
      (== n k) 0.0
      :else (- (- (ln (inc n)))
               (Beta/logBeta (inc (- n k)) (inc k))))))

;; Square and cubic
(defn sq
  "Computes `x^2` (`x*x`).

  Parameters:

  - `x` (double): value to square.

  Returns the result as a double.

  See also [[pow2]] (identical implementation), [[cb]], [[pow]]."
  {:inline (fn [x] `(let [x# (double ~x)] (* x# x#)))
   :inline-arities #{1}}
  ^double [^double x] (* x x))

(defn pow2
  "Computes `x^2` (`x*x`). Identical implementation to [[sq]].

  Parameters:

  - `x` (double): value to square.

  Returns the result as a double.

  See also [[sq]], [[pow3]], [[pow]]."
  {:inline (fn [x] `(let [x# (double ~x)] (* x# x#)))
   :inline-arities #{1}}
  ^double [^double x] (* x x))

(defn cb
  "Computes `x^3` (`x*x*x`).

  Parameters:

  - `x` (double): value to cube.

  Returns the result as a double.

  See also [[pow3]] (identical implementation), [[sq]], [[pow]]."
  {:inline (fn [x] `(let [x# (double ~x)] (* x# x# x#)))
   :inline-arities #{1}}
  ^double [^double x] (* x x x))

(defn pow3
  "Computes `x^3` (`x*x*x`). Identical implementation to [[cb]].

  Parameters:

  - `x` (double): value to cube.

  Returns the result as a double.

  See also [[cb]], [[pow2]], [[pow]]."
  {:inline (fn [x] `(let [x# (double ~x)] (* x# x# x#)))
   :inline-arities #{1}}
  ^double [^double x] (* x x x))

(defn pow10
  "Computes `x^10`.

  Parameters:

  - `x` (double): value to raise to the 10th power.

  Returns the result as a double, computed as `((x^3)^3)*x` to minimize the number of multiplications.

  See also [[pow]], [[pow3]]."
  {:inline (fn [x] `(let [x# (double ~x)
                         v# (* x# x# x#)]
                     (* v# v# v# x#)))
   :inline-arities #{1}}
  ^double [^double x] (let [v (* x x x)] (* v v v x)))

(defn safe-sqrt
  "Computes `sqrt(x)`, returning `0.0` instead of `##NaN` for negative `x`.

  Parameters:

  - `value` (double): value to take the square root of.

  Returns the result as a double, `0.0` when `value` is negative.

  See also [[sqrt]], [[qsqrt]]."
  ^double [^double value]
  (if (neg? value) 0.0 (FastMath/sqrt value)))

(defn qsqrt
  "Approximated [[sqrt]] using binary operations with error `1.0E-2`."
  {:inline (fn [x] `(. FastMath (sqrtQuick (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (sqrtQuick x)))

(defn rqsqrt
  "Reciprocal of [[qsqrt]]. Quick and less accurate."
  {:inline (fn [x] `(. FastMath (invSqrtQuick (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (invSqrtQuick x)))

(defn hypot
  "Calculates the Euclidean norm (distance from the origin) for 2 or 3 arguments.

  It uses a numerically stable algorithm to avoid intermediate overflow or underflow.

  See also [[hypot-sqrt]] for the direct calculation."
  {:inline (fn ([x y] `(. FastMath (hypot (double ~x) (double ~y))))
             ([x y z] `(. FastMath (hypot (double ~x) (double ~y) (double ~z)))))
   :inline-arities #{2 3}}
  (^double [^double x ^double y]
   (FastMath/hypot x y))
  (^double [^double x ^double y ^double z]
   (FastMath/hypot x y z)))

(defn hypot-sqrt
  "Calculates the Euclidean norm (distance from the origin) using a direct sqrt of sum of squares.

  Note: This method can be less numerically stable than [[hypot]] for inputs with vastly different magnitudes."
  {:inline (fn ([x y] `(. FastMath (sqrt (+ (sq ~x) (sq ~y)))))
             ([x y z] `(. FastMath (sqrt (+ (sq ~x) (sq ~y) (sq ~z))))))
   :inline-arities #{2 3}}
  (^double [^double x ^double y]
   (FastMath/sqrt (+ (* x x) (* y y))))
  (^double [^double x ^double y ^double z]
   (FastMath/sqrt (+ (* x x) (* y y) (* z z)))))

;; distance
(defn dist
  "Computes the Euclidean distance between two 2D points `(x1,y1)` and `(x2,y2)`.

  Parameters:

  - `[x1 y1]`, `[x2 y2]` (pairs of doubles): the two points, for the 2-arity form.
  - `x1`, `y1`, `x2`, `y2` (doubles): the two points' coordinates directly, for the 4-arity form.

  Returns the distance as a double.

  See also [[qdist]] (fast, less accurate version), [[hypot-sqrt]]."
  {:inline (fn [x1 y1 x2 y2] `(hypot-sqrt (- ~x2 ~x1) (- ~y2 ~y1)))
   :inline-arities #{4}}
  (^double [[^double x1 ^double y1] [^double x2 ^double y2]] (dist x1 y1 x2 y2))
  (^double [^double x1 ^double y1 ^double x2 ^double y2]
   (FastMath/sqrt (+ (sq (- x2 x1)) (sq (- y2 y1))))))

(defn qdist
  "Computes the Euclidean distance between two 2D points, using [[qsqrt]] instead of `sqrt` for a faster, less accurate result.

  Parameters:

  - `[x1 y1]`, `[x2 y2]` (pairs of doubles): the two points, for the 2-arity form.
  - `x1`, `y1`, `x2`, `y2` (doubles): the two points' coordinates directly, for the 4-arity form.

  Returns the approximate distance as a double.

  See also [[dist]], [[qsqrt]]."
  {:inline (fn [x1 y1 x2 y2] `(. FastMath (sqrtQuick (+ (sq (- ~x2 ~x1)) (sq (- ~y2 ~y1))))))
   :inline-arities #{4}}
  (^double [[^double x1 ^double y1] [^double x2 ^double y2]] (qdist x1 y1 x2 y2))
  (^double [^double x1 ^double y1 ^double x2 ^double y2]
   (FastMath/sqrtQuick (+ (sq (- x2 x1)) (sq (- y2 y1))))))

;; Rounding functions
(defn floor
  "Rounds `x` towards negative infinity.

  Parameters:

  - `x` (double): the number to round.
  - `scale` (double, optional): if given, `x` is divided by `scale`, floored, then multiplied back by `scale` -- i.e. rounds to the nearest multiple of `scale` towards negative infinity.

  Returns the rounded value as a double.

  See also [[ceil]], [[round]], [[qfloor]]."
  {:inline (fn ([x] `(. FastMath (floor (double ~x))))
             ([x scale] `(* (. FastMath (floor (double (/ ~x ~scale)))) ~scale)))
   :inline-arities #{1 2}}
  (^double [^double x] (FastMath/floor x))
  (^double [^double x ^double scale] (* (FastMath/floor (/ x scale)) scale)))

(defn ceil
  "Rounds `x` towards positive infinity.

  Parameters:

  - `x` (double): the number to round.
  - `scale` (double, optional): if given, `x` is divided by `scale`, ceiled, then multiplied back by `scale` -- i.e. rounds to the nearest multiple of `scale` towards positive infinity.

  Returns the rounded value as a double.

  See also [[floor]], [[round]], [[qceil]]."
  {:inline (fn ([x] `(. FastMath (ceil (double ~x))))
             ([x scale] `(* (. FastMath (ceil (double (/ ~x ~scale)))) ~scale)))
   :inline-arities #{1 2}}
  (^double [^double x] (FastMath/ceil x))
  (^double [^double x ^double scale] (* (FastMath/ceil (/ x scale)) scale)))

(defn round
  "Round to a `long` value. See: [[rint]], [[qround]]."
  {:inline (fn [x] `(. FastMath (round (double ~x))))
   :inline-arities #{1}} 
  ^long [^double x] (FastMath/round x))

(defn rint
  "Rounds `x` to the nearest integer, returned as a double, using round-half-to-even.

  Parameters:

  - `x` (double): the number to round.
  - `scale` (double, optional): if given, `x` is divided by `scale`, rounded, then multiplied back by `scale` -- i.e. rounds to the nearest multiple of `scale`, ties rounding to even.

  Returns the rounded value as a double.

  See also [[round]], [[round-even]], [[qround]]."
  {:inline (fn ([x] `(. FastMath (rint (double ~x))))
             ([x scale] `(* (. FastMath (rint (double (/ ~x ~scale)))) ~scale)))
   :inline-arities #{1 2}}
  (^double [^double x] (FastMath/rint x))
  (^double [^double x ^double scale] (* (FastMath/rint (/ x scale)) scale)))

(defn round-even
  "Round evenly, IEEE / IEC rounding"
  {:inline (fn [x] `(. FastMath (roundEven (double ~x))))
   :inline-arities #{1}} 
  ^long [^double x] (FastMath/roundEven x))

(defn qfloor
  "Fast version of [[floor]]. Returns `long`."
  {:inline (fn [x] `(. PrimitiveMath (fastFloor (double ~x))))
   :inline-arities #{1}} 
  ^long [^double x] (PrimitiveMath/fastFloor x))

(defn qceil
  "Fast version of [[ceil]]. Returns `long`."
  {:inline (fn [x] `(. PrimitiveMath (fastCeil (double ~x))))
   :inline-arities #{1}} 
  ^long [^double x] (PrimitiveMath/fastCeil x))

(defn qround
  "Fast version of [[round]]. Returns `long`"
  {:inline (fn [x] `(. PrimitiveMath (fastRound (double ~x))))
   :inline-arities #{1}} 
  ^long [^double x] (PrimitiveMath/fastRound x))

(defn remainder
  "From `FastMath` doc: returns dividend - divisor * n,
  where n is the mathematical integer closest to dividend/divisor. Returned value in `[-|divisor|/2,|divisor|/2]`"
  {:inline (fn [dividend divisor] `(. FastMath (remainder (double ~dividend) (double ~divisor))))
   :inline-arities #{2}}
  ^double [^double dividend ^double divisor]
  (. FastMath (remainder dividend divisor)))

(defn abs
  "Absolute value."
  {:inline (fn [x] `(. Math (abs ~x)))
   :inline-arities #{1}}
  ^double [^double x] (Math/abs x))

(defn iabs
  "Absolute value, `long` version. See [[abs]]."
  {:inline (fn [x] `(let [m# (>> ~x 63)] (bit-xor (+ m# ~x) m#)))
   :inline-arities #{1}
   :deprecated "Use long-abs."}
  ^long [^long x] (let [m (>> x 63)] (bit-xor (long-add m x) m)))

(defn long-abs
  "Absolute value of a `long`, computed with a branch-free bitwise trick.

  Parameters:

  - `x` (long): the input value.

  Returns `x` if non-negative, `-x` otherwise, as a `long`. For `x` equal to `Long/MIN_VALUE`, returns `x` unchanged (its true absolute value overflows the `long` range), matching `Math/abs`'s behavior for `long`.

  See also [[abs]]."
  {:inline (fn [x] `(let [m# (>> ~x 63)] (bit-xor (+ m# ~x) m#)))
   :inline-arities #{1}}
  ^long [^long x] (let [m (>> x 63)] (bit-xor (long-add m x) m)))

(defn trunc
  "Truncates the fractional part of `v`, rounding towards zero.

  Parameters:

  - `v` (double): the number to truncate.

  Returns the integer part of `v` as a double, preserving sign: [[ceil]] for negative `v`, [[floor]] otherwise.

  See also [[itrunc]]."
  ^double [^double v] (if (neg? v) (ceil v) (floor v)))

(defn itrunc
  "Truncates the fractional part of `v`, rounding towards zero.

  Parameters:

  - `v` (double): the number to truncate.

  Returns the integer part of `v` as a `long`, preserving sign: [[qceil]] for negative `v`, [[qfloor]] otherwise.

  See also [[trunc]]."
  ^long [^double v] (if (neg? v) (qceil v) (qfloor v)))

;; return approximate value
(defn approx
  "Rounds `v` to a given number of decimal places, half-up.

  Parameters:

  - `v` (double): the number to round.
  - `digits` (long, optional, default `2`): number of decimal places to keep. `0` rounds to the nearest integer; negative values round to the nearest multiple of a power of ten (e.g. `-1` rounds to the nearest 10).

  Returns the rounded value as a double.

  See also [[approx-eq]], [[delta-eq]]."
  (^double [^double v] (Precision/round v (int 2)))
  (^double [^double v ^long digits] (Precision/round v (int digits))))

(defn approx-eq
  "Checks whether `a` and `b` round to the same value at a given number of decimal places.

  Parameters:

  - `a`, `b` (doubles): the two numbers to compare.
  - `digits` (long, optional, default `2`): number of decimal places used for rounding (see [[approx]]).

  Returns `true` if `(== a b)`, or if [[approx]] of `a` and `b` are equal at the given number of digits. Returns `false` otherwise.

  This equality check can be inaccurate near a rounding boundary: `1.004999` and `1.005001` differ by only `0.000002` but round to `1.0` and `1.01` respectively at 2 digits, so `approx-eq` reports them as unequal -- prefer [[delta-eq]] for a tolerance-based comparison.

  See also [[approx]], [[delta-eq]]."
  ([^double a ^double b] (or (== (approx a) (approx b)) (== a b)))
  ([^double a ^double b ^long digits] (or (== (approx a digits)
                                              (approx b digits))
                                          (== a b))))

(defn delta-eq
  "Checks if two floating-point numbers `a` and `b` are approximately equal within given tolerances.

  The check returns true if `abs(a - b)` is less than a combined tolerance, or if `(== a b)`.

  - 2-arity `(delta-eq a b)`: Uses a default absolute tolerance of `1.0e-6`.
  - 3-arity `(delta-eq a b abs-tol)`: Uses the provided absolute tolerance `abs-tol`.
  - 4-arity `(delta-eq a b abs-tol rel-tol)`: Uses both absolute and relative tolerances. The combined tolerance is `max(abs-tol, rel-tol * max(abs(a), abs(b)))`.

  This function is useful for comparing floating-point numbers where exact equality checks (`==`) may fail due to precision issues."
  ([^double a ^double b] (delta-eq a b 1.0e-6))
  ([^double a ^double b ^double accuracy]
   (or (< (Math/abs (- a b)) accuracy) (== a b)))
  ([^double a ^double b ^double abs-tol ^double rel-tol]
   (or (< (Math/abs (- a b)) (max abs-tol (* rel-tol (max (Math/abs a) (Math/abs b))))) (== a b))))

(def ^{:doc "Alias for [[approx-eq]]"} approx= approx-eq)
(def ^{:doc "Alias for [[delta-eq]]"} delta= delta-eq)

(defn near-zero?
  "Checks whether `x` is close enough to zero to be treated as zero.

  Parameters:

  - `x` (double): the value to check.
  - `abs-tol` (double, optional, default `1.0e-6`): absolute tolerance.
  - `rel-tol` (double, optional, default `0.0`): relative tolerance, scaled by `(abs x)` itself.

  Returns `true` if `(abs x)` is less than `(max abs-tol (* rel-tol (abs x)))`, `false` otherwise.

  Since the relative-tolerance term is scaled by `x`'s own magnitude, it can only widen the threshold beyond `abs-tol` when `rel-tol` is `1.0` or greater -- for any smaller (i.e. any realistic) `rel-tol`, the 3-arity form behaves identically to the 2-arity form and `rel-tol` has no effect. This mirrors the standard convention that a relative tolerance compared against a target of zero contributes nothing.

  See also [[delta-eq]]."
  ([^double x] (near-zero? x 1.0e-6))
  ([^double x ^double abs-tol] (< (Math/abs x) abs-tol))
  ([^double x ^double abs-tol ^double rel-tol]
   (let [ax (Math/abs x)] (< ax (max abs-tol (* rel-tol ax)) ))))

(defn frac
  "Unsigned fractional part of `v`.

  Parameters:

  - `v` (double): the input value.

  Returns `(abs (- v (long v)))`, i.e. the magnitude of the part of `v` remaining after truncation towards zero. Always in the range `[0.0, 1.0)`.

  See also [[sfrac]] for the signed version."
  ^double [^double v] (Math/abs (- v (unchecked-long v))))

(defn sfrac
  "Signed fractional part of `v`.

  Parameters:

  - `v` (double): the input value.

  Returns `(- v (trunc v))`, i.e. the part of `v` remaining after truncation towards zero, keeping `v`'s sign. Always in the range `(-1.0, 1.0)`.

  See also [[frac]] for the unsigned version."
  ^double [^double v] (- v (trunc v)))

;; Find power of 2 exponent for double number where  
;; \\(2^(n-1)\leq x\leq 2^n\\)  
;; where n-1 is result of `low-2-exp` and n is result of `high-2-exp`
;; `(low-2-exp TWO_PI) => 2` \\(2^2\eq 4\leq 6.28\\)  
;; `(high-2-exp TWO_PI) => 3` \\(6.28\leq 2^3\eq 8\\)
(defn low-2-exp
  "Finds the greatest integer `n` such that `2^n` does not exceed the absolute value of `x`.

  Equivalent to the floor of the base-2 logarithm of `|x|`. Together with [[high-2-exp]] it brackets `|x|` between two consecutive powers of two.

  Parameters:

  - `x` (double): Input value.

  Returns the exponent `n` as a long.

  See also [[high-2-exp]], [[low-exp]]."
  ^long [^double x] (-> x Math/abs log2 floor unchecked-long))

(defn high-2-exp
  "Finds the smallest integer `n` such that `2^n` is not smaller than the absolute value of `x`.

  Equivalent to the ceiling of the base-2 logarithm of `|x|`. Together with [[low-2-exp]] it brackets `|x|` between two consecutive powers of two.

  Parameters:

  - `x` (double): Input value.

  Returns the exponent `n` as a long.

  See also [[low-2-exp]], [[high-exp]]."
  ^long [^double x] (-> x Math/abs log2 ceil unchecked-long))

(defn low-exp
  "Finds the greatest integer `n` such that `b^n` does not exceed the absolute value of `x`.

  Equivalent to the floor of the base-`b` logarithm of `|x|`. Together with [[high-exp]] it brackets `|x|` between two consecutive powers of `b`.

  Parameters:

  - `b` (double): Base, must be positive and not equal to `1.0`.
  - `x` (double): Input value.

  Returns the exponent `n` as a long.

  See also [[high-exp]], [[low-2-exp]]."
  ^long [^double b ^double x] (->> x Math/abs (logb b) floor unchecked-long))

(defn high-exp
  "Finds the smallest integer `n` such that `b^n` is not smaller than the absolute value of `x`.

  Equivalent to the ceiling of the base-`b` logarithm of `|x|`. Together with [[low-exp]] it brackets `|x|` between two consecutive powers of `b`.

  Parameters:

  - `b` (double): Base, must be positive and not equal to `1.0`.
  - `x` (double): Input value.

  Returns the exponent `n` as a long.

  See also [[low-exp]], [[high-2-exp]]."
  ^long [^double b ^double x] (->> x Math/abs (logb b) ceil unchecked-long))

(defn power-of-two?
  "Checks if `v` is a power of two, v=2^p for some p. Only for positive values."
  [^long v]
  (and (pos? v) (zero? (bit-and v (long-dec v)))))

(defn round-up-pow2
  "Rounds a positive `long` integer up to the smallest power of 2 greater than or equal to the input value."
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
  "Returns the closest representable double value greater than `v`, moving toward positive infinity.

  This steps `v` by one ulp (unit in the last place). The optional `delta` argument repeats the step `delta` times, moving further away from `v`.

  Parameters:

  - `v` (double): Starting value.
  - `delta` (long): Number of steps to take, must be non-negative. Defaults to `1`.

  Returns the resulting double value.

  Returns `NaN` when `v` is `NaN`. Returns `##Inf` when `v` is `##Inf`.

  See also [[prev-double]]."
  {:inline (fn [v] `(. FastMath (nextUp (double ~v))))
   :inline-arities #{1}}
  (^double [^double v]
   (FastMath/nextUp v))
  (^double [^double v ^long delta]
   (nth (iterate next-double v) delta)))

(defn prev-double
  "Returns the closest representable double value smaller than `v`, moving toward negative infinity.

  This steps `v` by one ulp (unit in the last place). The optional `delta` argument repeats the step `delta` times, moving further away from `v`.

  Parameters:

  - `v` (double): Starting value.
  - `delta` (long): Number of steps to take, must be non-negative. Defaults to `1`.

  Returns the resulting double value.

  Returns `NaN` when `v` is `NaN`. Returns `##-Inf` when `v` is `##-Inf`.

  See also [[next-double]]."
  {:inline (fn [v] `(. FastMath (nextDown (double ~v))))
   :inline-arities #{1}}
  (^double [^double v]
   (FastMath/nextDown v))
  (^double [^double v ^long delta]
   (nth (iterate prev-double v) delta)))

(defn double-high-bits
  "Returns high word from double as bits"
  {:inline (fn [v] `(bit-and (>>> (Double/doubleToRawLongBits (double ~v)) 32) 0xffffffff))
   :inline-arities #{1}}
  ^long [^double v]
  (bit-and (>>> (Double/doubleToRawLongBits v) 32) 0xffffffff))

(defn double-low-bits
  "Returns low word from double as bits"
  {:inline (fn [v] `(bit-and (Double/doubleToRawLongBits (double ~v)) 0xffffffff))
   :inline-arities #{1}}
  ^long [^double v]
  (bit-and (Double/doubleToRawLongBits v) 0xffffffff))

(defn double-bits
  "Returns double as 64-bits (long)"
  {:inline (fn [v] `(. Double (doubleToRawLongBits (double ~v))))
   :inline-arities #{1}}
  ^long [^double v]
  (Double/doubleToRawLongBits v))

(defn bits->double
  "Convert 64 bits to double"
  {:inline (fn [v] `(. Double (longBitsToDouble (double ~v))))
   :inline-arities #{1}}
  ^double [^long v]
  (Double/longBitsToDouble v))

(defn double-exponent
  "Extract exponent information from double"
  {:inline (fn [v] `(. FastMath (getExponent (double ~v))))
   :inline-arities #{1}}
  ^long [^double v]
  (FastMath/getExponent v))

(defn double-significand
  "Extract significand from double"
  {:inline (fn [v] `(bit-and (Double/doubleToRawLongBits (double ~v)) 4503599627370495))
   :inline-arities #{1}}
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

(defn leading-zero-bits
  "Leading zero bits"
  {:inline (fn [x] `(. Long (numberOfLeadingZeros (long ~x))))
   :inline-arities #{1}}
  ^long [^long v]
  (Long/numberOfLeadingZeros v))

(defn trailing-zero-bits
  "Trailing zero bits"
  {:inline (fn [x] `(. Long (numberOfTrailingZeros (long ~x))))
   :inline-arities #{1}}
  ^long [^long v]
  (Long/numberOfTrailingZeros v))

(defn most-significant-bit
  "Returns the most significant bit position. Can be treated as floor(log2(v)) for positive, integer `v`.

  Returns `-1` for `0` and `63` for negative numbers."
  {:inline (fn [x] `(PrimitiveMath/subtract 63 (. Long (numberOfLeadingZeros (long ~x)))))
   :inline-arities #{1}}
  ^long [^long v]
  (- 63 (Long/numberOfLeadingZeros v)))

(defn least-significant-bit
  "Returns the least significant bit position.

  Returns `-1` for `0`"
  ^long [^long v]
  (if (zero? v)
    -1
    (Long/numberOfTrailingZeros v)))

(defn ulp
  "Unit in the Last Place, distance between next value larger than `x` and `x`"
  {:inline (fn [x] `(. FastMath (ulp (double ~x))))
   :inline-arities #{1}}
  ^double [^double x]  (FastMath/ulp x))

;; More constants

(def ^{:const true :tag 'double :doc "Value of 0x1.fffffffffffffp-1d = 0.(9)"}
  double-one-minus-epsilon (Double/parseDouble "0x1.fffffffffffffp-1d"))

;; \\(\sqrt{2}\\)
(def ^{:const true :tag 'double :doc "Value of $\\sqrt{2}$"} SQRT2 (sqrt 2.0))
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\sqrt{2}}{2}$"} SQRT2_2 (* 0.5 SQRT2))

;; \\(\sqrt{3}\\)
(def ^{:const true :tag 'double :doc "Value of $\\sqrt{3}$"} SQRT3 (sqrt 3.0))
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\sqrt{3}}{2}$"} SQRT3_2 (* 0.5 (sqrt 3.0)))
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\sqrt{3}}{3}$"} SQRT3_3 (/ (sqrt 3.0) 3.0))
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\sqrt{3}}{4}$"} SQRT3_4 (/ (sqrt 3.0) 4.0))

;; \\(\sqrt{5}\\)
(def ^{:const true :tag 'double :doc "Value of $\\sqrt{5}$"} SQRT5 (sqrt 5.0))

;; \\(\sqrt{\pi}\\)
(def ^{:const true :tag 'double :doc "Value of $\\sqrt{\\pi}$"} SQRTPI (sqrt PI))
(def ^{:const true :tag 'double :doc "Value of $\\sqrt{2\\pi}$"} SQRT2PI (sqrt TWO_PI))
(def ^{:const true :tag 'double :doc "Value of $\\sqrt{\\frac{1}{2}\\pi}$"} SQRT_HALFPI (sqrt HALF_PI))

;; 
(def ^{:const true :tag 'double :doc "Golden ratio $\\phi$"} PHI (* (inc SQRT5) 0.5))
(def ^{:const true :tag 'double :doc "Silver ratio $\\delta_S$"} SILVER (inc SQRT2))

;; math.h predefined constants names
(def ^{:const true :tag 'double :doc "Value of $\\mathrm{e}$"} M_E E)
(def ^{:const true :tag 'double :doc "Value of $\\log_{2}{e}$"} M_LOG2E LOG2E)
(def ^{:const true :tag 'double :doc "Value of $\\log_{10}{e}$"} M_LOG10E LOG10E)
(def ^{:const true :tag 'double :doc "Value of $\\ln{2}$"} M_LN2 LN2)
(def ^{:const true :tag 'double :doc "Value of $\\ln{10}$"} M_LN10 LN10)
(def ^{:const true :tag 'double :doc "Value of $\\pi$"} M_PI PI)
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\pi}{2}$"} M_PI_2 HALF_PI)
(def ^{:const true :tag 'double :doc "Value of $\\frac{\\pi}{4}$"} M_PI_4 QUARTER_PI)
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{\\pi}$"} M_1_PI (/ PI))
(def ^{:const true :tag 'double :doc "Value of $\\frac{2}{\\pi}$"} M_2_PI (/ 2.0 PI))
(def ^{:const true :tag 'double :doc "Value of $\\frac{2}{\\sqrt\\pi}$"} M_2_SQRTPI (/ 2.0 SQRTPI))
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{\\sqrt{2\\pi}}$"} INV_SQRT2PI (/ 1.0 SQRT2PI))
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{\\sqrt\\pi}$"} INV_SQRTPI (/ 1.0 SQRTPI))
(def ^{:const true :tag 'double :doc "Value of $\\sqrt{\\frac{2}{\\pi}}$"} SQRT_2_PI (sqrt M_2_PI))
(def ^{:const true :tag 'double :doc "Value of $\\sqrt{2}$"} M_SQRT2 SQRT2)
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{\\sqrt{2}}$"} M_SQRT1_2 (/ SQRT2))
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{\\sqrt{2}}$"} INV_SQRT_2 M_SQRT1_2)
(def ^{:const true :tag 'double :doc "Value of $2\\pi$"} M_TWOPI TWO_PI)
(def ^{:const true :tag 'double :doc "Value of $\\frac{3\\pi}{4}$"} M_3PI_4 (* PI 0.75))
(def ^{:const true :tag 'double :doc "Value of $\\sqrt\\pi$"} M_SQRT_PI SQRTPI)
(def ^{:const true :tag 'double :doc "Value of $\\sqrt{3}$"} M_SQRT3 SQRT3)
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{\\ln{10}}$"} M_IVLN10 (/ LN10))
(def ^{:const true :tag 'double :doc "Value of $\\ln{2}$"} M_LOG2_E LN2)
(def ^{:const true :tag 'double :doc "Value of $\\frac{1}{\\ln{2}}$"} M_INVLN2 (/ LN2))

(defmacro constrain
  "Clamp `value` to the range `[mn,mx]`."
  [value mn mx]
  `(max (min ~value ~mx) ~mn))

(defn norm
  "Linearly maps `v` from one range to another.

  With two arities: the 3-argument form normalizes `v` from the range `[start,stop]` to `[0,1]`. The 5-argument form maps `v` from the range `[start1,stop1]` to the range `[start2,stop2]`. Both forms are equivalent to a linear interpolation and extrapolate when `v` lies outside the source range.

  Parameters:

  - `v` (double): Value to map.
  - `start`, `stop` (doubles): Source range for the 3-arity normalization to `[0,1]`.
  - `start1`, `stop1` (doubles): Source range for the 5-arity mapping.
  - `start2`, `stop2` (doubles): Target range for the 5-arity mapping.

  Returns the mapped value as a double.

  See also [[mnorm]] (macro version), [[make-norm]] (returns a reusable mapping function), [[constrain]] (clamps a value to a range)."
  {:inline (fn
             ([v start stop] `(PrimitiveMath/norm (double ~v) (double ~start) (double ~stop)))
             ([v start1 stop1 start2 stop2] `(PrimitiveMath/norm (double ~v)
                                                                 (double ~start1) (double ~stop1)
                                                                 (double ~start2) (double ~stop2))))
   :inline-arities #{3 5}}
  (^double [^double v ^double start ^double stop] ;; norm
   (PrimitiveMath/norm v start stop))
  ([v start1 stop1 start2 stop2] ;; map
   (PrimitiveMath/norm (double v) (double start1) (double stop1) (double start2) (double stop2))))

(defmacro mnorm
  "Macro version of [[norm]]."
  ([v start stop]
   `(PrimitiveMath/norm (double ~v) (double ~start) (double ~stop)))
  ([v start1 stop1 start2 stop2]
   `(PrimitiveMath/norm (double ~v) (double ~start1) (double ~stop1) (double ~start2) (double ~stop2))))

(defn make-norm
  "Creates a reusable function that linearly maps values from a fixed source range.

  This is a partially applied version of [[norm]], useful when the same source range (and optionally the same target range) is reused for many values, avoiding repeated range arguments.

  Parameters:

  - `start`, `stop` (doubles): Fixed source range.
  - `dstart`, `dstop` (doubles): Fixed target range. When omitted, the returned function accepts them on every call.

  Returns a function of a `double` value `v` to a `double`. When `dstart` and `dstop` are not provided here, the returned function has arity `[v dstart dstop]`; when they are provided, the returned function has arity `[v]`.

  See also [[norm]], [[mnorm]]."
  ([^double start ^double stop]
   (fn ^double [^double v ^double dstart ^double dstop]
     (PrimitiveMath/norm v start stop dstart dstop)))
  ([^double start ^double stop ^double dstart ^double dstop]
   (fn ^double [^double v]
     (PrimitiveMath/norm v start stop dstart dstop))))

(defn cnorm
  "Constrained version of norm. Result of [[norm]] is applied to [[constrain]] to `[0,1]` or `[start2,stop2]` ranges."
  {:inline (fn
             ([v start stop]
              `(constrain (PrimitiveMath/norm (double ~v) (double ~start) (double ~stop)) 0.0 1.0))
             ([v start1 stop1 start2 stop2]
              `(let [st2# (double ~start2)
                     sp2# (double ~stop2)]
                 (constrain (PrimitiveMath/norm (double ~v)
                                                (double ~start1) (double ~stop1)
                                                st2# sp2#) st2# sp2#))))
   :inline-arities #{3 5}}
  ([v start1 stop1 start2 stop2]
   (constrain (PrimitiveMath/norm v start1 stop1 start2 stop2) (double start2) (double stop2)))
  (^double [v ^double start ^double stop]
   (constrain (PrimitiveMath/norm v start stop) 0.0 1.0)))

;;; Interpolation functions

;; Linear interpolation between `start` and `stop`.
(defn lerp
  "Linear interpolation between `start` and `stop` for amount `t`. See also [[mlerp]], [[cos-interpolation]], [[quad-interpolation]] or [[smooth-interpolation]]."
  {:inline (fn [start stop t] `(let [s# (double ~start)]
                                (+ s# (* (double ~t) (- (double ~stop) s#)))))
   :inline-arities #{3}}
  ^double [^double start ^double stop ^double t]
  (+ start (* t (- stop start))))

(defmacro mlerp
  "[[lerp]] as macro. For inline code. See also [[lerp]], [[cos-interpolation]], [[quad-interpolation]] or [[smooth-interpolation]]."
  [start stop t]
  `(+ (double ~start) (* (double ~t) (- (double ~stop) (double ~start)))))

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
  "Wraps `value` cyclically into the range `[start,stop)`, similar to `openFrameworks`' `ofWrap`.

  When `value` lies outside the range, it is shifted by an integer multiple of the range width until it falls back inside, effectively behaving like a floating-point modulo operation over the range. This is useful for wrapping angles, cyclic coordinates, or any periodic quantity. `start` and `stop` do not need to be ordered, the range is normalized internally.

  Parameters:

  - `[start stop]` (sequence of two doubles), `value` (double): Range provided as a pair, plus the value to wrap.
  - `start`, `stop`, `value` (doubles): Range boundaries and the value to wrap, given directly.

  Returns the wrapped value as a double, always within `[start,stop)`. Returns `stop` when `start` equals `stop`.

  See also [[constrain]] (clamping instead of wrapping), [[norm]]."
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
        ^double den (reduce + eaxs)]
    (reduce + (map (fn [^double x ^double eax]
                     (/ (* x eax) den)) xs eaxs))))

;; https://iquilezles.org/articles/smin/

(defmacro ^:private smooth-max-kernel [g a b k] `(- ~b (* ~k (~g (/ (- ~b ~a) ~k)))))

(defn- sm-sigmoid  ^double [^double x] (/ x (- 1.0 (exp2 (- x)))))

(defn- sm-circular ^double [^double x]
  (cond
    (> x 1.0) x
    (< x -1.0) 0.0
    :else (inc (* 0.5 (- x (sqrt (- 2.0 (* x x))))))))

(defn- sm-quadratic ^double [^double x]
  (cond
    (> x 1.0) x
    (< x -1.0) 0.0
    :else (* 0.25 (inc (* x (+ 2.0 x))))))

(defn- sm-cubic ^double [^double x]
  (cond
    (> x 1.0) x
    (< x -1.0) 0.0
    :else (* SIXTH (inc (- (* 3.0 x (inc x))
                           (Math/abs (* x x x)))))))

(defn- sm-quartic ^double [^double x]
  (cond
    (> x 1.0) x
    (< x -1.0) 0.0
    :else (let [x+ (inc x)]
            (* 0.0625 x+ x+ (- 3.0 (* x (- x 2.0)))))))


(defn smooth-max
  "Returns a smooth approximation of the maximum value over a sequence `xs`.

  Unlike a hard `max`, this function is differentiable everywhere and is
  controlled by the sharpness parameter `alpha`. As `alpha` increases toward
  infinity, the result converges to the true maximum. Negative `alpha` values
  produce a smooth minimum approximation (except for `:p-norm`).

  Parameters:

  - `xs` (sequence of numbers): Input values.
  - `alpha` (double): Sharpness parameter. Defaults to `2.0`. Larger positive
    values approximate the true maximum more closely; negative values approximate
    the minimum (except `:p-norm`).
  - `family` (keyword): Smoothing family to use. Defaults to `:lse`. Available families:
    - `:lse` - LogSumExp: numerically stable, globally smooth over all elements.
    - `:boltzmann` - Boltzmann operator (weighted average); best suited for small `alpha` values.
    - `:mellowmax` - like `:lse` but mean-normalised, independent of sequence length.
    - `:smu` - smooth maximum unit; pairwise reduction, epsilon `= 1/|alpha|` controls the rounding radius.
    - `:p-norm` - p-norm-based smooth absolute maximum; always returns non-negative values; does not act as smooth minimum for negative `alpha`.
    - `:sigmoid`, `:circular`, `:quadratic`, `:cubic`, `:quartic`, `:exponential` - pairwise reductions based on Inigo Quilez smooth minimum kernels for SDFs.

  Returns the smooth maximum as a double.

  See also [[smooth-max-kernel]]."
  (^double [xs] (smooth-max xs 2.0))
  (^double [xs ^double alpha] (smooth-max xs alpha :lse))
  (^double [xs ^double alpha family]
   (case family
     :lse (/ (logsumexp (scale-xs xs alpha)) alpha)
     :boltzmann (smooth-max-boltzmann xs alpha)
     :mellowmax (/ (- (logsumexp (scale-xs xs alpha))
                      (log (count xs))) alpha)
     :p-norm (pow (reduce + (map (fn [^double x]
                                   (pow (Math/abs x) alpha)) xs)) (/ alpha))
     :smu (let [epsilon (Math/abs (/ alpha))]
            (reduce (if (pos? alpha)
                      (fn [^double a ^double b]
                        (* 0.5 (+ a b (sqrt (+ (sq (- a b))
                                               epsilon)))))
                      (fn [^double a ^double b]
                        (* 0.5 (+ a b (- (sqrt (+ (sq (- a b))
                                                  epsilon))))))) xs))
     :sigmoid (let [k (* M_LN2 (- (/ alpha)))]
                (reduce (fn [^double a ^double b]
                          (smooth-max-kernel sm-sigmoid a b k)) xs))
     :circular (let [k (* 3.414213562373096 (- (/ alpha)))]
                 (reduce (fn [^double a ^double b]
                           (smooth-max-kernel sm-circular a b k)) xs))
     :quadratic (let [k (* 4.0 (- (/ alpha)))]
                  (reduce (fn [^double a ^double b]
                            (smooth-max-kernel sm-quadratic a b k)) xs))
     :cubic (let [k (* 6.0 (- (/ alpha)))]
              (reduce (fn [^double a ^double b]
                        (smooth-max-kernel sm-cubic a b k)) xs))
     :quartic (let [k (* 5.333333333333333 (- (/ alpha)))]
                (reduce (fn [^double a ^double b]
                          (smooth-max-kernel sm-quartic a b k)) xs))
     :exponential (reduce (fn [^double a ^double b]
                            (/ (log2 (+ (exp2 (* a alpha)) (exp2 (* b alpha)))) alpha)) xs))))

;;

(defn nan?
  "Check if a number is a NaN"
  {:inline (fn [v] `(Double/isNaN (double ~v))) :inline-arities #{1}}
  [^double v]
  (Double/isNaN v))

(defn inf?
  "Check if a number is an infinite (positive or negative)."
  {:inline (fn [v] `(Double/isInfinite (double ~v))) :inline-arities #{1}}
  [^double v]
  (Double/isInfinite v))

(defn pos-inf?
  "Check if a number is positively infinite."
  {:inline (fn [v] `(== (double ~v) ##Inf)) :inline-arities #{1}}
  [^double v]
  (== v ##Inf))

(defn neg-inf?
  "Check if a number is negatively infinite."
  {:inline (fn [v] `(== (double ~v) ##-Inf)) :inline-arities #{1}}
  [^double v]
  (== v ##-Inf))

(defn invalid-double?
  "Check if a number is not finite double (NaN or ±Inf)."
  {:inline (fn [v] `(bool-not (Double/isFinite (double ~v)))) :inline-arities #{1}}
  [^double v]
  (bool-not (Double/isFinite v)))

(defn valid-double?
  "Check if a number is finite double."
  {:inline (fn [v] `(Double/isFinite (double ~v))) :inline-arities #{1}}
  [^double v]
  (Double/isFinite v))

(defn between?
  "Check if given number is within the range [x,y]."
  {:inline (fn [x y v] `(<= (double ~x) (double ~v) (double ~y)))
   :inline-arities #{3}}
  ([[^double x ^double y] ^double v] (<= x v y))
  ([^double x ^double y ^double v] (<= x v y)))

(defn between-?
  "Check if given number is within the range (x,y]."
  {:inline (fn [x y v] `(let [v# (double ~v)]
                         (and (< (double ~x) v#) (<= v# (double ~y)))))
   :inline-arities #{3}}
  ([[^double x ^double y] ^double v] (and (< x v) (<= v y)))
  ([^double x ^double y ^double v] (and (< x v) (<= v y))))

;;

(defn absolute-error
  "Absolute error between two values"
  (^double [^double v ^double v-approx]
   (abs (- v v-approx))))

(defn relative-error
  "Relative error between two values"
  (^double [^double v ^double v-approx]
   (abs (/ (- v v-approx) v))))

;; intervals

(defn slice-range 
  "Generates a sequence of `cnt` evenly spaced double values covering a range.

  The range can be given explicitly, derived from a collection of data, or defaults to `[0.0, 1.0]` when only `cnt` is provided.

  Parameters:

  - `[cnt]` (long): Number of points; range defaults to `[0.0, 1.0]`.
  - `[data cnt]` (sequence, long): Range is `[(min data), (max data)]`. Non-finite values (`NaN`, infinities) are removed from `data` before computing the range.
  - `[start end cnt]` (double, double, long): Explicit inclusive range `[start, end]`.

  The resulting sequence is inclusive, i.e. it starts at `start` (or the derived minimum) and ends at `end` (or the derived maximum). When `cnt` is `1`, a single value equal to the midpoint of the range is returned instead.

  Returns a sequence of `cnt` doubles, or an empty sequence when `cnt` is `0`.

  See also [[make-norm]]."
  ([data ^long cnt]
   (let [d (sort (remove invalid-double? data))]
     (slice-range (first d) (last d) cnt)))
  ([^double start ^double end ^long cnt] (if (= cnt 1)
                                           (list (+ start (* 0.5 (- end start))))
                                           (map (make-norm 0.0 (dec cnt) start end) (range cnt))))
  ([^long cnt] (slice-range 0.0 1.0 cnt)))

(defn cut
  "Divides a numerical range into `breaks` equally spaced, half-open intervals.

  The function first generates `breaks + 1` equally spaced boundary points across the range using [[slice-range]], then pairs consecutive points into intervals. Each interval is closed on the right, i.e. `(lower, upper]`, so that every value in the range falls into exactly one interval. The lower bound of the very first interval is nudged one ulp downwards (via [[prev-double]]) so that the exact minimum of the range is included too.

  Parameters:

  - `[data breaks]` (sequence, long): The range is `[(min data), (max data)]`. Non-finite values (`NaN`, infinities) are removed from `data` before computing the range.
  - `[x1 x2 breaks]` (double, double, long): Explicit range `[x1, x2]`.
  - `breaks` (long): Desired number of intervals, must be positive.

  Returns a sequence of `breaks` 2-element vectors `[lower-bound upper-bound]`.

  See also [[slice-range]], [[co-intervals]]."
  ([data ^long breaks]
   (let [d (sort (remove invalid-double? data))]
     (cut (first d) (last d) breaks)))
  ([^double x1 ^double x2 ^long breaks]
   (let [[[^double start end] & r] (->> (slice-range x1 x2 (long-inc breaks))
                                        (partition 2 1))]
     (conj r (list (prev-double start) end)))))

(defn co-intervals
  "Divides `data` into `number` overlapping intervals, each containing a similar count of values.

  Unlike [[cut]], which produces equally spaced intervals, this function determines interval boundaries from the sorted data itself so that every interval covers roughly the same number of observations, with consecutive intervals overlapping by a given proportion. This replicates the behaviour of R's `co.intervals()` function. Non-finite values (`NaN`, infinities) are removed from `data` before processing.

  Parameters:

  - `data` (sequence of numbers): Values to partition.
  - `number` (long): Desired number of intervals. Defaults to `6`.
  - `overlap` (double): Desired overlap proportion between consecutive intervals, in range `[0.0, 1.0]`. Defaults to `0.5`.

  Returns a sequence of intervals, each represented as a 2-element vector `[lower-bound upper-bound]`. The resulting number of intervals may be smaller than `number` when the data does not support that many distinct overlapping groups.

  See also [[cut]], [[group-by-intervals]]."
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
         eps (* 0.5 (double (if (seq diffs) (reduce min diffs) 0.0)))]
     (for [[[^double px ^double cx] [^double py ^double cy]] (map vector
                                                                  (partition 2 1 (conj x1 (dec ^double (first x1))))
                                                                  (partition 2 1 (conj xr (dec ^double (first xr)))))
           :when (or (pos? (- cx px))
                     (pos? (- cy py)))]
       [(- cx eps) (+ cy eps)]))))

(defn group-by-intervals
  "Groups values from a sequence `coll` into specified intervals.

  The function partitions the values in `coll` based on which interval they fall into.
  Each interval is a 2-element vector `[lower upper]`. Values are included in an interval if they are strictly greater than the `lower` bound and less than or equal to the `upper` bound, using [[between-?]].

  Arguments:
  - `intervals`: A sequence of 2-element vectors representing the intervals `[lower upper]`.
  - `coll`: A collection of numerical values to be grouped.

  If `intervals` are not provided, the function first calculates overlapping intervals using [[co-intervals]] from the values in `coll`, and then groups the values into these generated intervals.

  Returns a map where keys are the interval vectors and values are sequences of the numbers from `coll` that fall within that interval."
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
                              (recur (>> (long-sub a b) 1) b)
                              (recur (>> (long-sub b a) 1) a))))

(defn gcd
  "Fast binary greatest common divisor (Stein's algorithm)"
  ^long [^long a ^long b]
  (gcd- (long-abs a) (long-abs b)))

(defn lcm
  "Fast binary least common multiplier."
  ^long [^long a ^long b]
  (if (> a b)
    (* b (/ a (gcd- (long-abs a) (long-abs b))))
    (* a (/ b (gcd- (long-abs a) (long-abs b))))))

;; arithmetic-geometric-mean

(defn agm
  "Calculates the arithmetic-geometric mean (agM) of two numbers `x` and `y`.

  The agM is computed iteratively: at each step the arithmetic mean and the geometric mean of the current pair replace `x` and `y`, and both sequences converge to the same limit, which is returned. This limit lies between `x` and `y`.

  Parameters:

  - `x`, `y` (doubles): The two numbers to average.
  - `abs-tol` (double): Absolute tolerance used to detect convergence, i.e. `x` and `y` are considered equal. Defaults to `1.0e-12`.
  - `max-iters` (long): Maximum number of iterations allowed before giving up. Defaults to `100`.

  Both `x` and `y` should be non-negative, since the geometric mean step involves a square root.

  Returns the arithmetic-geometric mean as a double.

  Throws an exception if convergence is not reached within `max-iters` iterations."
  (^double [^double x ^double y] (agm x y 1.0e-12))
  (^double [^double x ^double y ^double abs-tol] (agm x y abs-tol 100))
  (^double [^double x ^double y ^double abs-tol ^long max-iters ]
   (if (zero? max-iters)
     (throw (ex-info "agM Convergence failed." {:x x :y y :diff (abs (- x y))}))
     (if (delta-eq x y abs-tol)
       (* 0.5 (+ x y))
       (recur (* 0.5 (+ x y))
              (sqrt (* x y))
              abs-tol
              (dec max-iters))))))

;;

(defn sample
  "Samples a function `f` by evaluating it at evenly spaced points within a numerical range.

  Generates `number-of-values` points in the specified range `[domain-min, domain-max]` (inclusive) and applies `f` to each point.

  Arguments:

  - `f`: The function to sample. Should accept a single `double` argument.
  - `number-of-values`: The total number of points to generate (a positive `long`).
  - `domain-min`: The lower bound of the sampling range (a `double`). Defaults to `0.0`.
  - `domain-max`: The upper bound of the sampling range (a `double`). Defaults to `1.0`.
  - `domain?`: A boolean flag. If `true`, returns pairs `[x, (f x)]`. If `false` (default), returns just `(f x)`.

  Arities:

  - `[f number-of-values]`: Samples `f` in `[0.0, 1.0]`. Returns `(f x)` values.
  - `[f number-of-values domain?]`: Samples `f` in `[0.0, 1.0]`. Returns `[x, (f x)]` pairs if `domain?` is true, otherwise `(f x)` values.
  - `[f domain-min domain-max number-of-values]`: Samples `f` in `[domain-min, domain-max]`. Returns `(f x)` values.
  - `[f domain-min domain-max number-of-values domain?]`: Samples `f` in `[domain-min, domain-max]`. Returns `[x, (f x)]` pairs if `domain?` is true, otherwise `(f x)` values.

  The points are generated linearly from `domain-min` to `domain-max`. If `number-of-values` is 1, it samples only the midpoint of the range.

  Returns a sequence of `double` values or vectors `[double, double]` depending on `domain?`."
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
  "Assigns ranks to values in a collection, handling ties according to a specified strategy.
  Ranks are 0-based indices indicating the position of each element in the sorted collection.

  Arguments:
  
  - `vs`: A collection of comparable values.
  - `ties`: The tie-breaking strategy (keyword, optional, default `:average`).
    Supported strategies:
    - `:average`: Assign the average rank to all tied values.
    - `:first`: Assign ranks based on their appearance order in the input.
    - `:last`: Assign ranks based on their appearance order in the input (reverse of `:first`).
    - `:random`: Assign random ranks to tied values.
    - `:min`: Assign the minimum rank to all tied values.
    - `:max`: Assign the maximum rank to all tied values.
    - `:dense`: Assign consecutive ranks without gaps (like `data.table::frank` in R).
  - `desc?`: If true, rank in descending order (boolean, optional, default `false`).

  Returns a sequence of rank values corresponding to the input elements."
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
                      (fn ^double [v] (/ ^double (reduce + (map first v)) (count v))))
             m (map (fn [[k v]] [k (tie-fn v)]) indexed-sorted-map)
             m (if (= ties :dense)
                 (map-indexed (fn [id [k _]]
                                [k id]) (sort-by second m))
                 m)]
         (map (into {} m) vs))))))

(def rank1 ^{:doc "[[rank]] with indexing statring from 1"}
  (comp (partial map inc) rank))

(defn order
  "Computes the permutation of indices that would sort the input collection `vs`.

  The result is a sequence of 0-based indices such that applying them to the original collection
  using `(map #(nth vs %) result)` yields a sorted sequence.

  Arguments:

  - `vs`: A collection of comparable values.
  - `decreasing?`: Optional boolean (default false). If true, the indices permute for a descending sort.

  Returns: A sequence of 0-based indices."
  ([vs] (order vs false))
  ([vs decreasing?]
   (->> (map-indexed vector vs)
        (sort-by second (if decreasing?
                          clojure.core/>
                          clojure.core/<))
        (map first))))

;;

(def double-array-type (Class/forName "[D"))
(def double-double-array-type (Class/forName "[[D"))

(def ^{:doc "Convert double array into sequence.

  Alias for `seq`."} double-array->seq seq)

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
  '[* + - / > < >= <= == abs rem quot mod bit-and-not bit-set bit-clear bit-test bit-flip bit-or bit-and bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd? bool-not << >> >>> not==])

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
  ([] (unuse-primitive-operators #{}))
  ([skip-set]
   (when (using-primitive-operators?)
     (doseq [v (remove skip-set vars-to-exclude)]
       (ns-unmap *ns* v))
     (refer-clojure :exclude (seq skip-set)))))

;;;;

(defn fast+
  {:inline (primitivemath-nary-inline 'add nil 0.0)
   :inline-arities >=0?
   :doc "Primitive and inlined `+` as a function"
   :deprecated "Use `+` instead"}
  (^double [] 0.0)
  (^double [^double a] a)
  (^double [^double a ^double b] (. PrimitiveMath (add a b)))
  ([a b & r] (reduce fast+ (. PrimitiveMath (add (double a) (double b))) r)))

(defn fast-
  {:inline (primitivemath-nary-inline 'subtract 'negate 0.0)
   :inline-arities >=0?
   :doc "Primitive and inlined `-` as a function"
   :deprecated "Use `-` instead"}
  (^double [] 0.0)
  (^double [^double a] (. PrimitiveMath (negate a)))
  (^double [^double a ^double b] (. PrimitiveMath (subtract a b)))
  ([a b & r] (reduce fast- (. PrimitiveMath (subtract (double a) (double b))) r)))

(defn fast*
  {:inline (primitivemath-nary-inline 'multiply nil 1.0)
   :inline-arities >=0?
   :doc "Primitive and inlined `*` as a function"
   :deprecated "Use `*` instead"}
  (^double [] 1.0)
  (^double [^double a] a)
  (^double [^double a ^double b] (. PrimitiveMath (multiply a b)))
  ([a b & r] (reduce fast* (. PrimitiveMath (multiply (double a) (double b))) r)))

(defn fast-div
  {:inline (primitivemath-nary-inline 'divide 'reciprocal 0.0)
   :inline-arities >=0?
   :doc "Primitive and inlined `/` as a function"
   :deprecated "Use `/` instead"}
  (^double [] 1.0)
  (^double [^double a] (. PrimitiveMath (reciprocal a)))
  (^double [^double a ^double b] (. PrimitiveMath (divide a b)))
  ([a b & r] (reduce fast-div (. PrimitiveMath (divide (double a) (double b))) r)))

(defn fast-max
  {:inline (primitivemath-nary-inline 'max)
   :inline-arities >=1?
   :doc "Primitive and inlined `max` as a function"
   :deprecated "Use `max` instead"}
  (^double [^double a] a)
  (^double [^double a ^double b] (. PrimitiveMath (max a b)))
  ([a b & r] (reduce fast-max (. PrimitiveMath (max (double a) (double b))) r)))

(defn fast-min
  {:inline (primitivemath-nary-inline 'min)
   :inline-arities >=1?
   :doc "Primitive and inlined `min` as a function"
   :deprecated "Use `min` instead"}
  (^double [^double a] a)
  (^double [^double a ^double b] (. PrimitiveMath (min a b)))
  ([a b & r] (reduce fast-min (. PrimitiveMath (min (double a) (double b))) r)))

(defn fast-identity
  {:inline (fn [x] `~x) :inline-arities #{1}
   :doc "Identity on double."
   :deprecated "Use `identity-double` instead"}
  ^double [^double a] a)


;;;;
