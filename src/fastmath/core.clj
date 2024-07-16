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
   :exclude [* + - / > < >= <= == rem quot mod bit-or bit-and bit-and-not bit-set bit-clear bit-test bit-flip bit-xor bit-not bit-shift-left bit-shift-right unsigned-bit-shift-right inc dec zero? neg? pos? min max even? odd? abs integer?])
  (:import [net.jafama FastMath]
           [fastmath.java PrimitiveMath]
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
   :doc "Primitive and inlined `+`."}
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
   :doc "Primitive and inlined `+`. Coerces arguments and returned values to longs."}
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
   :doc "Primitive and inlined `-`."}
  (^double [^double a] (. PrimitiveMath (negate a)))
  (^double [^double a ^double b] (. PrimitiveMath (subtract a b)))
  (^double [^double a ^double b ^double c] (. PrimitiveMath (subtract (. PrimitiveMath (subtract a b)) c)))
  (^double [^double a ^double b ^double c ^double d]
   (. PrimitiveMath (subtract (. PrimitiveMath (subtract (. PrimitiveMath (subtract a b)) c)) d)))
  ([a b c d & r] (reduce - (- (double a) (double b) (double c) (double d)) r)))

(defn long-sub
  {:inline (primitivemath-nary-inline-long 'subtract 'negate)
   :inline-arities >=1?
   :doc "Primitive and inlined `-`. Coerces arguments and returned values to longs."}
  (^long [^long a] (. PrimitiveMath (negate a)))
  (^long [^long a ^long b] (. PrimitiveMath (subtract a b)))
  (^long [^long a ^long b ^long c] (. PrimitiveMath (subtract (. PrimitiveMath (subtract a b)) c)))
  (^long [^long a ^long b ^long c ^long d]
   (. PrimitiveMath (subtract (. PrimitiveMath (subtract (. PrimitiveMath (subtract a b)) c)) d)))
  ([a b c d & r] (reduce long-sub (long-sub a b c d) r)))

(defn *
  {:inline (primitivemath-nary-inline 'multiply nil 1.0)
   :inline-arities >=0?
   :doc "Primitive and inlined `*`."}
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
   :doc "Primitive and inlined `*`. Coerces arguments and returned values to longs."}
  (^long [] 1)
  (^long [^long a] a)
  (^long [^long a ^long b] (. PrimitiveMath (multiply a b)))
  (^long [^long a ^long b ^long c] (. PrimitiveMath (add (. PrimitiveMath (multiply a b)) c)))
  (^long [^long a ^long b ^long c ^long d]
   (. PrimitiveMath (multiply (. PrimitiveMath (multiply (. PrimitiveMath (multiply a b)) c)) d)))
  ([a b c d & r] (reduce long-mult (long-mult a b c d) r)))

(defn /
  {:inline (primitivemath-nary-inline 'divide 'reciprocal)
   :inline-arities >=1?
   :doc "Primitive and inlined `/`."}
  (^double [^double a] (. PrimitiveMath (reciprocal a)))
  (^double [^double a ^double b] (. PrimitiveMath (divide a b)))
  (^double [^double a ^double b ^double c] (. PrimitiveMath (divide (. PrimitiveMath (divide a b)) c)))
  (^double [^double a ^double b ^double c ^double d]
   (. PrimitiveMath (divide (. PrimitiveMath (divide (. PrimitiveMath (divide a b)) c)) d)))
  ([a b c d & r] (reduce / (/ (double a) (double b) (double c) (double d)) r)))

(defn long-div
  {:inline (primitivemath-nary-inline-long 'subtract 'reciprocal)
   :inline-arities >=1?
   :doc "Primitive and inlined `/`. Coerces to arguments and returned values to longs."}
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
   :doc "Primitive and inlined `min`."}
  (^double [^double a] a)
  (^double [^double a ^double b] (. PrimitiveMath (min a b)))
  (^double [^double a ^double b ^double c] (. PrimitiveMath (min (. PrimitiveMath (min a b)) c)))
  (^double [^double a ^double b ^double c ^double d]
   (. PrimitiveMath (min (. PrimitiveMath (min (. PrimitiveMath (min a b)) c)) d)))
  ([a b c d & r] (reduce min (min (double a) (double b) (double c) (double d)) r)))

(defn long-min
  {:inline (primitivemath-nary-inline-long 'min)
   :inline-arities >=1?
   :doc "Primitive and inlined `min`. Coerces arguments and returned values to longs."}
  (^long [^long a] a)
  (^long [^long a ^long b] (. PrimitiveMath (min a b)))
  (^long [^long a ^long b ^long c] (. PrimitiveMath (min (. PrimitiveMath (min a b)) c)))
  (^long [^long a ^long b ^long c ^long d]
   (. PrimitiveMath (min (. PrimitiveMath (min (. PrimitiveMath (min a b)) c)) d)))
  ([a b c d & r] (reduce long-min (long-min a b c d) r)))

(defn max
  {:inline (primitivemath-nary-inline 'max)
   :inline-arities >=1?
   :doc "Primitive and inlined `max`."}
  (^double [^double a] a)
  (^double [^double a ^double b] (. PrimitiveMath (max a b)))
  (^double [^double a ^double b ^double c] (. PrimitiveMath (max (. PrimitiveMath (max a b)) c)))
  (^double [^double a ^double b ^double c ^double d]
   (. PrimitiveMath (max (. PrimitiveMath (max (. PrimitiveMath (max a b)) c)) d)))
  ([a b c d & r] (reduce max (max (double a) (double b) (double c) (double d)) r)))

(defn long-max
  {:inline (primitivemath-nary-inline-long 'max)
   :inline-arities >=1?
   :doc "Primitive and inlined `max`. Coerces arguments and returned values to longs."}
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
  "Primitive math equality function."
  {:inline (primitivemath-nary-inline-predicate 'eq)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (eq a b)))
  ([a b & r]
   (boolean (and (== (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (eq x y)) (reduced false) y)) r)))))
(defn eq 
  "Primitive math equality function."
  {:inline (primitivemath-nary-inline-predicate 'eq)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (eq  a b)))
  ([a b & r]
   (boolean (and (eq (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (eq x y)) (reduced false) y)) r)))))

(defn < 
  "Primitive math less-then function."
  {:inline (primitivemath-nary-inline-predicate 'lt)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (lt a b)))
  ([a b & r]
   (boolean (and (< (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (lt x y)) (reduced false) y)) r)))))

(defn > 
  "Primitive math greater-than function."
  {:inline (primitivemath-nary-inline-predicate 'gt)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (gt a b)))
  ([a b & r]
   (boolean (and (> (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (gt x y)) (reduced false) y)) r)))))

(defn <= 
  "Primitive math less-and-equal function."
  {:inline (primitivemath-nary-inline-predicate 'lte)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (lte a b)))
  ([a b & r]
   (boolean (and (<= (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (lte x y)) (reduced false) y)) r)))))

(defn >= 
  "Primitive math greater-and-equal function."
  {:inline (primitivemath-nary-inline-predicate 'gte)
   :inline-arities >=1?}
  ([_] true)
  ([^double a ^double b] (. PrimitiveMath (gte a b)))
  ([a b & r]
   (boolean (and (>= (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (gte x y)) (reduced false) y)) r)))))

(defn not==
  "Not equality. For more than two arguments, pairwise not equality is checked.

  `(not== 1 2 1)` === `(and (not= 1 2) (not= 2 1))`"
  {:inline (fn ([_] false)
             ([a b] `(. PrimitiveMath (neq ~a ~b)))
             ([a b & r] `(and (. PrimitiveMath (neq ~a ~b))
                              ~@(map (fn [[x y]] `(. PrimitiveMath (neq ~x ~y)))
                                     (partition 2 1 r)))))
   :inline-arities >=1?}
  ([_] false)
  ([^double a ^double b] (. PrimitiveMath (neq a b)))
  ([a b & r]
   (boolean (and (not== (double a) (double b))
                 (reduce (fn [^double x ^double y]
                           (if-not (. PrimitiveMath (neq x y)) (reduced false) y)) r)))))

;;;;;;;;;;;;;;

(defn bit-and
  "x ∧ y - bitwise AND"
  {:inline (primitivemath-nary-inline-long 'bitAnd)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitAnd x y)))
  ([x y & r] (reduce bit-and (. PrimitiveMath (bitAnd x y)) r)))

(defn bit-nand
  "~(x ∧ y) - bitwise NAND"
  {:inline (primitivemath-nary-inline-long 'bitNand)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitNand x y)))
  ([x y & r] (reduce bit-nand (. PrimitiveMath (bitNand x y)) r)))

(defn bit-and-not
  "x ∧ ~y - bitwise AND (with complement second argument)"
  {:inline (primitivemath-nary-inline-long 'bitAndNot)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitAndNot x y)))
  ([x y & r] (reduce bit-and-not (. PrimitiveMath (bitAndNot x y)) r)))

(defn bit-or
  "x ∨ y - bitwise OR"
  {:inline (primitivemath-nary-inline-long 'bitOr)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitOr x y)))
  ([x y & r] (reduce bit-or (. PrimitiveMath (bitOr x y)) r)))

(defn bit-nor
  "~(x ∨ y) - bitwise NOR"
  {:inline (primitivemath-nary-inline-long 'bitNor)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitNor x y)))
  ([x y & r] (reduce bit-nor (. PrimitiveMath (bitNor x y)) r)))

(defn bit-xor
  "x⊕y - bitwise XOR"
  {:inline (primitivemath-nary-inline-long 'bitXor)
   :inline-arities >=1?}
  (^long [^long x] x)
  (^long [^long x ^long y] (. PrimitiveMath (bitXor x y)))
  ([x y & r] (reduce bit-xor (. PrimitiveMath (bitXor x y)) r)))

(defn bit-xnor
  "~(x⊕y) - bitwise XNOR"
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
  "Primitive boolean xor"
  {:inline (primitivemath-nary-inline 'xor)
   :inline-arities >=2?}
  ([x y] (. PrimitiveMath (xor (boolean x) (boolean y))))
  ([x y & r] (reduce bool-xor (. PrimitiveMath (xor (boolean x) (boolean y))) r)))

(defn xor
  "Primitive boolean xor"
  {:inline (primitivemath-nary-inline 'xor)
   :inline-arities >=2?}
  ([x y] (. PrimitiveMath (xor (boolean x) (boolean y))))
  ([x y & r] (reduce xor (. PrimitiveMath (xor (boolean x) (boolean y))) r)))

;;;;

(defn negative-zero?
  "Check if zero is negative, ie. -0.0"
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
   :doc "Identity on double."}
  ^long [^long a] a)

(defn integer?
  "Check if given real number is an integer."
  {:inline (fn [v] `(== ~v (FastMath/rint ~v))) :inline-arities #{1}}
  [^double v]
  (== v (FastMath/rint v)))

;; macros for polynomials

(defn- ->fma
  [x y z]
  (if (< jvm-version 9)
    `(+ ~z (* ~x ~y))
    `(Math/fma (double ~x) (double ~y) (double ~z))))

(defmacro ^:private fma-macro
  [x y z]
  (if (< jvm-version 9)
    `(+ ~z (* ~x ~y))
    `(Math/fma ~x ~y ~z)))

(defn muladd
  "`(x y z)` -> `(+ z (* x y))` or `Math/fma` for java 9+"
  {:inline (fn [x y z] (->fma x y z))
   :inline-arities #{3}}
  ^double [^double x ^double y ^double z]
  (fma-macro x y z))

(defn fma
  "`(x y z)` -> `(+ z (* x y))` or `Math/fma` for java 9+"
  {:inline (fn [x y z] (->fma x y z))
   :inline-arities #{3}}
  ^double [^double x ^double y ^double z]
  (fma-macro x y z))

(defn negmuladd
  "`(x y z)` -> `(+ z (* -x y))` or `Math/fma` for java 9+"
  {:inline (fn [x y z] (->fma x y z))
   :inline-arities #{3}}
  ^double [^double x ^double y ^double z]
  (fma-macro (- x) y z))

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

(defonce ^{:const true :tag 'double :doc "ULP(1)/2. Half of the smallest difference between 1.0 and next possible double floating point number."}
  MACHINE-EPSILON (* 0.5 (double (loop [d (double 1.0)]
                                   (if (not== 1.0 (+ 1.0 (* d 0.5)))
                                     (recur (* d 0.5))
                                     d)))))

(def ^{:const true :tag 'double :doc "5*ULP(1)"}
  MACHINE-EPSILON10 (* 10.0 MACHINE-EPSILON))

(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{1}{3}\\\\)"} THIRD      0.333333333333333333333333)
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{1}{3}\\\\)"} ONE_THIRD  0.333333333333333333333333)
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{2}{3}\\\\)"} TWO_THIRD  0.666666666666666666666666)
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{2}{3}\\\\)"} TWO_THIRDS 0.666666666666666666666666)
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{1}{6}\\\\)"} SIXTH      0.166666666666666666666666)
(def ^{:const true :tag 'double :doc "Value of \\\\(\\frac{1}{6}\\\\)"} ONE_SIXTH  0.166666666666666666666666)

(defn sin
  "sin(x)"
  {:inline (fn [x] `(. FastMath (sin (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (sin x)))

(defn sinpi
  "sin(pi*x)"
  {:inline (fn [x] `(. FastMath (sin (* PI (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (sin (* PI x))))

(defn cos
  "cos(x)"
  {:inline (fn [x] `(. FastMath (cos (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (cos x)))

(defn cospi
  "cos(pi*x)"
  {:inline (fn [x] `(. FastMath (cos (* PI (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (cos (* PI x))))

(defn tan
  "tan(x)"
  {:inline (fn [x] `(. FastMath (tan (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (tan x)))

(defn tanpi
  "tan(pi*x)"
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
  "cot(x)"
  {:inline (fn [x] `(/ (tan (double ~x))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/tan x))))

(defn cotpi
  "cot(pi*x)"
  {:inline (fn [x] `(/ (tan (* PI (double ~x)))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/tan (* PI x)))))

(defn sec
  "sec(x)"
  {:inline (fn [x] `(/ (cos (double ~x))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/cos x))))

(defn secpi
  "sec(pi*x)"
  {:inline (fn [x] `(/ (cos (* PI (double ~x)))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/cos (* PI x)))))

(defn csc
  "csc(x)"
  {:inline (fn [x] `(/ (sin (double ~x))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/sin x))))

(defn cscpi
  "csc(pi*x)"
  {:inline (fn [x] `(/ (sin (* PI (double ~x)))))
   :inline-arities #{1}}
  (^double [^double x] (/ (FastMath/sin (* PI x)))))

;; Additional cyclometric functions

(defn acot
  "acot(x)"
  {:inline (fn [x] `(- HALF_PI (atan (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (- HALF_PI (FastMath/atan x)))

(defn asec
  "asec(x)"
  {:inline (fn [x] `(acos (/ ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acos (/ 1.0 x)))

(defn acsc
  "acsc(x)"
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
  "Hyperbolic cotangent"
  {:inline (fn [x] `(/ (tanh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (/ (FastMath/tanh x)))

(defn sech
  "Hyperbolic secant"
  {:inline (fn [x] `(/ (cosh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (/ (FastMath/cosh x)))

(defn csch
  "Hyperbilic cosecant"
  {:inline (fn [x] `(/ (sinh (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (/ (FastMath/sinh x)))

;; Additional inverse hyperbolic functions
(defn acoth
  "Area hyperbolic cotangent"
  {:inline (fn [x] `(atanh (/ ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/atanh (/ x)))

(defn asech
  "Area hyperbolic secant"
  {:inline (fn [x] `(acosh (/ ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acosh (/ x)))

(defn acsch
  "Area hyperbolic cosecant"
  {:inline (fn [x] `(asinh (/ ~x)))
   :inline-arities #{1}}
  ^double [^double v] (FastMath/asinh (/ v)))

;; historical

(defn crd
  "Chord"
  {:inline (fn [x] `(* 2.0 (sin (* 0.5 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (* 2.0 (FastMath/sin (* 0.5 x))))

(defn acrd
  "Inverse chord"
  {:inline (fn [x] `(* 2.0 (asin (* 0.5 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (* 2.0 (FastMath/asin (* 0.5 x))))

(defn versin
  "Versine"
  {:inline (fn [x] `(- 1.0 (cos (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (- 1.0 (FastMath/cos x)))

(defn coversin
  "Coversine"
  {:inline (fn [x] `(- 1.0 (sin (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (- 1.0 (FastMath/sin x)))

(defn vercos
  "Vercosine"
  {:inline (fn [x] `(inc (cos (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (inc (FastMath/cos x)))

(defn covercos
  "Covercosine"
  {:inline (fn [x] `(inc (sin (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (inc (FastMath/sin x)))

(defn aversin
  "Arc versine"
  {:inline (fn [x] `(acos (- 1.0 ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acos (- 1.0 x)))

(defn acoversin
  "Arc coversine"
  {:inline (fn [x] `(asin (- 1.0 ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/asin (- 1.0 x)))

(defn avercos
  "Arc vecosine"
  {:inline (fn [x] `(acos (dec ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acos (dec x)))

(defn acovercos
  "Arc covercosine"
  {:inline (fn [x] `(asin (dec ~x)))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/asin (dec x)))

(defn haversin
  "Haversine formula for value or lattitude and longitude pairs."
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
  "Hacoversine"
  {:inline (fn [x] `(* 0.5 (- 1.0 (sin (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x] (* 0.5 (- 1.0 (FastMath/sin x))))

(defn havercos
  "Havercosine"
  {:inline (fn [x] `(* 0.5 (inc (cos (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x] (* 0.5 (inc (FastMath/cos x))))

(defn hacovercos
  "Hacovercosine"
  {:inline (fn [x] `(* 0.5 (inc (sin (double ~x)))))
   :inline-arities #{1}}
  ^double [^double x] (* 0.5 (inc (FastMath/sin x))))

(defn ahaversin
  "Arc haversine"
  {:inline (fn [x] `(acos (- 1.0 (* 2.0 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acos (- 1.0 (* 2.0 x))))

(defn ahacoversin
  "Arc hacoversine"
  {:inline (fn [x] `(asin (- 1.0 (* 2.0 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/asin (- 1.0 (* 2.0 x))))

(defn ahavercos
  "Arc havecosine"
  {:inline (fn [x] `(acos (dec (* 2.0 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/acos (dec (* 2.0 x))))

(defn ahacovercos
  "Arc hacovercosine"
  {:inline (fn [x] `(asin (dec (* 2.0 ~x))))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/asin (dec (* 2.0 x))))

(defn exsec
  "Exsecant"
  {:inline (fn [x] `(dec (sec ~x)))
   :inline-arities #{1}}
  ^double [^double x] (dec (sec x)))

(defn excsc
  "Excosecant"
  {:inline (fn [x] `(dec (csc ~x)))
   :inline-arities #{1}}
  ^double [^double x] (dec (csc x)))

(defn aexsec
  "Arc exsecant"
  {:inline (fn [x] `(asec (inc ~x)))
   :inline-arities #{1}}
  ^double [^double x] (asec (inc x)))

(defn aexcsc
  "Arc excosecant"
  {:inline (fn [x] `(acsc (inc ~x)))
   :inline-arities #{1}}
  ^double [^double x] (acsc (inc x)))

(defn haversine-dist
  "Haversine distance `d` for `r=1`"
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

(defn log
  "log(x)=ln(x) or logarithm with given `base`."
  {:inline (fn ([x] `(. FastMath (log (double ~x))))
             ([base x] `(/ (. FastMath (log (double ~x))) (. FastMath (log (double ~base))))))
   :inline-arities #{1 2}}
  (^double [^double x] (. FastMath (log x)))
  (^double [^double base ^double x] (/ (. FastMath (log x)) (. FastMath (log base)))))

(defn ln
  "log(x)=ln(x)"
  {:inline (fn [x] `(. FastMath (log (double ~x))))}
  ^double [^double x] (. FastMath (log x)))

(defn log10
  "log_10(x)"
  {:inline (fn [x] `(. FastMath (log10 (double ~x))))}
  ^double [^double x] (. FastMath (log10 x)))

(defn log1p
  "log(1+x) for small x"
  {:inline (fn [x] `(. FastMath (log1p (double ~x))))}
  ^double [^double x] (. FastMath (log1p x)))

(defn expm1
  "exp(x)-1 for small x"
  {:inline (fn [x] `(. FastMath (expm1 (double ~x))))}
  ^double [^double x] (. FastMath (expm1 x)))

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
  "log(1-exp(x)), x<0"
  ^double [^double x]
  (if (< x LOG_HALF)
    (FastMath/log1p (- (FastMath/exp x)))
    (FastMath/log (- (FastMath/expm1 x)))))

(defn log2mexp
  "log(2-exp(x))"
  {:inline (fn [x] `(FastMath/log1p (- (FastMath/expm1 (double ~x)))))
   :inline-arities #{1}}
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
  {:inline (fn [x] `(FastMath/log (- (FastMath/log1p (- (double ~x))))))
   :inline-arities #{1}}
  ^double [^double x]
  (FastMath/log (- (FastMath/log1p (- x)))))

(defn loglog
  "-log(-log(x))"
  {:inline (fn [x] `(- (FastMath/log (- (FastMath/log ~x)))))
   :inline-arities #{1}}
  ^double [^double x]
  (- (FastMath/log (- (FastMath/log x)))))

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
  {:inline (fn [x] `(- (FastMath/expm1 (- (FastMath/exp (double ~x))))))
   :inline-arities #{1}}
  ^double [^double x]
  (- (FastMath/expm1 (- (FastMath/exp x)))))

(defn expexp
  "exp(-exp(-x))"
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
  "\\\\(\\sqrt{x}\\\\)"
  {:inline (fn [x] `(. FastMath (sqrt (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (. FastMath (sqrt x)))

(defn cbrt
  "\\\\(\\sqrt[3]{x}\\\\)"
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
(def ^{:const true :tag 'double :doc "\\\\(\\frac{180}{\\pi}\\\\)"} rad-in-deg (/ 180.0 PI))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{\\pi}{180}\\\\)"} deg-in-rad (/ PI 180.0))

(defn radians
  "Convert degrees into radians."
  {:inline (fn [deg] `(* deg-in-rad ~deg))
   :inline-arities #{1}}
  ^double [^double deg] (* deg-in-rad deg))

(defn degrees
  "Convert radians into degrees."
  {:inline (fn [rad] `(* rad-in-deg ~rad))
   :inline-arities #{1}}
  ^double [^double rad] (* rad-in-deg rad))

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
  {:inline (fn [x] `(/ (inc (FastMath/exp (- ~x)))))
   :inline-arities #{1}}
  ^double [^double x]
  (/ (inc (FastMath/exp (- x)))))

(def ^{:doc "Alias for [[sigmoid]]"} logistic sigmoid)

(defn logit
  "Logit function"
  {:inline (fn [x] `(FastMath/log (/ ~x (- 1.0 ~x))))
   :inline-arities #{1}}
  ^double [^double x]
  (FastMath/log (/ x (- 1.0 x))))

(defn log2
  "Logarithm with base 2.

  \\\\(\\ln_2{x}\\\\)"
  {:inline (fn [x] `(* (FastMath/log (double ~x)) INV_LN2))
   :inline-arities #{1}}
  ^double [^double x]
  (* (FastMath/log x) INV_LN2))

;; \\(\log_b x\\)
(defn logb
  "Logarithm with base `b`.

  \\\\(\\ln_b{x}\\\\)"
  {:inline (fn [b x] `(/ (FastMath/log (double ~x)) (FastMath/log (double ~b))))
   :inline-arities #{2}}
  ^double [^double b ^double x]
  (/ (FastMath/log x) (FastMath/log b)))

(defn logcosh
  "log(cosh(x))"
  {:inline (fn [x] `(let [absx# (FastMath/abs (double ~x))]
                     (- (+ absx# (log1pexp (* -2.0 absx#))) LN2)))
   :inline-arities #{1}}
  ^double [^double x]
  (let [absx (FastMath/abs x)]
    (- (+ absx (log1pexp (* -2.0 absx))) LN2)))

;; \\(\log_2 e\\)
(def ^{:const true :tag 'double :doc "\\\\(\\log_{2}{e}\\\\)"} LOG2E (log2 E))

;; \\(\log_{10} e\\)
(def ^{:const true :tag 'double :doc "\\\\(\\log_{10}{e}\\\\)"} LOG10E (log10 E))

;; Powers (normal, quick)

;; use Math here due to some fastmath innacuracies

(defn pow
  "Power of a number"
  {:inline (fn [x exponent] `(. Math (pow (double ~x) (double ~exponent))))
   :inline-arities #{2}}
  ^double [^double x ^double exponent] (. Math (pow x exponent)))

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
    (exp (Gamma/logGamma (double (inc n))))))

(defn inv-factorial "Inverse of factorial, 1/n!"
  ^double [^long n]
  (if (< n 21)
    (/ 1.0 (long (factorial20-table n)))
    (exp (- (Gamma/logGamma (double (inc n)))))))

(defn log-factorial
  "Log factorial, alias to log-gamma"
  {:inline (fn [x] `(Gamma/logGamma (double (inc (long ~x)))))
   :inline-arities #{1}}
  ^double [^long x] (Gamma/logGamma (double (inc x))))

(defn falling-factorial-int
  "Falling (descending) factorial for integer n."
  ^double [^long n ^double x]
  (if (not-neg? n)
    (loop [i (long 0)
           v 1.0]
      (if (== i n) v
          (recur (inc i) (* v (- x i)))))
    (/ (falling-factorial-int (- n) (- x n)))))

(defn falling-factorial
  "Falling (descending) factorial."
  ^double [^double n ^double x]
  (if (integer? n)
    (falling-factorial-int n x)
    (let [x+ (inc x)]
      (/ (Gamma/gamma x+)
         (Gamma/gamma (- x+ n))))))

(defn rising-factorial-int
  "Rising (Pochhammer) factorial for integer n."
  ^double [^long n ^double x]
  (if (not-neg? n)
    (loop [i (long 0)
           v 1.0]
      (if (== i n) v
          (recur (inc i) (* v (+ x i)))))
    (/ (rising-factorial-int (- n) (+ x n)))))

(defn rising-factorial
  "Rising (Pochhammer) factorial."
  ^double [^double n ^double x]
  (if (integer? n)
    (rising-factorial-int n x)
    (/ (Gamma/gamma (+ x n))
       (Gamma/gamma x))))

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
                    (Beta/logBeta (inc (- n k)) (inc k)))))))

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
               (Beta/logBeta (inc (- n k)) (inc k))))))

;; Square and cubic
(defn sq "Same as [[pow2]]. \\\\(x^2\\\\)"
  {:inline (fn [x] `(let [x# (double ~x)] (* x# x#)))
   :inline-arities #{1}}
  ^double [^double x] (* x x))

(defn pow2 "Same as [[sq]]. \\\\(x^2\\\\)"
  {:inline (fn [x] `(let [x# (double ~x)] (* x# x#)))
   :inline-arities #{1}}
  ^double [^double x] (* x x))

(defn pow3 "\\\\(x^3\\\\)"
  {:inline (fn [x] `(let [x# (double ~x)] (* x# x# x#)))
   :inline-arities #{1}}
  ^double [^double x](* x x x))

(defn cb "\\\\(x^3\\\\)"
  {:inline (fn [x] `(let [x# (double ~x)] (* x# x# x#)))
   :inline-arities #{1}}
  ^double [^double x] (* x x x))

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
  "Hypot.
  See also [[hypot-sqrt]]."
  {:inline (fn ([x y] `(. FastMath (hypot (double ~x) (double ~y))))
             ([x y z] `(. FastMath (hypot (double ~x) (double ~y) (double ~z)))))
   :inline-arities #{2 3}}
  (^double [^double x ^double y]
   (FastMath/hypot x y))
  (^double [^double x ^double y ^double z]
   (FastMath/hypot x y z)))

(defn hypot-sqrt
  "Hypot, sqrt version: \\\\(\\sqrt{x^2+y^2}\\\\) or \\\\(\\sqrt{x^2+y^2+z^2}\\\\)."
  {:inline (fn ([x y] `(. FastMath (sqrt (+ (sq ~x) (sq ~y)))))
             ([x y z] `(. FastMath (sqrt (+ (sq ~x) (sq ~y) (sq ~z))))))
   :inline-arities #{2 3}}
  (^double [^double x ^double y]
   (FastMath/sqrt (+ (* x x) (* y y))))
  (^double [^double x ^double y ^double z]
   (FastMath/sqrt (+ (* x x) (* y y) (* z z)))))

;; distance
(defn dist
  "Euclidean distance between points `(x1,y1)` and `(x2,y2)`. See [[fastmath.vector]] namespace to see other metrics which work on vectors."
  {:inline (fn [x1 y1 x2 y2] `(hypot-sqrt (- ~x2 ~x1) (- ~y2 ~y1)))
   :inline-arities #{4}}
  (^double [[^double x1 ^double y1] [^double x2 ^double y2]] (dist x1 y1 x2 y2))
  (^double [^double x1 ^double y1 ^double x2 ^double y2]
   (FastMath/sqrt (+ (sq (- x2 x1)) (sq (- y2 y1))))))

(defn qdist
  "Quick version of Euclidean distance between points. [[qsqrt]] is used instead of [[sqrt]]."
  {:inline (fn [x1 y1 x2 y2] `(. FastMath (sqrtQuick (+ (sq (- ~x2 ~x1)) (sq (- ~y2 ~y1))))))
   :inline-arities #{4}}
  (^double [[^double x1 ^double y1] [^double x2 ^double y2]] (qdist x1 y1 x2 y2))
  (^double [^double x1 ^double y1 ^double x2 ^double y2]
   (FastMath/sqrtQuick (+ (sq (- x2 x1)) (sq (- y2 y1))))))

;; Rounding functions
(defn floor
  "\\\\(\\lfloor x \\rfloor\\\\). See: [[qfloor]].

  Rounding is done to a multiply of scale value (when provided)."
  {:inline (fn ([x] `(. FastMath (floor (double ~x))))
             ([x scale] `(* (. FastMath (floor (double (/ ~x ~scale)))) ~scale)))
   :inline-arities #{1 2}}
  (^double [^double x] (FastMath/floor x))
  (^double [^double x ^double scale] (* (FastMath/floor (/ x scale)) scale)))

(defn ceil
  "\\\\(\\lceil x \\rceil\\\\). See: [[qceil]].

  Rounding is done to a multiply of scale value (when provided)."
  {:inline (fn ([x] `(. FastMath (ceil (double ~x))))
             ([x scale] `(* (. FastMath (ceil (double (/ ~x ~scale)))) ~scale)))
   :inline-arities #{1 2}}
  (^double [^double x] (FastMath/ceil x))
  (^double [^double x ^double scale] (* (FastMath/ceil (/ x scale)) scale)))

(defn round "Round to `long`. See: [[rint]], [[qround]]."
  {:inline (fn [x] `(. FastMath (round (double ~x))))
   :inline-arities #{1}} 
  ^long [^double x] (FastMath/round x))

(defn rint
  "Round to `double`. See [[round]], [[qround]].

  Rounding is done to a multiply of scale value (when provided)."
  {:inline (fn ([x] `(. FastMath (rint (double ~x))))
             ([x scale] `(* (. FastMath (rint (double (/ ~x ~scale)))) ~scale)))
   :inline-arities #{1 2}}
  (^double [^double x] (FastMath/rint x))
  (^double [^double x ^double scale] (* (FastMath/rint (/ x scale)) scale)))

(defn round-even
  "Round evenly (like in round in R), IEEE / IEC rounding"
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
  "\\\\(|x|\\\\) - `double` version. See [[iabs]]."
  {:inline (fn [x] `(. FastMath (abs (double ~x))))
   :inline-arities #{1}}
  ^double [^double x] (FastMath/abs x))

(defn iabs
  "\\\\(|x|\\\\) - `long` version. See [[abs]]."
  {:inline (fn [x] `(. Math (abs (long ~x))))
   :inline-arities #{1}}
  ^long [^long x] (Math/abs x))

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
  ([^double a ^double b] (or (== (approx a) (approx b)) (== a b)))
  ([^double a ^double b ^long digits] (or (== (approx a digits)
                                              (approx b digits))
                                          (== a b))))

(defn delta-eq
  "Checks equality for given absolute accuracy (default `1.0e-6`).

  Version with 4-arity accepts absolute and relative accuracy."
  ([^double a ^double b] (delta-eq a b 1.0e-6))
  ([^double a ^double b ^double accuracy]
   (or (< (abs (- a b)) accuracy) (== a b)))
  ([^double a ^double b ^double abs-tol ^double rel-tol]
   (or (< (abs (- a b)) (max abs-tol (* rel-tol (max (abs a) (abs b))))) (== a b))))

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
  ^long [^double x] (-> x abs log2 floor unchecked-long))

(defn high-2-exp
  "Find lowest exponent (power of 2) which is greater or equal `x`. See [[low-2-exp]]."
  ^long [^double x] (-> x abs log2 ceil unchecked-long))

(defn low-exp
  "Find greatest exponent for base `b` which is lower or equal `x`. See also [[high-exp]]."
  ^long [^double b ^double x] (->> x abs (logb b) floor unchecked-long))

(defn high-exp
  "Find lowest exponent for base `b` which is higher or equal`x`. See also [[low-exp]]."
  ^long [^double b ^double x] (->> x abs (logb b) ceil unchecked-long))

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
  {:inline (fn [v] `(. FastMath (nextUp (double ~v))))
   :inline-arities #{1}}
  (^double [^double v]
   (FastMath/nextUp v))
  (^double [^double v ^long delta]
   (nth (iterate next-double v) delta)))

(defn prev-double
  "Next double value. Optional value `delta` sets step amount."
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

(defn ulp
  "Unit in the Last Place, distance between next value larger than `x` and `x`"
  {:inline (fn [x] `(. FastMath (ulp (double ~x))))
   :inline-arities #{1}}
  ^double [^double x]  (FastMath/ulp x))

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
(def ^{:const true :tag 'double :doc "\\\\(\\frac{1}{\\sqrt{2\\pi}}\\\\)"} INV_SQRT2PI (/ 1.0 SQRT2PI))
(def ^{:const true :tag 'double :doc "\\\\(\\frac{1}{\\sqrt\\pi}\\\\)"} INV_SQRTPI (/ 1.0 SQRTPI))
(def ^{:const true :tag 'double :doc "\\\\(\\sqrt{\\frac{2}{\\pi}}\\\\)"} SQRT_2_PI (sqrt M_2_PI))
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
  {:inline (fn [value] `(let [v# (double ~value)]
                         (if (pos? v#) 1.0 (if (neg? v#) -1.0 0.0))))
   :inline-arities #{1}}
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
  {:inline (fn [value] `(if (neg? (double ~value)) -1.0 1.0))
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

(defmacro constrain
  "Clamp `value` to the range `[mn,mx]`."
  [value mn mx]
  `(max (min ~value ~mx) ~mn))

(defn norm
  "Normalize `v` from the range `[start,stop]` to the range `[0,1]` or map `v` from the range `[start1,stop1]` to the range `[start2,stop2]`. See also [[make-norm]]."
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
  "Make [[norm]] function for given range. Resulting function accepts `double` value (with optional target `[dstart,dstop]` range) and returns `double`."
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
        ^double den (reduce + eaxs)]
    (reduce + (map (fn [^double x ^double eax]
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
     :p-norm (pow (reduce + (map (fn [^double x]
                                   (pow (abs x) alpha)) xs)) (/ alpha))
     :smu (let [epsilon (/ alpha)]
            (reduce (fn [^double a ^double b]
                      (* 0.5 (+ a b (sqrt (+ (sq (- a b))
                                             epsilon))))) xs)))))

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
         eps (* 0.5 (double (if (seq diffs) (reduce min diffs) 0.0)))]
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
                      (fn ^double [v] (/ ^double (reduce + (map first v)) (count v))))
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
