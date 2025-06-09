### You'll be writing a documentation (book chapters and function docstrings) for `fastmath.stats` namespace. 

### Please write a short paragraph for `clay/complex_quaternions.clj` book. Do not copy content of the `clay/complex_quaternions.clj` file, write only asked chapter. Put text as a Clojure comment block. Put actual Clojure code without comment block and markdown wrapper.  Follow the convention of other chapters.

- Write a paragraph for the "## Other" chapter.
- Describe all mentioned functions, group them in points and describe them in brief
- When possible add latex formulas
- Illustrate with examples
- Illustrate with charts if applicable

Write only asked chapter. Refer to previous chapters as examples of chapter structure.




```clojure
;; ## Other

;; This section contains a collection of miscellaneous utility functions that don't fit neatly into the preceding categories, including integer number theory functions, boolean operations, data structure conversions, and error calculation utilities.

;; ::: {.callout-tip title="Defined functions"}
;; * `gcd`, `lcm`
;; * `bool-not`, `bool-xor`, `xor`
;; * `identity-long`, `identity-double`
;; * `relative-error`, `absolute-error`
;; * `seq->double-array`, `seq->double-double-array`, `double-array->seq`, `double-double-array->seq`
;; * `use-primitive-operators`, `unuse-primitive-operators`
;; :::

;; ### Number Theory

;; Functions for calculating the greatest common divisor (GCD) and least common multiple (LCM) of integers.

;; * `gcd(a, b)`: Computes the greatest common divisor of two `long` integers `a` and `b`. Uses the binary GCD algorithm (Stein's algorithm), which is typically faster than the Euclidean algorithm for integers.
;; $$ \operatorname{gcd}(a, b) = \max \{ d \in \mathbb{Z}^+ : d \mid a \land d \mid b \} $$
;; * `lcm(a, b)`: Computes the least common multiple of two `long` integers `a` and `b`.
;; $$ \operatorname{lcm}(a, b) = \frac{|a \cdot b|}{\operatorname{gcd}(a, b)} $$

(utls/examples-note
  (m/gcd 12 18)
  (m/gcd 35 49)
  (m/gcd 17 23)
  (m/lcm 12 18)
  (m/lcm 35 49)
  (m/lcm 17 23))

;; ### Boolean and Identity Utilities

;; Basic boolean operations and identity functions for primitive types.

;; * `bool-not(x)`: Primitive boolean NOT. Returns `true` if `x` is logically false, `false` otherwise.
;; * `bool-xor(x, y)` / `xor(x, y)`: Primitive boolean XOR (exclusive OR). Returns `true` if exactly one of `x` or `y` is logically true. Accepts multiple arguments, chaining the XOR operation.
;; * `identity-long(x)`: Returns its `long` argument unchanged. Useful for type hinting or guaranteeing a `long` primitive.
;; * `identity-double(x)`: Returns its `double` argument unchanged. Useful for type hinting or guaranteeing a `double` primitive.

(utls/examples-note
  (m/bool-not true)
  (m/bool-not false)
  (m/xor true false)
  (m/xor true true)
  (m/xor true false true)
  (m/identity-long 5)
  (m/identity-double 5.0))

;; ### Error Calculation

;; Functions to compute the difference between a value and its approximation.

;; * `absolute-error(v, v-approx)`: Computes the absolute difference between a value `v` and its approximation `v-approx`.
;; $$ \operatorname{absolute-error}(v, v_{\text{approx}}) = |v - v_{\text{approx}}| $$
;; * `relative-error(v, v-approx)`: Computes the relative difference between a value `v` and its approximation `v-approx`.
;; $$ \operatorname{relative-error}(v, v_{\text{approx}}) = \left| \frac{v - v_{\text{approx}}}{v} \right| $$

(utls/examples-note
  (m/absolute-error 10.0 10.01)
  (m/relative-error 10.0 10.01)
  (m/absolute-error 1.0e-6 1.01e-6)
  (m/relative-error 1.0e-6 1.01e-6))

;; ### Array and Sequence Conversion

;; Utilities for converting between Clojure sequences and primitive Java arrays (`double[]` and `double[][]`). These are often necessary for interoperation with Java libraries or for performance-critical operations.

;; * `seq->double-array(vs)`: Converts a sequence `vs` into a `double[]` array. If `vs` is already a `double[]`, it is returned directly. If `vs` is a single number, returns a `double[]` of size 1 containing that number.
;; * `double-array->seq(res)`: Converts a `double[]` array `res` into a sequence. This is an alias for Clojure's built-in `seq` function which works correctly for Java arrays.
;; * `seq->double-double-array(vss)`: Converts a sequence of sequences `vss` into a `double[][]` array. If `vss` is already a `double[][]`, it is returned directly.
;; * `double-double-array->seq(res)`: Converts a `double[][]` array `res` into a sequence of sequences.

(utls/examples-note
  (m/seq->double-array [1 2 3])
  (m/double-array->seq (double-array [1.0 2.0 3.0]))
  (m/seq->double-double-array [[1 2] [3 4]])
  (m/double-double-array->seq (into-array (map double-array [[1 2] [3 4]]))))

;; ### Primitive Operators Toggle

;; Macros to replace or restore core Clojure math functions with `fastmath.core` primitive-specialized versions. This can improve performance but should be used carefully as the `fastmath.core` versions have specific type and return value behaviors.

;; * `use-primitive-operators`: Replaces a select set of `clojure.core` functions (like `+`, `-`, `*`, `/`, comparison operators, bitwise operators, `inc`, `dec`, predicates) with their primitive-specialized `fastmath.core` macro equivalents within the current namespace.
;; * `unuse-primitive-operators`: Reverts the changes made by `use-primitive-operators`, restoring the original `clojure.core` functions in the current namespace. Recommended practice is to call this at the end of any namespace that calls `use-primitive-operators` (especially important in Clojure 1.12+).

;; Example usage (typically done at the top or bottom of a namespace):

```clojure
;; (ns my-namespace (:require [fastmath.core :as m]))
;; (m/use-primitive-operators)
;; ... your code using primitive math ...
;; (m/unuse-primitive-operators) ;; or at the end of the file
```

;; ## Constants

(utls/gen-constants 'fastmath.core #{'jvm-version})

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.core)
```

### 

<!-- Local Variables: -->
<!-- gptel-model: gemini-2.5-flash-preview-04-17 -->
<!-- gptel--backend-name: "Gemini" -->
<!-- gptel--bounds: ((response (752 6084))) -->
<!-- End: -->
