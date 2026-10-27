;; # Namespace scope
;;
;; Functions to manipulate Vec2 as Complex numbers.
;; Implementation based on Apache Commons Math, https://ece.uwaterloo.ca/~dwharder/C++/CQOST/src/Complex.cpp
;; and, for the inverse trigonometric and hyperbolic functions, on the algorithms of
;; Hull, Fairgrieve and Tang (1997) as implemented in the FreeBSD msun library (catrig.c)
;;
;; The functions `asinh`, `acos`, `atanh` and their helpers (`hard-work`, `f-helper`, `sum-squares`,
;; `real-part-reciprocal`, the constants and the branch conditions) are based on FreeBSD
;; `lib/msun/src/catrig.c`, which is distributed under the following licence:
;;
;; Copyright (c) 2012 Stephen Montgomery-Smith <stephen@FreeBSD.ORG>
;; All rights reserved.
;;
;; Redistribution and use in source and binary forms, with or without
;; modification, are permitted provided that the following conditions
;; are met:
;; 1. Redistributions of source code must retain the above copyright
;;    notice, this list of conditions and the following disclaimer.
;; 2. Redistributions in binary form must reproduce the above copyright
;;    notice, this list of conditions and the following disclaimer in the
;;    documentation and/or other materials provided with the distribution.
;;
;; THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
;; ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
;; IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
;; ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
;; FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
;; DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
;; OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
;; HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
;; LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
;; OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
;; SUCH DAMAGE.

(ns fastmath.complex
  "Complex numbers operations.

  Complex numbers are represented using the [[Vec2]] type defined in the `fastmath.vector` namespace.
  A complex number $z = a + bi$ corresponds to the vector `(Vec2. a b)`.

  To create a complex number, use [[complex]], [[vec2]], or convert from sequences/arrays using functions like [[fastmath.vector/seq->vec2]].

  The implementation correctly handles floating-point special values such as `##Inf`, `##NaN`, and distinguishes between `+0.0` and `-0.0` where necessary. Branch cuts follow C99 (Annex G) and the signed zero of the imaginary part selects the side of the cut: the result for `x+0.0i` is the limit from above, for `x-0.0i` from below.

  This namespace provides standard complex number operations like arithmetic (+, -, *, /), exponentiation (pow, exp, log), trigonometric functions (sin, cos, tan, etc.), and their inverses, as well as utility functions like `abs`, `arg`, `conjugate`, etc.

  The functions are accurate over the whole range of doubles: they do not overflow or underflow in intermediate results and do not lose the small component near zero, near the unit circle and near the branch points. The exception are `asec`, `acsc`, `acot`, `asech`, `acsch` and `acoth`, computed from the rounded reciprocal `1/z`: close to their branch points they are less accurate (see the functions)."
  (:refer-clojure :exclude [abs zero?])
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2 Vec4]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^{:doc "The imaginary unit, z=0+1i."} I (Vec2. 0.0 1.0))
(def ^{:doc "Minus the imaginary unit, z=0-1i. Alias of [[-I]]."} I- (Vec2. 0.0 -1.0))
(def ^{:doc "Minus the imaginary unit, z=0-1i."} -I (Vec2. 0.0 -1.0))
(def ^{:doc "The real unit, z=1+0i."} ONE (Vec2. 1.0 0.0))
(def ^{:doc "The real number two, z=2+0i."} TWO (Vec2. 2.0 0.0))
(def ^{:doc "The complex zero, z=0+0i."} ZERO (Vec2. 0.0 0.0))
(def ^{:doc "The real number pi as a complex number, z=pi+0i."} PI (Vec2. m/PI 0.0))

;; constants for the algorithms below

(def ^{:private true :const true :tag 'double} eps m/MACHINE-EPSILON2)
(def ^{:private true :const true :tag 'double} recip-eps 4.503599627370496E15)
(def ^{:private true :const true :tag 'double} sqrt-6-eps-4 9.125060374972143E-9) ;; sqrt(6 eps) / 4
(def ^{:private true :const true :tag 'double} sqrt-3-eps-2 1.2904784139758925E-8) ;; sqrt(3 eps) / 2
(def ^{:private true :const true :tag 'double} four-sqrt-min 5.966672584960166E-154) ;; 4 sqrt(Double/MIN_NORMAL)
(def ^{:private true :const true :tag 'double} sqrt-min 1.4916681462400413E-154)
(def ^{:private true :const true :tag 'double} a-crossover 10.0)
(def ^{:private true :const true :tag 'double} b-crossover 0.6417)

(defn complex
  "Creates a complex number represented as a [[Vec2]].

  Takes optional real and imaginary parts. Without arguments returns `0+0i`, with one argument the imaginary part is `0.0`.

  Examples:

  (complex 1 2) ;; => #vec2 [1.0, 2.0]
  (complex 3)   ;; => #vec2 [3.0, 0.0]
  (complex)     ;; => #vec2 [0.0, 0.0]"
  ([^double a ^double b] (Vec2. a b))
  ([^double a] (Vec2. a 0.0))
  ([] ZERO))

(defn ensure-complex
  "Converts `v` to a complex number.

  A number (including a ratio) becomes a complex number with zero imaginary part, an Apache Commons Math `Complex` is converted to [[Vec2]] and any other value is returned unchanged (a [[Vec2]] stays itself)."
  [v]
  (cond
    (number? v) (Vec2. (double v) 0.0)

    (instance? org.apache.commons.math3.complex.Complex v)
    (Vec2. (.getReal ^org.apache.commons.math3.complex.Complex v)
           (.getImaginary ^org.apache.commons.math3.complex.Complex v))

    :else v))

(defn re "Real part of `z` as a double." ^double [^Vec2 z] (.x z))
(defn im "Imaginary part of `z` as a double." ^double [^Vec2 z] (.y z))

(defn real?
  "Returns `true` when the imaginary part of `z` is zero (`+0.0` or `-0.0`)."
  [^Vec2 z] (m/zero? (.y z)))

(defn imaginary?
  "Returns `true` when the real part of `z` is zero (`+0.0` or `-0.0`), a pure imaginary number (zero included)."
  [^Vec2 z] (m/zero? (.x z)))

(defn zero?
  "Returns `true` when both parts of `z` are zero (the sign of a zero does not matter)."
  [^Vec2 z] (v/is-zero? z))

(defn inf?
  "Returns `true` when any part of `z` is infinite (also when the other part is `NaN`)."
  [^Vec2 z] (or (m/inf? (.x z))
                (m/inf? (.y z))))

(defn nan?
  "Returns `true` when any part of `z` is `NaN`."
  [^Vec2 z] (or (m/nan? (.x z))
                (m/nan? (.y z))))

(defn invalid?
  "Returns `true` when `z` has an infinite or `NaN` part. See also [[valid?]]."
  [z] (or (inf? z) (nan? z)))

(defn valid?
  "Returns `true` when both parts of `z` are finite numbers. See also [[invalid?]]."
  [z] (not (invalid? z)))

(defn delta-eq
  "Compares complex numbers `q1` and `q2` with the given absolute `accuracy` (`1.0e-6` by default) in each part.

  The comparison is absolute, so use a relative check for large magnitudes. It returns `false` for infinite parts, because `Inf - Inf` is `NaN`.

  See also [[fastmath.vector/delta-eq]]."
  ([q1 q2]
   (v/delta-eq q1 q2))
  ([q1 q2 ^double accuracy]
   (v/delta-eq q1 q2 accuracy)))

(defn csgn
  "Complex signum function.

  Takes `z` or the parts `re` and `im`. Returns `0.0` for the zero. For any other number returns the sign of the real part (`-1.0` or `1.0`), not the unit complex number.

  See also [[arg]], [[normalize]]."
  (^double [^double re ^double im]
   (if (and (m/zero? re)
            (m/zero? im)) 0.0 (m/sgn re)))
  (^double [^Vec2 z]
   (if (zero? z) 0.0 (m/sgn (.x z)))))

(defn flip
  "Exchanges the real and imaginary parts: `a+bi` becomes `b+ai`. It is `i * conjugate(z)`."
  ^Vec2 [^Vec2 z] (Vec2. (.y z) (.x z)))

(defn abs
  "Calculates the magnitude (absolute value) of the complex number `z`.

  The magnitude is computed without intermediate overflow or underflow, so it is correct for the whole range of doubles (for example `1e200+1e200i` or `1e-200+1e-200i`). It is `Infinity` when any part is infinite, also when the other one is `NaN`.

  See also [[norm]] (the squared magnitude)."
  ^double [^Vec2 z]
  (m/hypot (.x z) (.y z)))

(defn norm
  "Calculates the squared magnitude (norm) of the complex number `z`, `a^2+b^2`.

  Also known as the squared modulus or squared Euclidean norm. Note that [[fastmath.quaternion/norm]] is the magnitude (not squared). The value overflows to `Infinity` when `a^2+b^2` is not representable.

  See also [[abs]], [[abs-sq]]."
  ^double [^Vec2 z]
  (if (inf? z)
    ##Inf
    (m/+ (m/sq (.x z)) (m/sq (.y z)))))

(defn abs-sq
  "Squared magnitude of `z`, alias of [[norm]]."
  ^double [^Vec2 z]
  (norm z))

(defn normalize
  "Normalizes `z` to the unit magnitude, `z / abs(z)`. The zero stays zero.

  See also [[abs]], [[csgn]]."
  ^Vec2 [^Vec2 z] (v/normalize z))

(defn arg
  "Argument (phase angle) of the complex number `z`, `atan2(im, re)`, in the range `[-pi, pi]`.

  The sign of a zero imaginary part is respected: `arg(-1-0i)` is `-pi`, `arg(-1+0i)` is `pi`.

  See also [[abs]], [[log]]."
  ^double [^Vec2 z]
  (m/atan2 (.y z) (.x z)))

(defn proj
  "Projection onto the Riemann sphere (C99 `cproj`).

  Any number with an infinite part is mapped to the single point at infinity, `Infinity` as the real part and a zero imaginary part with the sign of the original imaginary part. Other numbers (also `NaN`) are returned unchanged.

  See also [[inf?]]."
  ^Vec2 [^Vec2 z]
  (if (inf? z)
    (Vec2. ##Inf (m/copy-sign 0.0 (.y z)))
    z))

(defn add
  "Sum of complex numbers.

  Without arguments returns `0`, with one argument the argument itself. More than two arguments are added from left to right."
  (^Vec2 [] ZERO)
  (^Vec2 [^Vec2 z] z)
  (^Vec2 [^Vec2 z1 ^Vec2 z2]
   (Vec2. (m/+ (.x z1) (.x z2))
          (m/+ (.y z1) (.y z2))))
  (^Vec2 [z1 z2 & more]
   (reduce add (add z1 z2) more)))

(defn adds
  "Adds the real number `v` to the complex number `z`."
  ^Vec2 [^Vec2 z ^double v]
  (Vec2. (m/+ (.x z) v) (.y z)))

(defn neg
  "Negation of the complex number, `-z`."
  ^Vec2 [^Vec2 z]
  (Vec2. (m/- (.x z)) (m/- (.y z))))

(defn sub
  "Difference of complex numbers, `z1 - z2 - ...`.

  With one argument returns the negation [[neg]]. More than two arguments are subtracted from the first one, from left to right."
  (^Vec2 [^Vec2 z] (neg z))
  (^Vec2 [^Vec2 z1 ^Vec2 z2]
   (Vec2. (m/- (.x z1) (.x z2))
          (m/- (.y z1) (.y z2))))
  (^Vec2 [z1 z2 & more]
   (reduce sub (sub z1 z2) more)))

(defn scale
  "Multiplication of the complex number `z` by the real number `v`."
  ^Vec2 [^Vec2 z ^double v]
  (Vec2. (m/* (.x z) v)
         (m/* (.y z) v)))

(defn mult
  "Product of complex numbers.

  Without arguments returns `1`, with one argument the argument itself. More than two arguments are multiplied from left to right.

  See also [[muladd]], [[scale]], [[mult-I]]."
  (^Vec2 [] ONE)
  (^Vec2 [^Vec2 z] z)
  (^Vec2 [^Vec2 z1 ^Vec2 z2]
   (let [a (.x z1)
         b (.y z1)
         c (.x z2)
         d (.y z2)]
     (Vec2. (m/- (m/* a c) (m/* b d))
            (m/+ (m/* a d) (m/* b c)))))
  (^Vec2 [z1 z2 & more]
   (reduce mult (mult z1 z2) more)))

(defn mult-I
  "Multiplication by `0+i`, `i*z`. Exact, a rotation by 90 degrees."
  ^Vec2 [^Vec2 z] (Vec2. (m/- (.y z)) (.x z)))

(defn mult-I-
  "Multiplication by `0-i`, `-i*z`. Exact, a rotation by -90 degrees."
  ^Vec2 [^Vec2 z] (Vec2. (.y z) (m/- (.x z))))

(defn muladd
  "Multiply-add of complex numbers, `z + x*y`.

  See also [[mult]], [[add]]."
  ^Vec2 [x y z]
  (add z (mult x y)))

(defn sq
  "Square of the complex number, `z^2`."
  ^Vec2 [^Vec2 z]
  (let [a (.x z)
        b (.y z)]
    (Vec2. (m/- (m/* a a) (m/* b b))
           (m/* 2.0 a b))))

(defn conjugate
  "Complex conjugate, `a-bi` for `a+bi`."
  ^Vec2 [^Vec2 z]
  (Vec2. (.x z) (m/- (.y z))))

(defn div
  "Division of two complex numbers, `z1 / z2`.

  Smith's algorithm: the result does not overflow or underflow in the middle when the quotient is representable (for example `1e200/1e200` or `1e-200/1e-200` with complex parts). Division by zero gives `NaN` parts.

  See also [[reciprocal]], [[mult]]."
  ^Vec2 [^Vec2 z1 ^Vec2 z2]
  (let [a (.x z1)
        b (.y z1)
        c (.x z2)
        d (.y z2)]
    (if (m/>= (m/abs c) (m/abs d))
      (let [r (m// d c)
            den (m/+ c (m/* d r))]
        (Vec2. (m// (m/+ a (m/* b r)) den)
               (m// (m/- b (m/* a r)) den)))
      (let [r (m// c d)
            den (m/+ d (m/* c r))]
        (Vec2. (m// (m/+ (m/* a r) b) den)
               (m// (m/- (m/* b r) a) den))))))

(defn reciprocal
  "Reciprocal, `1/z`.

  Computed without intermediate overflow or underflow, correct for the whole range of doubles. The reciprocal of zero has infinite parts with the signs of the zeros (`1/(+0+0i)` is `Inf-Inf*i`), the reciprocal of an infinite number is zero.

  See also [[div]]."
  ^Vec2 [^Vec2 z]
  (let [a (.x z)
        b (.y z)]
    (cond
      (zero? z) (Vec2. (m/copy-sign ##Inf a)
                       (m/- (m/copy-sign ##Inf b)))
      (inf? z) (Vec2. (m/copy-sign 0.0 a)
                      (m/- (m/copy-sign 0.0 b)))
      (nan? z) (Vec2. ##NaN ##NaN)
      :else (if (m/>= (m/abs a) (m/abs b))
              (let [r (m// b a)
                    den (m/+ a (m/* b r))]
                (Vec2. (m// 1.0 den) (m// (m/- r) den)))
              (let [r (m// a b)
                    den (m/+ b (m/* a r))]
                (Vec2. (m// r den) (m// -1.0 den)))))))

(defn sqrt
  "Principal square root of the complex number, `sqrt(z)`.

  The real part of the result is not negative; the imaginary part has the sign of the imaginary part of `z` (the sign of a zero imaginary part included). The root is computed without cancellation and without intermediate overflow or underflow, so it is accurate to a few ulps over the whole range of doubles (for example `1e200+1e200i`, `1e-200+1e-200i`, `-1e-10+1e-30i`). The zero is returned as `0` with the imaginary part of the argument. Special values (C99): an infinite imaginary part gives `Infinity+Infinity*i` (with the sign of the imaginary part) whatever the real part is, also `NaN`.

  See also [[sqrt1z]], [[cbrt]], [[nth-root]], [[abs]]."
  ^Vec2 [^Vec2 z]
  (cond
    (zero? z) (Vec2. 0.0 (.y z))
    (m/inf? (.y z)) (Vec2. ##Inf (.y z)) ;; C99: whatever the real part is, also NaN
    :else
    (let [big? (m/> (m/max (m/abs (.x z)) (m/abs (.y z))) 1.0e300)
          ;; a quarter of a huge number keeps |a| + hypot(a, b) finite
          k (if big? 0.25 1.0)
          back (if big? 2.0 1.0)
          a (m/* k (.x z))
          b (m/* k (.y z))
          t (m/sqrt (m/* 0.5 (m/+ (m/abs a) (m/hypot a b))))]
      (if (m/>= a 0.0)
        (Vec2. (m/* back t) (m/* back (m// b (m/* 2.0 t))))
        (Vec2. (m/* back (m// (m/abs b) (m/* 2.0 t))) (m/* back (m/copy-sign t b)))))))

(defn sqrt1z
  "Computes `sqrt(1-z^2)`, the companion of [[asin]] and [[acos]] (`sin(acos z)`)."
  ^Vec2 [z]
  (->> z
       (mult z)
       (sub ONE)
       (sqrt)))

;; exp

(defn exp
  "Complex exponential, `exp(z) = exp(a) (cos(b) + i sin(b))`.

  Entire function, no branch cut. For a real argument returns a real number with the sign of the zero imaginary part kept. The result does not overflow too early: `exp(710+i)` is finite in the parts which are representable.

  Special values (C99): an infinite real part with a non finite imaginary part gives `Infinity+NaN*i` for `+Infinity` and a zero for `-Infinity`.

  See also [[log]], [[expm1]], [[cis]], [[pow]]."
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/exp x) y)
      (m/zero? x) (Vec2. (m/cos y) (m/sin y))
      ;; C99: an infinite real part and a non finite imaginary part
      (and (m/inf? x) (or (m/inf? y) (m/nan? y))) (if (m/pos? x)
                                                    (Vec2. ##Inf ##NaN)
                                                    (Vec2. -0.0 (m/copy-sign 0.0 y)))
      (m/> x 700.0) (let [h (m/exp (m/* 0.5 x))]
                      (Vec2. (m/* h (m/cos y) h)
                             (m/* h (m/sin y) h)))
      :else (let [e (m/exp x)]
              (Vec2. (m/* e (m/cos y))
                     (m/* e (m/sin y)))))))

(defn cis
  "Returns `cos(theta) + i sin(theta)`, a point on the unit circle for the real angle `theta`.

  See also [[exp]], [[arg]]."
  ^Vec2 [^double theta]
  (Vec2. (m/cos theta) (m/sin theta)))

(defn expm1
  "Computes `exp(z) - 1` accurately for `z` near zero.

  The real part is computed as `expm1(a) cos(b) - 2 sin^2(b/2)`, so it does not cancel.

  See also [[exp]], [[log1p]]."
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/expm1 x) y)
      (m/> x 700.0) (let [w (exp z)] (Vec2. (m/dec (.x w)) (.y w)))
      :else (let [s (m/sin (m/* 0.5 y))
                  w (exp z)]
              (Vec2. (m/- (m/* (m/expm1 x) (m/cos y)) (m/* 2.0 s s))
                     (.y w))))))

;; log

(defn- log-abs
  "log of the magnitude, log|a+bi|, without overflow and underflow, accurate near the unit circle."
  ^double [^double a ^double b]
  (cond
    (or (m/inf? a) (m/inf? b)) ##Inf
    (or (m/nan? a) (m/nan? b)) ##NaN
    :else
    (let [x (m/abs a)
          y (m/abs b)
          hi (m/max x y)
          lo (m/min x y)]
      (if (m/zero? hi)
        ##-Inf
        (let [p (m/* hi hi)
              q (m/* lo lo)
              mag2 (m/+ p q )]
          (if (and (m/>= mag2 0.5) (m/<= mag2 2.0))
            ;; hi^2 + lo^2 - 1 as a double-double, exact enough for any distance to the unit circle
            (let [pe (m/fma hi hi (m/- p))
                  qe (m/fma lo lo (m/- q))
                  s1 (m/dec p)
                  ;; p - 1 is exact for p in [0.5, 2], two-sum otherwise
                  bb1 (m/- s1 p)
                  e1 (m/+ (m/- p (m/- s1 bb1)) (m/- -1.0 bb1))
                  s2 (m/+ s1 q)
                  bb2 (m/- s2 s1)
                  e2 (m/+ (m/- s1 (m/- s2 bb2)) (m/- q bb2))]
              (m/* 0.5 (m/log1p (m/+ s2 e1 e2 pe qe))))
            (let [r (m// lo hi)]
              (m/+ (m/log hi) (m/* 0.5 (m/log1p (m/* r r)))))))))))

(defn log
  "Principal value of the natural logarithm, `log(z) = log|z| + i arg(z)`.

  The imaginary part is in the range `[-pi, pi]`. The branch cut is the negative real axis, the sign of a zero imaginary part selects the side. The real part is computed without overflow or underflow (`1e200+1e200i` gives `461.8`) and without cancellation near the unit circle. `log(0)` is `-Infinity` with the argument of the signed zeros (`0`, `-0` or `+-pi`).

  See also [[log1p]], [[logb]], [[exp]], [[arg]]."
  ^Vec2 [^Vec2 z]
  (let [a (.x z)
        b (.y z)]
    (Vec2. (log-abs a b) (m/atan2 b a))))

(defn log1p
  "Computes `log(1+z)` accurately for `z` near zero.

  Principal value, the branch cut is `z` real and less than -1.

  See also [[log]], [[expm1]]."
  ^Vec2 [^Vec2 z]
  (let [a (.x z)
        b (.y z)]
    (if (m/< (m/max (m/abs a) (m/abs b)) 0.5)
      (Vec2. (m/* 0.5 (m/log1p (m/fma b b (m/* a (m/+ 2.0 a)))))
             (m/atan2 b (m/+ 1.0 a)))
      (log (Vec2. (m/+ 1.0 a) b)))))

(defn logb
  "Logarithm of `z` with the complex base `b`, `log(z)/log(b)`.

  Both arguments can be numbers. See also [[log]]."
  ^Vec2 [z b]
  (div (log (ensure-complex z)) (log (ensure-complex b))))

;; power

(defn- pow-int
  "Integer power by repeated squaring, `n` is not zero."
  ^Vec2 [^Vec2 z ^long n]
  (let [a (.x z)
        b (.y z)]
    (if (m/zero? b)
      ;; real base: no NaN from 0*Inf in the squares, the sign of the zero imaginary part is computed
      (let [s (m/* (m/copy-sign 1.0 b)
                   (if (m/< n 0) -1.0 1.0)
                   (if (and (m/neg? a) (m/even? n)) -1.0 1.0))]
        (Vec2. (m/fpow a n) (m/copy-sign 0.0 s)))
      (loop [r nil
             base (if (m/neg? n) (reciprocal z) z)
             e (m/abs n)]
        (if (m/zero? e)
          r
          (let [r (if (m/odd? e)
                    (if (nil? r) base (mult r base))
                    r)]
            (if (m/one? e)
              r
              (recur r (sq base) (m/bit-shift-right e 1)))))))))

(defn pow
  "Complex power `z1^z2`, the principal value `exp(z2 log(z1))`.

  Integer real exponents with absolute value up to 32 are computed by repeated multiplication, so they are exact for exactly representable results (`(1+2i)^3 = -11-2i`). Other exponents use the logarithm and the exponential, which are free of overflow in the intermediate results.

  Special cases:

  * any non zero base with a zero exponent gives `1`;
  * zero base: `0` when the real part of the exponent is positive, `Infinity+NaN*i` (the complex infinity) when it is negative, `NaN+NaN*i` when it is zero (also for `0^0`).

  See also [[exp]], [[log]], [[sq]], [[nth-root]]."
  ^Vec2 [^Vec2 z1 ^Vec2 z2]
  (let [zr (.x z2)
        zi (.y z2)]
    (cond
      (zero? z1) (cond
                   (m/pos? zr) ZERO
                   (m/neg? zr) (Vec2. ##Inf ##NaN)
                   :else (Vec2. ##NaN ##NaN))
      (zero? z2) ONE
      (and (m/zero? zi)
           (m/== zr (m/rint zr))
           (m/<= (m/abs zr) 32.0)) (pow-int z1 (long zr))
      (and (m/zero? zi)
           (m/zero? (.y z1))
           (m/pos? (.x z1))) (Vec2. (m/pow (.x z1) zr) (m/* (.y z1) zr))
      :else (let [l (log z1)
                  lr (.x l)
                  li (.y l)]
              (exp (Vec2. (m/- (m/* zr lr) (m/* zi li))
                          (m/+ (m/* zr li) (m/* zi lr))))))))

(defn nth-root
  "Principal `n`-th root of `z`, `z^(1/n)` for the integer `n` (positive or negative).

  The root has the argument `arg(z)/n`, the other roots are obtained by multiplying by the roots of unity. The root of zero is zero for `n>0` and `Infinity+NaN*i` for `n<0`, the root of an infinite number is infinite for `n>0` and zero for `n<0`.

  Throws `ex-info` when `n` is zero.

  See also [[sqrt]], [[cbrt]], [[pow]]."
  ^Vec2 [^Vec2 z ^long n]
  (when (m/zero? n)
    (throw (ex-info "n-th root requires a non zero n" {:option :n :value n :expected "not 0"})))
  (cond
    (zero? z) (if (m/pos? n) (Vec2. 0.0 (.y z)) (Vec2. ##Inf ##NaN))
    (inf? z) (if (m/pos? n) (Vec2. ##Inf (m/copy-sign 0.0 (.y z))) ZERO)
    (nan? z) (Vec2. ##NaN ##NaN)
    :else (let [r (m/exp (m// (log-abs (.x z) (.y z)) n))
                t (m// (m/atan2 (.y z) (.x z)) n)]
            (Vec2. (m/* r (m/cos t)) (m/* r (m/sin t))))))

(defn cbrt
  "Principal cube root of the complex number `z`.

  For a non negative real number the result is the real cube root (exact for perfect cubes). The branch cut is the negative real axis.

  See also [[nth-root]], [[sqrt]]."
  ^Vec2 [^Vec2 z]
  (let [h (abs z)]
    (cond
      (and (m/zero? (.y z)) (m/not-neg? (.x z))) (Vec2. (m/cbrt (.x z)) (.y z))
      (or (zero? z) (m/inf? h) (m/nan? h)) (nth-root z 3)
      :else (let [r (m/cbrt h)
                  t (m// (m/atan2 (.y z) (.x z)) 3.0)]
              (Vec2. (m/* r (m/cos t)) (m/* r (m/sin t)))))))

;;;;;; trig

(defn- scaled-exp
  "0.5 exp(|t|) a, where the exponential is split in two factors: no spurious overflow for |t| up to about 1400."
  ^double [^double a ^double t]
  (let [h (m/exp (m/* 0.5 (m/abs t)))]
    (m/* 0.5 h a h)))

(defn sin
  "Complex sine, `sin(a+bi) = sin(a) cosh(b) + i cos(a) sinh(b)`.

  Entire function. A part overflows only when it is not representable (the other part can be finite).

  Special values (C99): an infinite real part with a zero imaginary part gives `NaN` and the zero imaginary part, with two infinite parts `NaN+-Infinity*i`.

  See also [[cos]], [[asin]], [[sinh]]."
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/sin x) (if (or (m/inf? x) (m/nan? x)) y (m/* y (m/cos x))))
      (m/zero? x) (Vec2. x (m/sinh y))
      (and (or (m/inf? x) (m/nan? x)) (m/inf? y)) (Vec2. ##NaN y)
      (m/> (m/abs y) 700.0) (Vec2. (scaled-exp (m/sin x) y)
                                   (m/* (m/copy-sign 1.0 y) (scaled-exp (m/cos x) y)))
      :else (Vec2. (m/* (m/sin x) (m/cosh y))
                   (m/* (m/cos x) (m/sinh y))))))

(defn cos
  "Complex cosine, `cos(a+bi) = cos(a) cosh(b) - i sin(a) sinh(b)`.

  Entire function. A part overflows only when it is not representable (the other part can be finite).

  Special values (C99): an infinite real part with a zero imaginary part gives `NaN` and a zero imaginary part, with an infinite imaginary part `Infinity+NaN*i`.

  See also [[sin]], [[acos]], [[cosh]]."
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/cos x) (if (or (m/inf? x) (m/nan? x))
                                     (m/copy-sign 0.0 (m/* -1.0 (if (m/nan? x) 1.0 (m/copy-sign 1.0 x)) y))
                                     (m/* -1.0 y (m/sin x))))
      (and (or (m/inf? x) (m/nan? x)) (m/inf? y)) (Vec2. ##Inf ##NaN)
      (m/zero? x) (Vec2. (m/cosh y)
                         (m/copy-sign 0.0 (m/* -1.0 (m/copy-sign 1.0 x) y)))
      (m/> (m/abs y) 700.0) (Vec2. (scaled-exp (m/cos x) y)
                                   (m/* -1.0 (m/copy-sign 1.0 y) (scaled-exp (m/sin x) y)))
      :else (Vec2. (m/* (m/cos x) (m/cosh y))
                   (m/* -1.0 (m/sin x) (m/sinh y))))))

(defn sinh
  "Complex hyperbolic sine, `sinh(a+bi) = sinh(a) cos(b) + i cosh(a) sin(b)`.

  Entire function.

  Special values (C99): `+-Infinity` real part with an infinite or `NaN` imaginary part gives `+-Infinity+NaN*i`, `+-0` with a `NaN` imaginary part gives `+-0+NaN*i`.

  See also [[cosh]], [[asinh]], [[sin]]."
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/sinh x) y)
      (m/zero? x) (Vec2. (let [c (m/cos y)]
                           (if (m/nan? c) x (m/copy-sign 0.0 (m/* (m/copy-sign 1.0 x) c))))
                         (m/sin y))
      (and (m/inf? x) (or (m/inf? y) (m/nan? y))) (Vec2. x ##NaN)
      (m/> (m/abs x) 700.0) (Vec2. (m/* (m/copy-sign 1.0 x) (scaled-exp (m/cos y) x))
                                   (scaled-exp (m/sin y) x))
      :else (Vec2. (m/* (m/sinh x) (m/cos y))
                   (m/* (m/cosh x) (m/sin y))))))

(defn cosh
  "Complex hyperbolic cosine, `cosh(a+bi) = cosh(a) cos(b) + i sinh(a) sin(b)`.

  Entire function. The zero parts have the signs as in C99.

  Special values (C99): an infinite real part with an infinite or `NaN` imaginary part gives `Infinity+NaN*i`.

  See also [[sinh]], [[acosh]], [[cos]]."
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/cosh x)
                         (m/copy-sign 0.0 (m/* (m/copy-sign 1.0 x) (m/copy-sign 1.0 y))))
      (m/zero? x) (Vec2. (m/cos y) (m/copy-sign 0.0 (m/* (m/copy-sign 1.0 x) (m/sin y))))
      (and (m/inf? x) (or (m/inf? y) (m/nan? y))) (Vec2. ##Inf ##NaN)
      (m/> (m/abs x) 700.0) (Vec2. (scaled-exp (m/cos y) x)
                                   (m/* (m/copy-sign 1.0 x) (scaled-exp (m/sin y) x)))
      :else (Vec2. (m/* (m/cosh x) (m/cos y))
                   (m/* (m/sinh x) (m/sin y))))))

(defn tanh
  "Complex hyperbolic tangent.

  Computed with the formulas of Hull, Fairgrieve and Tang, which do not overflow for large `|a|` (for `|a|>=22` the real part is `+-1` up to rounding) and for `b` near `pi/2 + k pi`.

  Special values (C99): an infinite real part with a non finite imaginary part gives `+-1` and a zero imaginary part.

  See also [[tan]], [[coth]], [[atanh]]."
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)]
    (cond
      (m/zero? y) (Vec2. (m/tanh x) y)
      (m/zero? x) (Vec2. x (m/tan y))
      (m/nan? x) (Vec2. ##NaN ##NaN)
      (m/>= (m/abs x) 22.0) (let [em (m/exp (m/- (m/abs x)))]
                              (Vec2. (m/copy-sign 1.0 x)
                                     ;; C99: an infinite real part and a non finite imaginary part give a zero (numpy, Julia agree)
                                     (if (and (m/inf? x) (or (m/inf? y) (m/nan? y)))
                                       (m/copy-sign 0.0 y)
                                       (m/* 4.0 (m/sin y) (m/cos y) em em))))
      :else (let [t (m/tan y)
                  beta (m/inc (m/* t t))
                  s (m/sinh x)
                  rho (m/sqrt (m/inc (m/* s s)))
                  denom (m/inc (m/* beta s s))]
              (Vec2. (m// (m/* beta rho s) denom)
                     (m// t denom))))))

(defn tan
  "Complex tangent, `tan(z) = -i tanh(iz)`.

  Does not overflow for large `|b|` (the imaginary part is `+-1` up to rounding). Poles at `pi/2 + k pi` give large values.

  See also [[tanh]], [[cot]], [[atan]]."
  ^Vec2 [^Vec2 z]
  (let [w (tanh (Vec2. (m/- (.y z)) (.x z)))]
    (Vec2. (.y w) (m/- (.x w)))))

(defn coth
  "Complex hyperbolic cotangent, `1/tanh(z)`.

  See also [[tanh]], [[acoth]]."
  ^Vec2 [^Vec2 z]
  (if (zero? z)
    (Vec2. (m/coth (.x z)) (.y z))
    (reciprocal (tanh z))))

(defn cot
  "Complex cotangent, `1/tan(z)`.

  See also [[tan]], [[acot]]."
  ^Vec2 [^Vec2 z]
  (if (zero? z)
    (Vec2. (m/cot (.x z)) (.y z))
    (reciprocal (tan z))))

(defn sec
  "Complex secant, `1/cos(z)`.

  See also [[cos]], [[asec]]."
  ^Vec2 [^Vec2 z]
  (reciprocal (cos z)))

(defn csc
  "Complex cosecant, `1/sin(z)`.

  See also [[sin]], [[acsc]]."
  ^Vec2 [^Vec2 z]
  (if (zero? z)
    (Vec2. (m/csc (.x z)) (m/- (.y z)))
    (reciprocal (sin z))))

(defn sech
  "Complex hyperbolic secant, `1/cosh(z)`.

  See also [[cosh]], [[asech]]."
  ^Vec2 [^Vec2 z]
  (reciprocal (cosh z)))

(defn csch
  "Complex hyperbolic cosecant, `1/sinh(z)`.

  See also [[sinh]], [[acsch]]."
  ^Vec2 [^Vec2 z]
  (if (zero? z)
    (Vec2. (m/csch (.x z)) (m/- (.y z)))
    (reciprocal (sinh z))))

;; inverse functions
;;
;; core: asinh, acos and atanh, the others are derived from them by the identities
;; asin(z) = -i asinh(iz), atan(z) = -i atanh(iz), acosh(z) = +-i acos(z)

(defn- signbit?
  [^double x]
  (m/neg? (m/copy-sign 1.0 x)))

(defn- nan-case
  "Result for a NaN argument of asinh: an infinite real part gives `(+-Inf, NaN)` with the sign of the real part (C99), an infinite imaginary part `(Inf, NaN)`."
  ^Vec2 [^double x ^double y]
  (cond
    (m/inf? x) (Vec2. x ##NaN)
    (m/inf? y) (Vec2. ##Inf ##NaN)
    :else (Vec2. ##NaN ##NaN)))

(defn- f-helper
  "(hypot(a, b) - b) / 2 without cancellation"
  ^double [^double a ^double b ^double h]
  (cond
    (m/neg? b) (m/* 0.5 (m/- h b))
    (m/zero? b) (m/* 0.5 a)
    :else (m// (m// (m/* a a) (m/+ h b)) 2.0)))

(def ^{:private true :const true :tag 'double} eps2-128 (m// (m/* eps eps) 128.0))
(def ^{:private true :const true :tag 'double} two-over-eps (m// 2.0 eps))
(def ^{:private true :const true :tag 'double} four-over-eps2 (m// 4.0 (m/* eps eps)))


(defn- hard-work
  "Core of asinh and acos for non negative `x`, `y` (Hull, Fairgrieve, Tang).

  With `A = (|z+i| + |z-i|)/2`, `B = y/A` returns `Vec4`:
  * `x` - `acosh(A)`,
  * `y` - `1.0` when `asin(B)` is accurate, `0.0` otherwise,
  * `z` - `B` or `sqrt(A^2-y^2)`,
  * `w` - the (scaled) `y`, used with `atan2` when `B` is not usable."
  ^Vec4 [^double x ^double y]
  (let [r (m/hypot x (m/inc y))
        s (m/hypot x (m/dec y))
        a (m/max 1.0 (m/* 0.5 (m/+ r s)))
        rx (if (m/< a a-crossover)
             (cond
               (and (m/one? y) (m/< x eps2-128)) (m/sqrt x)
               (m/>= x (m/* eps (m/abs (m/dec y)))) (let [am1 (m/+ (f-helper x (m/+ 1.0 y) r)
                                                                   (f-helper x (m/- 1.0 y) s))]
                                                      (m/log1p (m/+ am1 (m/sqrt (m/* am1 (m/inc a))))))
               (m/< y 1.0) (m// x (m/sqrt (m/* (m/- 1.0 y) (m/+ 1.0 y))))
               :else (m/log1p (m/+ (m/dec y) (m/sqrt (m/* (m/dec y) (m/inc y))))))
             (m/log (m/+ a (m/sqrt (m/dec (m/* a a))))))]
    (if (m/< y four-sqrt-min)
      (Vec4. rx 0.0 (m/* a two-over-eps) (m/* y two-over-eps))
      (let [b (m// y a)]
        (if (m/> b b-crossover)
          (cond
            (and (m/one? y) (m/< x eps2-128))
            (Vec4. rx 0.0 (m/* (m/sqrt x) (m/sqrt (m/* 0.5 (m/+ a y)))) y)

            (m/>= x (m/* eps (m/abs (m/dec y))))
            (let [amy (m/+ (f-helper x (m/inc y) r)
                           (f-helper x (m/dec y) s))]
              (Vec4. rx 0.0 (m/sqrt (m/* amy (m/+ a y))) y))

            (m/> y 1.0)
            (Vec4. rx 0.0 (m// (m/* x four-over-eps2 y) (m/sqrt (m/* (m/inc y) (m/dec y)))) (m/* y four-over-eps2))

            :else
            (Vec4. rx 0.0 (m/sqrt (m/* (m/- 1.0 y) (m/+ 1.0 y))) y))
          (Vec4. rx 1.0 b y))))))

(defn asinh
  "Inverse hyperbolic sine, the principal value.

  The real part has the sign of `a`, the imaginary part is in `[-pi/2, pi/2]` and has the sign of `b`. Branch cuts are the imaginary axis for `|b|>1`. The algorithm (Hull, Fairgrieve, Tang) is accurate near zero, near the branch points `+-i` and for arguments up to the maximum double.

  Special values (C99): `asinh(NaN+0i)` keeps the zero imaginary part, `+-Infinity` with a `NaN` imaginary part gives `+-Infinity+NaN*i`.

  See also [[sinh]], [[asin]], [[acosh]], [[atanh]]."
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)
        ax (m/abs x)
        ay (m/abs y)]
    (cond
      ;; C99 G.6.2.2: casinh(NaN + 0i) is NaN + 0i; casin, cacos and catanh give NaN + NaN i (numpy and Julia agree)
      (or (m/nan? x) (m/nan? y)) (if (m/zero? y)
                                   (Vec2. ##NaN y)
                                   (nan-case x y))
      (or (m/> ax recip-eps) (m/> ay recip-eps))
      (Vec2. (m/copy-sign (m/+ (log-abs x y) m/LN2) x)
             (m/copy-sign (m/atan2 ay ax) y))
      (or (and (m/zero? x) (m/zero? y))
          (and (m/< ax sqrt-6-eps-4) (m/< ay sqrt-6-eps-4))) z
      :else (let [^Vec4 w (hard-work ax ay)
                  rx (.x w)
                  ry (if (m/> (.y w) 0.5)
                       (m/asin (.z w))
                       (m/atan2 (.w w) (.z w)))]
              (Vec2. (m/copy-sign rx x) (m/copy-sign ry y))))))

(defn asin
  "Inverse sine, the principal value, `asin(z) = -i asinh(iz)`.

  The real part is in `[-pi/2, pi/2]`. Branch cuts are the real axis for `|a|>1`; the sign of a zero imaginary part selects the side (`asin(2+0i)` is `pi/2 + 1.317i`, `asin(2-0i)` is `pi/2 - 1.317i`). Accurate near zero, near `+-1` and for large arguments.

  See also [[sin]], [[acos]], [[asinh]], [[sqrt1z]]."
  ^Vec2 [^Vec2 z]
  (let [w (asinh (Vec2. (.y z) (.x z)))]
    (Vec2. (.y w) (.x w))))

(defn acos
  "Inverse cosine, the principal value.

  The real part is in `[0, pi]`, the imaginary part has the sign opposite to `b`. Branch cuts are the real axis for `|a|>1`; the sign of a zero imaginary part selects the side. Accurate near `1`, `-1` and for large arguments (Hull, Fairgrieve, Tang).

  See also [[cos]], [[asin]], [[acosh]]."
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)
        ax (m/abs x)
        ay (m/abs y)]
    (cond
      (or (m/nan? x) (m/nan? y)) (cond
                                   (and (m/zero? x) (m/nan? y)) (Vec2. m/HALF_PI ##NaN) ;; C99: cacos(+-0 + iNaN)
                                   (m/inf? y) (Vec2. ##NaN (m/- y))
                                   (m/inf? x) (Vec2. ##NaN (m/copy-sign ##Inf y))
                                   :else (Vec2. ##NaN ##NaN))
      (or (m/> ax recip-eps) (m/> ay recip-eps))
      (let [ry (m/+ (log-abs x y) m/LN2)]
        (Vec2. (m/abs (m/atan2 y x))
               (if (signbit? y) ry (m/- ry))))
      (and (m/one? x) (m/zero? y)) (Vec2. 0.0 (m/- y))
      (and (m/< ax sqrt-6-eps-4) (m/< ay sqrt-6-eps-4)) (Vec2. (m/- m/HALF_PI x) (m/- y))
      :else (let [^Vec4 w (hard-work ay ax)
                  ry (.x w)
                  usable? (m/> (.y w) 0.5)
                  neg-x? (signbit? x)
                  rx (if usable?
                       (if neg-x? (m/acos (m/- (.z w))) (m/acos (.z w)))
                       (m/atan2 (.z w) (if neg-x? (m/- (.w w)) (.w w))))]
              (Vec2. rx (if (signbit? y) ry (m/- ry)))))))

(defn acosh
  "Inverse hyperbolic cosine, the principal value, `acosh(z) = +-i acos(z)`.

  The real part is not negative, the imaginary part is in `[-pi, pi]` and has the sign of `b`. The branch cut is the real axis for `a<1`; the sign of a zero imaginary part selects the side.

  See also [[cosh]], [[acos]], [[asinh]]."
  ^Vec2 [^Vec2 z]
  (let [w (acos z)
        rx (.x w)
        ry (.y w)]
    (cond
      (and (m/nan? rx) (m/nan? ry)) (Vec2. ry rx)
      (m/nan? rx) (Vec2. (m/abs ry) rx)
      (m/nan? ry) (Vec2. ry ry)
      :else (Vec2. (m/abs ry) (m/copy-sign rx (.y z))))))

(defn- sum-squares
  ^double [^double x ^double y]
  (if (m/< y sqrt-min)
    (m/* x x)
    (m/+ (m/* x x) (m/* y y))))

(defn- real-part-reciprocal
  "Re(1/(x+iy)) = x/(x^2+y^2) without overflow."
  ^double [^double x ^double y]
  (cond
    (m/inf? x) (m// x)
    (m/inf? y) (m/copy-sign 0.0 x)
    :else
    (let [mx (m/max (m/abs x) (m/abs y))
          xs (m// x mx)
          ys (m// y mx)]
      (m// (m// xs (m/+ (m/* xs xs) (m/* ys ys))) mx))))

(defn atanh
  "Inverse hyperbolic tangent, the principal value.

  The imaginary part is in `[-pi/2, pi/2]`. Branch cuts are the real axis for `|a|>1`; the sign of a zero imaginary part selects the side. `atanh(+-1)` is `+-Infinity` with a zero imaginary part. Accurate near zero, near the branch points `+-1` and for large arguments.

  Special values (C99): an infinite imaginary part gives `+-0` with the sign of the real part and `+-pi/2` (finite real part); `atanh(+-1)` is `+-Infinity`.

  See also [[tanh]], [[atan]], [[asinh]]."
  ^Vec2 [^Vec2 z]
  (let [x (.x z)
        y (.y z)
        ax (m/abs x)
        ay (m/abs y)]
    (cond
      (and (m/zero? y) (m/<= ax 1.0)) (Vec2. (if (m/one? ax)
                                               (m/copy-sign ##Inf x)
                                               (m/* 0.5 (m/log1p (m// (m/* 2.0 x) (m/- 1.0 x)))))
                                             y)
      (m/zero? x) (Vec2. x (m/atan y))
      (or (m/nan? x) (m/nan? y)) (if (or (m/inf? x) (m/inf? y))
                                   (Vec2. (m/copy-sign 0.0 x) (m/copy-sign m/HALF_PI y))
                                   (Vec2. ##NaN ##NaN))
      (or (m/> ax recip-eps) (m/> ay recip-eps)) (Vec2. (real-part-reciprocal x y)
                                                        (m/copy-sign m/HALF_PI y))
      (and (m/< ax sqrt-3-eps-2) (m/< ay sqrt-3-eps-2)) z
      :else (let [rx (if (and (m/one? ax) (m/< ay eps))
                       (m/* 0.5 (m/- m/LN2 (m/log ay)))
                       (m/* 0.25 (m/log1p (m// (m/* 4.0 ax) (sum-squares (m/dec ax) ay)))))
                  ry (cond
                       (m/one? ax) (m/* 0.5 (m/atan2 2.0 (m/- ay)))
                       (m/< ay eps) (m/* 0.5 (m/atan2 (m/* 2.0 ay) (m/* (m/- 1.0 ax) (m/+ 1.0 ax))))
                       :else (m/* 0.5 (m/atan2 (m/* 2.0 ay) (m/- (m/* (m/- 1.0 ax) (m/+ 1.0 ax)) (m/* ay ay)))))]
              (Vec2. (m/copy-sign rx x) (m/copy-sign ry y))))))

(defn atan
  "Inverse tangent, the principal value, `atan(z) = -i atanh(iz)`.

  The real part is in `[-pi/2, pi/2]`. Branch cuts are the imaginary axis for `|b|>1`; the sign of a zero real part selects the side. `atan(+-i)` has an infinite imaginary part. Accurate near zero, near `+-i` and for large arguments.

  Special values (C99): an infinite real part gives `+-pi/2` with a zero imaginary part of the sign of the imaginary part.

  See also [[tan]], [[atanh]], [[acot]]."
  ^Vec2 [^Vec2 z]
  (let [w (atanh (Vec2. (.y z) (.x z)))]
    (Vec2. (.y w) (.x w))))

(defn- inverse-argument
  "1/z for the inverse functions of the reciprocal. At zero it is the limit along the real axis, an infinite real part with the sign of the real part of `z` and a zero imaginary part (the side of the cut is selected by the sign of the zero imaginary part of `z`)."
  ^Vec2 [^Vec2 z]
  (if (zero? z)
    (Vec2. (m/copy-sign ##Inf (.x z)) (m/- (m/copy-sign 0.0 (.y z))))
    (reciprocal z)))

(defn asec
  "Inverse secant, the principal value, `asec(z) = acos(1/z)`.

  Branch cuts follow [[acos]]. Computed from `1/z`, which is rounded: close to the branch points `1` the relative error grows (about `eps/distance`; 8e-10 at the distance 1.4e-8). Away from them it is a few ulps. At zero the value is the limit along the real axis: `asec(0)` is `0+Infinity*i` (`pi` for `-0`), the sign of the infinite part follows the sign of the zero imaginary part.

  See also [[sec]], [[acos]]."
  ^Vec2 [^Vec2 z]
  (acos (inverse-argument z)))

(defn acsc
  "Inverse cosecant, the principal value, `acsc(z) = asin(1/z)`.

  Branch cuts follow [[asin]]. Computed from `1/z`, which is rounded: close to the branch points `+-1` the relative error grows (about `eps/sqrt(distance)`; 3e-11 at the distance 3e-12). Away from them it is a few ulps. At zero the value is the limit along the real axis: `acsc(0)` is `pi/2-Infinity*i`.

  See also [[csc]], [[asin]]."
  ^Vec2 [^Vec2 z]
  (asin (inverse-argument z)))

(defn acot
  "Inverse cotangent, the principal value, `acot(z) = atan(1/z)`.

  The real part is in `[-pi/2, pi/2]` and it is odd, `acot(-x) = -acot(x)` (the other common definition, `pi/2 - atan(z)`, differs on the negative real axis and on the cut). Branch cut: the segment of the imaginary axis between `-i` and `i`. Computed from `1/z`, which is rounded: close to the branch points `+-i` the relative error grows (about `eps/distance`; 3e-7 at the distance 2e-11). Away from them it is a few ulps.

  See also [[cot]], [[atan]]."
  ^Vec2 [^Vec2 z]
  (atan (inverse-argument z)))

(defn asech
  "Inverse hyperbolic secant, the principal value, `asech(z) = acosh(1/z)`.

  The real part is not negative. Branch cuts are the negative real axis and the real axis for `a>1`. Computed from `1/z`, which is rounded: close to the branch points `1` the relative error grows (about `eps/distance`; 8e-10 at the distance 1.4e-8). Away from them it is a few ulps. `asech(0)` is `Infinity-0i` (`Infinity-pi*i` for `-0`; the sign of the imaginary part is the opposite of the sign of the zero imaginary part of `z`).

  See also [[sech]], [[acosh]]."
  ^Vec2 [^Vec2 z]
  (acosh (inverse-argument z)))

(defn acsch
  "Inverse hyperbolic cosecant, the principal value, `acsch(z) = asinh(1/z)`.

  Branch cut: the segment of the imaginary axis between `-i` and `i`. Computed from `1/z`, which is rounded: close to the branch points `+-i` the relative error grows (about `eps/sqrt(distance)`; 1e-11 at the distance 2e-11). Away from them it is a few ulps. `acsch(0)` is `Infinity-0i` (`-Infinity-0i` for `-0`; the sign of the zero is the opposite of the sign of the zero imaginary part of `z`).

  See also [[csch]], [[asinh]]."
  ^Vec2 [^Vec2 z]
  (asinh (inverse-argument z)))

(defn acoth
  "Inverse hyperbolic cotangent, the principal value, `acoth(z) = atanh(1/z)`.

  Branch cut: the segment of the real axis between `-1` and `1`. `acoth(+-1)` is `+-Infinity`. Computed from `1/z`, which is rounded: close to the branch points `+-1` the relative error grows (about `eps/distance`; 1e-6 at the distance 3e-12). Away from them it is a few ulps.

  See also [[coth]], [[atanh]]."
  ^Vec2 [^Vec2 z]
  (atanh (inverse-argument z)))

