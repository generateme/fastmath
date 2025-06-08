(ns complex-quaternions
  (:require [fastmath.complex :as cplx]
            [fastmath.quaternion :as quat]
            [fastmath.vector :as v]
            [fastmath.dev.codox :as codox]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.clay :as utls]
            [scicloj.kindly.v4.kind :as kind]
            [fastmath.core :as m]
            [fastmath.matrix :as mat]))

;; # Complex numbers and quaternions {.unnumbered}

;; ## Complex numbers

;; ::: {.callout-tip title="Defined functions"}
;; * `complex`, `ensure-complex`
;; :::


;; Complex numbers extend the concept of real numbers by including an imaginary unit $i$, where $i^2 = -1$.

;; In `fastmath`, complex numbers are represented as 2-dimensional vectors using the `fastmath.vector/Vec2` type, where the x-component is the real part and the y-component is the imaginary part.

;; You can create complex numbers using the dedicated `complex` function, the general 2D vector constructor `fastmath.vector/vec2`, or by converting sequences/arrays using functions like `fastmath.vector/seq->vec2`. `ensure-complex` accepts single argument which can be a real or complex number. 

(utls/examples-note
  (cplx/complex)
  (cplx/complex 1.0)
  (cplx/complex 1.0 -2.0)
  (v/vec2 1.0 -2.0)
  (v/seq->vec2 [1.0 -2.0])
  (cplx/ensure-complex 2.0)
  (cplx/ensure-complex (cplx/complex 1.0 -2.0)))

;; The implementation handles floating-point subtleties, including special values (`##Inf`, `##NaN`) and signed zero (`+0.0`, `-0.0`), for robust calculations across a wide range of inputs.

;; This section covers constants, predicates, and a comprehensive set of operations including basic arithmetic, exponentiation, logarithms, trigonometric functions, and their inverses.

;; Plots are based on domain coloring in HSB space, where
;;
;; * hue - represents argument
;; * saturation/brightness - represents fractional part of logarithm of magnitude

(gg/->image (gg/complex-function identity {:title "f(z)=z"}))

;; We will use these numbers to illustrate operations:

(def z1 (cplx/complex 0.5 -1.5))
(def z2 (cplx/complex -0.5 2.5))
(def z3 (cplx/complex 1.5 0.5))

(utls/examples-note z1 z2 z3)

;; ### Constants

;; List of defined complex constants.

(kind/table
 {:column-names ["constant" "z" "value"]
  :row-vectors [['I (kind/md "$0+i$") (str cplx/I)]
                ['-I (kind/md "$0-i$") (str cplx/-I)]
                ['I- (kind/md "$0-i$") (str cplx/I-)]
                ['ZERO (kind/md "$0+0i$") (str cplx/ZERO)]
                ['ONE (kind/md "$1+0i$") (str cplx/ONE)]
                ['TWO (kind/md "$2+0i$") (str cplx/TWO)]
                ['PI (kind/md "$\\pi+0i$") (str cplx/PI)]]})

;; ### Basic operations

;; This section provides a collection of functions for performing fundamental operations on complex numbers.
;; These operations include extracting components, calculating magnitude and argument, finding conjugates,
;; performing arithmetic (addition, subtraction, multiplication, division), negation, reciprocals, and square roots.

;; ::: {.callout-tip title="Defined functions"}
;; * `re`, `im`
;; * `delta-eq`, `csgn`, `flip`
;; * `abs`, `norm`, `normalize`, `conjugate`
;; * `arg`
;; * `add`, `adds`, `sub`, `neg`
;; * `scale`, `mult`, `mult-I`, `mult-I-`, `muladd`, `sq`
;; * `div`, `reciprocal`
;; * `sqrt`, `sqrt1z`
;; :::
;;
;; *   **Components**
;;
;;     *   `re`: Returns the real part of a complex number $z$. For $z = a + bi$, returns $a$.
;;     *   `im`: Returns the imaginary part of a complex number $z$. For $z = a + bi$, returns $b$.

(utls/examples-note
  (cplx/re z1)
  (cplx/im z1))

;; *   **Properties and Comparisons**
;;
;;     *   `abs`: Calculates the magnitude or modulus of a complex number $|z|$. For $z = a + bi$: $$|z| = \sqrt{a^2 + b^2}$$ This corresponds to the length of the vector `(Vec2. a b)`.
;;     *   `norm`: Calculates the squared magnitude $|z|^2$. For $z = a + bi$: $$|z|^2 = a^2 + b^2$$ This avoids the square root calculation and is often used when comparing magnitudes.
;;     *   `arg`: Computes the argument (phase angle) of a complex number $z$. This is the angle $\phi$ such that $z = |z| (\cos \phi + i \sin \phi)$, and is typically in the range $(-\pi, \pi]$. Calculated as: $$\operatorname{arg}(z)=\operatorname{atan2}(\operatorname{im}(z), \operatorname{re}(z))$$
;;     *   `normalize`: Sets magnitude of a complex number $z$ to 1.0: $$\frac{z}{|z|}$$
;;     *   `conjugate`: Returns the complex conjugate of $z$, denoted $\bar{z}$. For $z = a + bi$: $$\bar{z} = a - bi$$ Geometrically, this is a reflection across the real axis.
;;     *   `flip`: Swaps the real and imaginary parts of a complex number. For $z = a + bi$, returns $b + ai$.
;;     *   `delta-eq`: Checks if two complex numbers are approximately equal within a given tolerance. This is essential for floating-point comparisons and equivalent to checking if $|z_1 - z_2| < \text{tolerance}$.
;;     *   `csgn`: Implements the complex signum function. Returns 0 for the zero complex number. For non-zero $z$, it returns the sign of the real part $\operatorname{sgn}(\operatorname{re}(z))$.

(utls/examples-note
  z1
  (cplx/abs z1)
  (cplx/norm z1)
  (cplx/arg z1)
  (cplx/normalize z1)
  (cplx/abs (cplx/normalize z1))
  (cplx/conjugate z1)
  (cplx/flip z1)
  z2
  (cplx/abs z2)
  (cplx/norm z2)
  (cplx/arg z2)
  (cplx/normalize z2)
  (cplx/abs (cplx/normalize z2))
  (cplx/conjugate z2)
  (cplx/flip z2)
  (cplx/csgn z1)
  (cplx/csgn z2)
  (cplx/csgn cplx/ZERO)
  (cplx/delta-eq (cplx/complex 0.001) cplx/ZERO)
  (cplx/delta-eq (cplx/complex 0.001) cplx/ZERO 0.01))

;; *   **Arithmetic Operations**
;;
;;     *   `add`: Computes the sum of two complex numbers, $z_1 + z_2$. If $z_1 = a + bi$ and $z_2 = c + di$: $$z_1 + z_2 = (a+c) + (b+d)i$$
;;     *   `adds`: Adds a real scalar $s$ to a complex number $z$. For $z = a + bi$: $$z + s = (a+s) + bi$$
;;     *   `sub`: Computes the difference of two complex numbers, $z_1 - z_2$. If $z_1 = a + bi$ and $z_2 = c + di$: $$z_1 - z_2 = (a-c) + (b-d)i$$
;;     *   `neg`: Returns the negation of a complex number, $-z$. For $z = a + bi$: $$-z = -a - bi$$
;;     *   `scale`: Multiplies a complex number $z$ by a real scalar $s$. For $z = a + bi$: $$z \cdot s = (as) + (bs)i$$
;;     *   `mult`: Computes the product of two complex numbers, $z_1 \cdot z_2$. If $z_1 = a + bi$ and $z_2 = c + di$: $$z_1 \cdot z_2 = (ac - bd) + (ad + bc)i$$
;;     *   `mult-I`: Multiplies a complex number $z$ by the imaginary unit $i$. For $z = a + bi$: $$z \cdot i = -b + ai$$
;;     *   `mult-I-`: Multiplies a complex number $z$ by $-i$. For $z = a + bi$: $$z \cdot (-i) = b - ai$$
;;     *   `muladd`: Computes the fused multiply-add operation $$(x \cdot y) + z$$ for complex numbers $x, y, z$.
;;     *   `sq`: Computes the square of a complex number, $z^2$. For $z = a + bi$: $$z^2 = (a^2 - b^2) + 2abi$$
;;     *   `div`: Computes the division of two complex numbers, $z_1 / z_2$. If $z_1 = a + bi$ and $z_2 = c + di$: $$z_1 / z_2 = \frac{z_1 \cdot \bar{z_2}}{|z_2|^2} = \frac{(ac+bd) + (bc-ad)i}{c^2+d^2}$$
;;     *   `reciprocal`: Computes the reciprocal of a complex number, $1/z$. For $z = a + bi$: $$1/z = \frac{\bar{z}}{|z|^2} = \frac{a - bi}{a^2+b^2}$$

(utls/examples-note
  (cplx/add z1 z2)
  (cplx/adds z1 -3.5)
  (cplx/sub z1 z2)
  (cplx/neg z1)
  (cplx/scale z1 5)
  (cplx/mult z1 z2)
  (cplx/mult-I z1)
  (cplx/mult-I- z1)
  (cplx/muladd z1 z2 z3)
  (cplx/sq z1)
  (cplx/div z1 z2)
  (cplx/reciprocal z1))

(kind/table
 [[(gg/->image (gg/complex-function cplx/sq {:title "f(z)=z^2"}))
   (gg/->image (gg/complex-function cplx/reciprocal {:title "f(z)=1/z"}))]])

;; *   **Roots**
;;
;;     *   `sqrt`: Computes the principal square root of a complex number: $$\sqrt{z} = \sqrt{\frac{|a| + x}{2}} + i \cdot \operatorname{sgn}(b) \cdot \sqrt{\frac{|z| - a}{2}}$$
;;     *   `sqrt1z`: Computes $\sqrt{1 - z^2}$ for a complex number $z$.
;;


(utls/examples-note
  (cplx/sqrt z1)
  (cplx/sqrt z2)
  (cplx/sqrt z3)
  (cplx/sqrt1z z1)
  (cplx/sqrt1z z2)
  (cplx/sqrt1z z3))

(kind/table
 [[(gg/->image (gg/complex-function cplx/sqrt {:title "f(z)=sqrt(z)"}))
   (gg/->image (gg/complex-function cplx/sqrt1z {:title "f(z)=sqrt(1-z^2)"}))]])


;; ### Predicates

;; Predicates are functions that return a boolean value, indicating whether a complex number satisfies a specific condition. They are useful for checking the nature or state of a complex number, such as whether it is real, imaginary, zero, infinite, or invalid.

;; ::: {.callout-tip title="Defined functions"}
;; * `real?`, `imaginary?`
;; * `zero?`
;; * `inf?`,`nan?`
;; * `invalid?`, `valid?`
;; :::
;;
;; *   `real?`: Checks if a complex number `z` has a zero imaginary part. $z = a + bi$ is real if $b = \operatorname{im}(z) = 0$.
;; *   `imaginary?`: Checks if a complex number `z` has a zero real part. $z = a + bi$ is pure imaginary if $a = \operatorname{re}(z) = 0$.
;; *   `zero?`: Checks if `z` is the complex number $0+0i$, which is true if both real and imaginary parts are zero.
;; *   `inf?`: Checks if either the real or imaginary part of `z` is positive or negative infinity (`##Inf` or `##-Inf`).
;; *   `nan?`: Checks if either the real or imaginary part of `z` is Not a Number (`##NaN`).
;; *   `invalid?`: Checks if `z` is either `inf?` or `nan?`. Represents values that are not finite numbers.
;; *   `valid?`: Checks if `z` is a finite, non-NaN number. This is the opposite of `invalid?`.

(utls/examples-note
  (cplx/real? cplx/ONE)
  (cplx/real? cplx/I)
  (cplx/imaginary? cplx/ONE)
  (cplx/imaginary? cplx/I)
  (cplx/zero? cplx/ZERO)
  (cplx/zero? cplx/ONE)
  (cplx/inf? (cplx/complex ##Inf 1))
  (cplx/nan? (cplx/complex 1 ##NaN))
  (cplx/invalid? (cplx/complex ##Inf 1))
  (cplx/invalid? (cplx/complex 1.0 2.0))
  (cplx/valid? (cplx/complex 1.0 2.0))
  (cplx/valid? (cplx/complex ##Inf 1)))


;; ### Power and logarithms

;; Complex numbers extend the real-valued exponential and logarithm functions. Complex exponentiation allows raising a complex base to a complex power, while the complex logarithm is the inverse operation of exponentiation. Unlike real logarithms, the complex logarithm is multi-valued, but the functions provided here compute the principal value.

;; ::: {.callout-tip title="Defined functions"}
;; * `pow`, `exp`
;; * `log`, `logb`
;; :::

;; *   `exp`: Computes the complex exponential $e^z$. For $z = x + iy$: $$e^z = e^{x+iy} = e^x (\cos y + i \sin y)$$
;; *   `log`: Computes the principal value of the complex natural logarithm $\log z$. This is given by: $$\log z = \ln|z| + i \arg(z)$$ where $\ln|z|$ is the real natural logarithm of the magnitude $|z|$, and $\arg(z)$ is the principal argument of $z$.
;; *   `logb`: Computes the logarithm of $z$ with a complex base $b$: $$\log_b z = \frac{\log z}{\log b}$$
;; *   `pow`: Computes the complex power $z_1^{z_2}$. This is defined as $$z_1^{z_2} = e^{z_2 \log z_1}$$ where $\log z_1$ is the principal value of the logarithm. Handles various edge cases, including $0^0$.

(utls/examples-note
  (cplx/exp z1)
  (cplx/log z1)
  (cplx/logb z1 z2)
  (cplx/pow z1 z2)
  (cplx/pow z3 z1)
  (cplx/pow cplx/ZERO cplx/ZERO)
  (cplx/pow cplx/ZERO (cplx/complex 1 1))
  (cplx/pow (cplx/complex 1 1) cplx/ZERO))

(kind/table
 [[(gg/->image (gg/complex-function cplx/exp {:title "f(z)=exp(z)"}))
   (gg/->image (gg/complex-function cplx/log {:title "f(z)=log(z)"}))]
  [(gg/->image (gg/complex-function #(cplx/pow z1 %) {:title "f(z)=z1^z"}))
   (gg/->image (gg/complex-function #(cplx/pow % z1) {:title "f(z)=z^z1"}))]
  [(gg/->image (gg/complex-function #(cplx/logb z1 %) {:title "f(z)=log_z(z1)"}))
   (gg/->image (gg/complex-function #(cplx/logb % z2) {:title "f(z)=log_z2(z)"}))]])

;; ### Trigonometric and Hyperbolic

;; Complex trigonometric and hyperbolic functions extend their real-valued counterparts to the complex plane. They are defined in terms of the complex exponential function:
;;
;; $$ \sin z = \frac{e^{iz} - e^{-iz}}{2i} \quad \cos z = \frac{e^{iz} + e^{-iz}}{2} $$
;; $$ \sinh z = \frac{e^{z} - e^{-z}}{2} \quad \cosh z = \frac{e^{z} + e^{-z}}{2} $$
;;
;; This leads to relationships between complex trigonometric and hyperbolic functions, such as $\sin(z) = -i \sinh(iz)$ and $\cos(z) = \cosh(iz)$. These functions are implemented using formulas that provide numerical stability, especially for large arguments.

;; ::: {.callout-tip title="Defined functions"}
;; * `sin`, `cos`, `tan`
;; * `sec`, `csc`, `cot`
;; * `sinh`, `cosh`, `tanh`
;; * `sech`, `csch`, `coth`
;; * `asin`, `acos`, `atan`
;; * `asec`, `acsc`, `acot`
;; * `asinh`, `acosh`, `atanh`
;; * `asech`, `acsch`, `acoth`
;; :::

;; *   **Standard Trigonometric Functions**: `sin`, `cos`, `tan`, `sec`, `csc`, `cot`
;;     *   These are the direct extensions of real trigonometric functions. For $z = x + iy$:
;;         $$ \sin z = \sin x \cosh y + i \cos x \sinh y $$
;;         $$ \cos z = \cos x \cosh y - i \sin x \sinh y $$
;;         $$ \tan z = \frac{\sin x \cos x + i \sinh y \cosh y}{\cos^2 x + \sinh^2 y} $$
;;     *   `sec(z)`, `csc(z)`, and `cot(z)` are computed as the reciprocals of `cos(z)`, `sin(z)`, and `tan(z)`, respectively.
;;
;; *   **Hyperbolic Trigonometric Functions**: `sinh`, `cosh`, `tanh`, `sech`, `csch`, `coth`
;;     *   These are also extensions of real hyperbolic functions. For $z = x + iy$:
;;         $$ \sinh z = \sinh x \cos y + i \cosh x \sin y $$
;;         $$ \cosh z = \cosh x \cos y + i \sinh x \sin y $$
;;         $$ \tanh z = \frac{\sinh x \cosh x + i \sin y \cos y}{\sinh^2 x + \cos^2 y} $$
;;     *   `sech(z)`, `csch(z)`, and `coth(z)` are computed as the reciprocals of `cosh(z)`, `sinh(z)`, and `tanh(z)`, respectively.
;;
;; *   **Inverse Trigonometric Functions**: `asin`, `acos`, `atan`, `asec`, `acsc`, `acot`
;;     *   These are the multi-valued inverse functions. The implementations return the principal values, which are related to complex logarithms. For example:
;;         $$ \arcsin z = -i \log(iz + \sqrt{1-z^2}) $$
;;         $$ \arccos z = -i \log(z + i\sqrt{1-z^2}) $$
;;         $$ \arctan z = \frac{i}{2} \log\left(\frac{1-iz}{1+iz}\right) $$
;;         $$ \operatorname{arcsec} z = \arccos(1/z) $$
;;         $$ \operatorname{arccsc} z = \arcsin(1/z) $$
;;         $$ \operatorname{arccot} z = \operatorname{atan}(1/z) $$
;;
;;
;; *   **Inverse Hyperbolic Trigonometric Functions**: `asinh`, `acosh`, `atanh`, `asech`, `acsch`, `acoth`
;;     *   These are the multi-valued inverse hyperbolic functions. The implementations return the principal values, which are also related to complex logarithms. For example:
;;         $$ \operatorname{arsinh} z = \log(z + \sqrt{z^2+1}) $$
;;         $$ \operatorname{arcosh} z = \log(z + \sqrt{z^2-1}) $$
;;         $$ \operatorname{artanh} z = \frac{1}{2} \log\left(\frac{1+z}{1-z}\right) $$
;;         $$ \operatorname{arsech} z = \operatorname{arcosh}(1/z) $$
;;         $$ \operatorname{arcsch} z = \operatorname{arsinh}(1/z) $$
;;         $$ \operatorname{arcoth} z = \operatorname{artanh}(1/z) $$

;;
;; These functions are continuous and analytic everywhere in the complex plane, except at certain points (poles) where the denominator is zero for tangent, cotangent, secant, cosecant, and at branch cuts for the inverse functions.
;;

(utls/examples-note
  (cplx/sin z1)
  (cplx/cos z1)
  (cplx/tan z1)
  (cplx/sec z1)
  (cplx/csc z1)
  (cplx/cot z1)
  (cplx/asin z1)
  (cplx/acos z1)
  (cplx/atan z1)
  (cplx/asec z1)
  (cplx/acsc z1)
  (cplx/acot z1)
  (cplx/sinh z1)
  (cplx/cosh z1)
  (cplx/tanh z1)
  (cplx/sech z1)
  (cplx/csch z1)
  (cplx/coth z1)
  (cplx/asinh z1)
  (cplx/acosh z1)
  (cplx/atanh z1)
  (cplx/asech z1)
  (cplx/acsch z1)
  (cplx/acoth z1))

(kind/table
 [["sin" "cos" "tan"]
  [(gg/->image (gg/complex-function cplx/sin {:title "f(z)=sin(z)"}))
   (gg/->image (gg/complex-function cplx/cos {:title "f(z)=cos(z)"}))
   (gg/->image (gg/complex-function cplx/tan {:title "f(z)=tan(z)"}))]
  ["sec" "csc" "cot"]
  [(gg/->image (gg/complex-function cplx/sec {:title "f(z)=sec(z)"}))
   (gg/->image (gg/complex-function cplx/csc {:title "f(z)=csc(z)"}))
   (gg/->image (gg/complex-function cplx/cot {:title "f(z)=cot(z)"}))]
  ["asin" "acos" "atan"]
  [(gg/->image (gg/complex-function cplx/asin {:title "f(z)=asin(z)"}))
   (gg/->image (gg/complex-function cplx/acos {:title "f(z)=acos(z)"}))
   (gg/->image (gg/complex-function cplx/atan {:title "f(z)=atan(z)"}))]
  ["asec" "acsc" "acot"]
  [(gg/->image (gg/complex-function cplx/asec {:title "f(z)=asec(z)"}))
   (gg/->image (gg/complex-function cplx/acsc {:title "f(z)=acsc(z)"}))
   (gg/->image (gg/complex-function cplx/acot {:title "f(z)=acot(z)"}))]
  ["sinh" "cosh" "tanh"]
  [(gg/->image (gg/complex-function cplx/sinh {:title "f(z)=sinh(z)"}))
   (gg/->image (gg/complex-function cplx/cosh {:title "f(z)=cosh(z)"}))
   (gg/->image (gg/complex-function cplx/tanh {:title "f(z)=tanh(z)"}))]
  ["sech" "csch" "coth"]
  [(gg/->image (gg/complex-function cplx/sech {:title "f(z)=sech(z)"}))
   (gg/->image (gg/complex-function cplx/csch {:title "f(z)=csch(z)"}))
   (gg/->image (gg/complex-function cplx/coth {:title "f(z)=coth(z)"}))]
  ["asinh" "acosh" "atanh"]
  [(gg/->image (gg/complex-function cplx/asinh {:title "f(z)=asinh(z)"}))
   (gg/->image (gg/complex-function cplx/acosh {:title "f(z)=acosh(z)"}))
   (gg/->image (gg/complex-function cplx/atanh {:title "f(z)=atanh(z)"}))]
  ["asech" "acsch" "acoth"]
  [(gg/->image (gg/complex-function cplx/asech {:title "f(z)=asech(z)"}))
   (gg/->image (gg/complex-function cplx/acsch {:title "f(z)=acsch(z)"}))
   (gg/->image (gg/complex-function cplx/acoth {:title "f(z)=acoth(z)"}))]])

;; ## Quaternions

;; ::: {.callout-tip title="Defined functions"}
;; * `quaternion`, `complex->quaternion`, `ensure-quaternion`
;; :::

;; Quaternions extend complex numbers with three imaginary units: $i$, $j$, and $k$, satisfying the relations $i^2 = j^2 = k^2 = ijk = -1$. They are particularly useful for representing 3D rotations concisely and avoiding gimbal lock.
;;
;; In `fastmath`, quaternions are represented as 4-dimensional vectors using the `fastmath.vector/Vec4` type. A quaternion $q = a + bi + cj + dk$ corresponds to the vector `(Vec4. a b c d)`, where $a$ is the scalar (real) part, and $(b, c, d)$ is the vector (imaginary) part.
;;
;; You can create quaternions using the dedicated `quaternion` function, or by converting numbers or complex numbers using `ensure-quaternion`.
;;
;; The implementation handles floating-point subtleties, including special values (`##Inf`, `##NaN`), for robust calculations across a wide range of inputs.
;;
;; This section covers constants, creation, accessors, predicates, basic operations, and advanced functions including power, logarithms, trigonometric functions, and functions specific to 3D rotations.
;;

(utls/examples-note
  (quat/quaternion 1.0 -2.0 3.0 -4.0)
  (quat/quaternion 5.0 [1.0 2.0 3.0])
  (quat/quaternion 7.0)
  (quat/complex->quaternion (cplx/complex 1.0 -2.0))
  (quat/ensure-quaternion 5.0)
  (quat/ensure-quaternion (cplx/complex 1.0 -2.0)))

;;
;; We will use these quaternions to illustrate operations:

(def q1 (quat/quaternion 0.5 -1.5 2.0 -0.5))
(def q2 (quat/quaternion -0.5 2.5 -1.0 1.5))
(def q3 (quat/quaternion 1.5 0.0 0.0 0.0)) 
(def q4 (quat/quaternion 0.0 1.0 -2.0 3.0))

(utls/examples-note q1 q2 q3 q4)

;; All of the trigonometric, hyperbolic, `exp`, `log` and `sqrt` functions are based on the following general formula. Let $r_q$ is the real part and $v_q$ is the vector (imaginary) part of quaternion $q$. Let $z_{in}$ is defined as follows:

;; $$z_{in} = r_q + |v_q|i$$

;; Then, if $f_z$ is given operation on complex number, corresponding quaternion version $f_q$ is:

;; $$f_q(q) = \Re(f_z(z_{in})) + \frac{\Im(f_z(z_{in}))}{|v_q|}v_q$$

;; ### Constants
;;
;; Predefined quaternion constants for common values.

(kind/table
 {:column-names ["constant" "q" "value"]
  :row-vectors [['ZERO (kind/md "$0+0i+0j+0k$") (str quat/ZERO)]
                ['ONE (kind/md "$1+0i+0j+0k$") (str quat/ONE)]
                ['I (kind/md "$0+i+0j+0k$") (str quat/I)]
                ['J (kind/md "$0+0i+j+0k$") (str quat/J)]
                ['K (kind/md "$0+0i+0j+k$") (str quat/K)]
                ['-I (kind/md "$0-i+0j+0k$") (str quat/-I)]
                ['-J (kind/md "$0+0i-j+0k$") (str quat/-J)]
                ['-K (kind/md "$0+0i+0j-k$") (str quat/-K)]]})

;; ### Basic operations
;;
;; Functions to create quaternions and access their scalar and vector components.
;;

;; ::: {.callout-tip title="Defined functions"}
;; * `scalar`, `re`, `vector`, `im-i`, `im-j`, `im-k`
;; * `delta-eq`, `qsgn`, `arg`, `norm`, `normalize`, `conjugate`, `reciprocal`
;; * `add`, `adds`, `sub`, `scale`, `mult`, `div`, `neg`, `sq`, `sqrt`
;; :::

;; * **Components**

;; *   `scalar`, `re`: Returns the scalar (real) part $a$ of a quaternion $a+bi+cj+dk$.
;; *   `im-i`, `im-j`, `im-k`: Return the coefficients of the $i$, $j$, and $k$ components, respectively.
;; *   `vector`: Returns the vector (imaginary) part $(b, c, d)$ as a `Vec3`.

(utls/examples-note
  (quat/scalar q1)
  (quat/re q1)
  (quat/vector q1)
  (quat/im-i q1)
  (quat/im-j q1)
  (quat/im-k q1))

;; * **Properties and Comparisons**

;; *   `arg`: Computes the argument $\theta$ such that $q = |q| (\cos \theta + \mathbf{v} \sin \theta)$, where $\mathbf{v}$ is a unit 3D vector. Calculated as: $$\operatorname{arg}(z)=\operatorname{atan2}(|v|, \operatorname{re}(q))$$
;; *   `norm`: Calculates the magnitude (Euclidean norm) of the quaternion: $$|q| = \sqrt{a^2+b^2+c^2+d^2}$$
;; *   `normalize`: Returns the unit quaternion: $$\hat{q} = q / |q|$$
;; *   `conjugate`: Returns the conjugate: $$\bar{q} = a - bi - cj - dk$$
;; *   `delta-eq`: Checks if two quaternions are approximately equal within a tolerance.
;; *   `qsgn`: Computes the complex signum, returning the sign of the scalar part or 0 for zero.

(utls/examples-note
  q1
  (quat/arg q1)
  (quat/norm q1)
  (quat/normalize q1)
  (quat/norm (quat/normalize q1))
  (quat/conjugate q1)
  (quat/qsgn q1)
  q2
  (quat/arg q2)
  (quat/norm q2)
  (quat/normalize q2)
  (quat/conjugate q2)
  (quat/qsgn q2)
  (quat/qsgn quat/ZERO)
  (quat/delta-eq (quat/quaternion 1.0000001) quat/ONE)
  (quat/delta-eq (quat/quaternion 1.0001) quat/ONE))

;; * **Arithemtic Operations**

;; Basic arithmetic extended to quaternions.
;;
;; *   `add`: Sum of two quaternions: $$q_1 + q_2 = (a_1+a_2) + (b_1+b_2)i + (c_1+c_2)j + (d_1+d_2)k$$
;; *   `adds`: Adds a real scalar $s$ to a quaternion: $q = a+bi+cj+dk$: $$q+s = (a+s)+bi+cj+dk$$
;; *   `sub`: Difference of two quaternions: $q_1 - q_2 = (a_1-a_2) + (b_1-b_2)i + (c_1-c_2)j + (d_1-d_2)k$$
;; *   `scale`: Multiplies a quaternion $q$ by a real scalar $s$: $$q \cdot s = (as) + (bs)i + (cs)j + (ds)k$$
;; *   `mult`: Quaternion multiplication. Note that quaternion multiplication is non-commutative ($q_1 q_2 \neq q_2 q_1$ in general).
;;     If $q_1 = a_1 + \mathbf{v}_1$ and $q_2 = a_2 + \mathbf{v}_2$, where $\mathbf{v}_1, \mathbf{v}_2$ are vector parts:
;;     $$ q_1 q_2 = (a_1 a_2 - \mathbf{v}_1 \cdot \mathbf{v}_2) + (a_1 \mathbf{v}_2 + a_2 \mathbf{v}_1 + \mathbf{v}_1 \times \mathbf{v}_2) $$
;; *   `div`: Quaternion division $$q_1 / q_2 = q_1 \cdot q_2^{-1}$$
;; *   `reciprocal`: Returns the reciprocal: $$q^{-1} = \bar{q} / |q|^2$$
;; *   `neg`: Negation: $$-q = -a - bi - cj - dk$$
;; *   `sq`: Square of a quaternion: $$q^2 = q \cdot q$$
;; *   `sqrt`: Computes the principal square root of a quaternion.

(utls/examples-note
  (quat/add q1 q2)
  (quat/adds q1 5.0)
  (quat/sub q1 q2)
  (quat/scale q1 3.0)
  (quat/mult q1 q2)
  (quat/mult q2 q1)
  (quat/div q1 q2)
  (quat/reciprocal q1)
  (quat/mult q1 (quat/reciprocal q1)) 
  (quat/neg q1)
  (quat/sq q1)
  (quat/sqrt q1))

;; ### Predicates

;; Functions to check the type or state of a quaternion and compute fundamental properties.

;; ::: {.callout-tip title="Defined functions"}
;; * `real?`, `imaginary?`
;; * `zero?`
;; * `inf?`, `nan?`
;; * `invalid?`, `valid?`
;; :::
;;
;; *   `real?`: Checks if the quaternion has a zero vector part ($b=c=d=0$).
;; *   `imaginary?`: Checks if the quaternion has a zero scalar part ($a=0$).
;; *   `zero?`: Checks if all components are zero ($a=b=c=d=0$).
;; *   `inf?`, `nan?`: Checks if any component is infinite or NaN, respectively.
;; *   `invalid?`, `valid?`: Checks if the quaternion is not a finite, non-NaN value.

(utls/examples-note
  (quat/real? q3)
  (quat/real? q1)
  (quat/imaginary? q4)
  (quat/imaginary? q1)
  (quat/zero? quat/ZERO)
  (quat/zero? q1)
  (quat/inf? (quat/quaternion ##Inf 1 2 3))
  (quat/nan? (quat/quaternion 1 ##NaN 2 3))
  (quat/invalid? (quat/quaternion 1 2 3 4))
  (quat/valid? (quat/quaternion 1 2 3 4)))

;; ### Power and Logarithms
;;
;; Extensions of exponential, logarithm, and power functions to quaternions.
;;
;; *   `exp`: Computes the quaternion exponential $e^q$.
;; *   `log`: Computes the principal value of the quaternion natural logarithm $\log q$.
;; *   `logb`: Computes the logarithm of $q$ with a quaternion base $b$: $\log_b q = (\log q) (\log b)^{-1}$.
;; *   `pow`: Computes the quaternion power $q_1^{q_2}$ defined as $e^{(\log q_1) q_2}$.

(utls/examples-note
  (quat/exp q1)
  (quat/log q1)
  (quat/logb q1 q2)
  (quat/pow q1 q2)
  (quat/pow q3 q1))

;; ### Trigonometric and Hyperbolic
;;
;; Extensions of standard and hyperbolic trigonometric functions to quaternions.

;; ::: {.callout-tip title="Defined functions"}
;; * `sin`, `cos`, `tan`
;; * `sec`, `csc`, `cot`
;; * `sinh`, `cosh`, `tanh`
;; * `sech`, `csch`, `coth`
;; * `asin`, `acos`, `atan`
;; * `asec`, `acsc`, `acot`
;; * `asinh`, `acosh`, `atanh`
;; * `asech`, `acsch`, `acoth`
;; :::

(utls/examples-note
  (quat/sin q1)
  (quat/cos q1)
  (quat/tan q1)
  (quat/sec q1)
  (quat/csc q1)
  (quat/cot q1)
  (quat/asin q1)
  (quat/acos q1)
  (quat/atan q1)
  (quat/asec q1)
  (quat/acsc q1)
  (quat/acot q1)
  (quat/sinh q1)
  (quat/cosh q1)
  (quat/tanh q1)
  (quat/sech q1)
  (quat/csch q1)
  (quat/coth q1)
  (quat/asinh q1)
  (quat/acosh q1)
  (quat/atanh q1)
  (quat/asech q1)
  (quat/acsch q1)
  (quat/acoth q1))

;; ### Rotation Functions
;;
;; Functions specifically for using quaternions to represent and manipulate 3D rotations.
;;

;; ::: {.callout-tip title="Defined functions"}
;; * `rotation-quaternion`, `rotate`
;; * `slerp`
;; * `to-euler`, `from-euler`
;; * `to-angles`, `from-angles`
;; * `to-rotation-matrix`, `from-rotation-matrix`
;; :::

;; *   `rotation-quaternion`: Creates a unit quaternion representing a rotation by a given `angle` around a `u` (Vec3) axis.
;; *   `rotate`: Rotates a 3D vector `in` (Vec3) using a quaternion `rotq`. Can also create the rotation quaternion directly from angle and axis.
;; *   `to-euler`, `from-euler`: Converts between a quaternion and ZYX (body 3-2-1) Euler angles [roll pitch yaw].
;; *   `to-angles`, `from-angles`: Converts between a quaternion and Tait-Bryan angles (z-y'-x'' intrinsic) [x y z].
;; *   `to-rotation-matrix`, `from-rotation-matrix`: Converts between a quaternion and a 3x3 rotation matrix.

(def rot-q (quat/rotation-quaternion m/HALF_PI (v/vec3 1 1 0)))

rot-q

(utls/examples-note
  (quat/rotate (v/vec3 1 0 0) rot-q)
  (quat/rotate (v/vec3 1 0 0) m/HALF_PI (v/vec3 1 1 0))
  (quat/from-euler m/HALF_PI 1 m/HALF_PI)
  (quat/to-euler (quat/from-euler m/HALF_PI 1 m/HALF_PI))
  (quat/from-angles m/HALF_PI 1 m/HALF_PI)
  (quat/to-angles (quat/from-angles m/HALF_PI 1 m/HALF_PI)))

;; For a rotation quaternion `rot-q` conversion between rotation matrix and back looks as follows

(quat/to-rotation-matrix rot-q)

(quat/from-rotation-matrix (mat/rows->mat3x3
                            [0.5 0.5 m/SQRT2_2]
                            [0.5 0.5 (- m/SQRT2_2)]
                            [(- m/SQRT2_2) m/SQRT2_2 0.0]))

;; #### Slerp

;; `slerp` function performs Spherical Linear Interpolation between two unit quaternions `q1` and `q2` based on a parameter `t` in [0, 1]. Useful for animating rotations.

(utls/examples
  (quat/slerp quat/ONE (quat/quaternion 0 0 0 1) 0.0)
  (quat/slerp quat/ONE (quat/quaternion 0 0 0 1) 0.5)
  (quat/slerp quat/ONE (quat/quaternion 0 0 0 1) 1.0))

;; Let's see how quaternion components changes for different `t`.

(def q-start (quat/quaternion -1 -2 -3 -4))
(def q-end (quat/quaternion 1 2 3 4))

(gg/->image (gg/functions [["re" (comp quat/re (partial quat/slerp q-start q-end))]
                           ["im i" (comp quat/im-i (partial quat/slerp q-start q-end))]
                           ["im j" (comp quat/im-j (partial quat/slerp q-start q-end))]
                           ["im k" (comp quat/im-k (partial quat/slerp q-start q-end))]]
                          {:title "SLERP between q-start and q-end"
                           :xlab "t"
                           :ylab "component value"
                           :legend-name "component"}))

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.complex)
(codox/make-public-fns-table-clay 'fastmath.quaternion)
