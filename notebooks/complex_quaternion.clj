^{:nextjournal.clerk/visibility :hide-ns
  :nextjournal.clerk/toc true}
(ns complex-quaternion
  {:clj-kondo/config '{:config-in-call {utils/table {:ignore [:unresolved-symbol]}
                                        utils/table2 {:ignore [:unresolved-symbol]}}}}
  (:require [fastmath.complex :as c]
            [fastmath.quaternion :as q]
            [nextjournal.clerk :as clerk]
            [utils :as u]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat]
            [fastmath.core :as m]))

;; # Complex numbers

#_:clj-kondo/ignore
(require '[fastmath.complex :as c])

;; Complex numbers $z = a + bi$ are represented internally as `fastmath.vector.Vec2` type and can be created with `complex` function.

(c/complex 1 2)

;; Implementation is based on the source code of [CQOST library](https://ece.uwaterloo.ca/~dwharder/C++/CQOST/src/Complex.cpp). Most of the functions distinguish $-0.0$ and $0.0$, it's visible when function is discontinuous at $0.0$.

(c/log (c/complex -2.0 -0.0))
(c/log (c/complex -2.0 0.0))

;; Plots are based on domain-coloring where (in HSB color space):
;;
;; * hue - represents argument
;; * brighness - represents $log_2$ of absolute value

(u/complex-graph identity)

;; ## Constants

;; Basic predefined values

^{::clerk/visibility :hide}
(clerk/table
 {:head ["symbol" "value"]
  :rows [['I (clerk/md "$i$")]
         ['-I (clerk/md "$-i$")]
         ['ONE (clerk/md "$1+0i$")]
         ['TWO (clerk/md "$2+0i$")]
         ['ZERO (clerk/md "$0+0i$")]]})

^{::clerk/visibility :hide}
(clerk/example
 c/I c/-I c/ONE c/TWO c/ZERO)

;; ## Functions

;; ### Creation

;; To create complex number call `complex` or `fastmath.vector/vec2`. When single argument is passed to `complex`, it creates complex representing real number. No args `complex` returns complex $0$ (zero). 

^{::clerk/visibility :hide}
(u/table2
 [[complex "Create complex number"]
  [fastmath.vector/vec2 "Create complex number"]])

^{::clerk/visibility :hide}
(clerk/example
 (c/complex -1 2)
 (c/complex 3)
 (c/complex)
 (v/vec2 -1 2)
 (class (c/complex 1 2)))

;; ### Fields

;; There are two functions `re` and `im` to access real and imaginary part of complex number. You can also use indexes: `0` for real part and `1` for imaginary.

^{::clerk/visibility :hide}
(u/table2
 [[re "Real part"]
  [im "Imaginary part"]])

(def complex-3-4i (c/complex 3 -4))

^{::clerk/visibility :hide}
(clerk/example
  (c/re complex-3-4i)
  (c/im complex-3-4i)
  (complex-3-4i 0)
  (complex-3-4i 1))

;; ### Predicates

^{::clerk/visibility :hide}
(u/table2
 [[real? "Is z a real number? ie. imaginary part = 0.0"]
  [imaginary? "Is z real part zero?"]
  [zero? "Is z is zero?"]
  [inf? "Is any part infinite?"]
  [nan? "Is any part NAN?"]])

^{::clerk/visibility :hide}
(clerk/example
  (c/real? c/ONE)
  (c/real? complex-3-4i)
  (c/imaginary? c/I)
  (c/imaginary? complex-3-4i)
  (c/zero? c/ZERO)
  (c/inf? (c/complex 0.0 ##Inf))
  (c/nan? (c/complex ##NaN 0.0)))

;; ### Basic operations

^{::clerk/visibility :hide}
(u/table2
 [[abs "Absolute value, magnitude"]
  [norm "Gauss norm of the complex number, abs squared"]
  [arg "Argument"]
  [conjugate "Complex conjugate"]
  [neg "Negation"]
  [add "Sum of two complex numbers"]
  [sub "Difference of two complex numbers"]
  [mult "Multiplication of two complex numbers"]
  [mult-I "Multiply by i"]
  [mult-I- "Multiply by -i"]
  [scale "Multiply by a scalar"]
  [div "Division of two complex numbers"]
  [reciprocal "Reciprocal of a complex number"]
  [flip "Exchange real and imaginary parts"]
  [csgn "Return 0 if zero complex number or sign of real part otherwise"]
  [delta-eq "Compare two complex numbers"]])

^{::clerk/visibility :hide}
(clerk/example
  (c/abs complex-3-4i)
  (c/norm complex-3-4i)
  (c/arg complex-3-4i)
  (c/conjugate complex-3-4i)
  (c/neg complex-3-4i)
  (c/add complex-3-4i complex-3-4i)
  (c/sub c/TWO complex-3-4i)
  (c/mult complex-3-4i complex-3-4i)
  (c/mult-I complex-3-4i)
  (c/mult-I- complex-3-4i)
  (c/scale complex-3-4i 2)
  (c/div c/TWO complex-3-4i)
  (c/reciprocal complex-3-4i)
  (c/flip complex-3-4i)
  (c/csgn complex-3-4i)
  (c/csgn c/ZERO))

;; To compare two complex numbers with some accuracy (default `1.0e-6`), call `delta-eq`.

^{::clerk/visibility :hide}
(clerk/example
  (c/delta-eq (c/complex 1.0001 0.0) c/ONE)
  (c/delta-eq (c/complex 1.0001 0.0) c/ONE 1.0e-3))

;; ### Trigonometric

;; #### Basic

^{::clerk/visibility :hide}
(u/table2
 [[sin (clerk/md "$\\sin(z)$")]
  [cos (clerk/md "$\\cos(z)$")]
  [tan (clerk/md "$\\tan(z)$")]
  [sec (clerk/md "$\\sec(z)$")]
  [csc (clerk/md "$\\csc(z)$")]
  [cot (clerk/md "$\\cot(z)$")]
  [asin (clerk/md "$\\arcsin(z)$")]
  [acos (clerk/md "$\\arccos(z)$")]
  [atan (clerk/md "$\\arctan(z)$")]
  [asec (clerk/md "$\\operatorname{arcsec}(z)$")]
  [acsc (clerk/md "$\\operatorname{arccsc}(z)$")]
  [acot (clerk/md "$\\operatorname{arccot}(z)$")]])

^{::clerk/visibility :hide}
(clerk/example
 (c/sin complex-3-4i)
 (c/cos complex-3-4i)
 (c/tan complex-3-4i)
 (c/sec complex-3-4i)
 (c/csc complex-3-4i)
 (c/cot complex-3-4i)
 (c/asin complex-3-4i)
 (c/acos complex-3-4i)
 (c/atan complex-3-4i)
 (c/asec complex-3-4i)
 (c/acsc complex-3-4i)
 (c/acot complex-3-4i))

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["sin" "cos" "tan"]
 [(u/complex-graph c/sin) (u/complex-graph c/cos) (u/complex-graph c/tan)]
 ["sec" "csc" "cot"]
 [(u/complex-graph c/sec) (u/complex-graph c/csc) (u/complex-graph c/cot)]
 ["asin" "acos" "atan"]
 [(u/complex-graph c/asin) (u/complex-graph c/acos) (u/complex-graph c/atan)]
 ["asec" "acsc" "acot"]
 [(u/complex-graph c/asec) (u/complex-graph c/acsc) (u/complex-graph c/acot)]]

;; #### Hyperbolic

^{::clerk/visibility :hide}
(u/table2
 [[sinh (clerk/md "$\\operatorname{sinh}(z)$")]
  [cosh (clerk/md "$\\operatorname{cosh}(z)$")]
  [tanh (clerk/md "$\\operatorname{tanh}(z)$")]
  [sech (clerk/md "$\\operatorname{sech}(z)$")]
  [csch (clerk/md "$\\operatorname{csch}(z)$")]
  [coth (clerk/md "$\\operatorname{coth}(z)$")]
  [asinh (clerk/md "$\\operatorname{arcsinh}(z)$")]
  [acosh (clerk/md "$\\operatorname{arccosh}(z)$")]
  [atanh (clerk/md "$\\operatorname{arctanh}(z)$")]
  [asech (clerk/md "$\\operatorname{arcsec}(z)$")]
  [acsch (clerk/md "$\\operatorname{arccsc}(z)$")]
  [acoth (clerk/md "$\\operatorname{arccot}(z)$")]])

^{::clerk/visibility :hide}
(clerk/example
 (c/sinh complex-3-4i)
 (c/cosh complex-3-4i)
 (c/tanh complex-3-4i)
 (c/sech complex-3-4i)
 (c/csch complex-3-4i)
 (c/coth complex-3-4i)
 (c/asinh complex-3-4i)
 (c/acosh complex-3-4i)
 (c/atanh complex-3-4i)
 (c/asech complex-3-4i)
 (c/acsch complex-3-4i)
 (c/acoth complex-3-4i))

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["sinh" "cosh" "tanh"]
 [(u/complex-graph c/sinh) (u/complex-graph c/cosh) (u/complex-graph c/tanh)]
 ["sech" "csch" "coth"]
 [(u/complex-graph c/sech) (u/complex-graph c/csch) (u/complex-graph c/coth)]
 ["asinh" "acosh" "atanh"]
 [(u/complex-graph c/asinh) (u/complex-graph c/acosh) (u/complex-graph c/atanh)]
 ["asech" "acsch" "acoth"]
 [(u/complex-graph c/asech) (u/complex-graph c/acsch) (u/complex-graph c/acoth)]]

;; ### Power and logarithm

^{::clerk/visibility :hide}
(u/table2
 [[sq (clerk/md "$z^2$")]
  [sqrt (clerk/md "$\\sqrt{z}$")]
  [sqrt1z (clerk/md "$\\sqrt{1-z^2}$")]
  [pow (clerk/md "$z^c$")]
  [exp (clerk/md "$e^z$")]
  [log (clerk/md "$\\operatorname{Log}(z)$, principal value")]
  [logb (clerk/md "$\\log_b(z)$")]])

^{::clerk/visibility :hide}
(clerk/example
 (c/sq complex-3-4i)
 (c/sqrt complex-3-4i)
 (c/sqrt1z complex-3-4i)
 (c/pow complex-3-4i complex-3-4i)
 (c/exp complex-3-4i)
 (c/log complex-3-4i)
 (c/logb complex-3-4i c/I))

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["sq" "sqrt" "sqrt1z"]
 [(u/complex-graph c/sq) (u/complex-graph c/sqrt) (u/complex-graph c/sqrt1z)]
 ["exp" "log" "z^(complex-3-4i)"]
 [(u/complex-graph c/exp) (u/complex-graph c/log) (u/complex-graph #(c/pow % complex-3-4i))]]

;; ## List of symbols

^{::clerk/visibility :hide ::clerk/no-cache true}
(u/make-public-fns-table 'fastmath.complex)

;; # Quaternions

#_:clj-kondo/ignore
(require '[fastmath.quaternion :as q])

;; Quaternions $z = a + bi + cj + dk$ are represented internally as `fastmath.vector.Vec4` type and can be created with `quaternion` function.

(q/quaternion 1 2 3 4)

;; ## Constants

;; Basic predefined values

^{::clerk/visibility :hide}
(clerk/table
 {:head ["symbol" "value"]
  :rows [['I (clerk/md "$0+1i+0j+0k$")]
         ['J (clerk/md "$0+0i+1j+0k$")]
         ['K (clerk/md "$0+0i+0j+1k$")]
         ['-I (clerk/md "$0-1i+0j+0k$")]
         ['-J (clerk/md "$0+0i-1j+0k$")]
         ['-K (clerk/md "$0+0i+0j-1k$")]
         ['ONE (clerk/md "$1+0i+0j+0k$")]
         ['ZERO (clerk/md "$0+0i+0j+0k$")]]})

^{::clerk/visibility :hide}
(clerk/example
  q/I q/J q/K q/-I q/-J q/-K q/ONE q/ZERO)

;; ## Functions

;; ### Creation

;; To create quaternion call `quaternion` or `fastmath.vector/vec4`. Arguments for `quaternion`:
;;
;; - 4 arguments - all four parts
;; - 2 arguments - scalar and vector parts
;; - 1 argument - only scalar part, returns quaternion representation of real number

^{::clerk/visibility :hide}
(u/table2
 [[quaternion "Create quaternion"]
  [fastmath.vector/vec4 "Create quaternion"]
  [complex->quaternion "Convert complex number into quaternion"]])

^{::clerk/visibility :hide}
(clerk/example
  (q/quaternion 1 2 3 4)
  (q/quaternion 1 [2 3 4])
  (q/quaternion 1)
  (fastmath.vector/vec4 1 2 3 4)
  (class (q/quaternion 1 2 3 4))
  (q/complex->quaternion (c/complex 1 2)))

;; ### Fields

^{::clerk/visibility :hide}
(u/table2
 [[scalar "Scalar / real part"]
  [re "Scalar / real part"]
  [vector "Vector part"]
  [im-i "i imaginary part"]
  [im-j "j imaginary part"]
  [im-k "k imaginary part"]])

(def quaternion-3-4i+2j-1k (q/quaternion 3 -4 2 -1))

^{::clerk/visibility :hide}
(clerk/example
  (q/scalar quaternion-3-4i+2j-1k)
  (q/re quaternion-3-4i+2j-1k)
  (q/vector quaternion-3-4i+2j-1k)
  (q/im-i quaternion-3-4i+2j-1k)
  (q/im-j quaternion-3-4i+2j-1k)
  (q/im-k quaternion-3-4i+2j-1k)
  (quaternion-3-4i+2j-1k 0)
  (quaternion-3-4i+2j-1k 1)
  (quaternion-3-4i+2j-1k 2)
  (quaternion-3-4i+2j-1k 3))

;; ### Predicates

^{::clerk/visibility :hide}
(u/table2
 [[real? "Is q a real number? ie. imaginary part = 0.0"]
  [imaginary? "Is q real part zero?"]
  [zero? "Is q is zero?"]
  [inf? "Is any part infinite?"]
  [nan? "Is any part NAN?"]])

^{::clerk/visibility :hide}
(clerk/example
  (q/real? q/ONE)
  (q/real? quaternion-3-4i+2j-1k)
  (q/imaginary? q/K)
  (q/imaginary? quaternion-3-4i+2j-1k)
  (q/zero? q/ZERO)
  (q/inf? (q/quaternion 0.0 0.0 0.0 ##Inf))
  (q/nan? (q/quaternion ##NaN 0.0 0.0 0.0)))

;; ### Basic operations

^{::clerk/visibility :hide}
(u/table2
 [[norm "Norm of the quaternion, length of the vector"]
  [arg "Argument of the quaternion, atan2(|vector(q)|,scalar(q))"]
  [normalize "Normalize quaternion"]
  [conjugate "Conjugate"]
  [neg "Negation"]
  [add "Sum of two quaternions"]
  [sub "Difference of two quaternions"]
  [mult "Multiplication of two quaternions"]
  [scale "Multiply by a scalar"]
  [div "Division of two quaternions"]
  [reciprocal "Reciprocal of the quaternion"]
  [qsgn "Return 0 for zero quaternion or sign of real part otherwise"]
  [delta-eq "Compare two quaternions"]
  [slerp "Interpolate quaternions, spherical linear interpolation"]])

^{::clerk/visibility :hide}
(clerk/example
  (q/norm quaternion-3-4i+2j-1k)
  (q/arg quaternion-3-4i+2j-1k)
  (q/conjugate quaternion-3-4i+2j-1k)
  (q/neg quaternion-3-4i+2j-1k)
  (q/add quaternion-3-4i+2j-1k quaternion-3-4i+2j-1k)
  (q/sub q/ONE quaternion-3-4i+2j-1k)
  (q/mult quaternion-3-4i+2j-1k quaternion-3-4i+2j-1k)
  (q/div (q/quaternion 0 1 2 3) quaternion-3-4i+2j-1k)
  (q/reciprocal quaternion-3-4i+2j-1k)
  (q/qsgn quaternion-3-4i+2j-1k)
  (q/qsgn q/ZERO)
  (q/slerp quaternion-3-4i+2j-1k q/K 0.001)
  (q/slerp quaternion-3-4i+2j-1k q/K 0.5)
  (q/slerp quaternion-3-4i+2j-1k q/K 0.999))

;; To compare two complex numbers with some accuracy (default `1.0e-6`), call `delta-eq`.

^{::clerk/visibility :hide}
(clerk/example
  (q/delta-eq (q/quaternion 1.0001 0.0 0.0 0.0001) q/ONE)
  (q/delta-eq (q/quaternion 1.0001 0.0 0.0 0.0001) q/ONE 1.0e-3))

;; ### Trigonometric

;; All trigonometric and hyperbolic (with their inversions) functions are derived from the complex operation where given function `f` is applied  for $z_{in} = r + |\vec{v}|i$ complex number, where $r$ is scalar and $\vec{v}$ is vector of the quaternion, and $|\vec{v}|$ is length of the vector $\vec{v}$. So, $z_{res} = f(z_{in})$. Then resulting quaternion has a form of $q = \Re(z_{res}) + \frac{\Im(z_{res})}{|\vec{v}|}\vec{v}$.

;; There results are the same as described in [SO](https://math.stackexchange.com/questions/1499095/how-to-calculate-sin-cos-tan-of-a-quaternion)

;; #### Basic

^{::clerk/visibility :hide}
(u/table2
 [[sin (clerk/md "$\\sin(q)$")]
  [cos (clerk/md "$\\cos(q)$")]
  [tan (clerk/md "$\\tan(q)$")]
  [sec (clerk/md "$\\sec(q)$")]
  [csc (clerk/md "$\\csc(q)$")]
  [cot (clerk/md "$\\cot(q)$")]
  [asin (clerk/md "$\\arcsin(q)$")]
  [acos (clerk/md "$\\arccos(q)$")]
  [atan (clerk/md "$\\arctan(q)$")]
  [asec (clerk/md "$\\operatorname{arcsec}(q)$")]
  [acsc (clerk/md "$\\operatorname{arccsc}(q)$")]
  [acot (clerk/md "$\\operatorname{arccot}(q)$")]])

^{::clerk/visibility :hide}
(clerk/example
  (q/sin  quaternion-3-4i+2j-1k)
  (q/cos  quaternion-3-4i+2j-1k)
  (q/tan  quaternion-3-4i+2j-1k)
  (q/sec  quaternion-3-4i+2j-1k)
  (q/csc  quaternion-3-4i+2j-1k)
  (q/cot  quaternion-3-4i+2j-1k)
  (q/asin quaternion-3-4i+2j-1k)
  (q/acos quaternion-3-4i+2j-1k)
  (q/atan quaternion-3-4i+2j-1k)
  (q/asec quaternion-3-4i+2j-1k)
  (q/acsc quaternion-3-4i+2j-1k)
  (q/acot quaternion-3-4i+2j-1k))

;; #### Hyperbolic

^{::clerk/visibility :hide}
(u/table2
 [[sinh (clerk/md "$\\operatorname{sinh}(q)$")]
  [cosh (clerk/md "$\\operatorname{cosh}(q)$")]
  [tanh (clerk/md "$\\operatorname{tanh}(q)$")]
  [sech (clerk/md "$\\operatorname{sech}(q)$")]
  [csch (clerk/md "$\\operatorname{csch}(q)$")]
  [coth (clerk/md "$\\operatorname{coth}(q)$")]
  [asinh (clerk/md "$\\operatorname{arcsinh}(q)$")]
  [acosh (clerk/md "$\\operatorname{arccosh}(q)$")]
  [atanh (clerk/md "$\\operatorname{arctanh}(q)$")]
  [asech (clerk/md "$\\operatorname{arcsec}(q)$")]
  [acsch (clerk/md "$\\operatorname{arccsc}(q)$")]
  [acoth (clerk/md "$\\operatorname{arccot}(q)$")]])

^{::clerk/visibility :hide}
(clerk/example
  (q/sinh  quaternion-3-4i+2j-1k)
  (q/cosh  quaternion-3-4i+2j-1k)
  (q/tanh  quaternion-3-4i+2j-1k)
  (q/sech  quaternion-3-4i+2j-1k)
  (q/csch  quaternion-3-4i+2j-1k)
  (q/coth  quaternion-3-4i+2j-1k)
  (q/asinh quaternion-3-4i+2j-1k)
  (q/acosh quaternion-3-4i+2j-1k)
  (q/atanh quaternion-3-4i+2j-1k)
  (q/asech quaternion-3-4i+2j-1k)
  (q/acsch quaternion-3-4i+2j-1k)
  (q/acoth quaternion-3-4i+2j-1k))

;; ### Power and logarithm

^{::clerk/visibility :hide}
(u/table2
 [[sq (clerk/md "$q^2$")]
  [sqrt (clerk/md "$\\sqrt{q}$")]
  [pow (clerk/md "$q^c$")]
  [exp (clerk/md "$e^q$")]
  [log (clerk/md "$\\operatorname{Log}(q)$, principal value")]
  [logb (clerk/md "$\\log_b(q)$")]])

^{::clerk/visibility :hide}
(clerk/example
  (q/sq quaternion-3-4i+2j-1k)
  (q/sqrt quaternion-3-4i+2j-1k)
  (q/pow quaternion-3-4i+2j-1k quaternion-3-4i+2j-1k)
  (q/exp quaternion-3-4i+2j-1k)
  (q/log quaternion-3-4i+2j-1k)
  (q/logb quaternion-3-4i+2j-1k q/K))

;; ### Rotations

;; 3d rotations can be treated in many ways, two options are possible:
;; * Euler angles, like presented in [wiki](https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles)
;; * Tait-Bryan angles, to achieve the same result as OpenGL rotations

^{::clerk/visibility :hide}
(u/table2
 [[rotation-quaternion "Create rotation quaternion around vector u and angle alpha"]
  [rotate "Rotate 3d `in` vector around axis `u` (the same as `fastmath.vector/axis-rotate`)"]
  [to-euler "Convert quaternion to Euler angles (ZYX, body 3-2-1), output: [roll, pitch, yaw]"]
  [from-euler "Convert Euler angles ZYX representation, [roll, pitch, yaw] to quaternion"]
  [to-angles "Convert quaternion to Tait–Bryan (z-y′-x\") angles, output: [x, y, z]"]
  [from-angles "Convert Tait-Bryan (z-y′-x\") angles to quaternion"]
  [to-rotation-matrix "Convert quaternion to a rotation matrix"]
  [from-rotation-matrix "Convert rotation matrix to a quaternion"]])

^{::clerk/visibility :hide}
(clerk/example
  (q/rotation-quaternion m/HALF_PI [1 1 0])
  (q/rotate [1 1 1] m/HALF_PI [1 1 0])
  (v/axis-rotate (v/vec3 1 1 1) m/HALF_PI (v/vec3 1 1 0))
  (q/to-euler (q/quaternion m/SQRT2_2 0 m/SQRT2_2 0))
  (q/from-euler 0 m/HALF_PI 0)
  (q/to-angles (q/quaternion m/SQRT2_2 0 m/SQRT2_2 0))
  (q/from-angles 0 m/HALF_PI 0)
  (q/rotate [1 1 1] (q/from-euler 0.1 0.2 0.3))
  (q/rotate [1 1 1] (q/from-angles 0.1 0.2 0.3))
  (v/rotate (v/vec3 1 1 1) 0.1 0.2 0.3)
  (q/from-angles 0.1 0.2 0.3)
  (q/to-rotation-matrix (q/from-angles 0.1 0.2 0.3))
  (q/from-rotation-matrix (mat/rotation-matrix-3d 0.1 0.2 0.3)))

;; ## List of symbols

^{::clerk/visibility :hide ::clerk/no-cache true}
(u/make-public-fns-table 'fastmath.quaternion)
